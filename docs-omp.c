#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

#define MIN_DOCS 500


typedef struct cabinet{
	int prev_num_docs;
	unsigned int num_docs;
	double *averages;
	double *new_averages;
} cabinet;

int num_subs, num_cabs, num_docs;
double **distance;
double **doc_subjects;			/* Matrix that maps the subjects to their documents */
cabinet *cabinets;		
int *doc_index;				/* Structure that keeps track of all documents */
int *modified;
omp_lock_t *cab_lock;			

/* Function that allocates and initializes the cabinets */
void create_cabinets(){
	int cab_i;

	cab_lock = (omp_lock_t *)malloc(sizeof(omp_lock_t) * num_cabs);
	cabinets = (cabinet*) malloc(sizeof(cabinet) * num_cabs);
	modified = (int*) calloc(num_cabs, sizeof(int));
	
	for(cab_i = 0; cab_i < num_cabs; cab_i++){
		cabinets[cab_i].num_docs = 0;
		cabinets[cab_i].prev_num_docs = 0;
		cabinets[cab_i].averages = (double*) calloc(num_subs, sizeof(double));
		cabinets[cab_i].new_averages = (double*) calloc(num_subs, sizeof(double));	
		modified[cab_i] = 1;
		omp_init_lock(&cab_lock[cab_i]);
	}
}

/* Function that reads information from the file and stores it in its 
   proper structures.                                                 */
void readAndStore(FILE *input_file){
	int doc_i;

	#pragma omp parallel for private(doc_i) if(num_docs > MIN_DOCS)
	for(doc_i = 0; doc_i < num_docs; doc_i++){
		int sub_i, doc_id, cab_id;
		char *save_ptr, *token, parameters[20000];
		
		fgets(parameters, sizeof(parameters), input_file);
		token = strtok_r(parameters, " ", &save_ptr);
		doc_id = atoi(token);
		cab_id = doc_id%num_cabs;
		doc_index[doc_id] = cab_id;
		
		for(sub_i = 0; sub_i < num_subs; sub_i++){
			token = strtok_r(NULL, " ", &save_ptr);
			doc_subjects[doc_id][sub_i] = atof(token);
		}
		
		omp_set_lock(&cab_lock[cab_id]);
		cabinets[cab_id].num_docs++;
		omp_unset_lock(&cab_lock[cab_id]);
	}
	fclose(input_file);
}

/* Function that initializes the averages for each cabinet */
void initializeAverages(){
	int cab_i, sub_i, doc_i, n_docs;

	for(doc_i = 0; doc_i < num_docs; doc_i++){
		cab_i = doc_index[doc_i];
		n_docs = cabinets[cab_i].num_docs;
		
		if(n_docs != 0){
			for(sub_i = 0; sub_i < num_subs; sub_i++)
				cabinets[cab_i].averages[sub_i] += (doc_subjects[doc_i][sub_i] / n_docs);		
		}
	}
		
}

/* Function that updates the cabinets */
void updateAverages(){
	int cab_i, sub_i;
	
	#pragma omp parallel for private (sub_i) schedule (static, 1) if(num_cabs >= omp_get_num_threads())
	for(cab_i = 0; cab_i < num_cabs; cab_i++){
		if(modified[cab_i]){
			int prev_num_docs = cabinets[cab_i].num_docs;
			int num_docs = prev_num_docs + cabinets[cab_i].prev_num_docs;
			for(sub_i = 0; sub_i < num_subs; sub_i++){
				if(num_docs == 0)
					cabinets[cab_i].averages[sub_i] = 0;
				else{				
					cabinets[cab_i].averages[sub_i] *= prev_num_docs;
					cabinets[cab_i].averages[sub_i] += cabinets[cab_i].new_averages[sub_i];
					cabinets[cab_i].averages[sub_i] /= num_docs;			
				}

				cabinets[cab_i].new_averages[sub_i] = 0;
			}
			cabinets[cab_i].num_docs = num_docs;
			cabinets[cab_i].prev_num_docs = 0;
		}
	}
}

/* Function that calculates the distance between a document and a cabinet */
double calculateDistance(double *subjects, double *averages){
	int sub_i;
	double current_distance = 0;
	
	for(sub_i = 0; sub_i < num_subs; sub_i++){
		double subtract = subjects[sub_i] - averages[sub_i];
		current_distance += (subtract * subtract);
	}

	return current_distance;
}

/* Function that calculates all the distances between documents and cabinets */
void calculateDistances(double** distance){
	int cab_i, doc_i;

	#pragma omp parallel for private(cab_i) if(num_docs > MIN_DOCS)
	for(doc_i = 0; doc_i < num_docs; doc_i++){
		for(cab_i = 0; cab_i < num_cabs; cab_i++)
			if(modified[cab_i])
				distance[doc_i][cab_i] = calculateDistance(doc_subjects[doc_i], cabinets[cab_i].averages);
	}
	memset(modified, 0 , sizeof(int)*num_cabs);
}

/* Function that finds the minimum distance between a document 
   and a cabinet.                                                   */
int findMinDistance(double* distances, int cabinet_id){
	int cab_i;
	double min_distance = distances[cabinet_id];
	
	for(cab_i = 0; cab_i < num_cabs; cab_i++){
		double current_distance = distances[cab_i];
		
		if(current_distance < min_distance){
			min_distance = current_distance;
			cabinet_id = cab_i;
		}
	}
	
	return cabinet_id;
}

/* Function that moves documents from one cabinet to another */
int changeDocuments(){
	int doc_i, moved_flag = 0;
	
	#pragma omp parallel for if(num_docs > MIN_DOCS)
	for(doc_i = 0; doc_i < num_docs; doc_i++){
		int id, cabs_i = doc_index[doc_i];
		
		id = findMinDistance(distance[doc_i], cabs_i);
		
		if(id != cabs_i){
			int sub_i;
			moved_flag = 1;
			
			omp_set_lock(&cab_lock[cabs_i]);
			for(sub_i = 0; sub_i < num_subs; sub_i++){
				cabinets[cabs_i].new_averages[sub_i] -=	doc_subjects[doc_i][sub_i];
			}
			cabinets[cabs_i].prev_num_docs--;
			modified[cabs_i] = 1;
			omp_unset_lock(&cab_lock[cabs_i]);
			
			omp_set_lock(&cab_lock[id]);
			for(sub_i = 0; sub_i < num_subs; sub_i++){
				cabinets[id].new_averages[sub_i] += doc_subjects[doc_i][sub_i];
			}
			cabinets[id].prev_num_docs++;
			modified[id] = 1;
			omp_unset_lock(&cab_lock[id]);
			
			doc_index[doc_i] = id;
		}
	}
	
	return moved_flag;
}

/* Function that writes the position of each document to a file */
void writeToFile(char *input_filename){
	int doc_i;
	FILE *output_file;
	
	char output_filename[500];
	memcpy( output_filename, input_filename, strlen(input_filename) - 3 );
	output_filename[strlen(input_filename) - 3] = '\0';

	strcat(output_filename, ".out");
	output_file = fopen(output_filename, "w+");

	for(doc_i = 0; doc_i < num_docs; doc_i++)
		fprintf(output_file, "%d %d\n", doc_i, doc_index[doc_i]);

	fclose(output_file);
}

/* Function that allocates a matrix of doubles */
double **allocateDoubleMatrix(int num_lines, int num_columns){
	int line_i;	
	double **matrix = (double**) malloc(sizeof(double*)*num_lines);
	
	for(line_i = 0; line_i < num_lines; line_i++)
		matrix[line_i] = calloc(num_columns, sizeof(double));

	return matrix;
}

/* Function that frees a matrix of doubles */
void freeDoubleMatrix(double **matrix, int num_lines){
	int line_i;	
	
	for(line_i = 0; line_i < num_lines; line_i++)
		free(matrix[line_i]);

	free(matrix);
}

/* Function that frees the allocated structures along the program */
void cleanup(){
	int cab_i;

	for(cab_i = 0; cab_i < num_cabs; cab_i++){
		free(cabinets[cab_i].averages);
		free(cabinets[cab_i].new_averages);
		omp_destroy_lock(&cab_lock[cab_i]);
	}

	free(cab_lock);
	free(doc_index);
	free(modified);
	free(cabinets);
}

int main(int argc, char *argv[]){
	
	
	double start = omp_get_wtime(), algorithm;
	int temp_cabs, moved_flag = 1;
	FILE *input_file;
	char *input_filename = (char*)malloc(sizeof(char)*(strlen(argv[1])+2));

	strcpy(input_filename, argv[1]);
	input_file = fopen(input_filename, "r");
	fscanf(input_file, "%d %d %d\n", &temp_cabs, &num_docs, &num_subs);

	if(argv[2] != NULL)
		num_cabs = atoi(argv[2]);
	else 
		num_cabs = temp_cabs;
	
	create_cabinets();
	doc_index = (int*) malloc(sizeof(int)*num_docs);
	distance = allocateDoubleMatrix(num_docs, num_cabs);
	doc_subjects = allocateDoubleMatrix(num_docs, num_subs);

	if(input_file != NULL)
		readAndStore(input_file);
	else 
		perror(input_filename);
		
	algorithm = omp_get_wtime();
	initializeAverages();

	while(moved_flag){
		
		updateAverages();
		calculateDistances(distance);
		moved_flag = changeDocuments();
	}	

	writeToFile(input_filename);
	cleanup();	
	freeDoubleMatrix(distance, num_docs);
	freeDoubleMatrix(doc_subjects, num_docs);
	free(input_filename);
	
	
	printf("Algorithm Time: %f \n", omp_get_wtime() - algorithm);
	printf("Elapsed Time: %f \n", omp_get_wtime() - start);
	
	return 0;
}

