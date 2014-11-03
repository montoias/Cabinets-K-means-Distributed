#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#define BUFFER_SIZE 20000
#define CHAR_BUFFER 20
#define FILENAME_BUFFER 500
#define ROOT 0
#define CHUNK_MSG 1
#define DOC_SUBS_RESULT 2
#define CHUNK_SIZE_MSG 3

/* Macro that yields 1 if a process has an extra doc */
#define HAS_EXTRA(proc, num_procs, num_docs) ((proc%num_procs >= num_procs - num_docs%num_procs) ? 1 : 0)

/* Global Variables */
int num_subs, num_cabs, num_docs, my_docs, num_procs, rank;
double *averages, *new_averages, *temp_new_averages;
double **distance, **doc_subjects;
int *doc_index, *cab_docs, *new_num_docs, *temp_cab_docs;
char *doc_chunk;

/* debug */
double initializeTime;

/* Function that sends three variables from process 0 to the other processes:
	num_cabs, num_docs, num_subs */
void shareInitializationValues(){
	int info[3];
	info[0] = num_cabs;
	info[1] = num_docs;
	info[2] = num_subs;
		
	MPI_Bcast(info, 3, MPI_INT, ROOT, MPI_COMM_WORLD);
}

/* Function that receives three variables from process:
	num_cabs, num_docs, num_subs */
void receiveInitializationValues(){
	int info[3];
	MPI_Bcast(info, 3, MPI_INT, ROOT, MPI_COMM_WORLD);
		
	num_cabs = info[0];
	num_docs = info[1];
	num_subs = info[2];
}

/* Function that allocates a matrix of doubles */
double **allocateDoubleMatrix(int num_lines, int num_columns){
	int line_i;	
	double **matrix = (double**) malloc(sizeof(double*)*num_lines);
	
	for(line_i = 0; line_i < num_lines; line_i++)
		matrix[line_i] = calloc(num_columns, sizeof(double));

	return matrix;
}

/* Function that allocates and initializes the data structures	*/
void initializeStructures(){
	averages = (double*) calloc(num_subs * num_cabs, sizeof(double));
	new_averages = (double*) calloc(num_subs * num_cabs, sizeof(double));
	new_num_docs = (int*) calloc(num_cabs, sizeof(int));
	
	my_docs = num_docs/num_procs + HAS_EXTRA(rank, num_procs, num_docs);
	doc_index = (int*) calloc(my_docs, sizeof(int));
	distance = allocateDoubleMatrix(my_docs, num_cabs);
	doc_subjects = allocateDoubleMatrix(my_docs, num_subs);
	
	if(rank == ROOT){
		/* NOTE: adding one more doc -> worst case while distributing docs */
		doc_chunk = (char*) malloc((my_docs+1)*num_subs*CHAR_BUFFER);		
		temp_cab_docs = (int*) calloc(num_cabs, sizeof(int));
		temp_new_averages = (double*) calloc(num_cabs * num_subs, sizeof(double));
		cab_docs = (int*) calloc(num_cabs , sizeof(int));
	} 
}

/* Function that reads from the file and sends the information to every
	process */
void sendFileChunks(FILE *input_file){
	int proc_i, doc_i, proc_docs, offset, param_len;
	int chunk = num_docs/num_procs;
	char parameters[BUFFER_SIZE];
	
	for(proc_i = 1; proc_i < num_procs+1; proc_i++){
		offset = 0;
		proc_docs = chunk + HAS_EXTRA(proc_i, num_procs, num_docs);
		for(doc_i = 0; doc_i < proc_docs; doc_i++){
			fgets(parameters, sizeof(parameters), input_file);
			param_len = strlen(parameters);
			memcpy(doc_chunk + offset, parameters, param_len);
			offset += param_len;
		}
		if(proc_i != num_procs){
			MPI_Send(&offset, 1, MPI_INT, proc_i, CHUNK_SIZE_MSG, MPI_COMM_WORLD);
			MPI_Send(doc_chunk, offset, MPI_CHAR, proc_i, CHUNK_MSG, MPI_COMM_WORLD);
		}
	}
	
	fclose(input_file);
}

/* Function that reads information received from process 0 and stores 
	it in its proper structures                                             */
void readAndStore(char *doc_chunk){
	int doc_i;
	char *line, *line_tok;
	
	line = strtok_r(doc_chunk, "\n", &line_tok);
	for(doc_i = 0; doc_i < my_docs; doc_i++){
		int sub_i, doc_id, cab_id, offset;
		char *token, *token_tok;
		
		token = strtok_r(line, " ", &token_tok);
		
		doc_id = atoi(token);
		cab_id = doc_id%num_cabs;
		doc_index[doc_i] = cab_id;
		offset = num_subs * cab_id;
		
		for(sub_i = 0; sub_i < num_subs; sub_i++){
			token = strtok_r(NULL, " ", &token_tok);
			doc_subjects[doc_i][sub_i] = atof(token);
			new_averages[offset++] += doc_subjects[doc_i][sub_i];
		}

		new_num_docs[cab_id]++;
		line = strtok_r(NULL, "\n", &line_tok);
	}
	
	free(doc_chunk);	
}

/* Function that updates the cabinets */
void updateAverages(){
	int cab_i, sub_i, prev_num_docs, new_cab_docs;
	
	MPI_Reduce(new_num_docs, temp_cab_docs, num_cabs, MPI_INT, MPI_SUM , ROOT, MPI_COMM_WORLD);
	MPI_Reduce(new_averages, temp_new_averages, num_cabs * num_subs, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	

	if(rank == ROOT){
		for(cab_i = 0; cab_i < num_cabs; cab_i++){	
			int offset = cab_i * num_subs;
			prev_num_docs = cab_docs[cab_i];
			new_cab_docs = prev_num_docs + temp_cab_docs[cab_i];
			
			for(sub_i = 0; sub_i < num_subs; sub_i++){
				if(new_cab_docs == 0)
					averages[offset++] = 0;
				else {				
					averages[offset] *= prev_num_docs;
					averages[offset] += temp_new_averages[offset];
					averages[offset++] /= new_cab_docs;
				}
			}
		
			cab_docs[cab_i] = new_cab_docs;
		}	
	}

	memset(new_num_docs, 0, sizeof(int) * num_cabs);
	memset(new_averages, 0, sizeof(double) * num_cabs * num_subs);

	MPI_Bcast(averages, num_cabs * num_subs, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);	
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
void calculateDistances(){
	int cab_i, doc_i;

	for(doc_i = 0; doc_i < my_docs; doc_i++){
		for(cab_i = 0; cab_i < num_cabs; cab_i++){
			distance[doc_i][cab_i] = calculateDistance(doc_subjects[doc_i], &averages[cab_i * num_subs]);
		}
	}
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
	int doc_i, moved_flag_aux = 0, moved_flag = 0;	

	for(doc_i = 0; doc_i < my_docs; doc_i++){
		int current_cab = doc_index[doc_i];
		int closest_cab = findMinDistance(distance[doc_i], current_cab);
		int cur_cab_offset = current_cab*num_subs;
		int clo_cab_offset = closest_cab*num_subs;
		
		if(current_cab != closest_cab){
			int sub_i;
			moved_flag_aux = 1;
			
			for(sub_i = 0; sub_i < num_subs; sub_i++){
				new_averages[cur_cab_offset++] -= doc_subjects[doc_i][sub_i];
				new_averages[clo_cab_offset++] += doc_subjects[doc_i][sub_i];
			}

			new_num_docs[current_cab]--;
			new_num_docs[closest_cab]++;
			doc_index[doc_i] = closest_cab;
		}
	}
	
	MPI_Allreduce(&moved_flag_aux, &moved_flag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

	return moved_flag;
}


/* Function that writes the position of each document to a file */
void writeToFile(char *input_filename){
	int doc_i, proc_i;
	FILE *output_file;
	MPI_Status status;
	int startDoc = 0;
	char output_filename[FILENAME_BUFFER];
	int chunk = num_docs/num_procs;
	int *temp_doc_index = malloc(sizeof(int) * (my_docs+1));
	
	memcpy(output_filename, input_filename, strlen(input_filename) - 3);
	output_filename[strlen(input_filename) - 3] = '\0';

	strcat(output_filename, ".out");
	output_file = fopen(output_filename, "w+");

	for(proc_i = 1; proc_i < num_procs; proc_i++){
		int proc_docs = chunk + HAS_EXTRA(proc_i, num_procs, num_docs);
		
		MPI_Recv(temp_doc_index, proc_docs , MPI_INT, proc_i, DOC_SUBS_RESULT, MPI_COMM_WORLD, &status);
		for(doc_i = 0; doc_i < proc_docs; doc_i++)
			fprintf(output_file, "%d %d\n", startDoc++, temp_doc_index[doc_i]);
			
	}
	
	for(doc_i = 0; doc_i < my_docs; doc_i++)
		fprintf(output_file, "%d %d\n", startDoc++, doc_index[doc_i]);

	free(temp_doc_index);
	fclose(output_file);
}

/* Function that frees a matrix of doubles */
void freeDoubleMatrix(double **matrix, int num_lines){
	int line_i;	
	
	for(line_i = 0; line_i < num_lines; line_i++)
		free(matrix[line_i]);

	free(matrix);
}

/* Function that frees the allocated structures along the execution of the program */
void cleanup(){
	if(rank == ROOT){		
		free(temp_cab_docs);
		free(temp_new_averages);
		free(cab_docs);
	} 
	
	freeDoubleMatrix(distance, my_docs);
	freeDoubleMatrix(doc_subjects, my_docs);
		
	free(averages);
	free(new_averages);
	free(new_num_docs);
	free(doc_index);
}

int main(int argc, char *argv[]){
	int temp_cabs;
	FILE *input_file;
	MPI_Status status;
	char *input_filename = (char*) malloc(sizeof(char)*(strlen(argv[1])+2));
	
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
	
	initializeTime = MPI_Wtime();

	if(rank == ROOT){
		strcpy(input_filename, argv[1]); 
		input_file = fopen(input_filename, "r");
		
		if(input_file == NULL){
			printf("ERROR reading file!");
			perror(input_filename);
			MPI_Finalize();
			return -1;
		}
		
		fscanf(input_file, "%d %d %d\n", &temp_cabs, &num_docs, &num_subs);
		num_cabs = (argv[2] != NULL) ? atoi(argv[2]) : temp_cabs;
			
		shareInitializationValues();
	}
	else 
		receiveInitializationValues();
		
	initializeStructures();

	if(rank == ROOT)
		sendFileChunks(input_file);
	else {
		int docChunkSize;
		MPI_Recv(&docChunkSize, 1, MPI_INT, ROOT, CHUNK_SIZE_MSG, MPI_COMM_WORLD, &status);
		doc_chunk = (char*) malloc(docChunkSize);	
		MPI_Recv(doc_chunk, docChunkSize, MPI_CHAR, ROOT, CHUNK_MSG, MPI_COMM_WORLD, &status);
	}	

	readAndStore(doc_chunk);	
	
 	printf("Initialization ended after :%f\n", MPI_Wtime() - initializeTime);
	
	while(1){
		updateAverages();
		
		calculateDistances();	
		
		if(!changeDocuments())
			break;
	}
		
	if(rank == ROOT)
		writeToFile(input_filename);
	else
		MPI_Send(doc_index, my_docs, MPI_INT, ROOT, DOC_SUBS_RESULT, MPI_COMM_WORLD);

	printf("Program ended after :%f\n", MPI_Wtime() - initializeTime);

	cleanup();	
	if(rank == ROOT)
		free(input_filename);
		
	MPI_Finalize();
	return 0;
}

