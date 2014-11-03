serial:
	gcc -ansi -pedantic -Wall ../Proj/docs-serial.c -o docs-serial
	mv docs-serial resultsSerial/

parallel:
	gcc -ansi -pedantic -Wall -posix -fopenmp ../Proj/docs-parallel.c -o docs-parallel
	mv docs-parallel resultsParallel/

debug-serial:
	gcc -ansi -pedantic -Wall -g docs-serial.c -o docs-serial

debug-parallel:
	gcc -ansi -pedantic -Wall -posix -fopenmp -g docs-parallel.c -o docs-serial

profile-parallel:
	'/home/paulo/ompp/bin/kinst-ompp' gcc -ansi -pedantic -posix -fopenmp docs-parallel.c -o docs-serial


clean:
	rm docs-serial
	rm tests/*d.out
