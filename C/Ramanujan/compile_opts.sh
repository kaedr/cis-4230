gcc -Wall -o noopt_pi pi.c -lgmp -lmpfr -lm
gcc -Wall -O1 -o 1opt_pi pi.c -lgmp -lmpfr -lm
gcc -Wall -O2 -o 2opt_pi pi.c -lgmp -lmpfr -lm
gcc -Wall -O3 -o 3opt_pi pi.c -lgmp -lmpfr -lm
