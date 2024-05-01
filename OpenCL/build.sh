gcc -o $1.exe -I/usr/local/cuda/include -L/usr/local/cuda/lib64/ -lm Timer.c $1.c -lOpenCL
