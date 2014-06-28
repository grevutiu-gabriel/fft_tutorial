#valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes 
CC=gcc

SOURCE = main.c

FLAGS=-Wall -lm -g -O3 -Wno-long-long -fno-builtin -lfftw3f -lfftw3
TARGET=fftw_demo

all:
	$(CC) $(SOURCE) $(FLAGS) -o $(TARGET)

clean:
