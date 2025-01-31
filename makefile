CC = gcc
FLAGS   = -std=c99 -c -fPIC -g -Wall -include lib/qhull/src/libqhull_r/qhull_ra.h
LFLAGS  = -L lib/qhull/lib -lqhull_r -lqhullstatic_r -lqhullstatic -lm
OUTPUT_SO  = bin/libgriddata.so
OUTPUT  = udi.exe

OBJECTS = \
	build/interpolation.o \
	build/delaunator.o \
	build/main.o

all: $(OBJECTS)
	mkdir -p bin
	mkdir -p build
	$(CC) -g $(OBJECTS) -o $(OUTPUT) $(LFLAGS)
	$(CC) -g $(OBJECTS) -shared -o $(OUTPUT_SO) $(LFLAGS)

build/main.o: src/main.c
	mkdir -p bin
	mkdir -p build
	$(CC) $(FLAGS) src/main.c -o build/main.o

build/interpolation.o: src/interpolation.c
	mkdir -p bin
	mkdir -p build
	$(CC) $(FLAGS) src/interpolation.c -o build/interpolation.o

build/delaunator.o: src/delaunator.c
	mkdir -p bin
	mkdir -p build
	$(CC) $(FLAGS) src/delaunator.c -o build/delaunator.o

clean:
	rm -Rf $(OUTPUT) $(OBJECTS) $(OUTPUT_SO) build



