CC = g++
CFLAGS = -I -Isrc  -Iinclude -Iextern/eigen
DEPS = src/Meshgen2d.cpp src/Meshutils.cpp src/GeoSpline.cpp
DEPS_O = Meshgen2d.o Meshutils.o GeoSpline.o
OBJ = test

make:
	$(CC) $(CFLAGS) -c -g  $(DEPS)

.PHONY: clean
clean:
	rm -rf *.o

.PHONY: opt
opt:
	$(CC) $(CFLAGS) -c -O3 $(DEPS)

.PHONY: run
run:
	$(CC) $(CFLAGS)  -g test2.cpp $(DEPS_O) -o $(OBJ)
	./$(OBJ)
	sudo cp test.vtk /mnt/c/Users/Jacob/Documents/test.vtk
