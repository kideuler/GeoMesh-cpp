CC = g++
CFLAGS = -Iinclude -I src/
DEPS = Meshgen2d.cpp MeshSmooth.cpp GeoSplines.cpp MeshBlossom.cpp Meshutils.cpp SurfaceMesh.cpp
DEPS_O = Meshgen2d.o MeshSmooth.o GeoSplines.o MeshBlossom.o Meshutils.o SurfaceMesh.o
OBJ = main

make:
	$(CC) $(CFLAGS) -c -g  $(DEPS)

.PHONY: clean
clean:
	rm -rf *.o

.PHONY: run
run:
	$(CC) $(CFLAGS)  -g main2.cpp $(DEPS_O) -o $(OBJ)
	./$(OBJ)
	sudo cp test.vtk /mnt/c/Users/Jacob/Documents/test.vtk
	sudo cp bunny.obj /mnt/c/Users/Jacob/Documents/bunny.obj
