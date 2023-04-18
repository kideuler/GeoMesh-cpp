CC = g++
CFLAGS = -I -Isrc  -Iinclude -Iextern/eigen
DEPS = src/*.cpp
DEPS_O = Meshgen2d.o Meshutils.o GeoSpline.o MeshSmooth.o MeshSurface.o FEM.o FacetColor.o kdTree.o
OBJ = test
OBJDIR = bin

make:
	mpic++ $(CFLAGS) -c -g $(DEPS)
.PHONY: clean
clean:
	rm -rf *.o *.vtk

.PHONY: debug
debug:
	$(CC) $(CFLAGS) -c -g $(DEPS)

.PHONY: multi
multi:
	mpic++ $(CFLAGS) -c -g $(DEPS)

# test 1
.PHONY: 1
1:
	$(CC) $(CFLAGS)  -g test.cpp $(DEPS_O) -o $(OBJ)
	./$(OBJ) 1
	sudo cp test1.vtk /mnt/c/Users/Jacob/Documents/test1.vtk
# test 2
.PHONY: 2
2:
	$(CC) $(CFLAGS)  -g test.cpp $(DEPS_O) -o $(OBJ)
	./$(OBJ) 2
	sudo cp test2.vtk /mnt/c/Users/Jacob/Documents/test2.vtk
# test 3
.PHONY: 3
3:
	$(CC) $(CFLAGS)  -g test.cpp $(DEPS_O) -o $(OBJ)
	./$(OBJ) 3
	sudo cp test3.vtk /mnt/c/Users/Jacob/Documents/test3.vtk
# test 4
.PHONY: 4
4:
	$(CC) $(CFLAGS)  -g test.cpp $(DEPS_O) -o $(OBJ)
	./$(OBJ) 4
	sudo cp test4.vtk /mnt/c/Users/Jacob/Documents/test4.vtk
# test 5
.PHONY: 5
5:
	$(CC) $(CFLAGS)  -g test.cpp $(DEPS_O) -o $(OBJ)
	./$(OBJ) 5
	sudo cp test5.vtk /mnt/c/Users/Jacob/Documents/test5.vtk
# test 6
.PHONY: 6
6:
	$(CC) $(CFLAGS)  -g test.cpp $(DEPS_O) -o $(OBJ)
	./$(OBJ) 6
	sudo cp test6.vtk /mnt/c/Users/Jacob/Documents/test6.vtk
	# test 6
.PHONY: 7
7:
	$(CC) $(CFLAGS)  -g test.cpp $(DEPS_O) -o $(OBJ)
	./$(OBJ) 7
.PHONY: 8
8:
	mpic++ $(CFLAGS)  -g test.cpp $(DEPS_O) -o $(OBJ)
	mpirun -np 4 ./$(OBJ) 8
