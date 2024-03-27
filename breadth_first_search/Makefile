all: default grade debug

default: main.cpp bfs.cpp
	g++ -I../ -std=c++17 -fopenmp -O3 -g -o bfs main.cpp bfs.cpp ../common/graph.cpp ref_bfs.a

debug: main.cpp bfs.cpp
	g++ -I../ -std=c++17 -fopenmp -O3 -pg -DVERBOSE -g -o bfs_debug main.cpp bfs.cpp ../common/graph.cpp ref_bfs.a

grade: grade.cpp bfs.cpp
	g++ -I../ -std=c++17 -fopenmp -O3 -g -o bfs_grader grade.cpp bfs.cpp ../common/graph.cpp ref_bfs.a
clean:
	rm -rf bfs_grader bfs  *~ *.*~
