# Google AI Agents Build Error Fixing Test

This repository is created for testing Google AI agents in fixing build errors of a C++ codebase.

## How to Reproduce the Errors

To reproduce the demonstrated build errors, simply execute the following command:

```bash
cd breadth_first_search && make
```
# Expected compile errors
``` bash
g++ -I../ -std=c++17 -fopenmp -O3 -g -o bfs main.cpp bfs.cpp ../common/graph.cpp ref_bfs.a
bfs.cpp: In function ‘void vertex_set_init(vertex_set, int)’:
bfs.cpp:25:9: error: base operand of ‘->’ has non-pointer type ‘vertex_set’
   25 |     list->max_vertices = count;
      |         ^~
bfs.cpp:26:9: error: base operand of ‘->’ has non-pointer type ‘vertex_set’
   26 |     list->vertices = (int *)malloc(sizeof(int) * list->max_vertices);
      |         ^~
bfs.cpp:26:54: error: base operand of ‘->’ has non-pointer type ‘vertex_set’
   26 |     list->vertices = (int *)malloc(sizeof(int) * list->max_vertices);
      |                                                      ^~
bfs.cpp:27:22: error: cannot convert ‘vertex_set’ to ‘vertex_set*’
   27 |     vertex_set_clear(list);
      |                      ^~~~
      |                      |
      |                      vertex_set
bfs.cpp:18:35: note:   initializing argument 1 of ‘void vertex_set_clear(vertex_set*)’
   18 | void vertex_set_clear(vertex_set *list)
      |                       ~~~~~~~~~~~~^~~~
bfs.cpp: In function ‘void bfs_top_down(Graph, solution*)’:
bfs.cpp:113:21: error: could not convert ‘& list1’ from ‘vertex_set*’ to ‘vertex_set’
  113 |     vertex_set_init(&list1, graph->num_nodes);
      |                     ^~~~~~
      |                     |
      |                     vertex_set*
bfs.cpp:114:21: error: could not convert ‘& list2’ from ‘vertex_set*’ to ‘vertex_set’
  114 |     vertex_set_init(&list2, graph->num_nodes);
      |                     ^~~~~~
      |                     |
      |                     vertex_set*
bfs.cpp:119:25: error: could not convert ‘(vset + ((sizetype)(((long unsigned int)i) * 16)))’ from ‘vertex_set*’ to ‘vertex_set’
  119 |         vertex_set_init(&vset[i], graph->num_nodes);
      |                         ^~~~~~~~
      |                         |
      |                         vertex_set*
bfs.cpp: In function ‘void bfs_hybrid(Graph, solution*)’:
bfs.cpp:236:13: error: request for member ‘distance’ in ‘sol’, which is of pointer type ‘solution*’ (maybe you meant to use ‘->’ ?)
  236 |         sol.distance[i] = NOT_VISITED_MARKER;
      |             ^~~~~~~~
bfs.cpp:240:21: error: could not convert ‘& list1’ from ‘vertex_set*’ to ‘vertex_set’
  240 |     vertex_set_init(&list1, graph->num_nodes);
      |                     ^~~~~~
      |                     |
      |                     vertex_set*
bfs.cpp:241:21: error: could not convert ‘& list2’ from ‘vertex_set*’ to ‘vertex_set’
  241 |     vertex_set_init(&list2, graph->num_nodes);
      |                     ^~~~~~
      |                     |
      |                     vertex_set*
bfs.cpp:246:25: error: could not convert ‘(vset + ((sizetype)(((long unsigned int)i) * 16)))’ from ‘vertex_set*’ to ‘vertex_set’
  246 |         vertex_set_init(&vset[i], graph->num_nodes);
      |                         ^~~~~~~~
      |                         |
      |                         vertex_set*
make: *** [Makefile:4: default] Error 1
```

## Hint
This will result in 11 compile errors, but actually, only 3 modifications are required to eliminate them.

## Environment Information
  * Ubuntu 22.04
  * g++ 11.4.0