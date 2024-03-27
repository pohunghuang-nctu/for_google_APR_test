#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>
#include <algorithm>
#include <cstring>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1
//#define VERBOSE 1

void vertex_set_clear(vertex_set *list)
{
    list->count = 0;
}

void vertex_set_init(vertex_set list, int count)
{
    list->max_vertices = count;
    list->vertices = (int *)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set *frontier,
    vertex_set *new_frontier,
    vertex_set *vset,
    bool *bf,
    int *distances,
    bool hybrid)
{
    const int *outgoing_starts = g->outgoing_starts;
    const Vertex *outgoing_edges = g->outgoing_edges;
    const int g_num_nodes = g->num_nodes;
    const int g_num_edges = g->num_edges;
    const int current_distance = distances[frontier->vertices[0]] + 1;
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        int *local_frontiers = vset[tid].vertices;
        int local_count = 0;
        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < frontier->count; i++)
        {

            int node = frontier->vertices[i];

            int start_edge = outgoing_starts[node];
            int end_edge = (node == g_num_nodes - 1)
                            ? g_num_edges
                            : outgoing_starts[node + 1];

            // attempt to add all neighbors to the new frontier
            //#pragma omp for schedule(static, 2)
            for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
            {
                int outgoing = outgoing_edges[neighbor];

                //if (__sync_bool_compare_and_swap(&distances[outgoing], NOT_VISITED_MARKER, current_distance))
                if (distances[outgoing] == NOT_VISITED_MARKER)
                //if (__sync_bool_compare_and_swap(&bf[outgoing], false, true))
                //if (!bf[outgoing])
                {
                    
                    // int insert;
                    // #pragma omp atomic capture
                    // {
                    //     insert = distances[outgoing]
                        distances[outgoing] = current_distance;
                        bf[outgoing] = true;
                    // }
                    // if (insert == NOT_VISITED_MARKER)
                    // {
                        local_frontiers[local_count++] = outgoing;
                        //if (hybrid)
                        
                    //}

                }
            }
        }
        int *current_new_frontier_ptr;
        #pragma omp critical
        {
            current_new_frontier_ptr = new_frontier->vertices + new_frontier->count;
            new_frontier->count += local_count;
        }
        memcpy(current_new_frontier_ptr, local_frontiers, sizeof(int) * local_count);
        //free(local_frontiers);
    }
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution *sol)
{
    //printf("Running top-down BFS\n");
    //fflush(stdout);
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);
    int max_threads = omp_get_max_threads();
    vertex_set *vset = (vertex_set*)malloc(sizeof(vertex_set) * max_threads);
    for (int i = 0; i < max_threads; i++)
    {
        vertex_set_init(&vset[i], graph->num_nodes);
    }
    vertex_set *frontier = &list1;
    vertex_set *new_frontier = &list2;
    // initialize all nodes to NOT_VISITED
    #pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    bool *bf = (bool *)malloc(sizeof(bool) * graph->num_nodes);
    std::fill(bf, bf + graph->num_nodes, false);
    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;
    bf[ROOT_NODE_ID] = true;

    while (frontier->count != 0)
    {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);
        top_down_step(graph, frontier, new_frontier, vset, bf, sol->distances, true );

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        std::swap(frontier, new_frontier);

    }
}

int bottom_up_step(Graph g, bool *frontier, bool *new_frontier, int iteration, int *distances) {
    int frontier_count = 0;
    int cur_dist = iteration + 1;
    const int *incoming_starts = g->incoming_starts;
    const Vertex *incoming_edges = g->incoming_edges;
    const int g_num_nodes = g->num_nodes;
    const int g_num_edges = g->num_edges;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 512) reduction(+:frontier_count)
        for (int v = 0; v < g->num_nodes; v++)
        {
            if (distances[v] == NOT_VISITED_MARKER)
            {
                int start_edge = incoming_starts[v];
                int end_edge = (v == g_num_nodes - 1)
                                ? g_num_edges
                                : incoming_starts[v + 1];
                for (int i = start_edge; i < end_edge; i++)
                {
                    int u = incoming_edges[i];
                    if (frontier[u])
                    {
                            distances[v] = cur_dist;
                            new_frontier[v] = true;
                            frontier_count += 1;
                            break;
                    }
                }
            }
        }
    }
    return frontier_count;
}

void bfs_bottom_up(Graph graph, solution *sol)
{
    //printf("Running bottom-up BFS\n");
    //fflush(stdout);
    #pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    //frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    bool *frontier = (bool *)malloc(sizeof(bool) * graph->num_nodes);
    std::fill(frontier, frontier + graph->num_nodes, false);
    bool *new_frontier = (bool *)malloc(sizeof(bool) * graph->num_nodes);
    std::fill(new_frontier, new_frontier + graph->num_nodes, false);
    frontier[ROOT_NODE_ID] = true;

    int frontier_count = 1;
    sol->distances[ROOT_NODE_ID] = 0;

    int iteration = 0;
    while (frontier_count > 0)
    { 
#ifdef VERBOSE
        //printf("frontier=%-10d\n", frontier->count);
        double start_time = CycleTimer::currentSeconds();
#endif        
        //vertex_set_clear(new_frontier);
        frontier_count = bottom_up_step(graph, frontier, new_frontier, iteration, sol->distances);
        // swap pointers
        std::swap(frontier, new_frontier);
        iteration++;
        
#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier_count, end_time - start_time);
#endif               
    }
}

void bfs_hybrid(Graph graph, solution *sol)
{
    //printf("Running Hybrid BFS\n");
    //fflush(stdout);
    // initialize all nodes to NOT_VISITED
    #pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
        sol.distance[i] = NOT_VISITED_MARKER;    
    // prepare for top-down
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);
    int max_threads = omp_get_max_threads();
    vertex_set *vset = (vertex_set*)malloc(sizeof(vertex_set) * max_threads);
    for (int i = 0; i < max_threads; i++)
    {
        vertex_set_init(&vset[i], graph->num_nodes);
    }
    vertex_set *frontier = &list1;
    vertex_set *new_frontier = &list2;
    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;

    // prepare for bottom-up
    bool *bf = (bool *)malloc(sizeof(bool) * graph->num_nodes);
    std::fill(bf, bf + graph->num_nodes, false);
    bool *new_bf = (bool *)malloc(sizeof(bool) * graph->num_nodes);
    std::fill(new_bf, new_bf + graph->num_nodes, false);
    bf[ROOT_NODE_ID] = true;
    int frontier_count = 1;

    // for both, we set the root node to be visited
    sol->distances[ROOT_NODE_ID] = 0;
    int iteration = 0;    

    // top-down first, until number of frontiers is less than 0.02 of total nodes
    //printf("Running top-down\n");
    //fflush(stdout);
    while (frontier->count != 0 && frontier->count < graph->num_nodes * 0.02)
    {
#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif
        vertex_set_clear(new_frontier);
        top_down_step(graph, frontier, new_frontier, vset, bf, sol->distances, true);
#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif
        std::swap(frontier, new_frontier);
        iteration++;
    }
    frontier_count = frontier->count;
    // bottom-up
    //printf("Running bottom-up\n");
    //fflush(stdout);
    while (frontier_count > 0)
    {
#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif
        frontier_count = bottom_up_step(graph, bf, new_bf, iteration, sol->distances);
#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif
        std::swap(bf, new_bf);
        iteration++;
    }
}
