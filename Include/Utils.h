#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <vector>

// Useful constants
constexpr double PI = 3.14159265358979323846;
const double EPS = 1e-3;  // small tolerance for floating-point comparison

// Data types
using cord_t      = double; // Type for x and y co-ordinates
using distance_t  = double; // Type for distance between two vertices 
using demand_t    = double; // Type for demand of a customer
using capacity_t  = double; // Type for capacity of the vehicles
using node_t      = int; // Type for unique id for a node 
using angle_t     = double;

// Useful things for error handling
extern std::ostream& ERROR_FILE;
void handle_error(const char*, int, std::ostream& out, const std::string, const bool);
#define HANDLE_ERROR(msg, exit_flag) handle_error(__FILE__, __LINE__, ERROR_FILE, msg, exit_flag)

class Point
{
    /*
    * 2-D cordinates & demand of a customer (or) depot
    */    
public:
    cord_t x;
    cord_t y;
    demand_t demand;

    Point();
    Point(cord_t, cord_t, demand_t);
};

class Min_Heap_Node 
{
public:
    node_t u; // vertex already present in MST
    node_t v; // vertex going to be present in MST
    distance_t weight; // weight of the edge

    Min_Heap_Node(node_t, node_t, distance_t);
    Min_Heap_Node();
    void reset_node(node_t, node_t, distance_t);
};

class Min_Heap 
{
private:
    void swap_nodes(Min_Heap_Node&, Min_Heap_Node&);
    void heapify_up(node_t);
    void heapify_down(node_t);
public:
    node_t   last_index;
    Min_Heap_Node*    vec;
    node_t*  vertex_to_index_map;

    Min_Heap(const int);
    void DecreaseKey(const Min_Heap_Node&);
    Min_Heap_Node pop();
    bool empty();
    ~Min_Heap();
};

class Unit_Vector_2D
{
    /*
    * Unit_Vector_2D: Class to maintain 2-dimensional vector
    */
private:
    cord_t x, y; // x i + y j vector in 2-D plane

    distance_t get_norm(
        const cord_t, 
        const cord_t) const;
    angle_t degree_to_radian(
        const angle_t) const;

public:
    Unit_Vector_2D(
        cord_t, 
        cord_t);
    Unit_Vector_2D(
        const Point&, 
        const Point&);
    Unit_Vector_2D(
        const Unit_Vector_2D&, 
        const angle_t);
    bool operator==(
        const Unit_Vector_2D&) const;
    bool is_in_between(
        const Unit_Vector_2D&, 
        const Unit_Vector_2D&) const;
    void print(
        std::ostream&) const;
    Unit_Vector_2D(cord_t x1, cord_t y1, cord_t x2, cord_t y2) 
    {
        x = x2 - x1;
        y = y2 - y1;

        auto norm = get_norm(x, y);
        x /= norm;
        y /= norm;
    }
    Unit_Vector_2D() 
    {
        x = 0;
        y = 0;
    }
};

class Shared_Adj 
{
    /*
    * Shared_Adj: Shared adjacency list for MST's of all buckets 
    *             in Bucket-Partitioned-MDS which supports thread safe read 
    *             and write to data structures. 
    *             This data structure is thread safe for Bucket-Partitioned-MDS 
    *             sequence of operations.
    */

    const int num_buckets;
    const int cvrp_size;
    const node_t depot;
public:
    std::vector <std::vector <node_t>> depot_neighbours;
    std::vector <std::vector <node_t>> shared_adj_list; 

    Shared_Adj(
        const int, 
        const int, 
        const node_t
    );

    void add_edge(
        const int, 
        const node_t, 
        const node_t
    );

    const std::vector <node_t>& get_depot_neighbours(
        const int) const;

    const std::vector <node_t>& get_customer_neighbours(
        const int) const;
};

#endif