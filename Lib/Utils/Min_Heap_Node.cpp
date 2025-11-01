#include "Utils.h"

Min_Heap_Node::Min_Heap_Node(
    node_t _u, 
    node_t _v, 
    distance_t _weight)
    : u(_u), v(_v), weight(_weight) {}

Min_Heap_Node::Min_Heap_Node()
    : u(0), v(0), weight(0) {}

void Min_Heap_Node::reset_node(
    node_t _u, 
    node_t _v, 
    distance_t _weight)
{
    /*
    * reset_node: Member function to set the values of data members after constructing the object
    * @param _u: The vertex already present in MST
    * @param _v: The vertex that is a candidate to be present in MST
    * @param _weight: The weight of edge (_u, _v)
    */

    this->u = _u;
    this->v = _v;
    this->weight = _weight;
}