#include "Utils.h"
#include <climits>

void Min_Heap::swap_nodes(Min_Heap_Node& left, Min_Heap_Node& right) 
{
    Min_Heap_Node storage = left; 
    left = right;
    right = storage;
}

Min_Heap::Min_Heap(const size_t num_nodes)
{
    this->last_index = num_nodes - 1;
    this->vec = new Min_Heap_Node[num_nodes];
    this->vertex_to_index_map = new node_t[num_nodes];

    for(node_t i = 0; i < num_nodes; i++)
    {
        this->vertex_to_index_map[i] = i;
        this->vec[i].reset_node(-1, i, INT_MAX);
    }
}

void Min_Heap::heapify_up(node_t index) 
{
    while(index > 0) 
    {
        node_t parent = (index - 1) >> 1;
        if(vec[parent].weight > vec[index].weight) 
        {
            swap_nodes(vec[parent], vec[index]);
            vertex_to_index_map[vec[index].v] = index;
            index = parent;
        } else break;
    }
    vertex_to_index_map[vec[index].v] = index;
}

void Min_Heap::DecreaseKey(const Min_Heap_Node& node) 
{
    auto vertex = node.v;
    auto index  = vertex_to_index_map[vertex];
    
    if(index > last_index) 
    {
        // The vertex is already is present in MST and that is why it is out of heap
        // So no need to reduce further
        return;
    }

    if(vec[index].weight <= node.weight) 
    {
        // Current weight <= new weight, so no need to reduce further
        return;
    }

    vec[index].weight = node.weight;
    vec[index].u = node.u;
    heapify_up(index);
}

void Min_Heap::heapify_down(node_t index) 
{
    node_t left, right, small_weight_index;

    while (index <= last_index) 
    {
        left = (index << 1) + 1;
        right = (left) + 1;
        small_weight_index = index;

        if (left <= last_index && vec[left].weight < vec[small_weight_index].weight)
        {
            small_weight_index = left;
        }

        if (right <= last_index && vec[right].weight < vec[small_weight_index].weight) 
        {
            small_weight_index = right;
        }

        if (small_weight_index != index) 
        {
            swap_nodes(vec[small_weight_index], vec[index]);
            vertex_to_index_map[vec[index].v] = index;
            index = small_weight_index;
        } 
        else 
        {
            vertex_to_index_map[vec[index].v] = index;
            break;
        }
    }
}

Min_Heap_Node Min_Heap::pop() 
{
    swap_nodes(vec[0], vec[last_index]);
    vertex_to_index_map[vec[last_index].v] = last_index;
    last_index--;

    heapify_down(0);
    return vec[last_index + 1];
}

bool Min_Heap::empty() 
{
    return last_index < 0;
}

Min_Heap::~Min_Heap() 
{
    delete[] vec;
    delete[] vertex_to_index_map;
}