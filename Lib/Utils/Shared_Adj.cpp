#include "Utils.h"

Shared_Adj::Shared_Adj(
    const int cvrp_size, 
    const int num_buckets,
    const node_t depot) : cvrp_size(cvrp_size), num_buckets(num_buckets), depot(depot)
{
    /*
    * Shared_Adj: Constructor for creating shared adjacency list for CVRP space
    * @param num_buckets: Number of buckets the CVRP space is partitioned to
    * @param cvrp_size: Number of nodes (customers + depot) in CVRP space
    */

    this->shared_adj_list.resize(cvrp_size);
    this->depot_neighbours.resize(num_buckets);
}

void Shared_Adj::add_edge(
    const int bucket_index, 
    const node_t u, 
    const node_t v
)
{
    if(u == this->depot)
    {
        this->depot_neighbours[bucket_index].push_back(v);
    }
    else
    {
        this->shared_adj_list[u].push_back(v);
    }

    return;
}

const std::vector <node_t>& Shared_Adj::get_depot_neighbours(
    const int bucket_id) const 
{
    return this->depot_neighbours[bucket_id];
}

const std::vector <node_t>& Shared_Adj::get_customer_neighbours(
    const int vertex) const 
{
    return this->shared_adj_list[vertex];
}