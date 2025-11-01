#ifndef BUCKET_PARTITIONED_MDS_H
#define BUCKET_PARTITIONED_MDS_H

#include "Utils.h"
#include <vector>
#include <string>
#include <cfloat>
#include <iostream>

namespace Bucket_Partitioned_MDS
{
    class CVRP
    {
        /*
        * CVRP: Class for maintaining CVRP 
        */
        
    private: 
        capacity_t capacity;
        size_t size;
        std::vector <Point> node;
        std::string type;
        const node_t depot = 0; // depot is always 0

    public:
        CVRP(
            std::istream&);
        distance_t get_distance(
            const node_t, 
            const node_t) const;
        const Point& operator[](
            const node_t) const;
        void print(
            std::ostream&) const;
        // Getters
        capacity_t get_capacity() const;
        size_t get_size() const;
        node_t get_depot() const;
    };

    class Solution 
    {
    private:
        double time_for_solving;
        double cost;
        std::vector <std::vector <node_t>> routes;

        distance_t get_total_cost_of_routes(
            const CVRP& cvrp);
    public:
        Solution(
            const double, 
            const double,
            const std::vector <std::vector <node_t>>&);
        bool verify(
            const CVRP&) const;
        void print(
            std::ostream&) const;
    };

    class Solver 
    {
        /*
        * Solver: Bucket Partitioned MDS solver for solving CVRP
        */

    private:
        const double alpha;
        const int lambda;
        const int rho;

        void get_bucket(
            const CVRP& cvrp, 
            const double start_angle, 
            const double end_angle, 
            std::vector <node_t>& bucket, 
            std::vector <node_t>& reverse_map) const;


        void construct_random_mst(
            const CVRP& cvrp,
            const std::vector <node_t>& bucket,
            std::vector <node_t>& depot_neighbours, 
            std::vector <std::vector <node_t>>& adj) const;

        distance_t get_route_distance(
            const CVRP&, 
            const std::vector <node_t>&) const;

        void tsp_approx(
            const CVRP&, 
            std::vector <node_t>&, 
            std::vector <node_t>&, 
            node_t) const;

        std::vector<std::vector<node_t>> process_tsp_approx(
            const CVRP&, 
            const std::vector<std::vector<node_t>>&) const;

        void tsp_2opt(
            const CVRP&, 
            std::vector <node_t>&,
            std::vector <node_t>&, 
            size_t) const;
        
        std::vector<std::vector<node_t>> process_2OPT(
            const CVRP&,
            const std::vector<std::vector<node_t>>&) const;
        
        void process_routes(
            const CVRP&, 
            std::vector <std::vector<int>>&, 
            distance_t&) const;

        void get_routes(
            const CVRP& cvrp, 
            const std::vector <node_t>& depot_neighbours, 
            const std::vector <std::vector <node_t>>& mst_adj,
            std::vector <std::vector <node_t>>& routes, 
            distance_t& cost, 
            const std::vector <node_t>& reverse_map, 
            const int num_nodes) const;

    public:
        Solver(
            const double, 
            const int, 
            const int);

        Solution solve(
            const CVRP&) const;
    };
}

#endif