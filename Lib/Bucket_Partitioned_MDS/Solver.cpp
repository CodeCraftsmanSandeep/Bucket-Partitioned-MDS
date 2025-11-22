#include "Utils.h"
#include "Bucket_Partitioned_MDS.h"
#include <chrono>
#include <stack>
#include <random>
#include <algorithm>
#include <climits>

namespace Bucket_Partitioned_MDS
{
    void Solver::create_buckets(
        const CVRP& cvrp, 
        std::vector <std::vector <node_t>>& buckets) const 
    {
        const int num_buckets = buckets.size();
        const int N = cvrp.size();
        if(num_buckets == 1)
        {
            // Push all the vertices into one bucket
            buckets[0].reserve(N);
            for(node_t u = 0; u < N; u++)
            {
                buckets[0].push_back(u);
            }
            return;
        }

        std::vector<Unit_Vector_2D> seperating_vectors(num_buckets + 1);
        Unit_Vector_2D xaxis(1, 0); 
        seperating_vectors[0] = seperating_vectors[num_buckets] = xaxis;
        for(int i = 1; i < num_buckets; i++)
        {
            // Create a vector by rotating x-axis by i * alpha degrees
            seperating_vectors[i] = Unit_Vector_2D(xaxis, i * alpha);
        }

        // Depot goes into all buckets
        auto depot = cvrp.depot();
        node_t u = depot;
        for(int i = 0; i < num_buckets; i++)
        {
            buckets[i].push_back(u); 
        }

        // Assigning bucket to each vertex
        for(u = 1; u < N; u++)
        {
            Unit_Vector_2D vec(cvrp[depot], cvrp[u]);
            for(int i = 0; i < num_buckets; i++)
            {
                if(vec.is_in_between(seperating_vectors[i], seperating_vectors[i + 1]))
                {
                    // Pushing node u in the bucket i
                    buckets[i].push_back(u);
                    break;
                }
            }
        }

        return;
    }

    void Solver::get_bucket(
        const CVRP&             cvrp, 
        const double            start_angle, 
        const double            end_angle, 
        std::vector <node_t>&   bucket) const
    {
        /*
        * get_bucket: Stores all the customers present in [start_angle, end_angle) region in bucket along with depot 
        * @param start_angle: Angle made with positive x axis to start, in anti clock wise direction
        * @param end_angle: Angle made with positive x axis to end, in anti clock wise direction
        * @param bucket: Bucket to store the customers which lie in [start_angle, end_angle) region along with depot
        * @return : Returns nothing
        */

        // Get size and depot of CVRP
        const int N = cvrp.size();
        const node_t depot = cvrp.depot();

        // Depot is present in every bucket
        bucket.push_back(depot); 

        // Finding partition vectors 
        Unit_Vector_2D xaxis(1, 0);        
        Unit_Vector_2D start_vector(xaxis, start_angle); // rotating xaxis by start_angle in anti clock wise direction
        Unit_Vector_2D end_vector(xaxis, end_angle); // rotating xaxis by end_angle in anti clock wise direction

        // Pushing all customers in [start_angle, end_angle) region into bucket
        for(node_t u = 1; u < N; u++)  
        {
            Unit_Vector_2D vec(cvrp[depot], cvrp[u]); // creating vector depot->u 
            if(vec.is_in_between(start_vector, end_vector)) // check if vec is in [stat_vector, end_vector) region
            {
                bucket.push_back(u); // push into bucket if present 
            }
        }

        return;
    }

    void Solver::construct_mst(
        const CVRP&                         cvrp,
        const std::vector <node_t>&         bucket,
        std::vector <std::vector <node_t>>& adj) const
    {
        /*
        * construct_mst: Helper function to construct random MST on the nodes in bucket
        * @param cvrp: CVRP instance which is currently considered to solve
        * @param bucket: Nodes in bucket on which MST is to be constructed
        * @param adj: Adjacency list of MST
        * @return: Returns nothing
        */

        // Create min heap
        const int num_nodes = bucket.size();
        Min_Heap min_heap(num_nodes);

        // Starting from depot
        const node_t depot = 0;
        const node_t depot_index = 0;
        node_t u = depot;
        node_t u_index = depot_index;             
        min_heap.DecreaseKey(Min_Heap_Node(-1, depot_index, 0));    
        Min_Heap_Node min_node = min_heap.pop(); 
        node_t v;
        node_t v_index;
        for(v_index = 1; v_index < num_nodes; v_index++) 
        {
            min_heap.DecreaseKey(Min_Heap_Node(u_index, v_index, cvrp.distance(bucket[u_index], bucket[v_index]))); 
        }

        // Adding edges to MST 
        while(!min_heap.empty()) 
        { 
            // Get the minimum weight edge 
            min_node = min_heap.pop();
            u_index = min_node.u; // Index of the node in the bucket
            v_index = min_node.v; // Index of the neighbour of bucket[u_index]

            // Add the edge to the graph (v_index is added to MST)
            adj[u_index].push_back(v_index);                                
            adj[v_index].push_back(u_index);       

            // Get the corresponding vertex in entire CVRP space                                      
            v = bucket[v_index];                                                       

            // Loop over all neighbours of v_index 
            for(node_t w_index = 0; w_index < num_nodes; w_index++) 
            {
                node_t w = bucket[w_index];
                min_heap.DecreaseKey(Min_Heap_Node(v_index, w_index, cvrp.distance(v, w)));
            }
        }

        return;
    }

    distance_t Solver::get_route_distance(
        const CVRP&                 cvrp,
        const std::vector <node_t>& route) const
    {
        /*
        * get_route_distance: Helper function to calculate the cost
        * @param cvrp: CVRP instance 
        * @param route: Route 
        * @return: Returns distance a vehicle need to travel to cover the route
        */

        node_t prev_node = cvrp.depot();
        distance_t cost = 0;

        for(auto& node: route)
        {
            cost += cvrp.distance(prev_node, node);
            prev_node = node;
        }
        cost += cvrp.distance(prev_node, cvrp.depot());

        return cost;
    }

    void Solver::tsp_approx(
        const CVRP&             cvrp, 
        std::vector<node_t>&    cities, 
        std::vector<node_t>&    tour, 
        const int               ncities) const 
    {
        /*
        * tsp_approx: Travelling Salesman Problem (TSP) approximation
        * @param cvrp: CVRP instance 
        * @param cities: Cities in TSP
        * @param tour: Tour
        * @param ncities: Number of cities
        * @return: Returns nothing
        */

        node_t i, j;
        node_t ClosePt = 0;
        distance_t CloseDist;

        std::copy(cities.begin(), cities.end() - 1, tour.begin() + 1);
        tour[0] = cities.back();


        for (i = 1; i < ncities; i++) 
        {
            distance_t ThisX = cvrp[tour[i - 1]].x;
            distance_t ThisY = cvrp[tour[i - 1]].y;
            CloseDist = DBL_MAX;

            for (j = ncities - 1;; j--) 
            {
                distance_t ThisDist = (cvrp[tour[j]].x - ThisX) * (cvrp[tour[j]].x - ThisX);
                if (ThisDist <= CloseDist) 
                {
                    ThisDist += (cvrp[tour[j]].y - ThisY) * (cvrp[tour[j]].y - ThisY);
                    if (ThisDist <= CloseDist) 
                    {
                        if (j < i)
                        {
                            break;
                        }
                        CloseDist = ThisDist;
                        ClosePt = j;
                    }
                }
            }

            /*swapping tour[i] and tour[ClosePt]*/
            auto temp = tour[i];
            tour[i] = tour[ClosePt];
            tour[ClosePt] = temp;
        }

        return;
    }

    std::vector<std::vector<node_t>> Solver::process_tsp_approx(
        const CVRP&                             cvrp, 
        const std::vector<std::vector<node_t>>& sol_routes) const 
    {
        std::vector<std::vector<node_t>> modifiedRoutes;
        const int num_routes = sol_routes.size();
        modifiedRoutes.reserve(sol_routes.size());

        for(int i = 0; i < num_routes; i++)
        {
            // Processing solRoutes[i]
            int sz = sol_routes[i].size();
            std::vector<node_t> cities = sol_routes[i];
            cities.push_back(cvrp.depot()); 
            std::vector<node_t> tour(sz + 1);
            tsp_approx(cvrp, cities, tour, sz + 1);

            // the first element of the tour is now the depot. So, ignore tour[0] and insert the rest into the vector.
            // modifiedRoutes.emplace_back(tour.begin() + 1, tour.begin() + 1 + sz);
            std::vector<node_t> curr_route; 
            for (int kk = 1; kk < sz + 1; ++kk) 
            { 
                curr_route.push_back(tour[kk]); 
            } 
            modifiedRoutes.push_back(curr_route);
        }
        return modifiedRoutes;
    }

    void Solver::tsp_2OPT(
        const CVRP&             cvrp, 
        std::vector <node_t>&   cities, 
        std::vector <node_t>&   tour, 
        int                     ncities) const
    {
        /* 
        * tsp_2OPT: Finds the 2opt solution, repeats until optimization is found
        * @param cities: Contains the orginal solution, it is updated during the course of the 2opt-scheme to contain the 2opt solution
        * @param tour: Auxilary array for the function
        * @param ncities: Number of cities
        * @return: Returns nothing
        */

        int improve = 0;

        while (improve < 2) 
        {
            double best_distance = cvrp.distance(cvrp.depot(), cities[0]);  

            for (int jj = 1; jj < ncities; ++jj) 
            {
                best_distance   += cvrp.distance(cities[jj - 1], cities[jj]);
            }

            best_distance       += cvrp.distance(cvrp.depot(), cities[ncities - 1]);
            // 1x 2x 3x 4 5
            //  1 2  3  4 5
            for (int i = 0; i < ncities - 1; i++) 
            {
                for (int k = i + 1; k < ncities; k++) 
                {
                    for (int c = 0; c < i; ++c) 
                    {
                        tour[c] = cities[c];
                    }

                    int dec  = 0;
                    for (int c = i; c < k + 1; ++c) 
                    {
                        tour[c] = cities[k - dec];
                        dec++;
                    }
                    for (int c = k + 1; c < ncities; ++c) 
                    {
                        tour[c] = cities[c];
                    }

                    double new_distance = cvrp.distance(cvrp.depot(), tour[0]);
                    for (int jj = 1; jj < ncities; ++jj) 
                    {
                        new_distance    += cvrp.distance(tour[jj - 1], tour[jj]);
                    }
                    new_distance        += cvrp.distance(cvrp.depot(), tour[ncities - 1]);

                    if (new_distance < best_distance) 
                    {
                        // Improvement found so reset
                        improve = 0;
                        for (int jj = 0; jj < ncities; jj++)
                        {
                            cities[jj]  = tour[jj];
                        }
                        best_distance   = new_distance;
                    }
                }
            }
            improve++;
        }

        return;
    }

    std::vector<std::vector<node_t>> Solver::process_2OPT(
        const CVRP&                             cvrp, 
        const std::vector<std::vector<node_t>>& routes) const
    {
        /*
        * postprocess_2OPT: Helper function to optimize cost 
        * @param cvrp: CVRP instance 
        * @param routes: Routes to process
        * @return: Returns processed routes
        */

        std::vector<std::vector<node_t>> processed_routes;
        int nroutes = routes.size();

        for (int i = 0; i < nroutes; ++i) 
        {
            int sz = routes[i].size();
            std::vector<node_t> cities = routes[i];
            std::vector<node_t> tour(sz);
            std::vector<node_t> curr_route;

            // For sz <= 1, the cost of the path cannot change. So not running tsp_2opt
            if (sz > 2)
            {
                tsp_2OPT(cvrp, cities, tour, sz);
            }                         
            
            for (int kk = 0; kk < sz; ++kk) 
            {
                curr_route.push_back(cities[kk]);
            }
            
            processed_routes.push_back(curr_route);
        }

        return processed_routes;
    }

    void Solver::process_routes(
        const CVRP& cvrp,
        std::vector<std::vector<node_t>>& routes, 
        distance_t& total_cost) const
    {
        /*
        * process_routes: Processing routes to get as minimal as possible
        * @param cvrp: CVRP instance 
        * @param routes: routes to process
        * @param total_cost: cost of routes 
        * @return: Returns nothing
        */

        auto processed_routes1  = process_tsp_approx(cvrp, routes);
        processed_routes1       = process_2OPT(cvrp, processed_routes1);
        auto processed_routes2  = process_2OPT(cvrp, routes);

        double min_cost = 0;
        for(int i = 0; i < routes.size(); i++)
        {
            // Choose as minimum cost route as possible
            auto route_cost     = get_route_distance(cvrp, routes[i]);
            auto route_cost1    = get_route_distance(cvrp, processed_routes1[i]);
            auto route_cost2    = get_route_distance(cvrp, processed_routes2[i]); 

            if(std::min(route_cost1, route_cost2) >= route_cost)
            {
                min_cost        += route_cost;
            }
            else if(std::min(route_cost1, route_cost) >= route_cost2)
            {
                routes[i]       = processed_routes2[i];
                min_cost        += route_cost2;
            }
            else
            {
                routes[i]       = processed_routes1[i];
                min_cost        += route_cost1;   
            }
        }

        total_cost = min_cost;

        return;
    }

    void Solver::get_routes(
        const CVRP&                                 cvrp, 
        const std::vector <std::vector <node_t>>&   mst_adj, 
        const std::vector <node_t>&                 bucket, 
        std::vector <std::vector <node_t>>&         routes, 
        distance_t&                                 cost) const
    {
        /*
        * get_routes: Helper function to construct randomized routes from MST
        * @param cvrp: CVRP instance
        * @param mst_adj: Adjacency list of MST 
        * @param bucket: Nodes in the bucket
        * @param routes: Routes constructed from random DFS order of MST
        * @param const: Total cost of the routes constructed
        * @return: Returns nothing
        */

        // Helping variables & storage
        const int num_nodes = bucket.size();
        if(num_nodes == 1) return;          
        const node_t depot          = cvrp.depot();
        const node_t depot_index    = 0; // index of depot in the bucket
        std::vector <bool> visited(num_nodes, false);

        // Random engine
        std::random_device rd;
        std::mt19937 rng(rd());

        // Data structures for route which is under construction
        std::vector <node_t> curr_route;

        // Stack data structure for DFS
        std::stack <std::pair <int, node_t*>> rec;

        // Starting DFS with depot
        node_t v                    = depot;
        node_t v_index              = depot_index;
        visited[v_index]            = true;                    
        int neigh_size              = mst_adj[v_index].size();
        node_t* neigh               = new node_t[neigh_size];
        std::copy(mst_adj[v_index].begin(), mst_adj[v_index].end(), neigh);
        std::shuffle(neigh, neigh + neigh_size, rng);
        rec.push({neigh_size - 1, neigh});
        int index;
        node_t prev_node            = v;
        capacity_t residue_capacity = cvrp.capacity();

        // DFS iterative
        while(!rec.empty()) 
        {
            index = rec.top().first; 
            neigh = rec.top().second;

            while(index >= 0) 
            {
                v_index = neigh[index];
                v       = bucket[v_index];
                
                if(!visited[v_index]) 
                {
                    visited[v_index] = true;
                    if(residue_capacity >= cvrp[v].demand) 
                    {
                        // Push vertex v into current route
                        curr_route.push_back(v);
                        cost                += cvrp.distance(prev_node, v);
                        residue_capacity    -= cvrp[v].demand;
                        prev_node           = v;                    
                    } 
                    else 
                    {
                        // End the current route
                        routes.push_back(std::move(curr_route));
                        cost                += cvrp.distance(prev_node, depot); 
                        prev_node           = depot; 
                        residue_capacity    = cvrp.capacity(); 
                
                        // Start a new route
                        curr_route.push_back(v);
                        residue_capacity    -= cvrp[v].demand;
                        cost                += cvrp.distance(prev_node, v);
                        prev_node           = v; 
                    }

                    // Pushing vertex v into stack for DFS
                    rec.top().first = index - 1; 
                    neigh_size      = mst_adj[v_index].size();
                    neigh           = new node_t[neigh_size];
                    std::copy(mst_adj[v_index].begin(), mst_adj[v_index].end(), neigh);
                    std::shuffle(neigh, neigh + neigh_size, rng);
                    rec.push({neigh_size - 1, neigh});
                    break;
                }
                index--;
            }

            if(index < 0) 
            {
                // Every neighbour of vertex is visited
                free(rec.top().second);
                rec.pop();
            }
        }

        // Pushing remaining nodes in route into list of routes
        if(!curr_route.empty()) 
        { 
            routes.push_back(std::move(curr_route));
            cost += cvrp.distance(prev_node, depot); 
        }

        return;
    }

    Solver::Solver(
        const double _alpha, 
        const int _rho) 
        : alpha(_alpha), rho(_rho) 
    {
        /*
        * Solver: Constructor for setting all paramters to run Bucket Paritioned MDS algo
        * @param _alpha: Alpha value
        * @param _rho: Rho value
        */

        return;
    }

    Solution Solver::solve(
        const CVRP& cvrp) const
    {
        /*
        * solve: Function to solve CVRP with Bucket Partitioned MDS algo using paramters alpha and rho
        * @param cvrp: CVRP to solve
        * @return: Solution after solving the CVRP problem using Bucket Partitioned MDS 
        */
        double maxMB_before_execution = get_peak_mb();

        // Start execution timer
        auto start = std::chrono::high_resolution_clock::now();

        // Data structures for storing solution
        distance_t final_cost   = 0.0;
        std::vector <std::vector<int>> final_routes;

        // Useful values for solving
        const int num_buckets = std::ceil(360.00 / alpha); 
        std::vector <std::vector<node_t>> buckets(num_buckets);
        create_buckets(cvrp, buckets); 

        // Paritioning the problem for exploitation
        #pragma omp parallel for 
        for(int bucket_id = 0; bucket_id < num_buckets; bucket_id++) 
        {
            const std::vector <node_t>& bucket = buckets[bucket_id];

            // Useful data structures for exploitation
            std::vector <std::vector <node_t>> low_cost_routes;
            distance_t low_cost      = DBL_MAX;

            // Get MST
            std::vector <std::vector <node_t>> mst_adj(bucket.size());
            construct_mst(cvrp, bucket, mst_adj);

            // Exploitation: getting different routes from different DFS orders of MST
            #pragma omp parallel for 
            for(int _ = 0; _ < rho; _++)
            {
                // Finding random DFS order of the MST
                std::vector <std::vector <node_t>> routes;
                distance_t cost = 0.0;
                get_routes(cvrp, mst_adj, bucket, routes, cost);

                #pragma omp critical 
                { 
                    // Updating the cost if lower cost routes are found
                    if(low_cost > cost) 
                    {
                        low_cost_routes = routes;
                        low_cost        = cost;
                    }
                }
            }

            if(!low_cost_routes.empty())
            {
                process_routes(cvrp, low_cost_routes, low_cost);   
                #pragma omp critical 
                {
                    for(auto& route: low_cost_routes) 
                    { 
                        final_routes.push_back(std::move(route)); 
                    } 
                    final_cost += low_cost;
                }
            }
        }

        auto end                    = std::chrono::high_resolution_clock::now();
        double maxMB_after_execution = get_peak_mb();

        double execution_time = std::chrono::duration<double>(end - start).count();
        double maxMB_difference = maxMB_after_execution - maxMB_before_execution;

        return Solution(execution_time, maxMB_difference, final_cost, final_routes);
    }
}