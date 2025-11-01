#include "Utils.h"
#include "Bucket_Partitioned_MDS.h"
#include <chrono>
#include <stack>
#include <random>
#include <algorithm>
#include <climits>

namespace Bucket_Partitioned_MDS
{
    void Solver::get_bucket(
        const CVRP& cvrp, 
        const double start_angle, 
        const double end_angle, 
        std::vector <node_t>& bucket, 
        std::vector <node_t>& reverse_map) const
    {
        /*
        * get_bucket: Stores all the customers present in [start_angle, end_angle) region in bucket along with depot 
        * @param start_angle: Angle made with positive x axis to start, in anti clock wise direction
        * @param end_angle: Angle made with positive x axis to end, in anti clock wise direction
        * @param bucket: Bucket to store the customers which lie in [start_angle, end_angle) region along with depot
        * @return : Returns nothing
        */

        // Get size and depot of CVRP
        const size_t N = cvrp.get_size();
        const node_t depot = cvrp.get_depot();

        // Depot is present in every bucket
        bucket.push_back(depot); 

        // Finding partition vectors 
        Unit_Vector_2D xaxis(1, 0);        
        Unit_Vector_2D start_vector(xaxis, start_angle); // rotating xaxis by start_angle in anti clock wise direction
        Unit_Vector_2D end_vector(xaxis, end_angle); // rotating xaxis by end_angle in anti clock wise direction

        // Pushing all customers in [start_angle, end_angle) region into bucket
        int index = 1;
        for(node_t u = 1; u < N; u++)  
        {
            Unit_Vector_2D vec(cvrp[depot], cvrp[u]); // creating vector depot->u 
            if(vec.is_in_between(start_vector, end_vector)) // check if vec is in [stat_vector, end_vector) region
            {
                bucket.push_back(u); // push into bucket if present 
                reverse_map[u] = index++;
            }
        }

        return;
    }

    void Solver::construct_random_mst(
        const CVRP& cvrp,
        const std::vector <node_t>& bucket,
        std::vector <node_t>& depot_neighbours, 
        std::vector <std::vector <node_t>>& adj) const
    {
        /*
        * construct_random_mst: Helper function to construct random MST on the nodes in bucket
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
        int completed = 1; 
        node_t v;
        node_t v_index;
        for(v_index = 1; v_index < num_nodes; v_index++) 
        {
            min_heap.DecreaseKey(Min_Heap_Node(u_index, v_index, cvrp.get_distance(bucket[u_index], bucket[v_index]))); 
        }

        // Adding edges to MST 
        while(!min_heap.empty()) 
        { 
            // Get the minimum weight edge 
            min_node = min_heap.pop();
            u_index = min_node.u; // Index of the node in the bucket
            v_index = min_node.v; // Index of the neighbour of bucket[u_index]

            // Get the corresponding vertex in entire CVRP space                                      
            u = bucket[u_index];
            v = bucket[v_index];

            // Add the edge to the graph
            if(u == 0) depot_neighbours.push_back(v);
            else adj[u].push_back(v);

            if(v == 0) depot_neighbours.push_back(u);
            else adj[v].push_back(u);

            // v_index is added to MST
            completed++;                     

            // Loop over all neighbours of v_index 
            for(node_t w_index = 0; w_index < num_nodes; w_index++) 
            {
                node_t w = bucket[w_index];
                min_heap.DecreaseKey(Min_Heap_Node(v_index, w_index, cvrp.get_distance(v, w)));
            }
        }

        // Checking whether all bucket nodes are added in MST (or) not
        if(completed != num_nodes) 
        {
            HANDLE_ERROR("Not all nodes are included in the MST! Completed: " + std::to_string(completed) + ", Expected: " + std::to_string(num_nodes), true);
        }

        return;
    }

    distance_t Solver::get_route_distance(
        const CVRP& cvrp,
        const std::vector <node_t>& route) const
    {
        /*
        * get_route_distance: Helper function to calculate the cost
        * @param cvrp: CVRP instance 
        * @param route: Route 
        * @return: Returns distance a vehicle need to travel to cover the route
        */

        node_t prev_node = cvrp.get_depot();
        distance_t cost = 0;

        for(auto& node: route) 
        {
            cost += cvrp.get_distance(prev_node, node);
            prev_node = node;
        }
        cost += cvrp.get_distance(prev_node, cvrp.get_depot());

        return cost;
    }

    void Solver::tsp_approx(
        const CVRP& cvrp, 
        std::vector<node_t>& cities, 
        std::vector<node_t>& tour, 
        node_t ncities) const 
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

        for (i = 1; i < ncities; i++)
        {
            tour[i] = cities[i - 1];
        }
        
        tour[0] = cities[ncities - 1];

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
        const CVRP &cvrp, 
        const std::vector<std::vector<node_t>> &solRoutes) const 
    {
        std::vector<std::vector<node_t>> modifiedRoutes;
        size_t nroutes = solRoutes.size();

        for (size_t i = 0; i < nroutes; ++i) 
        {
            // Processing solRoutes[i]
            size_t sz = solRoutes[i].size();
            std::vector<node_t> cities(sz + 1);
            std::vector<node_t> tour(sz + 1);

            for (size_t j = 0; j < sz; ++j)
            {
                cities[j] = solRoutes[i][j];
            }

            cities[sz] = 0;  // the last node is the depot.
            this->tsp_approx(cvrp, cities, tour, sz + 1);

            // the first element of the tour is now the depot. So, ignore tour[0] and insert the rest into the vector.
            std::vector<node_t> curr_route;
            for (size_t kk = 1; kk < sz + 1; ++kk) 
            {
                curr_route.push_back(tour[kk]);
            }

            modifiedRoutes.push_back(curr_route);
        }
        return modifiedRoutes;
    }

    void Solver::tsp_2opt(
        const CVRP& cvrp, 
        std::vector <node_t>& cities, 
        std::vector <node_t>& tour, 
        size_t ncities) const
    {
        /* 
        * tsp_2opt: Finds the 2opt solution, repeats until optimization is found
        * @param cities: Contains the orginal solution, it is updated during the course of the 2opt-scheme to contain the 2opt solution
        * @param tour: Auxilary array for the function
        * @param ncities: Number of cities
        * @return: Returns nothing
        */

        size_t improve = 0;

        while (improve < 2) 
        {
            double best_distance = 0.0;
            best_distance += cvrp.get_distance(cvrp.get_depot(), cities[0]);  

            for (size_t jj = 1; jj < ncities; ++jj) 
            {
                best_distance += cvrp.get_distance(cities[jj - 1], cities[jj]);
            }

            best_distance += cvrp.get_distance(cvrp.get_depot(), cities[ncities - 1]);
            // 1x 2x 3x 4 5
            //  1 2  3  4 5
            for (size_t i = 0; i < ncities - 1; i++) 
            {
                for (size_t k = i + 1; k < ncities; k++) 
                {
                    for (size_t c = 0; c < i; ++c) 
                    {
                        tour[c] = cities[c];
                    }

                    size_t dec = 0;
                    for (size_t c = i; c < k + 1; ++c) 
                    {
                        tour[c] = cities[k - dec];
                        dec++;
                    }
                    for (size_t c = k + 1; c < ncities; ++c) 
                    {
                        tour[c] = cities[c];
                    }

                    double new_distance = cvrp.get_distance(cvrp.get_depot(), tour[0]);
                    for (size_t jj = 1; jj < ncities; ++jj) 
                    {
                        new_distance += cvrp.get_distance(tour[jj - 1], tour[jj]);
                    }
                    new_distance += cvrp.get_distance(cvrp.get_depot(), tour[ncities - 1]);

                    if (new_distance < best_distance) 
                    {
                        // Improvement found so reset
                        improve = 0;
                        for (size_t jj = 0; jj < ncities; jj++)
                        {
                            cities[jj] = tour[jj];
                        }
                        best_distance = new_distance;
                    }
                }
            }
            improve++;
        }

        return;
    }

    std::vector<std::vector<node_t>> Solver::process_2OPT(
        const CVRP& cvrp, 
        const std::vector<std::vector<node_t>>& routes) const
    {
        /*
        * postprocess_2OPT: Helper function to optimize cost 
        * @param cvrp: CVRP instance 
        * @param routes: Routes to process
        * @return: Returns processed routes
        */

        std::vector<std::vector<node_t>> processed_routes;
        size_t nroutes = routes.size();

        for (size_t i = 0; i < nroutes; ++i) 
        {
            size_t sz = routes[i].size();
            std::vector<node_t> cities = routes[i];
            std::vector<node_t> tour(sz);
            std::vector<node_t> curr_route;

            // For sz <= 1, the cost of the path cannot change. So not running tsp_2opt
            if (sz > 2)
            {
                this->tsp_2opt(cvrp, cities, tour, sz);
            }                         
            
            for (size_t kk = 0; kk < sz; ++kk) 
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

        auto processed_routes1 = process_tsp_approx(cvrp, routes);
        processed_routes1 = process_2OPT(cvrp, processed_routes1);
        auto processed_routes2 = process_2OPT(cvrp, routes);

        double min_cost = 0;
        for(int i = 0; i < routes.size(); i++)
        {
            // Choose as minimum cost route as possible
            auto route_cost = get_route_distance(cvrp, routes[i]);
            auto route_cost1 = get_route_distance(cvrp, processed_routes1[i]);
            auto route_cost2 = get_route_distance(cvrp, processed_routes2[i]); 

            if(std::min(route_cost, route_cost1) > route_cost2)
            {
                routes[i] = processed_routes2[i];
                min_cost += route_cost2;
            }
            else if(std::min(route_cost, route_cost2) > route_cost1)
            {
                routes[i] = processed_routes1[i];
                min_cost += route_cost1;
            }
            else
            {
                min_cost += route_cost;
            }
        }

        total_cost = min_cost;

        return;
    }
    
    void Solver::get_routes(
        const CVRP& cvrp, 
        const std::vector <node_t>& depot_neighbours, 
        const std::vector <std::vector <node_t>>& mst_adj,
        std::vector <std::vector <node_t>>& routes, 
        distance_t& cost, 
        const std::vector <node_t>& reverse_map, 
        const int num_nodes) const
    {
        /*
        * get_routes: Helper function to construct randomized routes from MST
        * @param cvrp: CVRP instance
        * @param mst_adj: Adjacency list of MST 
        * @param routes: Routes constructed from random DFS order of MST
        * @param cost: Total cost of the routes constructed
        * @return: Returns nothing
        */

        // Helping variables & storage      
        const node_t depot = 0;
        const node_t depot_index = 0;
        std::vector <bool> visited(num_nodes, false);

        // Random engine
        std::random_device rd;
        std::mt19937 rng(rd());

        // Data structures for route which is under construction
        std::vector <node_t> curr_route;
        distance_t curr_route_cost = 0.0;

        // Stack data structure for DFS (rec = recursive stack)
        std::stack <std::pair <size_t, node_t*>> rec;

        // Starting DFS with depot
        node_t v = depot;
        node_t v_index = depot_index;
        visited[v_index] = true;
        int covered = 1;
        size_t neigh_size = depot_neighbours.size();
        node_t* neigh = new node_t[neigh_size];
        std::copy(depot_neighbours.begin(), depot_neighbours.end(), neigh);
        std::shuffle(neigh, neigh + neigh_size, rng);
        rec.push({neigh_size - 1, neigh});
        int index;
        node_t prev_node = v;
        capacity_t residue_capacity = cvrp.get_capacity();

        // DFS iterative
        while(!rec.empty()) 
        {
            index = rec.top().first;
            neigh = rec.top().second;
            
            while(index >= 0) 
            {
                v = neigh[index];
                v_index = reverse_map[v];

                if(!visited[v_index]) 
                {
                    visited[v_index] = true;
                    if(residue_capacity >= cvrp[v].demand) 
                    {
                        // Push vertex v into current route
                        curr_route.push_back(v);
                        curr_route_cost += cvrp.get_distance(prev_node, v);
                        residue_capacity -= cvrp[v].demand;
                        prev_node = v;                    
                    } 
                    else 
                    {
                        // End the current route
                        covered += curr_route.size(); 
                        routes.push_back(std::move(curr_route));
                        curr_route_cost += cvrp.get_distance(prev_node, depot); 
                        cost += curr_route_cost;
                        curr_route.clear();
                        prev_node = depot; 
                        curr_route_cost = 0.0; 
                        residue_capacity = cvrp.get_capacity(); 
                
                        // Start a new route
                        curr_route.push_back(v);
                        residue_capacity -= cvrp[v].demand;
                        curr_route_cost += cvrp.get_distance(prev_node, v);
                        prev_node = v; 
                    }

                    // Pushing vertex v into stack for DFS
                    rec.top().first = index - 1; 
                    const std::vector <node_t>& vertex_neighbours = mst_adj[v];
                    neigh_size = vertex_neighbours.size();
                    neigh = new node_t[neigh_size];
                    std::copy(vertex_neighbours.begin(), vertex_neighbours.end(), neigh);
                    std::shuffle(neigh, neigh + neigh_size, rng);
                    rec.push({neigh_size - 1, neigh});
                    break;
                }
                index--;
            } 

            if(index < 0) 
            {
                // Every neighbour of vertex is visited
                delete[] neigh;
                rec.pop();
            }
        }

        // Pushing remaining nodes in route into list of routes
        if(!curr_route.empty()) 
        { 
            covered += curr_route.size(); 
            routes.push_back(std::move(curr_route));
            curr_route_cost += cvrp.get_distance(prev_node, depot); 
            cost += curr_route_cost; 
        }

        // // Making sure that every node in bucket is pushed into route
        // if(covered != num_nodes)
        // {
        //     HANDLE_ERROR("Not all nodes in bucket are pushed into routes", true);
        // }
        return;
    }

    Solver::Solver(
        const double _alpha, 
        const int _lambda, 
        const int _rho) 
        : alpha(_alpha), lambda(_lambda), rho(_rho) 
    {
        /*
        * Solver: Constructor for setting all paramters to run Bucket Paritioned MDS algo
        * @param _alpha: Alpha value
        * @param _lambda: Lambda value
        * @param _rho: Rho value
        */

        return;
    }

    Solution Solver::solve(
        const CVRP& cvrp) const
    {
        /*
        * solve: Function to solve CVRP with Bucket Partitioned MDS algo using paramters alpha, lambda and rho
        * @param cvrp: CVRP to solve
        * @return: Solution after solving the CVRP problem using Bucket Partitioned MDS 
        */

        // Start execution timer
        auto start = std::chrono::high_resolution_clock::now();

        // Data structures for storing solution
        distance_t final_cost = 0.0;
        std::vector <std::vector<int>> final_routes;

        // Useful values for solving
        const size_t N = cvrp.get_size();
        const node_t depot = cvrp.get_depot();
        const int num_buckets = std::ceil(360.00 / this->alpha);

        // Thread-safe shared adj list for Bucket-Partitioned-MDS steps
        std::vector <std::vector <node_t>> mst_adj(N);
        std::vector <node_t> reverse_map(N, 0);

        // Partitioning the CVRP space
        #pragma omp parallel for 
        for(int bucket_id = 0; bucket_id < num_buckets; bucket_id++) 
        {
            // Finding the vertices in this partition
            const angle_t start_angle = (bucket_id) * this->alpha;
            const angle_t end_angle = std::min(360.00, (bucket_id + 1) * this->alpha);
            std::vector <node_t> bucket;
            get_bucket(cvrp, start_angle, end_angle, bucket, reverse_map);

            // Useful data structures for exploration & exploitaion
            std::vector <std::vector <node_t>> low_cost_routes;
            distance_t low_cost = INT_MAX;

            // Get MST
            std::vector <node_t> depot_neighbours;
            this->construct_random_mst(cvrp, bucket, depot_neighbours, mst_adj);

            // Exploitation: getting different routes from different DFS orders of MST
            #pragma omp parallel for
            for(int _ = 1; _ <= rho; _++)
            {
                // Finding random DFS order of the MST
                std::vector <std::vector <node_t>> routes;
                distance_t cost = 0.0;
                this->get_routes(cvrp, depot_neighbours, mst_adj, routes, cost, reverse_map, bucket.size());

                #pragma omp critical 
                {
                    // Updating the cost if lower cost is found
                    if(low_cost > cost) 
                    {
                        low_cost_routes = routes;
                        low_cost = cost;
                    }
                }
            }

            if(!low_cost_routes.empty())
            {
                this->process_routes(cvrp, low_cost_routes, low_cost);   
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
        
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed_time = std::chrono::duration<double>(end - start).count();

        return Solution(elapsed_time, final_cost, final_routes);
    }
}