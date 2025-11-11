#include "Bucket_Partitioned_MDS.h"
#include <vector>
#include <cmath>

namespace Bucket_Partitioned_MDS
{
    Solution::Solution(
        const double _time_for_solving, 
        const double _cost, 
        const std::vector <std::vector <node_t>>& _routes) 
    {
        /*
        * Solution: Constructor for creating the object
        * @param _time_for_solving: Time consumed for solving the problem
        * @param _cost: Cost required by the vehicles to cover the _routes
        * @param _routes: Routes that should be traversed by the vehicles
        */

        time_for_solving = _time_for_solving;
        cost = _cost;
        routes = _routes;
    }

    bool Solution::verify(
        const CVRP& cvrp) const
    {
        /*
        * verify: Verifies the correctness of the solution
        *   1. For every route, the vehicle capacity constraint is respected
        *   2. All customers appear in exactly one route
        *   3. The cost returned by the solver is same as actual cost of the routes
        * @param cvrp: CVRP object to solve
        * @return: Returns true if the above verification steps are satisfied else false
        */

        std::vector <int> appear_count(cvrp.size(), 0);
        distance_t computed_cost = 0;
        for(const auto& route: routes)
        {
            distance_t curr_route_cost = cvrp.distance(cvrp.depot(), route[0]);
            demand_t route_demand = 0;
            for(int j = 1; j < route.size(); ++j)
            {
                curr_route_cost += cvrp.distance(route[j - 1], route[j]);
                route_demand += cvrp[route[j-1]].demand;
                appear_count[route[j-1]]++;
            }
            curr_route_cost += cvrp.distance(route.back(), cvrp.depot());
            route_demand += cvrp[route.back()].demand;
            appear_count[route.back()]++;

            computed_cost += curr_route_cost;
            if(route_demand > cvrp.capacity())
            {
                HANDLE_ERROR("Vehicle constriant is voilated in the routes", false);
                return false;
            }
        }

        for(node_t customer = 1; customer < cvrp.size(); customer++)
        {
            if(appear_count[customer] != 1) 
            {
                HANDLE_ERROR("Customer " + std::to_string(customer) + " is appearing " + std::to_string(appear_count[customer]) + " but expected to appear only once", false);
                return false;
            }
        }

        if(std::fabs(this->cost - computed_cost) > EPS) 
        {
            HANDLE_ERROR("The computed cost (" + std::to_string(computed_cost) + ") is not equal to solution cost (" + std::to_string(this->cost) + ")", false);
            return false;
        }
        return true;
    }

    void Solution::print(
        std::ostream& output) const 
    {
        /*
        * print: Function to print the execution time, cost and routes 
        * @param output: Output stream the contents should be printed
        * @return: Returns nothing
        */
        
        output << "Execution time for solving (sec): " << this->time_for_solving << "\n";
        output << "Cost: " << this->cost << "\n";

        for(int i = 0; i < this->routes.size(); i++)
        {
            output << "Route #" << i + 1 << ":";
            for(auto& node: this->routes[i])
            {
                output << " " << node;
            }
            output << "\n";
        }
    }
}