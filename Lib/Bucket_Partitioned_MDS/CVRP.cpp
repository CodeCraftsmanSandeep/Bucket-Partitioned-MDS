#include "Bucket_Partitioned_MDS.h"
#include "Utils.h"

#include <stdexcept>
#include <iomanip>
#include <cmath>

namespace Bucket_Partitioned_MDS
{
    CVRP::CVRP(
        std::istream& input)
    {
        /* 
        * CVRP: Constructor for CVRP class which reads and stores CVRP
        */
 
        std::string line;
        
        // Ignoring first 3 lines
        for (int i = 0; i < 3; ++i)
        {
            getline(input, line);
        }

        // Dimension
        getline(input, line);
        _size = stof(line.substr(line.find(":") + 2));

        // Distance type
        getline(input, line);
        type = line.substr(line.find(":") + 2);

        // Capacity
        getline(input, line);
        _capacity = stof(line.substr(line.find(":") + 2));

        // Ignoring "NORD_COORD_SECTION" line
        getline(input, line);

        // Allocate nodes
        node.resize(_size);

        // Parsing node co-ordinates
        for (int i = 0; i < _size; ++i) {
            getline(input, line);

            std::stringstream iss(line);
            int id;
            std::string xStr, yStr;

            iss >> id >> xStr >> yStr;
            node[i].x = stof(xStr);
            node[i].y = stof(yStr);
        }

        // Skip "DEMAND_SECTION" line
        getline(input, line);

        // Parsing demands
        for (int i = 0; i < _size; ++i) {
            getline(input, line);
            std::stringstream iss(line);

            int id;
            std::string dStr;
            iss >> id >> dStr;

            node[i].demand = stof(dStr);
        }
    }
    
    distance_t CVRP::distance(
        const node_t u, 
        const node_t v) const
    {
        /*
        * distance: Function to get distance between two nodes
        * @param u: Node u
        * @param v: Node v
        * @return: Distance between u and v
        */

        if (u < 0 || u >= _size)
        {
            HANDLE_ERROR("Index out of bounds: " + std::to_string(u) + " must be >= 0 and < " + std::to_string(_size), true);
        }

        if (v < 0 || v >= _size)
        {
            HANDLE_ERROR("Index out of bounds: " + std::to_string(v) + " must be >= 0 and < " + std::to_string(_size), true);
        }

        if (u == v) 
        {
            return 0.0; // Distance to itself is zero
        }

        // Euclidian distance calculating on the fly
        return std::sqrt((node[u].x - node[v].x) * (node[u].x - node[v].x) + 
            (node[u].y - node[v].y) * (node[u].y - node[v].y));
    }

    const Point& CVRP::operator[](
        const node_t index) const 
    {
        /*
        * operator[]: Function to access the Point at the index 
        * @param index: Index of the node to access
        */

        if (index < 0 || index >= _size)
        {
            HANDLE_ERROR("Index out of bounds: " + std::to_string(index) + " must be >= 0 and < " + std::to_string(_size), true);
        }
        return this->node[index]; 
    }

    void CVRP::print (
        std::ostream& output) const
    {
        /* 
        * print: Function to print CVRP
        * @param output: Output stream to print 
        * @return: Returns nothing
        */

        output << "SIZE: " << _size << "\n";
        output << "Capacity: " << _capacity << "\n";
        output << std::setw(6) << "NODE" 
                << std::setw(10) << "X"  
                << std::setw(10) << "Y" 
                << std::setw(10) << "DEMAND" << "\n";

        for (int i = 0; i < _size; ++i) 
        {
            output << std::setw(6) << i
                    << std::setw(10) << node[i].x 
                    << std::setw(10) << node[i].y 
                    << std::setw(10) << node[i].demand 
                    << "\n";
        }

        return;
    }

    capacity_t CVRP::capacity() const 
    {
        /*
        * capacity: Getter function to return capacity
        * @return: Returns capacity of the vehicle
        */

        return _capacity;
    }

    int CVRP::size() const 
    {
        /*
        * size: Getter function to return size
        * @return: Returns size of the CVRP
        */

        return _size;
    }

    node_t CVRP::depot() const 
    {
        /*
        * depot: Getter function to return depot
        * @return: Returns the depot 
        */

        return _depot;
    }
}