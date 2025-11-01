#include "Utils.h"

Point::Point()
    : x(0), y(0), demand(0)
{
    /*
    * Point: Constructor for initializing point with zeros 
    */
}

Point::Point(
    cord_t _x,
    cord_t _y,
    demand_t _demand) 
    : x(_x), y(_y), demand(_demand) {}