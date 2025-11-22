
#include "Utils.h"
#include <iostream>
#include <cmath>

distance_t Unit_Vector_2D::get_norm(
    const cord_t x, 
    const cord_t y) const
{
    /*
    * get_norm: Helper function for returning L2 norm
    * @param x: Length of vector in x direction
    * @param y: Length of vector in y direction
    * @return: Returns norm
    */

    return sqrt(x * x + y * y);
}

angle_t Unit_Vector_2D::degree_to_radian(
    const angle_t theta_deg) const 
{
    /*
    * degree_to_radian: Helper function for converting angle in degrees to radian
    * @param theta_deg: Theta in degrees
    * @return: Returns theta in radians
    */

    return theta_deg * PI / 180.0;
}

Unit_Vector_2D::Unit_Vector_2D(
    cord_t _x, 
    cord_t _y)
{
    /*
    * Unit_Vector_2D: Constructor for constructing unit vector with given vector co-ordinates
    * @param _x: Length of vector in x direction
    * @param _y: Length of vector in y direction
    * @return: Returns nothing
    */ 

    auto norm_val = get_norm(_x, _y);
    this->x = _x / norm_val;
    this->y = _y / norm_val;

    return;
}

Unit_Vector_2D::Unit_Vector_2D(
    const Point& point1, 
    const Point& point2)
{
    /*
    * Unit_Vector_2D: Constructor for constructing unit vector along p1 ~~> p2
    * @param point1: Tail (initial point) of vector
    * @param point2: Head (terminal point) of vector
    * @return: Returns nothing
    */

    auto _x = point2.x - point1.x;
    auto _y = point2.y - point1.y;
    auto norm_val = get_norm(_x, _y);
    this->x = _x / norm_val;
    this->y = _y / norm_val;

    return;
}

Unit_Vector_2D::Unit_Vector_2D(
    const Unit_Vector_2D& v, 
    const angle_t theta_deg) 
{
    /*
    * Unit_Vector_2D: Constructor for constructing vector by rotating vector `v` by angle `theta` in anti clock wise direction
    * @param v: Unit vector v
    * @param theta_deg: Angle to rotate v in anti clock wise direction
    */
        
    auto theta_rad = this->degree_to_radian(theta_deg);
    this->x = v.x * std::cos(theta_rad) - v.y * std::sin(theta_rad);
    this->y = v.x * std::sin(theta_rad) + v.y * std::cos(theta_rad);

    return;
}

bool Unit_Vector_2D::operator==(
    const Unit_Vector_2D& other) const 
{
    /*
    * operator==: Defining operator to check equality with vector
    * @param other: Unit vector `other` to check equality
    * @return: Returns true if `other` vector is equal to this vector
    */
    
    return (std::fabs(this->x - other.x) < EPS) && (std::fabs(this->y - other.y) < EPS);
}


bool Unit_Vector_2D::is_in_between(
    const Unit_Vector_2D& vec1, 
    const Unit_Vector_2D& vec2) const 
{
    /*
    * is_in_between: Function for checking whether this vector is in R = [vec1, vec2) region 
        i.e vec1 swiped in anti clock wise direction until vec2
    * @param vec1: Starting vector of the intereset region R
    * @param vec2: Ending vector of the interest region R
    * @return: Returns true if this vector is inside R
    */
    
    if (vec1 == vec2) 
    {
        HANDLE_ERROR("The angle between starting & ending vectors of a partition cannot be zero!", true);
    }
    
    if (vec1 == *this) 
    {
        return true; // vector is same as vec1
    }

    if (vec2 == *this) 
    {
        return false; // vector is same as vec2
    }

    // Compute cross products
    cord_t cross12 = vec1.x * vec2.y - vec1.y * vec2.x;     // vec1 × vec2
    cord_t cross1p = vec1.x * this->y - vec1.y * this->x;   // vec1 × vecp (vecp is `this` vector)
    cord_t crossp2 = this->x * vec2.y - this->y * vec2.x;   // vecp × vec2 (vecp is `this` vector)

    if(cross12 >= 0) 
    {
        // Angle between 1 and 2 is <= 180 degress CCW
        return cross1p >= 0 && crossp2 >= 0;
    }
    else
    {
        // Angle between 1 and 2 is > 180 degrees CCW
        return !(cross1p < 0 && crossp2 < 0);
    }
}

void Unit_Vector_2D::print(
    std::ostream& output) const
{ 
    /*
    * print: Function to print `this` unit vector
    * @param output: output stream where the vector is to be printed
    * @return: Returns nothing
    */

    output << "Unit Vector: (" << this->x << ", " << this->y << ")" << "\n";

    return;
}