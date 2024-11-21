#include "Objects.hpp"
#include <cmath>
#include <iostream>

#define PI M_PI


//Uniform_Magnetic_Field
Uniform_Magnetic_Field::Uniform_Magnetic_Field(Vector const &direction, long double const &magnitude){ //passing by reference to reduce amount of memory used
    this->magnetic_field_vector = Vector::adjust_mag(direction, magnitude);
    this->name = "UMF";
}
// Uniform_Magnetic_Field() : direction(Vector(0,0,0)), magnitude(0) {} //default constructor
Vector Uniform_Magnetic_Field::magnetic_field_strength(Vector const &position){
    return this->magnetic_field_vector;
}
Vector Uniform_Magnetic_Field::electric_field_strength(Vector const &position){
    return Vector(0,0,0);
}



//Uniform_Electric_Field
Uniform_Electric_Field::Uniform_Electric_Field(Vector const &direction, long double const &magnitude){ //passing by reference to reduce amount of memory used
    this->electric_field_vector = Vector::adjust_mag(direction, magnitude);
    this->name = "UEF";
}
// Uniform_Magnetic_Field() : direction(Vector(0,0,0)), magnitude(0) {} //default constructor
Vector Uniform_Electric_Field::magnetic_field_strength(Vector const &position){
    return Vector(0,0,0);
}
Vector Uniform_Electric_Field::electric_field_strength(Vector const &position){
    return this->electric_field_vector;
}


//point charge
Static_Point_Charge::Static_Point_Charge(Vector const &origin, long double charge){
    this->charge = charge;
    this->origin = origin;
    this->name = "SPC";
    this->k = 8.99e9;
}
Vector Static_Point_Charge::magnetic_field_strength(Vector const &position){
    return Vector(0,0,0);
}
Vector Static_Point_Charge::electric_field_strength(Vector const &position){

    // Vector OP = Vector::sub(position, this->origin); //OLD
    Vector OP = position - this->origin; //defining a vector from the origin to the point in space
    long double r = Vector::norm(OP); //compute distance from point charge
    //equation is E = kq/r^2. 
    long double magnitude = (this->k * this->charge)/std::pow(r, 2);
    OP = Vector::adjust_mag(OP, magnitude);

    return OP;
}

//Wire

Wire::Wire(Vector const &origin, Vector const &wire_direction, long double i_wire){
    this->i_wire = i_wire;
    this->wire_direction = wire_direction;
    this->origin = origin;
    this->name = "Wire";
    this->mu_0 = mu_0;
}
Wire::Wire() : origin(Vector(0,0,0)), wire_direction(Vector(0,0,0)), i_wire(0) {}
//return the field vector at a given position
Vector Wire::magnetic_field_strength(Vector const &position){  
    //first, we need to figure out the distance from the wire, r.
    //we know that the wire is a line with direction vector wire_direction and that the particle is at a point described by the position vector.
    //we use the perp formula
    //consider O the origin of wire and P the position of particle

    // Vector OP = Vector::sub(position, this->origin); //OLD
    Vector OP = position - this->origin;
    // Vector perp = Vector::sub(OP, Vector::sc_mult(this->wire_direction, (Vector::dot(OP,this->wire_direction))/(Vector::dot(wire_direction, wire_direction)))); //OLD
    Vector perp = OP - (((OP % this->wire_direction)/(wire_direction % wire_direction)) * this->wire_direction);
    long double radius = Vector::norm(perp); // dist from the wire to point

    //now, recall B = k/r
    long double k_wire = this->mu_0*i_wire/(2*PI);  //shorthand for the B equation: B = (mu/2PI)*(I/r) = k/r

    long double field_strength = k_wire/radius; // calculate field strength magnitude
    // Vector field_dir = Vector::cross(this->wire_direction, perp);  //OLD
    Vector field_dir = this->wire_direction * perp; //determine the field direction vector
    long double scaling_constant = field_strength/Vector::norm(field_dir);  //determine the scaling constant required to get the right magnitude

    // Vector field_vec = Vector::sc_mult(field_dir, scaling_constant); //OLD
    Vector field_vec = scaling_constant * field_dir; // compute final field vector
    return field_vec;
}
Vector Wire::electric_field_strength(Vector const &position){
    return Vector(0,0,0);
}
void Wire::save_to_file(std::ostream &file,  Vector const &particle_init_velocity, long double t, long double dt){ 
    //we are going to use the initial velocity of the particle to scale the line.
    Vector point_on_wire = this->origin;
    long double scaling_constant = Vector::norm(particle_init_velocity);
    for(long i{0}; i<t/dt;i++){
        //now we save the point on the wire
        Vector::save_to_file(file, point_on_wire, CSV_F);
        std::cout<<"wire: "<<i*100/(t/dt)<<'%'<<" \n";
        // Vector::print(point_on_wire, CSV_F);
        // point_on_wire = Vector::add(point_on_wire, Vector::sc_mult(this->wire_direction, 0.00002)); //OLD
        point_on_wire = point_on_wire + (0.00002 * this->wire_direction);
    }
}

