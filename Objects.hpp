// Objects
#include "Vector.hpp"

#ifndef OBJECTS_H
#define OBJECTS_H

class Object;
class Wire;
class Uniform_Electric_Field;
class Uniform_Magnetic_Field;
class Static_Point_Charge;

typedef class Object{
public:
    std::string name;
    virtual Vector magnetic_field_strength(Vector const &position) = 0;
    virtual Vector electric_field_strength(Vector const &position) = 0;
} Object;

//various objects that can be simulated
typedef class Uniform_Magnetic_Field : public Object {
    Vector magnetic_field_vector;
public:
    Uniform_Magnetic_Field(Vector const &direction, long double const &magnitude);
    // Uniform_Magnetic_Field();
    Vector magnetic_field_strength(Vector const &position);
    Vector electric_field_strength(Vector const &position);
} Uniform_Magnetic_Field;

typedef class Uniform_Electric_Field : public Object {
    Vector electric_field_vector;
public:
    Uniform_Electric_Field(Vector const &direction, long double const &magnitude);
    // Uniform_Electric_Field() : direction(Vector(0,0,0)), magnitude(0) {} //default constructor
    Vector magnetic_field_strength(Vector const &position);
    Vector electric_field_strength(Vector const &position);
} Uniform_Electric_Field;

typedef class Static_Point_Charge : public Object {
    long double charge;
    Vector origin;
    long double k;
public:
    Static_Point_Charge(Vector const &origin, long double charge);
    Vector magnetic_field_strength(Vector const &position);
    Vector electric_field_strength(Vector const &position);
} Static_Point_Charge;

typedef class Wire : public Object{
    Vector wire_direction; 
    Vector origin;
    std::string name;
    long double mu_0; // permeability of free space
    long double i_wire; // current through the wire

public:
    Wire(Vector const &origin, Vector const &wire_direction, long double i_wire);
    Wire();
    //return the field vector at a given position
    Vector magnetic_field_strength(Vector const &position);
    Vector electric_field_strength(Vector const &position);
    void save_to_file(std::ostream &file,  Vector const &particle_init_velocity, long double t, long double dt);
} Wire;

#endif