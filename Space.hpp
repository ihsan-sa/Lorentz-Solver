#ifndef SPACE_H
#define SPACE_H

#include <vector>
#include "Particle.hpp"

class Particle;

//Enums used to define function/simulation options

class Space;

typedef class Space : public Object{
    std::vector<Object*> objects;
    std::vector<Particle*> particles;
public:
    void add_object(Object &object);
    void add_object(Particle &particle);
    Vector magnetic_field_strength(Vector const &position);
    Vector electric_field_strength(Vector const &position);
    Vector electric_field_strength(Vector const &position, std::string name);
    void delete_space();
    void simulate(long double t, long double dt, Simulation_Type opt);
    void print_objects();

} Space;



#endif