#include "Vector.hpp"
#include "Objects.hpp"
#include "Particle.hpp"
#include "Space.hpp"
#include <iostream>
#include <cmath>



Particle::Particle(Vector const &init_pos, Vector const &init_vel, long double charge, long double mass, std::string name){
    this->position = init_pos;
    this->velocity = init_vel;
    this->mass = mass;
    this->charge = charge;
    this->next_position = init_pos; //initialize next pos with something random
    this->name = name;
    std::string file_name = name + ".csv";
    this->Data = std::ofstream(file_name);
    this->Data<<"x,y,z"<<std::endl;
}
Vector Particle::magnetic_field_strength(Vector const &position){
    return Vector(0,0,0);
}
Vector Particle::electric_field_strength(Vector const &position){

    // Vector OP = Vector::sub(position, this->position); //OLD
    Vector OP = position  - this->position; //defining a vector from the origin to the point in space
    long double r = Vector::norm(OP); //compute distance from point charge
    //equation is E = kq/r^2. 
    long double magnitude = (this->k * this->charge)/std::pow(r, 2);
    OP = Vector::adjust_mag(OP, magnitude);

    return OP;
}
void Particle::save_position(){
    Vector::save_to_file(this->Data, this->position, CSV_F);
}


//Space class for defining a simulation space


//definition for Particle class

Vector Particle::lorentz_force(Space &space, Vector const &position, Vector const &velocity, Lorentz_Calculation_Opt opt){

    // //OLD WAY
    // Vector f_e = Vector::sc_mult(space.electric_field_strength(position, this->name), this->charge);
    // Vector f_b = Vector::cross(Vector::sc_mult(velocity, this->charge), space.magnetic_field_strength(position));
    // Vector net_force = Vector::add(f_e, f_b);

    Vector f_e = this->charge * space.electric_field_strength(position, this->name);
    Vector f_b = (this->charge * velocity) * space.magnetic_field_strength(position);
    Vector net_force = f_e + f_b;


    if(opt == NO_MASS) return net_force;
    else if(opt == WITH_MASS) return Vector::sc_mult(net_force, 1/this->mass);
    else return Vector(0,0,0); //FIX ERROR HANDLING
}
//Function to compute position using kinematic equations
void Particle::compute_position(Space &space, long double dt){
    Vector acceleration = this->lorentz_force(space, this->position, this->velocity, WITH_MASS);
    // Vector::print(force, 1);
    // Vector::print(acceleration, 1);

    // Vector dS = Vector::add(Vector::sc_mult(this->velocity, dt), Vector::sc_mult(acceleration, 0.5*pow(dt,2))); //OLD
    Vector dS = (dt * this->velocity) + (0.5*pow(dt,2) * acceleration);
    // Vector::print(dS, 1);
    // this->next_position = Vector::add(this->position, dS); //OLD
    this->next_position = this->position + dS;
    // Vector::print(this->velocity, 1);
    // this->velocity = Vector::add(this->velocity, Vector::sc_mult(acceleration, dt)); //OLD
    this->velocity = this->velocity + (dt * acceleration);
    // Vector::print(this->velocity, 1);
}

//Function to compute position using RK4 method
void Particle::compute_position_RK4_HYBRID(Space &space, long double dt){
    Vector v1 = this->lorentz_force(space, this->position, this->velocity, WITH_MASS);
    //OLD
    // Vector v2 = this->lorentz_force(space, this->position, Vector::add(this->velocity, Vector::sc_mult(v1, dt/2)), WITH_MASS);
    // Vector v3 = this->lorentz_force(space, this->position, Vector::add(this->velocity, Vector::sc_mult(v2, dt/2)), WITH_MASS);
    // Vector v4 = this->lorentz_force(space, this->position, Vector::add(this->velocity, Vector::sc_mult(v3, dt)), WITH_MASS);

    Vector v2 = this->lorentz_force(space, this->position, (this->velocity + (v1 * (dt/2))), WITH_MASS);
    Vector v3 = this->lorentz_force(space, this->position, (this->velocity + (v2 * (dt/2))), WITH_MASS);
    Vector v4 = this->lorentz_force(space, this->position, (this->velocity + (v3 * dt)), WITH_MASS);

    //OLD
    // this->velocity = Vector::add(this->velocity, Vector::sc_mult(Vector::add(Vector::add(v1, Vector::sc_mult(v2, 2)),Vector::add(Vector::sc_mult(v3,2),v4)), dt/6));
    this->velocity = this->velocity + ((v1 + (2 * v2) + (2 * v3) + v4) * (dt/6));

    // Vector s1 = Vector::add(this->position, Vector::sc_mult(this->velocity, dt)); //OLD
    Vector s1 = this->position + (dt * this->velocity);

    this->next_position = s1;
}   
Vector Particle::update_position(){
    this->position = this->next_position;
    return this->position;
}
void Particle::simulate(Space &space, long double t, long double dt, std::ofstream &Data, Simulation_Type sim_type){
    //repeat code so that it only needs to compare sim type once
    if (sim_type == KINEMATIC){
        for(long i{0};i<t/dt;i++){
            this->compute_position(space, dt);
            Vector pos = update_position();
            std::cout<<""<<i*100/(t/dt)<<'%'<<"\n";
            // Vector::print(pos, CSV_F);
            Vector::save_to_file(Data, pos, CSV_F);
        }
    }
    if (sim_type == RK4_HYBRID){
        for(long i{0};i<t/dt;i++){
            this->compute_position_RK4_HYBRID(space, dt);
            Vector pos = update_position();
            std::cout<<""<<i*100/(t/dt)<<'%'<<"\n";
            // Vector::print(pos, CSV_F);
            Vector::save_to_file(Data, pos, CSV_F);
        }
    }
}