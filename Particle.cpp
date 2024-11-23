#include "Vector.hpp"
#include "Objects.hpp"
#include "Particle.hpp"
#include "Space.hpp"
#include <iostream>
#include <cmath>



Particle::Particle(Vector const &init_pos, Vector const &init_vel, long double charge, long double mass, std::string name){
    position_ = init_pos;
    velocity_ = init_vel;
    mass_ = mass;
    charge_ = charge;
    next_position_ = init_pos; //initialize next pos with something random
    this->name = name;
    std::string file_name = name + ".csv";
    Data_ = std::ofstream(file_name);
    Data_<<"x,y,z"<<std::endl;
    approx_nbr_ = 1;

    //set past positions to the correct values
    past_positions_[0] = init_pos;
    past_velocities_[0] = init_vel;
    for(int i{1}; i<3;i++){
        past_positions_ [i] = Vector(0,0,0);
        past_velocities_[i] = Vector(0,0,0);
    }
}
Vector Particle::magnetic_field_strength(Vector const &position){
    return Vector(0,0,0);
}
Vector Particle::electric_field_strength(Vector const &position){

    // Vector OP = Vector::sub(position, position_); //OLD
    Vector OP = position  - position_; //defining a vector from the origin to the point in space
    long double r = Vector::norm(OP); //compute distance from point charge
    //equation is E = kq/r^2. 
    long double magnitude = (k_ * charge_)/std::pow(r, 2);
    OP = Vector::adjust_mag(OP, magnitude);

    return OP;
}
void Particle::save_position(){
    Vector::save_to_file(Data_, position_, CSV_F);
}


//Space class for defining a simulation space


//definition for Particle class

Vector Particle::lorentz_force(Space &space, Vector const &position, Vector const &velocity, Lorentz_Calculation_Opt opt){

    // //OLD WAY
    // Vector f_e = Vector::sc_mult(space.electric_field_strength(position, this->name), charge_);
    // Vector f_b = Vector::cross(Vector::sc_mult(velocity, charge_), space.magnetic_field_strength(position));
    // Vector net_force = Vector::add(f_e, f_b);

    Vector f_e = charge_ * space.electric_field_strength(position, this->name);
    Vector f_b = (charge_ * velocity) * space.magnetic_field_strength(position);
    Vector net_force = f_e + f_b;


    if(opt == NO_MASS) return net_force;
    else if(opt == WITH_MASS) return Vector::sc_mult(net_force, 1/mass_);
    else return Vector(0,0,0); //FIX ERROR HANDLING
}
//Function to compute position using kinematic equations
void Particle::compute_position(Space &space, long double dt){
    Vector acceleration = this->lorentz_force(space, position_, velocity_, WITH_MASS);
    // Vector::print(force, 1);
    // Vector::print(acceleration, 1);

    // Vector dS = Vector::add(Vector::sc_mult(velocity_, dt), Vector::sc_mult(acceleration, 0.5*pow(dt,2))); //OLD
    Vector dS = (dt * velocity_) + (0.5*pow(dt,2) * acceleration);
    // Vector::print(dS, 1);
    // next_position_ = Vector::add(position_, dS); //OLD
    next_position_ = position_ + dS;
    // Vector::print(velocity_, 1);
    // velocity_ = Vector::add(velocity_, Vector::sc_mult(acceleration, dt)); //OLD
    velocity_ = velocity_ + (dt * acceleration);
    // Vector::print(velocity_, 1);
}

//Function to compute position using RK4 method
Vector Particle::compute_velocity(Space  &space, long double dt){
    Vector a1 = this->lorentz_force(space, position_, velocity_, WITH_MASS);
    Vector a2 = this->lorentz_force(space, position_, (velocity_ + (a1 * (dt/2))), WITH_MASS);
    Vector a3 = this->lorentz_force(space, position_, (velocity_ + (a2 * (dt/2))), WITH_MASS);
    Vector a4 = this->lorentz_force(space, position_, (velocity_ + (a3 * dt)), WITH_MASS);

    Vector computed_velocity = velocity_ + ((a1 + (2 * a2) + (2 * a3) + a4) * (dt/6));
    return computed_velocity;
}
Vector Particle::compute_pos_update(Space  &space, long double dt){
    Vector a1 = this->lorentz_force(space, position_, velocity_, WITH_MASS);
    Vector a2 = this->lorentz_force(space, position_, (velocity_ + (a1 * (dt/2))), WITH_MASS);
    Vector a3 = this->lorentz_force(space, position_, (velocity_ + (a2 * (dt/2))), WITH_MASS);
    Vector a4 = this->lorentz_force(space, position_, (velocity_ + (a3 * dt)), WITH_MASS);

    Vector s1 = velocity_;
    Vector s2 = velocity_ + a1*(dt/2);
    Vector s3 = velocity_ + a2*(dt/2);
    Vector s4 = velocity_ + a3*(dt);

    Vector computed_position = position_ + (dt/6)*(s1 + (2*s2) + (2*s3) + s4);

    return computed_position;
}
void Particle::compute_position_RK4_HYBRID(Space &space, long double dt){
    
    Vector velocity_1 = compute_velocity(space, dt); //computing velocity at time t_1

    std::cout<<"v_i: "<<past_velocities_[0]<<" v_i-1: "<<past_velocities_[1]<<" v_i-2: "<<past_velocities_[2]<<std::endl;
    std::cout<<"s_i: "<<past_positions_[0]<<" s_i-1: "<<past_positions_[1]<<" s_i-2: "<<past_positions_[2]<<std::endl;



    //shift past velocities
    past_velocities_[2] = past_velocities_[1]; 
    past_velocities_[1] = past_velocities_[0];
    past_velocities_[0] = velocity_1;

    Vector position_1;
    //now, we compute the new position using Serge's algorithm

    if(approx_nbr_ == 1){
        //this means that we now have to compute the position for the first time
        Vector v_0 = velocity_; //init velocity
        Vector v_1 = velocity_1; //velocity at time t_1
        Vector v_2 = compute_velocity(space, dt*2); //velocity at time t_2
        Vector v_3 = compute_velocity(space, dt*3); //velocity at time t_3

        position_1 = position_ + ((dt/24) * ((9*v_0) + (19*v_1) - (5*v_2) + v_3)); //position at time t_0

    }else{

        position_1 = past_positions_[2] + ((dt/3)*(past_velocities_[2] + (4*past_velocities_[1]) + velocity_1));

    }

    //set past positions
    past_positions_[2] = past_positions_[1];
    past_positions_[1] = past_positions_[0];
    past_positions_[0] = position_1;

    approx_nbr_ ++;

    //

    // Vector v1 = velocity_;
    // Vector v2 = compute_velocity(space, dt/2);
    // Vector v3 = compute_velocity(space, dt/2);
    // Vector v4 = compute_velocity(space, dt);

    // Vector s1 = position_ + (dt/6)*(v1 + (2*v2) + (2*v3) + v4); //use this for the revised RK4 approx of position
    // Vector s1 = position_ + dt*velocity_1;    //use this for the old approx
    // Vector s1 = position_ + dt*velocity_;  //use this for an approx using the current val of velocity.
    // Vector s1 = compute_pos_update(space, dt); //use this for an approx using another "tandem" RK4

    velocity_ = velocity_1;
    next_position_ = position_1;
}   
Vector Particle::update_position(){
    position_ = next_position_;
    return position_;
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