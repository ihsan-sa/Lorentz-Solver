
#include "Vector.hpp"
#include "Objects.hpp"
#include "Space.hpp"
#include "Particle.hpp"
#include <chrono>
#include <iostream>


void Space::add_object(Object &object){ //const or not?
    objects.push_back(&object);
}
void Space::add_object(Particle &particle){ //const or not?
    objects.push_back(&particle);
    particles.push_back(&particle);
}
Vector Space::magnetic_field_strength(Vector const &position){
    Vector magnetic_field_strength(0,0,0);
    for(Object* object : objects){
        // magnetic_field_strength = Vector::add(magnetic_field_strength, object->magnetic_field_strength(position)); //OLD
        magnetic_field_strength = magnetic_field_strength + object->magnetic_field_strength(position);
    }
    return magnetic_field_strength;
}
Vector Space::electric_field_strength(Vector const &position){
    Vector electric_field_strength(0,0,0);
    for(Object* object : objects){
        // electric_field_strength = Vector::add(electric_field_strength, object->electric_field_strength(position)); //OLD
        electric_field_strength = electric_field_strength + object->electric_field_strength(position);
    }
    return electric_field_strength;
}
Vector Space::electric_field_strength(Vector const &position, std::string name){
    Vector electric_field_strength(0,0,0);
    for(Object* object : objects){
        if(object->name != name){
            // electric_field_strength = Vector::add(electric_field_strength, object->electric_field_strength(position)); //OLD
            electric_field_strength = electric_field_strength + object->electric_field_strength(position);

        }
    }
    return electric_field_strength;
}
void Space::delete_space(){
    for(Object* object : objects){
        delete object;
    }
    for(Particle* particle : particles){
        delete particle;
    }
}
void Space::simulate(long double t, long double dt, Simulation_Type opt){
    auto start = std::chrono::high_resolution_clock::now();
    if(opt == RK4_HYBRID){
        for(long i{0};i<t/dt;i++){
            for(Particle* particle : particles){
                particle->compute_position_RK4_HYBRID(*this, dt);
            }
            for(Particle* particle : particles){
                particle->update_position();
                particle->save_position();
            }
            std::cout<<""<<i*100/(t/dt)<<'%'<<"\n";
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout<<"Elapsed time: "<<std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()*10e-7<<"s"<<std::endl;

}
void Space::print_objects(){
    for(Object * obj : objects){
        std::cout<<obj->name<<std::endl;
    }
}

