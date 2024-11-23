#include "Vector.hpp"
#include "Objects.hpp"
#include "Space.hpp"
#include "Particle.hpp"
#include <iostream>

void sim1_UMF(){
    //Start by opening the data file and initializing it

    std::ofstream Data("data.csv");
    
    Data<<"x,y,z"<<std::endl;

    //creating a space
    Space space1;

    //creating objects
    Vector m1_direction(0,0,1);
    long double m1_strength = 1;
    Uniform_Magnetic_Field m1(m1_direction, m1_strength);

    //add object to space
    space1.add_object(m1);

    //create a particle
    Vector p1_init_pos(0,0,0);
    Vector p1_init_vel(1,1,1);
    long double p1_mass = 1;
    long double p1_charge = 1;
    Particle p1(p1_init_pos, p1_init_vel, p1_charge, p1_mass, "p1");

    //set simulation parameters
    long double t = 20;
    long double dt = 0.1;

    //run simulation
    p1.simulate(space1, t, dt, Data, RK4_HYBRID);

}

//test with both uniform magnetic and electric fields

void sim2_UMF_UEF(){

    //creating a space
    Space space1;

    //creating objects
    Vector m1_direction(0,0,1);
    long double m1_strength = 1;
    Uniform_Magnetic_Field m1(m1_direction, m1_strength);

    Vector e1_direction(1,0,0);
    long double e1_strength = 0.1;
    Uniform_Electric_Field e1(e1_direction, e1_strength);

    //add object to space
    space1.add_object(m1);
    space1.add_object(e1);

    //create a particle
    Vector p1_init_pos(0,0,0);
    Vector p1_init_vel(1,1,1);
    long double p1_mass = 1;
    long double p1_charge = 1;
    Particle p1(p1_init_pos, p1_init_vel, p1_charge, p1_mass, "p1");

    space1.add_object(p1);

    //set simulation parameters
    long double t = 20;
    long double dt = 0.1;

    //run simulation
    space1.simulate(t, dt, RK4_HYBRID);
    space1.print_objects();

}

//test with uniform mag field, uniform elec field and point charge
void sim3_UMF_UEF_PC(){
    //Start by opening the data file and initializing it
    std::ofstream Data("data.csv");
    Data<<"x,y,z"<<std::endl;
    //creating a space
    Space space1;

    //creating objects
    Vector m1_direction(0,0,1);
    long double m1_strength = 0;
    Uniform_Magnetic_Field m1(m1_direction, m1_strength);

    Vector e1_direction(1,0,0);
    long double e1_strength = 0.0;
    Uniform_Electric_Field e1(e1_direction, e1_strength);

    Static_Point_Charge e2(Vector(0,-600,100),0.0005);
    Static_Point_Charge e3(Vector(0,-400,750),-0.0005);
    Static_Point_Charge e4(Vector(-3000,0,0),0.001);

    //add object to space
    space1.add_object(m1);
    space1.add_object(e1);
    space1.add_object(e2);
    space1.add_object(e3);
    space1.add_object(e4);

    //create a particle
    Vector p1_init_pos(0,0,0);
    Vector p1_init_vel(0,0,0);
    long double p1_mass = 1;
    long double p1_charge = -1;
    Particle p1(p1_init_pos, p1_init_vel, p1_charge, p1_mass, "p1");

    space1.add_object(p1);
    //set simulation parameters
    long double t = 500;
    long double dt = 0.01;

    //run simulation
    space1.simulate(t,dt,RK4_HYBRID);
}

//particle orbiting a static charge

void sim4(){
//Start by opening the data file and initializing it


    //creating a space
    Space space1;

    //creating objects
    Static_Point_Charge e1(Vector(0,0,0), -10e-4);
    Static_Point_Charge e2(Vector(-3,4,3.5), -10e-4);


    //add object to space
    space1.add_object(e1);
    space1.add_object(e2);

    //create a particle
    Vector p1_init_pos(1,0,0);
    Vector p1_init_vel(0,1,1);
    long double p1_mass = 100;
    long double p1_charge = 10e-6;
    Particle p1(p1_init_pos, p1_init_vel, p1_charge, p1_mass, "p1");

    //set simulation parameters
    long double t = 20;
    long double dt = 0.0001;

    space1.add_object(p1);
    //run simulation
    space1.simulate(t, dt, RK4_HYBRID);


}

//

void primitive_sim(){

    //creating the magnetic field
    long double current = 1000;
    Vector wire_origin(0,0,-0.2);
    Vector wire_direction(0,0,1);

    Wire m1(wire_origin, wire_direction, current); //create wire magnetic field 

    //examples from the Gianocoli textbook

    // example 1:
    Vector::print(m1.magnetic_field_strength(Vector(0.1,0,0)));
    std::cout<<Vector::norm(m1.magnetic_field_strength(Vector(0.1,0,0)))<<std::endl;

    // other example:
    std::cout<<std::endl;
    Wire m2(Vector(0,0,0), Vector(0,0,1), 5);
    Wire m3(Vector(0.1,0,0), Vector(0,0,-1), 7);

    Vector sum = Vector::add(m2.magnetic_field_strength(Vector(0.05,0,0)), m3.magnetic_field_strength(Vector(0.05,0,0)));
    Vector::print(sum);
    std::cout<<Vector::norm(sum)<<std::endl;
}

//sim1 updated for simulation in space
void sim5_UMF(){
    //Start by opening the data file and initializing it

    std::ofstream Data("data.csv");
    
    Data<<"x,y,z"<<std::endl;

    //creating a space
    Space space1;

    //creating objects
    Vector m1_direction(0,0,1);
    long double m1_strength = 1;
    Uniform_Magnetic_Field m1(m1_direction, m1_strength);

    //add object to space
    space1.add_object(m1);

    //create a particle
    Vector p1_init_pos(0,0,0);
    Vector p1_init_vel(1,1,1);
    long double p1_mass = 1;
    long double p1_charge = 1;
    Particle p1(p1_init_pos, p1_init_vel, p1_charge, p1_mass, "p1");

    space1.add_object(p1);

    //set simulation parameters
    long double t = 20;
    long double dt = 0.1;

    space1.print_objects();
    //run simulation
    space1.simulate(t, dt, RK4_HYBRID);
}


void sim_3_particle(){

    Particle p1(Vector(1,0,0),Vector(0,0,0), 0.00005, 1, "p1");
    Particle p2(Vector(-1,0,0),Vector(0,0,0), -0.00005, 1, "p2");
    Particle p3(Vector(0,0,-1),Vector(0,0,0), -0.0001, 1, "p3");
    Particle p4(Vector(0,0,0),Vector(0,0,0), 0.000001, 1, "p4");


    Space space1;
    space1.add_object(p1);
    space1.add_object(p2);
    space1.add_object(p3);
    space1.add_object(p4);

    long double t = 3;
    long double dt = 0.0001;
    space1.simulate(t, dt, RK4_HYBRID);
}
