#include "Vector.hpp"

#include <iostream> //for std I/O
#include <cmath> //for basic math (eg. pow())
#include <vector> //for creating the dynamic storage of objects in the space
#include <fstream> //for writing to files
#include <chrono> //for timing
#include <string>
#include <cassert>

#define PI 3.141592
/*

Ihsan S. Novemeber 2024

*/

//Enums used to define function/simulation options


typedef enum Lorentz_Calculation_Opt{
    NO_MASS,
    WITH_MASS,
}Lorentz_Calculation_Opt;

typedef enum Simulation_Type{
    KINEMATIC, 
    RK4_HYBRID
} Simulation_Type;





//std Object class
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
    Uniform_Magnetic_Field(Vector const &direction, long double const &magnitude){ //passing by reference to reduce amount of memory used
        this->magnetic_field_vector = Vector::adjust_mag(direction, magnitude);
        this->name = "UMF";
    }
    // Uniform_Magnetic_Field() : direction(Vector(0,0,0)), magnitude(0) {} //default constructor
    Vector magnetic_field_strength(Vector const &position){
        return this->magnetic_field_vector;
    }
    Vector electric_field_strength(Vector const &position){
        return Vector(0,0,0);
    }
} Uniform_Magnetic_Field;

typedef class Uniform_Electric_Field : public Object {
    Vector electric_field_vector;
public:
    Uniform_Electric_Field(Vector const &direction, long double const &magnitude){ //passing by reference to reduce amount of memory used
        this->electric_field_vector = Vector::adjust_mag(direction, magnitude);
        this->name = "UEF";
    }
    // Uniform_Magnetic_Field() : direction(Vector(0,0,0)), magnitude(0) {} //default constructor
    Vector magnetic_field_strength(Vector const &position){
        return Vector(0,0,0);
    }
    Vector electric_field_strength(Vector const &position){
        return this->electric_field_vector;
    }
} Uniform_Electric_Field;

typedef class Static_Point_Charge : public Object {
    long double charge;
    Vector origin;
    long double k = 8.99e9;
public:
    Static_Point_Charge(Vector const &origin, long double charge){
        this->charge = charge;
        this->origin = origin;
        this->name = "SPC";
    }
    Vector magnetic_field_strength(Vector const &position){
        return Vector(0,0,0);
    }
    Vector electric_field_strength(Vector const &position){

        // Vector OP = Vector::sub(position, this->origin); //OLD
        Vector OP = position - this->origin; //defining a vector from the origin to the point in space
        long double r = Vector::norm(OP); //compute distance from point charge
        //equation is E = kq/r^2. 
        long double magnitude = (this->k * this->charge)/std::pow(r, 2);
        OP = Vector::adjust_mag(OP, magnitude);

        return OP;
    }
} Static_Point_Charge;

typedef class Wire : public Object{
    Vector wire_direction; 
    Vector origin;
    std::string name;
    long double mu_0 = 4*PI*pow(10,-7); // permeability of free space
    long double i_wire; // current through the wire

public:
    Wire(Vector const &origin, Vector const &wire_direction, long double i_wire){
        this->i_wire = i_wire;
        this->wire_direction = wire_direction;
        this->origin = origin;
        this->name = "Wire";
    }
    Wire() : origin(Vector(0,0,0)), wire_direction(Vector(0,0,0)), i_wire(0) {}
    //return the field vector at a given position
    Vector magnetic_field_strength(Vector const &position){  
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
        long double k_wire = mu_0*i_wire/(2*PI);  //shorthand for the B equation: B = (mu/2PI)*(I/r) = k/r

        long double field_strength = k_wire/radius; // calculate field strength magnitude
        // Vector field_dir = Vector::cross(this->wire_direction, perp);  //OLD
        Vector field_dir = this->wire_direction * perp; //determine the field direction vector
        long double scaling_constant = field_strength/Vector::norm(field_dir);  //determine the scaling constant required to get the right magnitude

        // Vector field_vec = Vector::sc_mult(field_dir, scaling_constant); //OLD
        Vector field_vec = scaling_constant * field_dir; // compute final field vector
        return field_vec;
    }
    Vector electric_field_strength(Vector const &position){
        return Vector(0,0,0);
    }
    void save_to_file(std::ostream &file,  Vector const &particle_init_velocity, long double t, long double dt){ 
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
} Wire;

class Space;
class Particle;

//Particle class for defining a particle
typedef class Particle : public Object{
    Vector position;
    Vector next_position;
    Vector velocity;
    long double charge;
    long double mass;
    long double k = 8.99e9;
    std::ofstream Data;
public:
    Particle(Vector const &init_pos, Vector const &init_vel, long double charge, long double mass, std::string name);
    Vector magnetic_field_strength(Vector const &position);
    Vector electric_field_strength(Vector const &position);
    Vector lorentz_force(Space &space, Vector const &position, Vector const &velocity, Lorentz_Calculation_Opt opt);
    void compute_position(Space &space, long double dt);
    void compute_position_RK4_HYBRID(Space &space, long double dt);
    Vector update_position();
    void simulate(Space &space, long double t, long double dt, std::ofstream &Data, Simulation_Type sim_type);
    void save_position();

}Particle;

//def for Particle class

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
typedef class Space : public Object{
    std::vector<Object*> objects;
    std::vector<Particle*> particles;
public:
    void add_object(Object &object){ //const or not?
        objects.push_back(&object);
    }
    void add_object(Particle &particle){ //const or not?
        objects.push_back(&particle);
        particles.push_back(&particle);
    }
    Vector magnetic_field_strength(Vector const &position){
        Vector magnetic_field_strength(0,0,0);
        for(Object* object : objects){
            // magnetic_field_strength = Vector::add(magnetic_field_strength, object->magnetic_field_strength(position)); //OLD
            magnetic_field_strength = magnetic_field_strength + object->magnetic_field_strength(position);
        }
        return magnetic_field_strength;
    }
    Vector electric_field_strength(Vector const &position){
        Vector electric_field_strength(0,0,0);
        for(Object* object : objects){
            // electric_field_strength = Vector::add(electric_field_strength, object->electric_field_strength(position)); //OLD
            electric_field_strength = electric_field_strength + object->electric_field_strength(position);
        }
        return electric_field_strength;
    }
    Vector electric_field_strength(Vector const &position, std::string name){
        Vector electric_field_strength(0,0,0);
        for(Object* object : objects){
            if(object->name != name){
                // electric_field_strength = Vector::add(electric_field_strength, object->electric_field_strength(position)); //OLD
                electric_field_strength = electric_field_strength + object->electric_field_strength(position);

            }
        }
        return electric_field_strength;
    }
    void delete_space(){
        for(Object* object : objects){
            delete object;
        }
        for(Particle* particle : particles){
            delete particle;
        }
    }
    void simulate(long double t, long double dt, Simulation_Type opt){
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
    void print_objects(){
        for(Object * obj : objects){
            std::cout<<obj->name<<std::endl;
        }
    }

} Space;

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

//initial test for a very simple uniform magnetic field

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
    long double dt = 0.01;

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


int main(){

    sim2_UMF_UEF();


    return 0;
}