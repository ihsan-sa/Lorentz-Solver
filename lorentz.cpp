#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>

#define PI 3.141592

typedef enum VFormat_Opt{
    NL, //new line
    NR, //no return
    CSV_F, //CSV Format
}VFormat_Opt;
typedef class Vector{
    long double x;
    long double y;
    long double z;

public:
    Vector(long double x, long double y, long double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Vector() : x(0), y(0), z(0){}
    static Vector cross(Vector v1, Vector v2){

        long double new_x = (v1.y * v2.z) - (v1.z * v2.y);
        long double new_y = (v1.z * v2.x) - (v1.x * v2.z);  
        long double new_z = (v1.x * v2.y) - (v1.y * v2.x);
        Vector result(new_x, new_y, new_z);

        return result;
    }
    static long double dot(Vector const &v1, Vector const &v2){
        return (v1.x*v2.x) + (v1.y*v2.y) + (v1.z*v2.z);
    }
    static Vector add(Vector const &v1, Vector const &v2){
        return Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
    }
    static Vector sub(Vector const &v1, Vector const &v2){
        return Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
    }
    static Vector sc_mult(Vector const &v1, long double k){
        return Vector(v1.x * k, v1.y * k, v1.z * k);
    }
    static void print(Vector const &v){
        std::cout<<"[ "<<v.x<<", "<<v.y<<", "<<v.z<<" ]"<<std::endl;
    }
    static void print(Vector const &v, VFormat_Opt opt ){
        if(opt == NL) std::cout<<"[ "<<v.x<<", "<<v.y<<", "<<v.z<<" ]"<<std::endl;
        else if(opt == NR) std::cout<<"[ "<<v.x<<", "<<v.y<<", "<<v.z<<" ]  ";
        else if(opt == CSV_F) std::cout<<v.x<<", "<<v.y<<", "<<v.z<<std::endl;
    }
    static void save_to_file(std::ostream &file, Vector const &v, VFormat_Opt opt){
        if(opt == NL) file<< "[ "<<v.x<<", "<<v.y<<", "<<v.z<<" ]"<<std::endl;
        else if(opt == NR) file<<"[ "<<v.x<<", "<<v.y<<", "<<v.z<<" ]  ";
        else if(opt == CSV_F) file<<v.x<<", "<<v.y<<", "<<v.z<<std::endl;
    }
    static long double norm(Vector const &v){
        return sqrt((v.x*v.x) + (v.y*v.y) + (v.z*v.z));
    }
} Vector;

typedef class Magnetic_field_generic{
public:
    virtual Vector field_vector(Vector const &position) = 0;
}Magnetic_field_generic;


typedef class Wire_Magnetic_Field : public Magnetic_field_generic{
    Vector wire_direction; 
    Vector origin;
    long double mu_0; // permeability of free space
    long double i_wire; // current through the wire

public:
    Wire_Magnetic_Field(Vector const &origin, Vector const &wire_direction, long double i_wire, long double mu_0){
        this->i_wire = i_wire;
        this->mu_0 = mu_0;
        this->wire_direction = wire_direction;
        this->origin = origin;
    }
    //return the field vector at a given position
    Vector field_vector(Vector const &position){  
        //first, we need to figure out the distance from the wire, r.
        //we know that the wire is a line with direction vector wire_direction and that the particle is at a point described by the position vector.
        //we use the perp formula
        //consider O the origin of wire and P the position of particle

        Vector OP = Vector::sub(position, this->origin);
        Vector perp = Vector::sub(OP, Vector::sc_mult(this->wire_direction, (Vector::dot(OP,this->wire_direction))/(Vector::dot(wire_direction, wire_direction))));
        long double radius = Vector::norm(perp); // dist from the wire to point

        //now, recall B = k/r
        long double k_wire = mu_0*i_wire/(2*PI);  //shorthand for the B equation: B = (mu/2PI)*(I/r) = k/r

        long double field_strength = k_wire/radius; // calculate field strength magnitude
        Vector field_dir = Vector::cross(this->wire_direction, perp);  //determine the field direction vector
        long double scaling_constant = field_strength/Vector::norm(field_dir);  //determine the scaling constant required to get the right magnitude

        Vector field_vec = Vector::sc_mult(field_dir, scaling_constant); // compute final field vector
        return field_vec;
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
            point_on_wire = Vector::add(point_on_wire, Vector::sc_mult(this->wire_direction, 0.00002));
        }
    }
} Wire_Magnetic_Field;


typedef class Particle{
    Vector position;
    Vector velocity;
    long double charge;
    long double mass;

public:
    Particle(Vector const &init_pos, Vector const &init_vel, long double charge, long double mass){
        this->position = init_pos;
        this->velocity = init_vel;
        this->mass = mass;
        this->charge = charge;
    }
    //for the moment, we just pass it one magnetic field
    Vector lorentz_force(Magnetic_field_generic &m){
        Vector field_vec  = m.field_vector(this->position);
        Vector force = Vector::cross(Vector::sc_mult(this->velocity,this->charge), field_vec);
        return force;
    } 
    //compute acceleration for a time t
    Vector compute_position(Magnetic_field_generic &m, long double dt){
        Vector force = this->lorentz_force(m);
        Vector acceleration = Vector::sc_mult(force, (1/this->mass));
        // Vector::print(force, 1);
        // Vector::print(acceleration, 1);

        Vector dS = Vector::add(Vector::sc_mult(this->velocity, dt), Vector::sc_mult(acceleration, 0.5*pow(dt,2)));
        // Vector::print(dS, 1);
        this->position = Vector::add(this->position, dS);
        // Vector::print(this->velocity, 1);
        this->velocity = Vector::add(this->velocity, Vector::sc_mult(acceleration, dt));
        // Vector::print(this->velocity, 1);
        return position;
    }
} Particle;

typedef class Uniform_Magnetic_Field : public Magnetic_field_generic{
    Vector magnetic_field;

public: 
    Uniform_Magnetic_Field(Vector const &magnetic_field){
        this->magnetic_field =  magnetic_field;
    }
    Vector field_vector(Vector const &position) {
        return this->magnetic_field;
    }

} Uniform_Magnetic_Field;


int main(){

    //simulation
    
    //initialize the csv
    std::ofstream Data("data.csv");
    Data<<"x, y,  z"<<std::endl;
    std::ofstream Wire_data("wire_data.csv");
    Wire_data<<"x, y,  z"<<std::endl;
    
    long double t = 3e-4;
    long double dt = 1e-8;

    //creating the magnetic field
    long double mu_0 = 4*PI*pow(10,-7);
    long double current = 1000;
    Vector wire_origin(0,0,-0.2);
    Vector wire_direction(0,0,1);


    Wire_Magnetic_Field m1(wire_origin, wire_direction, current, mu_0); //create wire magnetic field 
    
    //creating a uniform magnetic field
    Vector uniform_field(0,10e-7,10e-4);
    Uniform_Magnetic_Field m_uniform(uniform_field);

    //creating the particle
    Vector velocity(-1e3/3, -3e3/3, 0);
    Vector position(0.5,0.5,0);
    long double q = 1.6*pow(10, -19);
    long double m_p = 1.672621898*pow(10, -27);
    Particle p1(position, velocity, q, m_p);

    //print the initial lorentz force
    Vector::print(p1.lorentz_force(m1));
    std::cout<<Vector::norm(p1.lorentz_force(m1))<<"\n"<<std::endl;


    //simulation

    auto start = std::chrono::high_resolution_clock::now();

    // m1.save_to_file(Wire_data, velocity, t, dt);

    auto middle = std::chrono::high_resolution_clock::now();

    for(long i{0};i<t/dt;i++){
        Vector pos = p1.compute_position(m_uniform, dt);
        std::cout<<"Soln: "<<i*100/(t/dt)<<'%'<<"\n";
        // Vector::print(pos, CSV_F);
        Vector::save_to_file(Data, pos, CSV_F);
    }

    auto end = std::chrono::high_resolution_clock::now();

    auto wire_time = std::chrono::duration_cast<std::chrono::microseconds>(middle - start);
    auto sol_time = std::chrono::duration_cast<std::chrono::microseconds>(end - middle);
    auto total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout<<std::endl<<"Wire computation: "<<wire_time.count()*10e-7<<"s Soln Time: "<<sol_time.count()*10e-7<<"s Total time: "<<total_time.count()*10e-7<<"s"<<std::endl;

    //example 1:
    // Vector::print(m1.field_vector(Vector(0.1,0,0)));
    // std::cout<<Vector::norm(m1.field_vector(Vector(0.1,0,0)))<<std::endl;


    //other example:
    // std::cout<<std::endl;
    // Wire_Magnetic_Field m2(Vector(0,0,0), Vector(0,0,1), 5, mu_0);
    // Wire_Magnetic_Field m3(Vector(0.1,0,0), Vector(0,0,-1), 7, mu_0);

    // Vector sum = Vector::add(m2.field_vector(Vector(0.05,0,0)), m3.field_vector(Vector(0.05,0,0)));
    // Vector::print(sum);
    // std::cout<<Vector::norm(sum)<<std::endl;

    Data.close();
    return 0;
}