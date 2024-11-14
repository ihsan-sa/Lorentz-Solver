#include <iostream>
#include <cmath>

#define PI 3.141592

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
    static long double dot(Vector v1, Vector v2){
        return (v1.x*v2.x) + (v1.y*v2.y) + (v1.z*v2.z);
    }
    static Vector add(Vector v1, Vector v2){
        return Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
    }
    static Vector sub(Vector v1, Vector v2){
        return Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
    }
    static Vector sc_mult(Vector v1, long double k){
        return Vector(v1.x * k, v1.y * k, v1.z * k);
    }
    static void print(Vector v){
        std::cout<<"[ "<<v.x<<", "<<v.y<<", "<<v.z<<" ]"<<std::endl;
    }
    static void print(Vector v, bool line){
        if(!line) std::cout<<"[ "<<v.x<<", "<<v.y<<", "<<v.z<<" ]"<<std::endl;
        else std::cout<<"[ "<<v.x<<", "<<v.y<<", "<<v.z<<" ]  ";
    }
    static long double norm(Vector v){
        return sqrt((v.x*v.x) + (v.y*v.y) + (v.z*v.z));
    }
} Vector;

typedef class Wire_Magnetic_Field{
    Vector wire_direction; 
    Vector origin;
    long double mu_0; // permeability of free space
    long double i_wire; // current through the wire

public:
    Wire_Magnetic_Field(Vector origin, Vector wire_direction, long double i_wire, long double mu_0){
        this->i_wire = i_wire;
        this->mu_0 = mu_0;
        this->wire_direction = wire_direction;
        this->origin = origin;
    }
    //return the field vector at a given position
    Vector field_vector(Vector position){  
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
} Wire_Magnetic_Field;


typedef class Particle{
    Vector position;
    Vector velocity;
    long double charge;
    long double mass;

public:
    Particle(Vector init_pos, Vector init_vel, long double charge, long double mass){
        this->position = init_pos;
        this->velocity = init_vel;
        this->mass = mass;
        this->charge = charge;
    }
    //for the moment, we just pass it one magnetic field
    Vector lorentz_force(Wire_Magnetic_Field m){
        Vector field_vec  = m.field_vector(this->position);
        Vector force = Vector::cross(Vector::sc_mult(this->velocity,this->charge), field_vec);
        return force;
    } 
    //compute acceleration for a time t
    Vector compute_position(Wire_Magnetic_Field m, long double dt){
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


int main(){

    //example 1:
    long double q = 1.6*pow(10, -19);
    long double m_p = 1.672621898*pow(10, -27);
    long double mu_0 = 4*PI*pow(10,-7);
    long double dt = 2./1000;

    Wire_Magnetic_Field m1(Vector(0,0,0), Vector(0,0,1), 25, mu_0);
    Vector::print(m1.field_vector(Vector(0.1,0,0)));
    std::cout<<Vector::norm(m1.field_vector(Vector(0.1,0,0)))<<std::endl;

    Particle p1(Vector(0.1,0,0), Vector(1,0,0), q,m_p);
    Vector::print(p1.lorentz_force(m1));
    std::cout<<Vector::norm(p1.lorentz_force(m1))<<"\n"<<std::endl;

    for(int i{0};i<1000;i++){
        Vector::print(p1.compute_position(m1, dt));
    }




    //other example:
    // std::cout<<std::endl;
    // Wire_Magnetic_Field m2(Vector(0,0,0), Vector(0,0,1), 5, mu_0);
    // Wire_Magnetic_Field m3(Vector(0.1,0,0), Vector(0,0,-1), 7, mu_0);

    // Vector sum = Vector::add(m2.field_vector(Vector(0.05,0,0)), m3.field_vector(Vector(0.05,0,0)));
    // Vector::print(sum);
    // std::cout<<Vector::norm(sum)<<std::endl;

    return 0;
}