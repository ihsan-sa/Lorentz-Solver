#include <iostream>
#include <cmath>
#include <vector>
#include <memory>
#include <fstream>
#include <chrono>

/*

Ihsan S. Novemeber 2024


*/

//Create Vector class
typedef enum VFormat_Opt{
    NL, //new line
    NR, //no return
    CSV_F, //CSV Format
}VFormat_Opt;

typedef enum Lorentz_Calculation_Opt{
    NO_MASS,
    WITH_MASS,
}Lorentz_Calculation_Opt;

typedef enum Simulation_Type{
    KINEMATIC, 
    RK4_HYBRID
} Simulation_Type;
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
    static Vector adjust_mag(Vector const &v, long double const &mag){
        long double scaling_const = mag/Vector::norm(v);
        return Vector::sc_mult(v, scaling_const);
    }
    static long double x_get(Vector const &v) {return v.x;}
    static long double y_get(Vector const &v) {return v.y;}
    static long double z_get(Vector const &v) {return v.z;}
} Vector;

typedef class Object{
public:
    virtual Vector magnetic_field_strength(Vector const &position) = 0;
    virtual Vector electric_field_strength(Vector const &position) = 0;
} Object;

typedef class Uniform_Magnetic_Field : public Object {
    Vector magnetic_field_vector;
public:
    Uniform_Magnetic_Field(Vector const &direction, long double const &magnitude){ //passing by reference to reduce amount of memory used
        this->magnetic_field_vector = Vector::adjust_mag(direction, magnitude);
    }
    // Uniform_Magnetic_Field() : direction(Vector(0,0,0)), magnitude(0) {} //default constructor
    Vector magnetic_field_strength(Vector const &position){
        return this->magnetic_field_vector;
    }
    Vector electric_field_strength(Vector const &position){
        return Vector(0,0,0);
    }
} Uniform_Magnetic_Field;


typedef class Space : public Object{
    std::vector<Object*> objects;
public:
    void add_object(Object &object){ //const or not?
        objects.push_back(&object);
    }
    Vector magnetic_field_strength(Vector const &position){
        Vector magnetic_field_strength(0,0,0);
        for(Object* object : objects){
            magnetic_field_strength = Vector::add(magnetic_field_strength, object->magnetic_field_strength(position));
        }
        return magnetic_field_strength;
    }
    Vector electric_field_strength(Vector const &position){
        Vector electric_field_strength(0,0,0);
        for(Object* object : objects){
            electric_field_strength = Vector::add(electric_field_strength, object->electric_field_strength(position));
        }
        return electric_field_strength;
    }
    void delete_space(){
        for(Object* object : objects){
            delete object;
        }
    }

} Space;

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
    Vector lorentz_force(Space &space, Vector const &position, Vector const &velocity, Lorentz_Calculation_Opt opt){
        Vector f_e = Vector::sc_mult(space.electric_field_strength(position), this->charge);
        Vector f_b = Vector::cross(Vector::sc_mult(velocity, this->charge), space.magnetic_field_strength(position));
        Vector net_force = Vector::add(f_e, f_b);
        if(opt == NO_MASS) return net_force;
        else if(opt == WITH_MASS) return Vector::sc_mult(net_force, 1/this->mass);
        else return Vector(0,0,0); //FIX ERROR HANDLING
    }
    //Function to compute position using kinematic equations
    Vector compute_position(Space &space, long double dt){
        Vector acceleration = this->lorentz_force(space, this->position, this->velocity, WITH_MASS);
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

    //Function to compute position using RK4 method
    Vector compute_position_RK4_HYBRID(Space &space, long double dt){
        Vector v1 = this->lorentz_force(space, this->position, this->velocity, WITH_MASS);
        Vector v2 = this->lorentz_force(space, this->position, Vector::add(this->velocity, Vector::sc_mult(v1, dt/2)), WITH_MASS);
        Vector v3 = this->lorentz_force(space, this->position, Vector::add(this->velocity, Vector::sc_mult(v2, dt/2)), WITH_MASS);
        Vector v4 = this->lorentz_force(space, this->position, Vector::add(this->velocity, Vector::sc_mult(v3, dt)), WITH_MASS);
        this->velocity = Vector::add(this->velocity, Vector::sc_mult(Vector::add(Vector::add(v1, Vector::sc_mult(v2, 2)),Vector::add(Vector::sc_mult(v3,2),v4)), dt/6));

        Vector s1 = Vector::add(this->position, Vector::sc_mult(this->velocity, dt));
        this->position = s1;
        return this->position;
    }   

    void simulate(Space &space, long double t, long double dt, std::ofstream &Data, Simulation_Type sim_type){
        //repeat code so that it only needs to compare sim type once
        if (sim_type == KINEMATIC){
            for(long i{0};i<t/dt;i++){
                Vector pos = this->compute_position(space, dt);
                std::cout<<""<<i*100/(t/dt)<<'%'<<"\n";
                // Vector::print(pos, CSV_F);
                Vector::save_to_file(Data, pos, CSV_F);
            }
        }
        if (sim_type == RK4_HYBRID){
            for(long i{0};i<t/dt;i++){
                Vector pos = this->compute_position_RK4_HYBRID(space, dt);
                std::cout<<""<<i*100/(t/dt)<<'%'<<"\n";
                // Vector::print(pos, CSV_F);
                Vector::save_to_file(Data, pos, CSV_F);
            }
        }
    }


} Particle;


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
    Particle p1(p1_init_pos, p1_init_vel, p1_charge, p1_mass);

    //set simulation parameters
    long double t = 20;
    long double dt = 0.1;

    //run simulation
    p1.simulate(space1, t, dt, Data, RK4_HYBRID);

}



int main(){

    sim1_UMF();
    return 0;
}