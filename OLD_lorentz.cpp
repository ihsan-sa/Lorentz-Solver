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
    static Vector adjust_mag(Vector &v, long double const mag){
        long double scaling_const = mag/Vector::norm(v);
        return Vector::sc_mult(v, scaling_const);
    }
    static long double x_get(Vector const &v) {return v.x;}
    static long double y_get(Vector const &v) {return v.y;}
    static long double z_get(Vector const &v) {return v.z;}
} Vector;

typedef class Field_Generic{
public:
    virtual Vector field_vector(Vector const &position) = 0;
}Field_Generic;



typedef class Wire_Magnetic_Field : public Field_Generic{
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
    Wire_Magnetic_Field() : mu_0(0), origin(Vector(0,0,0)), wire_direction(Vector(0,0,0)), i_wire(0) {}
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
typedef class Uniform_Magnetic_Field : public Field_Generic{
    Vector magnetic_field;

public: 
    Uniform_Magnetic_Field(Vector const &magnetic_field){
        this->magnetic_field =  magnetic_field;
    }
    Uniform_Magnetic_Field() : magnetic_field(Vector(0,0,0)){}
    Vector field_vector(Vector const &position) {
        return this->magnetic_field;
    }

} Uniform_Magnetic_Field;

typedef class Uniform_Electric_Field : public Field_Generic{
    Vector electric_field;
public:
    Uniform_Electric_Field(Vector const &electric_field){
        this->electric_field = electric_field;
    }
    Uniform_Electric_Field() : electric_field(Vector(0,0,0)) {}
    Vector field_vector(Vector const &position){
        return this->electric_field;
    }

}Uniform_Electric_Field;

typedef class Radial_Electric_Field : public Field_Generic{
    Vector origin;
    long double k = 8.99e9;
    long double charge;
public:
    Radial_Electric_Field(Vector const &origin, long double charge){
        this->origin = origin;
        this->charge = charge;
    }
    Radial_Electric_Field() : origin(Vector(0,0,0)), k(8.99e9), charge(0){}
    Vector field_vector(Vector const &position){
        Vector OP = Vector::sub(position, this->origin);
        long double r = Vector::norm(OP);
        long double magnitude = (this->k * this->charge)/pow(r, 2);

        OP = Vector::adjust_mag(OP, magnitude);
        return OP;
    }
} Radial_Electric_Field;

//We will now define a simulation space which contains a number of fields.
typedef class Space{
    //add dynamic mem later
    int n_WMF = 0; //number of wire magnetic fields
    int n_UMF = 0;
    int n_UEF = 0;
    int n_REF = 0;
    Wire_Magnetic_Field wire_magnetic_fields[10];
    Uniform_Magnetic_Field uniform_magnetic_fields[10];
    Uniform_Electric_Field uniform_electric_fields[10];
    Radial_Electric_Field radial_electric_fields[10];

public:
    Space(){
        this->n_REF = 0;
        this->n_UEF = 0;
        this->n_UMF = 0;
        this->n_WMF = 0;
    }
    //adding fields 
    void add_WMF(Wire_Magnetic_Field const &m){
        this->wire_magnetic_fields[n_WMF] = m;
        n_WMF ++;
    }
    void add_UMF(Uniform_Magnetic_Field const &m){
    this->uniform_magnetic_fields[n_UMF] = m;
    n_UMF ++;
}
    void add_UEF(Uniform_Electric_Field const &m){
    this->uniform_electric_fields[n_UEF] = m;
    n_UEF ++;
}
    void add_REF(Radial_Electric_Field const &m){
        this->radial_electric_fields[n_REF] = m;
        n_REF ++;
    }

    Vector electric_field(Vector const &position){
        Vector net_electric_field(0,0,0);
        //loop through the elctric fields and add them up
        for(int i{0}; i<n_UEF;i++){
            net_electric_field = Vector::add(net_electric_field, uniform_electric_fields[i].field_vector(position));
        }
        for(int i{0}; i<n_REF;i++){
            net_electric_field = Vector::add(net_electric_field, radial_electric_fields[i].field_vector(position));
        }
        return net_electric_field;
    }
    Vector magnetic_field(Vector const &position){
        Vector net_magnetic_field(0,0,0);
        //loop through the elctric fields and add them up
        for(int i{0}; i<n_UMF;i++){
            net_magnetic_field = Vector::add(net_magnetic_field, uniform_magnetic_fields[i].field_vector(position));
        }
        for(int i{0}; i<n_WMF;i++){
            net_magnetic_field = Vector::add(net_magnetic_field, wire_magnetic_fields[i].field_vector(position));
        }
        return net_magnetic_field;
    }

}Space;

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
    Vector lorentz_force(Field_Generic &m){
        Vector field_vec  = m.field_vector(this->position);
        Vector force = Vector::cross(Vector::sc_mult(this->velocity,this->charge), field_vec);
        return force;
    } 
    //compute acceleration for a time t
    Vector compute_position(Field_Generic &m, long double dt){
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
    //temporarily overload the function for testing 
    Vector lorentz_force(Space &space){
        // Vector force = Vector::add(Vector::sc_mult(space.electric_field(this->position), this->charge), Vector::cross(Vector::sc_mult(this->velocity,this->charge), space.magnetic_field(this->position)));
        Vector F_e = Vector::sc_mult(space.electric_field(this->position), this->charge);
        Vector F_b = Vector::cross(Vector::sc_mult(this->velocity,this->charge), space.magnetic_field(this->position));
        Vector force = Vector::add(F_e,F_b);
        return force;
    } 
    //other version used for RK4
    Vector lorentz_force(Space &space, Vector const &pos){
        // Vector force = Vector::add(Vector::sc_mult(space.electric_field(this->position), this->charge), Vector::cross(Vector::sc_mult(this->velocity,this->charge), space.magnetic_field(this->position)));
        Vector F_e = Vector::sc_mult(space.electric_field(pos), this->charge);
        Vector F_b = Vector::cross(Vector::sc_mult(this->velocity,this->charge), space.magnetic_field(pos));
        Vector force = Vector::add(F_e,F_b);
        return force;
    } 
    //compute acceleration for a time t
    Vector compute_position(Space &space, long double dt){
        Vector force = this->lorentz_force(space);
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
    Vector compute_positionRK4(Space &space, long double dt){

        //velocity calculations

        Vector f1 = Vector::sc_mult(this->lorentz_force(space), 1/this->mass);
        long double m1 = Vector::x_get(f1);  //x component
        long double k1 = Vector::y_get(f1);  //y component
        long double l1 = Vector::z_get(f1);  //z component

        Vector f2 = Vector::sc_mult(this->lorentz_force(space, Vector::add(this->position, Vector::sc_mult(Vector(m1, k1, l1), dt/2))), 1/this->mass);
        long double m2 = Vector::x_get(f2); 
        long double k2 = Vector::y_get(f2); 
        long double l2 = Vector::z_get(f2); 

        Vector f3 = Vector::sc_mult(this->lorentz_force(space, Vector::add(this->position, Vector::sc_mult(Vector(m2, k2, l2), dt/2))), 1/this->mass);
        long double m3 = Vector::x_get(f3); 
        long double k3 = Vector::y_get(f3); 
        long double l3 = Vector::z_get(f3);

        Vector f4 = Vector::sc_mult(this->lorentz_force(space, Vector::add(this->position, Vector::sc_mult(Vector(m3, k3, l3), dt))), 1/this->mass);
        long double m4 = Vector::x_get(f4); 
        long double k4 = Vector::y_get(f4); 
        long double l4 = Vector::z_get(f4);

        // Vector new_velocity = Vector::add(this->velocity, );






        // Vector force = this->lorentz_force(space);
        // Vector acceleration = Vector::sc_mult(force, (1/this->mass));
        // // Vector::print(force, 1);
        // // Vector::print(acceleration, 1);

        // Vector dS = Vector::add(Vector::sc_mult(this->velocity, dt), Vector::sc_mult(acceleration, 0.5*pow(dt,2)));
        // // Vector::print(dS, 1);
        // this->position = Vector::add(this->position, dS);
        // // Vector::print(this->velocity, 1);
        // this->velocity = Vector::add(this->velocity, Vector::sc_mult(acceleration, dt));
        // // Vector::print(this->velocity, 1);
        return Vector(0,0,0); //CHANGE THIS!
    }

} Particle;


void sim2_wirefield(){
    
    //initialize the csv
    std::ofstream Data("data.csv");
    Data<<"x, y,  z"<<std::endl;
    std::ofstream Wire_data("wire_data.csv");
    Wire_data<<"x, y,  z"<<std::endl;
    
    //setting the simulation time
    long double t = 27;
    long double dt = 0.001;

    //creating the magnetic field
    long double mu_0 = 4*PI*pow(10,-7);
    long double current = 1000;
    Vector wire_origin(0,0,-0.2);
    Vector wire_direction(0,0,1);

    Wire_Magnetic_Field m1(wire_origin, wire_direction, current, mu_0); //create wire magnetic field 

    //creating the particle
    Vector velocity(8.3, 5.2, 1.2);
    Vector position(0.0,0.0,0);
    long double q = 1.6*pow(10, -19);
    long double m_p = 1.672621898*pow(10, -27);
    Particle p1(position, velocity, 1, 0.475);


    //print the initial lorentz force

    // Vector::print(p1.lorentz_force(m1));
    // std::cout<<Vector::norm(p1.lorentz_force(m1))<<"\n"<<std::endl;


    //simulation

    //save wire info the file
    auto start = std::chrono::high_resolution_clock::now();
    m1.save_to_file(Wire_data, velocity, t, dt);

    //run sim
    auto middle = std::chrono::high_resolution_clock::now();
    for(long i{0};i<t/dt;i++){
        Vector pos = p1.compute_position(m1, dt);
        std::cout<<"Soln: "<<i*100/(t/dt)<<'%'<<"\n";
        // Vector::print(pos, CSV_F);
        Vector::save_to_file(Data, pos, CSV_F);
    }

    auto end = std::chrono::high_resolution_clock::now();

    //compute run time
    auto wire_time = std::chrono::duration_cast<std::chrono::microseconds>(middle - start);
    auto sol_time = std::chrono::duration_cast<std::chrono::microseconds>(end - middle);
    auto total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout<<std::endl<<"Wire computation: "<<wire_time.count()*10e-7<<"s Soln Time: "<<sol_time.count()*10e-7<<"s Total time: "<<total_time.count()*10e-7<<"s"<<std::endl;

    //close files
    Data.close();
    Wire_data.close();
}

void sim1_ufield(){

    //compare to https://www.geogebra.org/m/xpRMzPgc
    
    //csv setup
    std::ofstream Data("data.csv");
    Data<<"x, y,  z"<<std::endl;
    
    //set simulation time
    long double t = 27;
    long double dt = 0.001;
    
    //creating a uniform magnetic field
    Vector uniform_field(0,0,0.55);
    Uniform_Magnetic_Field m_uniform(uniform_field);

    //creating the first particle
    Vector velocity(8.3, 5.2, 1.2);
    Vector position(0.0,0.0,0);
    long double q = 1.6*pow(10, -19); //charge of a proton
    long double m_p = 1.672621898*pow(10, -27); //mass of a proton
    Particle p1(position, velocity, 1, 0.475);

    //creating the second particle
    Vector velocity2(-10, 0, 1.3);
    Vector position2(75,0,0);
    Particle p2(position2, velocity2, -1, 0.475);


    //simulation

    //first particle
    auto start_particle1 = std::chrono::high_resolution_clock::now();
    for(long i{0};i<t/dt;i++){
        Vector pos = p1.compute_position(m_uniform, dt);
        std::cout<<"Particle 1: "<<i*100/(t/dt)<<'%'<<"\n";
        // Vector::print(pos, CSV_F);
        Vector::save_to_file(Data, pos, CSV_F);
    }

    //second particle
    auto start_particle2 = std::chrono::high_resolution_clock::now();
    for(long i{0};i<t/dt;i++){
        Vector pos = p2.compute_position(m_uniform, dt);
        std::cout<<"Particle 2: "<<i*100/(t/dt)<<'%'<<"\n";
        // Vector::print(pos, CSV_F);
        Vector::save_to_file(Data, pos, CSV_F);
    }

    auto end = std::chrono::high_resolution_clock::now();

    //computing elapsed time
    auto time_particle1 = std::chrono::duration_cast<std::chrono::microseconds>(start_particle2 - start_particle1);
    auto time_particle2 = std::chrono::duration_cast<std::chrono::microseconds>(end - start_particle2);
    auto total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start_particle1);
    std::cout<<std::endl<<"Particle 1 Time: "<<time_particle1.count()*10e-7<<"s Particle 2 Time: "<<time_particle2.count()*10e-7<<"s Total time elapsed: "<<total_time.count()*10e-7<<"s"<<std::endl;

    //close file
    Data.close();
}

void primitize_sim(){

    //creating the magnetic field
    long double mu_0 = 4*PI*pow(10,-7);
    long double current = 1000;
    Vector wire_origin(0,0,-0.2);
    Vector wire_direction(0,0,1);

    Wire_Magnetic_Field m1(wire_origin, wire_direction, current, mu_0); //create wire magnetic field 

    //examples from the Gianocoli textbook

    // example 1:
    Vector::print(m1.field_vector(Vector(0.1,0,0)));
    std::cout<<Vector::norm(m1.field_vector(Vector(0.1,0,0)))<<std::endl;

    // other example:
    std::cout<<std::endl;
    Wire_Magnetic_Field m2(Vector(0,0,0), Vector(0,0,1), 5, mu_0);
    Wire_Magnetic_Field m3(Vector(0.1,0,0), Vector(0,0,-1), 7, mu_0);

    Vector sum = Vector::add(m2.field_vector(Vector(0.05,0,0)), m3.field_vector(Vector(0.05,0,0)));
    Vector::print(sum);
    std::cout<<Vector::norm(sum)<<std::endl;
}

void sim3_ufield(){
    
    //csv setup
    std::ofstream Data("data.csv");
    Data<<"x, y,  z"<<std::endl;
    
    //set simulation time
    long double t = 27;
    long double dt = 0.001;
    
    //creating a uniform magnetic field
    Vector uniform_field(1,0,0);
    Uniform_Magnetic_Field m_uniform(uniform_field);

    //creating the first particle
    Vector velocity(1,0,10);
    Vector position(0.0,0.0,0);
    long double q = 1.6*pow(10, -19); //charge of a proton
    long double m_p = 1.672621898*pow(10, -27); //mass of a proton
    Particle p1(position, velocity, -1, 0.475);



    //simulation

    //first particle
    auto start_particle1 = std::chrono::high_resolution_clock::now();
    for(long i{0};i<t/dt;i++){
        Vector pos = p1.compute_position(m_uniform, dt);
        std::cout<<i*100/(t/dt)<<'%'<<"\n";
        // Vector::print(pos, CSV_F);
        Vector::save_to_file(Data, pos, CSV_F);
    }

    
    auto end = std::chrono::high_resolution_clock::now();

    //computing elapsed time
    auto total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start_particle1);
    
    std::cout<<std::endl<<"Total time elapsed: "<<total_time.count()*10e-7<<"s"<<std::endl;

    //close file
    Data.close();
}

void sim4_ufields(){
    
    //csv setup
    std::ofstream Data("data.csv");
    Data<<"x, y,  z"<<std::endl;
    
    //set simulation time
    long double t = 27;
    long double dt = 0.001;
    
    //creating a uniform magnetic field
    Vector uniform_field_1(1,0,0);
    Vector uniform_field_2(0,1.5,0);
    Vector uniform_field_3(0,0,2);
    Uniform_Magnetic_Field m_uniform1(uniform_field_1);
    Uniform_Magnetic_Field m_uniform2(uniform_field_2);
    Uniform_Magnetic_Field m_uniform3(uniform_field_3);


    //creating the space;
    Space space1;

    space1.add_UMF(m_uniform1);
    space1.add_UMF(m_uniform2);
    space1.add_UMF(m_uniform3);

    //creating the first particle
    Vector velocity(1,0,10);
    Vector position(0.0,0.0,0);
    long double q = 1.6*pow(10, -19); //charge of a proton
    long double m_p = 1.672621898*pow(10, -27); //mass of a proton
    Particle p1(position, velocity, -1, 0.475);



    //simulation

    //first particle
    auto start_particle1 = std::chrono::high_resolution_clock::now();
    for(long i{0};i<t/dt;i++){
        Vector pos = p1.compute_position(space1, dt);
        std::cout<<i*100/(t/dt)<<'%'<<"\n";
        // Vector::print(pos, CSV_F);
        Vector::save_to_file(Data, pos, CSV_F);
    }

    
    auto end = std::chrono::high_resolution_clock::now();

    //computing elapsed time
    auto total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start_particle1);
    
    std::cout<<std::endl<<"Total time elapsed: "<<total_time.count()*10e-7<<"s"<<std::endl;

    //close file
    Data.close();
}


void sim6_pointcharge(){
    
    //csv setup
    std::ofstream Data("data.csv");
    Data<<"x, y,  z"<<std::endl;
    
    //set simulation time
    long double t = 100;
    long double dt = 10e-5;
    
    //creating a uniform magnetic field
    Uniform_Electric_Field e1(Vector(0,0,0));
    Radial_Electric_Field e2(Vector(0,-600,100),0.0005);
    Radial_Electric_Field e3(Vector(0,-400,750),-0.0005);
    Radial_Electric_Field e4(Vector(-3000,0,0),0.001);
    Data<<"0,0,0"<<std::endl;
    Data<<"0,4,4"<<std::endl;
    long double mu_0 = 4*PI*pow(10,-7);
    Uniform_Magnetic_Field m1(Vector(0,0,0.1));
    Wire_Magnetic_Field m2(Vector(-4000, 3000,0),Vector(1,2,-2), 10000, mu_0);

    


    //creating the space;
    Space space1;
    space1.add_UEF(e1);
    space1.add_REF(e2);
    space1.add_REF(e3);
    space1.add_REF(e4);
    space1.add_UMF(m1);
    // space1.add_WMF(m2);
   

    //creating the first particle
    Vector velocity(0,0,0);
    Vector position(0,0,0);
    // Data<<"1,0,0"<<std::endl;
    long double q = 1.6*pow(10, -19); //charge of a proton
    long double m_p = 1.672621898*pow(10, -27); //mass of a proton
    Particle p1(position, velocity, -1, 1);



    //simulation

    //first particle
    auto start_particle1 = std::chrono::high_resolution_clock::now();
    for(long i{0};i<t/dt;i++){
        Vector pos = p1.compute_position(space1, dt);
        std::cout<<i*100/(t/dt)<<'%'<<"\n";
        // Vector::print(pos, CSV_F);
        Vector::save_to_file(Data, pos, CSV_F);
    }

    
    auto end = std::chrono::high_resolution_clock::now();

    //computing elapsed time
    auto total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start_particle1);
    
    std::cout<<std::endl<<"Total time elapsed: "<<total_time.count()*10e-7<<"s"<<std::endl;

    //close file
    Data.close();
}

void sim5_ufield(){
    
    //csv setup
    std::ofstream Data("data.csv");
    Data<<"x, y,  z"<<std::endl;
    
    //set simulation time
    long double t = 27;
    long double dt = 0.001;
    
    //creating a uniform magnetic field
    Vector uniform_field(1,0,0);
    Uniform_Magnetic_Field m_uniform(uniform_field);

    //creating the first particle
    Vector velocity(1,0,10);
    Vector position(0.0,0.0,0);
    long double q = 1.6*pow(10, -19); //charge of a proton
    long double m_p = 1.672621898*pow(10, -27); //mass of a proton
    Particle p1(position, velocity, -1, 0.475);



    //simulation

    //first particle
    auto start_particle1 = std::chrono::high_resolution_clock::now();
    for(long i{0};i<t/dt;i++){
        Vector pos = p1.compute_position(m_uniform, dt);
        std::cout<<i*100/(t/dt)<<'%'<<"\n";
        // Vector::print(pos, CSV_F);
        Vector::save_to_file(Data, pos, CSV_F);
    }

    
    auto end = std::chrono::high_resolution_clock::now();

    //computing elapsed time
    auto total_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start_particle1);
    
    std::cout<<std::endl<<"Total time elapsed: "<<total_time.count()*10e-7<<"s"<<std::endl;

    //close file
    Data.close();
}

void sim7_plot_field(){
    //plotting the electric field due to a point charge

    std::ofstream Data("field.csv");
    Data<<"x,y,z,magnitude,"<<std::endl;

    Radial_Electric_Field e1(Vector(0,0,0), 1);
    std::cout<<Vector::norm(e1.field_vector(Vector(1,0,0)))<<std::endl;
    std::cout<<Vector::norm(e1.field_vector(Vector(2,0,0)))<<std::endl;

    long double small_bound = -0.2;
    long double big_bound = -0.001;
    long double spacing = 0.001;

    for(long double x{small_bound}; x<big_bound; x+=spacing){
        for(long double y{small_bound}; y<big_bound; y+=spacing){
            for(long double z{0}; z<spacing; z+=spacing){
                std::cout<<x<<", "<<y<<", "<<z<<", "<<Vector::norm(e1.field_vector(Vector(x,y,z)))<<std::endl;
                Data<<x<<","<<y<<","<<z<<","<<Vector::norm(e1.field_vector(Vector(x,y,z)))<<","<<std::endl;
            }
        }
    }


}

int main(){
    std::cout<<"OLD SIM"<<std::endl;
    primitize_sim();
    std::cout<<"OLD SIM"<<std::endl;

    return 0;
}