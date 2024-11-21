#ifndef PARTICLE_H
#define PARTICLE_H


typedef enum Lorentz_Calculation_Opt{
    NO_MASS,
    WITH_MASS,
}Lorentz_Calculation_Opt;

typedef enum Simulation_Type{
    KINEMATIC, 
    RK4_HYBRID
} Simulation_Type;

class Space;
class Particle;

typedef class Particle : public Object{
    Vector position_;
    Vector next_position_;
    Vector velocity_;
    long double charge_;
    long double mass_;
    long double k_ = 8.99e9;
    std::ofstream Data_;
public:
    Particle(Vector const &init_pos, Vector const &init_vel, long double charge, long double mass, std::string name);
    Vector magnetic_field_strength(Vector const &position);
    Vector electric_field_strength(Vector const &position);
    Vector lorentz_force(Space &space, Vector const &position, Vector const &velocity, Lorentz_Calculation_Opt opt);
    void compute_position(Space &space, long double dt);
    void compute_position_RK4_HYBRID(Space &space, long double dt);
    Vector compute_velocity(Space &space, long double dt);
    Vector update_position();
    void simulate(Space &space, long double t, long double dt, std::ofstream &Data, Simulation_Type sim_type);
    void save_position();

} Particle;

#endif