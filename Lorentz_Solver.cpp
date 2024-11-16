#include <iostream>
#include <cmath>

/*

Ihsan S. Novemeber 2024

 --- OBJECTS: 

A simulation space contains objects (and later charged particles):
    - Wires (magnetic field)
        - Current 
        - Origin
        - Direction
    - Static point charges (electric field)
        - Charge
        - Origin
    - Uniform electric fields
        - Scaled field vector (pass in direction and magnitude)
    - Uniform magnetic fields
        - Scaled field vector (pass in direction and magnitude)
    - Magnetic dipoles (natural ie magnet) 
        - Origin (center of dipole)

        not user defined:
        - Dipole moment

    - Loops of wire (Magnetic dipole)
        - Current
        - Origin (center of circle)
        - Radius

        not user defined:
        - Dipole moment

All objects inherit an Object class. They have a methods named:
    - magnetic_field_strength(...)
    - electric_field_strength(...)


 --- SPACES:

Inside of the space, there is a dynamic array which holds all items type Object. 
Charged particles will also be added to spaces and allowed to interact later.

**Currently chsrged particles can be added to a space. 
The method simulate(...) can be called to simulate the path of the particle.

Spaces have methods:
    - magnetic_field_strength(...)
    - electric field_strength(...)
which return the sum of all the magnetic and electric fields respectively of the objects within it at a given position.

 --- PARTICLES:

Particles have a:   
    - mass
    - charge
    - velocity vector
    - position vector

Particles also have a method simulate(...) which computes the path the particle takes within a given space. 

 --- CONSTRUCTING AND RUNNING A SIMULATION

 1. Create a Space
 2. Create Objects and add them to the Space
 3. Create a Particle
 4. Set simulation parameters including:
    - Time step (dt)
    - Duration (time)
 5. Either add Particle to the space and call [space_name].simulate(...) OR call [particle_name].simulate([space_name], ...)
 

*/

typedef class Space{
    

} Space;




int main(){

    return 0;
}