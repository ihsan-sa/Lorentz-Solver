# Lorentz Solver

This is a Lorentz Solver which allows for the simulation of particles moving through a space with various objects within it. 

 --- OBJECTS: 

A simulation space contains objects with specified attributes:
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
        - Dipole moment 
    - Loops of wire (Magnetic dipole)
        - Current
        - Origin (center of circle)
        - Radius
        - Dipole moment (NOT user defined)

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
 

Below is an image of a simulation of two charged particles in a uniform magnetic field in the z-direction.

![Screenshot 2024-11-14 145923](https://github.com/user-attachments/assets/6f2a0060-970b-4c92-913d-30f6483274c4)
