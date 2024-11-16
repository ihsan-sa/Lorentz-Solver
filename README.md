# Lorentz Solver

This is a Lorentz Solver which simulates particles moving through a space with various objects within it. 

## Files in this repository

Use `Lorentz_Solver.cpp` as it is the most capable solver.

Read instructions below, then compile and run the program.

Then run the `parse_simulation_data.m` MATLAB script to parse and plot the data.

## Overview of the code structure

### **Objects**

A simulation _space_ contains _objects_ with specified attributes:
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


All _objects_ inherit an _Object_ class. They have a methods named:
- `magnetic_field_strength(...)`
- `electric_field_strength(...)`

_Spaces_ can also contain other _spaces_ since _Spaces_ also inherit the _Object_ class.

### **Spaces**

Inside of the _space_, there is a dynamic array which holds all items type _Object*_. 
Charged particles will also be added to spaces and allowed to interact later.

_Spaces_ have methods:
- `magnetic_field_strength(...)`
- `electric field_strength(...)`
  
which return the sum of all the magnetic and electric fields respectively of the objects within it at a given position.

### **Particles**

_Particles_ have a:   
- mass
- charge
- velocity vector
- position vector

_Particles_ also have a method `simulate(...)` which computes the path the particle takes within a given space. 

### Vectors

Vectors are defined in a _Vector_ class.

To create a vector, use `Vector [vector_name](x, y, z)` or `Vector(x, y, z)`, where `x`, `y` and `z` are of type `long double`.

## Constructing and running a simulation

**1. Initialize the data collection file**

Use `std::ofstream [file_name]("data.csv");` and then

`[file_name]<<"x,y,z"<<std::endl;` to write the column headers.

**2. Create a _Space_**
  
Create a _Space_ by using `Space [space_name];`

**3. Create _Objects_ and add them to the Space**

Create an _Object_ using `[Object_Type] [object_name](...);`.

Here are the available objects and their definitions:

- `Wire [wire_name](Vector const &origin, Vector const &wire_direction, long double i_wire);`
- `Static_Point_Charge [point_charge_name](Vector const &origin, long double charge);`
- `Uniform_Electric_Field [field_name](Vector const &direction, long double const &magnitude);`
- `Uniform_Magnetic_Field [field_name](Vector const &direction, long double const &magnitude);`
 
**4. Create a Particle**

Create a particle as follows:

`Particle [particle_name](Vector const &init_pos, Vector const &init_vel, long double charge, long double mass);`

**5. Set simulation parameters including:**

Set simulation time `long double t` and step size `long double dt`.

**6. Simulation**

run `simulate(Space &space, long double t, long double dt, std::ofstream &Data, Simulation_Type sim_type)`, where `Data` is your output file and `sim_type` is either `KINEMATIC` or `RK4_HYBRID`.

**7. Graphing**

Run the matlab script to graph the data in the csv output file. You may need to add the output file to your path. 

If you add a _Wire_ to your simulation, you may want to plot it. In that case, plot the _wire_data.csv_ as well. You may need to adjust the scaling factors in the `save_to_file` method in _Wire_.

## Example simulations

There are a few example simulations in the form of functions in the `Lorentz_Solver.cpp` file.

Simulation output when `sim1_UMF()` function is run and data is plotted.
![image](https://github.com/user-attachments/assets/db68fe48-ac94-45b9-a07d-4f6c6d0b815e)

Simulation output when `sim2_UMF_UEF()` function is run and data is plotted.
![image](https://github.com/user-attachments/assets/1daee448-6eff-4d93-8f61-c63c3924351d)

Simulation output when `sim4()` function is run and data is plotted.
![image](https://github.com/user-attachments/assets/e46e35d8-797a-41cd-af5b-585258268f1b)



Below is an image of a simulation of two charged particles in a uniform magnetic field in the z-direction.

![Screenshot 2024-11-14 145923](https://github.com/user-attachments/assets/6f2a0060-970b-4c92-913d-30f6483274c4)
