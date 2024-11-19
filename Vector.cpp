#include "Vector.hpp"

#include <iostream>
#include <cassert>
#include <cmath>



Vector::Vector(long double x, long double y, long double z){
    this->x = x;
    this->y = y;
    this->z = z;
}
Vector::Vector() : x(0), y(0), z(0){}
Vector Vector::cross(Vector v1, Vector v2){

    long double new_x = (v1.y * v2.z) - (v1.z * v2.y);
    long double new_y = (v1.z * v2.x) - (v1.x * v2.z);  
    long double new_z = (v1.x * v2.y) - (v1.y * v2.x);
    Vector result(new_x, new_y, new_z);

    return result;
}
long double Vector::dot(Vector const &v1, Vector const &v2){
    return (v1.x*v2.x) + (v1.y*v2.y) + (v1.z*v2.z);
}
Vector Vector::add(Vector const &v1, Vector const &v2){
    return Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}
Vector Vector::sub(Vector const &v1, Vector const &v2){
    return Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}
Vector Vector::sc_mult(Vector const &v1, long double k){
    return Vector(v1.x * k, v1.y * k, v1.z * k);
}
void Vector::print(Vector const &v){
    std::cout<<"[ "<<v.x<<", "<<v.y<<", "<<v.z<<" ]"<<std::endl;
}
void Vector::print(Vector const &v, VFormat_Opt opt ){
    if(opt == NL) std::cout<<"[ "<<v.x<<", "<<v.y<<", "<<v.z<<" ]"<<std::endl;
    else if(opt == NR) std::cout<<"[ "<<v.x<<", "<<v.y<<", "<<v.z<<" ]  ";
    else if(opt == CSV_F) std::cout<<v.x<<", "<<v.y<<", "<<v.z<<std::endl;
}
void Vector::save_to_file(std::ostream &file, Vector const &v, VFormat_Opt opt){
    if(opt == NL) file<< "[ "<<v.x<<", "<<v.y<<", "<<v.z<<" ]"<<std::endl;
    else if(opt == NR) file<<"[ "<<v.x<<", "<<v.y<<", "<<v.z<<" ]  ";
    else if(opt == CSV_F) file<<v.x<<", "<<v.y<<", "<<v.z<<std::endl;
}
long double Vector::norm(Vector const &v){
    return sqrt((v.x*v.x) + (v.y*v.y) + (v.z*v.z));
}
Vector Vector::adjust_mag(Vector const &v, long double const &mag){
    long double norm = Vector::norm(v);
    if(norm == 0){
        std::cout<<"ERR norm == 0"<<std::endl;
        assert(norm != 0);
    }
    long double scaling_const = mag/norm;
    return Vector::sc_mult(v, scaling_const);
}
long double Vector::x_get(Vector const &v) {return v.x;}
long double Vector::y_get(Vector const &v) {return v.y;}
long double Vector::z_get(Vector const &v) {return v.z;}

Vector Vector::operator+(Vector const &v2){ //addition
    return Vector(this->x + v2.x, this->y + v2.y, this->z + v2.z);
}
Vector Vector::operator+(Vector &v2){ //addition
    return Vector(this->x + v2.x, this->y + v2.y, this->z + v2.z);
}
Vector Vector::operator-(Vector const &v2){  //subtraction
    return Vector(this->x - v2.x, this->y - v2.y, this->z - v2.z);
}
Vector Vector::operator-(Vector &v2){  //subtraction
    return Vector(this->x - v2.x, this->y - v2.y, this->z - v2.z);
}
Vector Vector::operator*(Vector v2){ //cross prod

    long double new_x = (this->y * v2.z) - (this->z * v2.y);
    long double new_y = (this->z * v2.x) - (this->x * v2.z);  
    long double new_z = (this->x * v2.y) - (this->y * v2.x);
    Vector result(new_x, new_y, new_z);

    return result;
}
long double Vector::operator%(Vector const &v2){ //dot prod
    return (this->x*v2.x) + (this->y*v2.y) + (this->z*v2.z);
}
Vector Vector::operator*(long double k){ //sc mult
    return Vector(this->x * k, this->y * k, this->z * k);
}


//define vector ops
Vector operator*(long double k, Vector const &v1){ //sc mult
    return Vector(Vector::x_get(v1) * k, Vector::y_get(v1) * k,Vector::z_get(v1) * k);
}
std::ostream &operator<<(std::ostream &output, Vector const &v){ //print
    output <<"[ "<<Vector::x_get(v)<<", "<<Vector::y_get(v)<<", "<<Vector::z_get(v)<<" ]";
    return output;
}

Vector operator-(Vector const &v1, Vector v2){  //subtraction
    return Vector(Vector::x_get(v1)- Vector::x_get(v2), Vector::y_get(v1)- Vector::y_get(v2), Vector::z_get(v1)- Vector::z_get(v2));
}
Vector operator+(Vector const &v1, Vector v2){  //subtraction
    return Vector(Vector::x_get(v1) + Vector::x_get(v2), Vector::y_get(v1) + Vector::y_get(v2), Vector::z_get(v1) + Vector::z_get(v2));
}

    

