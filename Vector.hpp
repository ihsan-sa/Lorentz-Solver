// Vector header file
#ifndef VECTOR_H
#define VECTOR_H

#include <fstream>

class Vector;

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
    Vector(long double x, long double y, long double z);
    Vector();
    static Vector cross(Vector v1, Vector v2);
    static long double dot(Vector const &v1, Vector const &v2);
    static Vector add(Vector const &v1, Vector const &v2);
    static Vector sub(Vector const &v1, Vector const &v2);
    static Vector sc_mult(Vector const &v1, long double k);
    static void print(Vector const &v);
    static void print(Vector const &v, VFormat_Opt opt );
    static void save_to_file(std::ostream &file, Vector const &v, VFormat_Opt opt);
    static long double norm(Vector const &v);
    static Vector adjust_mag(Vector const &v, long double const &mag);
    static long double x_get(Vector const &v);
    static long double y_get(Vector const &v);
    static long double z_get(Vector const &v);

    Vector operator+(Vector const &v2);
    Vector operator+(Vector &v2);
    Vector operator-(Vector const &v2);
    Vector operator-(Vector &v2);
    Vector operator*(Vector v2);
    
    long double operator%(Vector const &v2);
    Vector operator*(long double k);
    
} Vector;

//define vector ops

Vector operator*(long double k, Vector const &v1);
std::ostream &operator<<(std::ostream &output, Vector const &v);
Vector operator-(Vector const &v1, Vector v2);
Vector operator+(Vector const &v1, Vector v2);

#endif