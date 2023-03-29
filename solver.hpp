#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <type_traits>
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <exception>

class Vector {
    double orts[3];
public:
    Vector() {
        orts[0] = 0;
        orts[1] = 0;
        orts[2] = 0;
    }
    Vector(double x, double y, double z) {
        orts[0] = x;
        orts[1] = y;
        orts[2] = z;
    }
    double operator[](int i) const {
        int k = i%3;
        return orts[k];
    }
    Vector operator-(const Vector& a) const {
        return Vector(orts[0] - a[0], orts[1] - a[1], orts[2] - a[2]);
    }
    Vector& operator()(double x, double y, double z) {
        orts[0] = x;
        orts[1] = y;
        orts[2] = z;
        return *this;
    }
    Vector& operator()(const Vector& a, const Vector& b) {
        for(int i = 0; i < 3; ++i)
            orts[i] = a[i + 1]*b[i + 2] - a[i + 2]*b[i + 1];
        return *this;
    }
    double operator*(const Vector& a) const {
        double tmp{};
        for(int i = 0; i < 3; ++i)
            tmp += orts[i]*a[i];
        return tmp;
    }
    Vector& operator += (const Vector& a) {
        orts[0] += a.orts[0];
        orts[1] += a.orts[1];
        orts[2] += a.orts[2];
        return *this;
    }
    Vector operator + (const Vector& a) {
        return Vector(orts[0] + a.orts[0], orts[1] + a.orts[1], orts[2] + a.orts[2]);
    }
    bool operator!=(const Vector& a) const {
        if(orts[0] == a[0] &&
            orts[1] == a[1] &&
            orts[2] == a[2]) {
            return true;
        } else return false;
    }
    double length() const {
        return std::sqrt(orts[0]*orts[0] + orts[1]*orts[1] + orts[2]*orts[2]);
    }
    double length(const Vector& a) const {
        return std::sqrt((orts[0] - a[0])*(orts[0] - a[0]) +
                         (orts[1] - a[1])*(orts[1] - a[1]) +
                         (orts[2] - a[2])*(orts[2] - a[2]));
    }
    Vector Normalize() const {
        return Vector(orts[0]/length(), orts[1]/length(), orts[2]/length());
    }
    static Vector Barycentric(const std::vector<Vector>& star) {
        double sX = .0;
        double sY = 0.;
        double sZ = 0.;
        double AllArea = 0.;
        double x1 = star[0].orts[0];
        double y1 = star[0].orts[1];
        double z1 = star[0].orts[2];
        for(size_t i = 2; i < star.size(); ++i) {
                double x2 = star[i - 1].orts[0];
                double y2 = star[i - 1].orts[1];
                double z2 = star[i - 1].orts[2];
                double x3 = star[i].orts[0];
                double y3 = star[i].orts[1];
                double z3 = star[i].orts[2];
                double dx1 = x3 - x1;
                double dy1 = y3 - y1;
                double dz1 = z3 - z1;
                double dx2 = x3 - x2;
                double dy2 = y3 - y2;
                double dz2 = z3 - z2;
                double cpx = dy1*dz2 - dz1*dy2;
                double cpy = dz1*dx2 - dx1*dz2;
                double cpz = dx1*dy2 - dy1*dx2;
                double area = sqrt(cpx*cpx + cpy*cpy + cpz*cpz)/2;
                sX += (x1 + x2 + x3)/3*area;
                sY += (y1 + y2 + y3)/3*area;
                sZ += (z1 + z2 + z3)/3*area;
                AllArea += area;
        }
        return Vector(sX/AllArea, sY/AllArea, sZ/AllArea);
    }

    friend std::ostream& operator << (std::ostream& os, const Vector& A) {
        os << "{" << A[0] << ", " << A[1] << ", " << A[2] << "}";
//                os << std::endl;
        return os;
    }
};

class Plane {
private:
    double a, b, c, d;
public:
    Plane(){
        a = 0.0;
        b = 0.0;
        c = 0.0;
        d = 0.0;
    }
    Plane(const Plane &P) {
        a = P.a;
        b = P.b;
        c = P.c;
        d = P.d;
    }
    Plane(double _a, double _b, double _c, double _d) {
        a = _a;
        b = _b;
        c = _c;
        d = _d;
    }
    Plane(const Vector &Normal, float _d) {
        a = Normal[0];
        b = Normal[1];
        c = Normal[2];
        d = _d;
    }
    double operator[](int i) const {
        if(i > 4 || i < 0) {
            throw std::runtime_error("Выпал из диапазона");
        }
        if(i == 0) return a;
        if(i == 1) return b;
        if(i == 2) return c;
        if(i == 3) return d;
    }
    Plane Normalize() {
        Plane Result;
        const double Distance = std::sqrt(a*a + b*b + c*c);
        Result.a = a/Distance;
        Result.b = b/Distance;
        Result.c = c/Distance;
        Result.d = d/Distance;
        return Result;
    }
    static Plane FromPointNormal(const Vector &Pt, const Vector &Normal) {
        Plane Result;
        Vector NormalizedNormal = Normal.Normalize();
        Result.a = NormalizedNormal[0];
        Result.b = NormalizedNormal[1];
        Result.c = NormalizedNormal[2];
        Result.d = -(Pt*NormalizedNormal);
        return Result;
    }
    Plane& FromPoints(const Vector &V0, const Vector &V1, const Vector &V2) {
        Vector Cros;
        Vector Normal = Cros(V1 - V0, V2 - V0);
        a = FromPointNormal(V0, Normal).a;
        b = FromPointNormal(V0, Normal).b;
        c = FromPointNormal(V0, Normal).c;
        d = FromPointNormal(V0, Normal).d;
        return *this;
    }
    static bool Intersection(const Plane &P1, const Plane &P2) {
        double Den = P1.a * P2.b - P1.b * P2.a;
        if(Den == 0.0) {
            return false;
        }
        Vector Cros;
        Vector pp1(P1.a, P1.b, P1.c);
        Vector pp2(P2.a, P2.b, P2.c);
        Vector D = Cros(pp1.Normalize(), pp2.Normalize());
        if(D.length() == 0.0) {
            return false;
        }

        return true;
    }
    double SignedDistance(const Vector &V0) const
    {
        return (a * V0[0] + b * V0[1] + c * V0[2] + d);
    }
    double get_z(double x, double y) const
    {
        return c != 0 ? -(a * x + b * y + d)/c : 0.0;
    }
    friend std::ostream& operator << (std::ostream& os, const Plane& A) {
        os << "{" << A.a << ", " << A.b << ", " << A.c << "; " << A.d << "}";
                os << std::endl;
        return os;
    }
};


#endif // SOLVER_HPP
