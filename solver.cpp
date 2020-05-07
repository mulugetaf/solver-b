
#include <complex>
#include <cmath>
#include <stdexcept>
#include "solver.hpp"

using namespace std;
using namespace solver;
/*RealVariable constructor*/
RealVariable::RealVariable()
{
    a = 0;
    b = 1;
    c = 0;
    n = 1;
}
/*ComplexVariable constructor*/
ComplexVariable::ComplexVariable()
{
    a = 0;
    b = 1;
    c = 0;
    realPart = ans.real();
    imaginaryPart = ans.imag();
}
/*LinearEquation
          return:x value of linear equation
          by using -(c/b) */
double solver::LinearEquation(solver::RealVariable &equation)
{
    double ans = -(equation.c / equation.b);
    if (equation.b == 0)
        throw std::invalid_argument("no solution");
    return ans;
    ;
};

/* 
        overload operator  ==          
        from right and left side
    */
RealVariable &solver::operator==(const solver::RealVariable &r, double x)
{
    RealVariable *ans = new RealVariable();
    ans->a = r.a;
    ans->b = r.b;
    ans->c = r.c;
    ans->n = r.n;

    if (x > 0)
    {
        ans->c -= x;
    }
    else
    {
        ans->c += -x;
    }
    if (ans->n == 1)
        ans->a = 0;
    return *ans;
};
RealVariable &solver::operator==(double x, const solver::RealVariable &r)
{
    return (r == x);
};
RealVariable &solver::operator==(const solver::RealVariable &r, int x)
{
    return (r == (double)x);
};
RealVariable &solver::operator==(int x, const solver::RealVariable &r)
{
    return (r == (double)x);
};
RealVariable &solver::operator==(const solver::RealVariable &r1, const solver::RealVariable &r2)
{
    RealVariable *ans = new RealVariable();
    ans->a = r1.a;
    ans->b = r1.b;
    ans->c = r1.c;
    ans->n = r1.n;
    if (ans->n == 1)
        ans->a = 0;
    if (r2.a > 0)
    {
        if (r2.n > 1)
            ans->a -= r2.a;
    }
    else
    {
        ans->a += -r2.a;
    }
    if (r2.b > 0)
    {
        ans->b -= r2.b;
    }
    else
    {
        ans->b += -r2.b;
    }

    if (r2.c > 0)
    {
        ans->c -= r2.c;
    }
    else
    {
        ans->c += -r2.c;
    }
    return *ans;
};
/* 
        overload operator  <=          
        from right and left side
        throw "unlegal operator" message exception
    */
RealVariable &solver::operator<=(const solver::RealVariable &d, double x)
{
    throw std::invalid_argument("unlegal operator");
};
RealVariable &solver::operator<=(double x, const solver::RealVariable &d)
{
    throw std::invalid_argument("unlegal operator");
};
RealVariable &solver::operator<=(const solver::RealVariable &d, int x)
{
    throw std::invalid_argument("unlegal operator");
};
RealVariable &solver::operator<=(int x, const solver::RealVariable &d)
{
    throw std::invalid_argument("unlegal operator");
};
RealVariable &solver::operator<=(const solver::RealVariable &r, const solver::RealVariable &d)
{
    throw std::invalid_argument("unlegal operator");
};
/* 
        overload operator  >=          
        from right and left side
        throw "unlegal operator" message exception
    */
RealVariable &solver::operator>=(const solver::RealVariable &d, double x)
{
    throw std::invalid_argument("unlegal operator");
};
RealVariable &solver::operator>=(double x, const solver::RealVariable &d)
{
    throw std::invalid_argument("unlegal operator");
};
RealVariable &solver::operator>=(const solver::RealVariable &d, int x)
{
    throw std::invalid_argument("unlegal operator");
};
RealVariable &solver::operator>=(int x, const solver::RealVariable &d)
{
    throw std::invalid_argument("unlegal operator");
};
RealVariable &solver::operator>=(const solver::RealVariable &r, const solver::RealVariable &d)
{
    throw std::invalid_argument("unlegal operator");
};
/* 
        overload operator RealVariablr + int         
        from right side
    */
RealVariable &solver::operator+(const solver::RealVariable &r, int x)
{
    RealVariable *ans = new RealVariable();
    ans->a = r.a;
    ans->b = r.b;
    ans->c = r.c + x;
    ans->n = r.n;
    return *ans;
};
RealVariable &solver::operator-(const solver::RealVariable &r, int x)
{
    RealVariable *ans = new RealVariable();
    ans->a = r.a;
    ans->b = r.b;
    ans->c = r.c - x;
    ans->n = r.n;
    return *ans;
};
RealVariable &solver::operator*(const solver::RealVariable &r, int x)
{
    RealVariable *ans = new RealVariable();

    ans->a = r.a * x;
    ans->b = r.b * x;
    ans->c = r.c * x;
    ans->n = r.n;
    return *ans;
};
RealVariable &solver::operator/(const solver::RealVariable &r, int x)
{
    RealVariable *ans = new RealVariable();
    if (x == 0)
        throw std::invalid_argument("devide by zero");
    ans->a = r.a / x;
    ans->b = r.b / x;
    ans->c = r.c / x;
    ans->n = r.n;
    return *ans;
};
RealVariable &solver::operator^(const solver::RealVariable &r, int x)
{
    RealVariable *ans = new RealVariable();
    /* x^n |   0 < n < 3    */
    if ((x < 1) || (x > 2))
    {
        throw std::invalid_argument("unlegal equation");
    }
    ans->a = 1;
    ans->b = 0;
    ans->c = r.c;
    ans->n = 2;
    return *ans;
};
/* 
        overload operator int + RealVariablr         
        from right side
    */
RealVariable &solver::operator+(int x, const solver::RealVariable &r)

{
    return (r + x);
};
RealVariable &solver::operator-(int x, const solver::RealVariable &r)

{
    return (r - x);
};
RealVariable &solver::operator*(int x, const solver::RealVariable &r)
{
    return (r * x);
};
RealVariable &solver::operator/(int x, const solver::RealVariable &r)

{
    return (r / x);
};
RealVariable &solver::operator^(int x, const solver::RealVariable &r)

{
    return (r ^ x);
};

/* 
        overload operator RealVariable + RealVariablr         
        from right side
    */
RealVariable &solver::operator+(const solver::RealVariable &r1, const solver::RealVariable &r2)
{
    RealVariable *ans = new RealVariable();
    ans->a = r1.a;
    ans->b = r1.b;
    ans->c = r1.c;
    ans->n = r1.n;
    if (ans->n == r2.n)
    {
        ans->a += r2.a;
    }
    ans->b += r2.b;
    ans->c += r2.c;
    return *ans;
};
RealVariable &solver::operator-(const solver::RealVariable &r1, const solver::RealVariable &r2)
{
    RealVariable *ans = new RealVariable();
    ans->a = r1.a;
    ans->b = r1.b;
    ans->c = r1.c;
    ans->n = r1.n;
    if (ans->n == r2.n)
    {
        ans->a -= r2.a;
    }
    ans->b -= r2.b;
    ans->c -= r2.c;
    return *ans;
};
RealVariable &solver::operator*(const solver::RealVariable &r1, const solver::RealVariable &r2)
{
    RealVariable *ans = new RealVariable();
    ans->a = r1.a;
    ans->b = r1.b;
    ans->c = r1.c;
    ans->n = r1.n;

    ans->n = +r2.n;

    ans->a *= r2.a;
    ans->a *= r2.b;
    ans->a *= r2.c;

    return *ans;
};
solver::RealVariable &solver::operator/(const solver::RealVariable &r1, const solver::RealVariable &r2)
{
    RealVariable *ans = new RealVariable();
    ans->a = r1.a;
    ans->b = r1.b;
    ans->c = r1.c;
    ans->n = r1.n;

    ans->n -= r2.n;
    ans->a /= r2.a;
    ans->b /= r2.b;
    ans->c /= r2.c;

    return *ans;
};

/* 
    ComplexVariable
        overload operator  ==          
        from right and left side
    */
ComplexVariable &solver::operator==(const solver::ComplexVariable &c, double x)
{
    ComplexVariable *ans = new ComplexVariable();
    ans->a = c.a;
    ans->b = c.b;
    ans->c = c.c;
    ans->n = c.n;

    if (x > 0)
    {
        ans->c -= x;
    }
    else
    {
        ans->c += -x;
    }
    if (ans->n == 1)
        ans->a = 0;
    return *ans;
};
ComplexVariable &solver::operator==(double x, const solver::ComplexVariable &c)
{
    return (c == x);
};
ComplexVariable &solver::operator==(const solver::ComplexVariable &c, int x)
{
    return (c == (double)x);
};
ComplexVariable &solver::operator==(int x, const solver::ComplexVariable &c)
{
    return (c == (double)x);
};
ComplexVariable &solver::operator==(const solver::ComplexVariable &c1, const solver::ComplexVariable &c2)
{
    ComplexVariable *ans = new ComplexVariable();
    ans->a = c1.a;
    ans->b = c1.b;
    ans->c = c1.c;
    ans->realPart = c1.realPart;
    ans->imaginaryPart = c1.imaginaryPart;
    ans->n = c1.n;
    if (ans->n == 1)
        ans->a = 0;
    if (c2.a > 0)
    {
        if (c2.n > 1)
            ans->a -= c2.a;
    }
    else
    {
        ans->a += -c2.a;
    }
    if (c2.b > 0)
    {
        ans->b -= c2.b;
    }
    else
    {
        ans->b += -c2.b;
    }
    if (c2.c > 0)
    {
        ans->c -= c2.c;
    }
    else
    {
        ans->c += -c2.c;
    }
    if (ans->realPart > 0)
    {
        ans->realPart -= c2.realPart;
    }
    else
    {
        ans->realPart += -c2.realPart;
    }
    if (ans->imaginaryPart > 0)
    {
        ans->imaginaryPart -= c2.imaginaryPart;
    }
    else
    {
        ans->imaginaryPart += -c2.imaginaryPart;
    }
    return *ans;
};
/* 
        overload operator  <=          
        from right and left side
        throw "unlegal operator" message exception
    */

ComplexVariable &solver::operator<=(double x, const solver::ComplexVariable &d)
{
    throw std::invalid_argument("unlegal operator");
};
ComplexVariable &solver::operator<=(const solver::ComplexVariable &d, int x)
{
    throw std::invalid_argument("unlegal operator");
};
ComplexVariable &solver::operator<=(int x, const solver::ComplexVariable &d)
{
    throw std::invalid_argument("unlegal operator");
};
ComplexVariable &solver::operator<=(solver::ComplexVariable d, solver::ComplexVariable s)
{
    throw std::invalid_argument("unlegal operator");
};
/* 
        overload operator  >=          
        from right and left side
        throw "unlegal operator" message exception
    */

ComplexVariable &solver::operator>=(double x, const solver::ComplexVariable &d)
{
    throw std::invalid_argument("unlegal operator");
};
ComplexVariable &solver::operator>=(const solver::ComplexVariable &d, int x)
{
    throw std::invalid_argument("unlegal operator");
};
ComplexVariable &solver::operator>=(int x, const solver::ComplexVariable &d)
{
    throw std::invalid_argument("unlegal operator");
};
ComplexVariable &solver::operator>=(solver::ComplexVariable d, solver::ComplexVariable s)
{
    throw std::invalid_argument("unlegal operator");
};

/* 
        overload operator ComplexVariable + int         
        from right side
    */
ComplexVariable &solver::operator+(const solver::ComplexVariable &r, int x)
{
    ComplexVariable *ans = new ComplexVariable();
    ans->a = r.a;
    ans->b = r.b;
    ans->c = r.c + x;
    ans->n = r.n;
    return *ans;
};
ComplexVariable &solver::operator-(const solver::ComplexVariable &r, int x)
{
    ComplexVariable *ans = new ComplexVariable();
    ans->a = r.a;
    ans->b = r.b;
    ans->c = r.c - x;
    ans->n = r.n;
    return *ans;
};
ComplexVariable &solver::operator*(const solver::ComplexVariable &r, int x)
{
    ComplexVariable *ans = new ComplexVariable();
    ans->a = r.a * x;
    ans->b = r.b * x;
    ans->c = r.c * x;
    ans->n = r.n;
    return *ans;
};
ComplexVariable &solver::operator/(const solver::ComplexVariable &r, int x)
{
    ComplexVariable *ans = new ComplexVariable();
    if (x == 0)
        throw std::invalid_argument("devide by zero");
    ans->a = r.a / x;
    ans->b = r.b / x;
    ans->c = r.c / x;
    ans->n = r.n;
    return *ans;
};
ComplexVariable &solver::operator^(const solver::ComplexVariable &r, int x)
{
    ComplexVariable *ans = new ComplexVariable();
    /* x^n |   0 < n < 3    */
    if ((x < 1) || (x > 2))
    {
        throw std::invalid_argument("unlegal equation");
    }
    ans->a = 1;
    ans->b = 0;
    ans->c = r.c;
    ans->n = 2;
    return *ans;
};
/* 
        overload operator int + ComplexVariable         
        from left side
    */
ComplexVariable &solver::operator+(int x, const solver::ComplexVariable &r)
{
    return (r + x);
};
ComplexVariable &solver::operator-(int x, const solver::ComplexVariable &r)
{
    return (r - x);
};
ComplexVariable &solver::operator*(int x, const solver::ComplexVariable &r)
{
    return (r * x);
};
ComplexVariable &solver::operator/(int x, const solver::ComplexVariable &r)
{
    return (r / x);
};
ComplexVariable &solver::operator^(int x, const solver::ComplexVariable &r)
{
    return (r ^ x);
};
/* 
        overload operator ComplexVariable + ComplexVariable         
        from right side
    */
ComplexVariable &solver::operator+(const solver::ComplexVariable &c1, const solver::ComplexVariable &c2)
{
    ComplexVariable *ans = new ComplexVariable();
    ans->a = c1.a;
    ans->b = c1.b;
    ans->c = c1.c;
    ans->n = c1.n;
    if (ans->n == c2.n)
    {
        ans->a += c2.a;
    }
    ans->b += c2.b;
    ans->c += c2.c;
    return *ans;
};
ComplexVariable &solver::operator-(const solver::ComplexVariable &c1, const solver::ComplexVariable &c2)
{
    ComplexVariable *ans = new ComplexVariable();
    ans->a = c1.a;
    ans->b = c1.b;
    ans->c = c1.c;
    ans->n = c1.n;
    if (ans->n == c2.n)
    {
        ans->a -= c2.a;
    }
    ans->b -= c2.b;
    ans->c -= c2.c;
    return *ans;
};
ComplexVariable &solver::operator*(const solver::ComplexVariable &c1, const solver::ComplexVariable &c2)
{
    ComplexVariable *ans = new ComplexVariable();
    ans->a = c1.a;
    ans->b = c1.b;
    ans->c = c1.c;
    ans->n = c1.n;

    ans->n = +c2.n;

    ans->a *= c2.a;
    ans->a *= c2.b;
    ans->a *= c2.c;

    return *ans;
};
ComplexVariable &solver::operator/(const solver::ComplexVariable &c1, const solver::ComplexVariable &c2)
{
    ComplexVariable *ans = new ComplexVariable();
    ans->a = c1.a;
    ans->b = c1.b;
    ans->c = c1.c;
    ans->n = c1.n;

    ans->n -= c2.n;

    ans->a *= c2.a;
    ans->a *= c2.b;
    ans->a *= c2.c;

    return *ans;
};
ComplexVariable &solver::operator^(const solver::ComplexVariable &r, const solver::ComplexVariable &t)
{
    ComplexVariable *ans = new ComplexVariable();
    return *ans;
};

/* 
        overload operator ComplexVariable + complex<double>         
        from right side
    */
ComplexVariable &solver::operator+(const solver::ComplexVariable &r, std::complex<double> t)
{

    ComplexVariable *ans = new ComplexVariable();
    ans->a = r.a;
    ans->b = r.b;
    ans->c = r.c;
    ans->n = r.n;
    ans->realPart += t.real();
    ans->imaginaryPart += t.imag();
    return *ans;
};
ComplexVariable &solver::operator-(const solver::ComplexVariable &r, std::complex<double> t)
{
    ComplexVariable *ans = new ComplexVariable();
    ans->a = r.a;
    ans->b = r.b;
    ans->c = r.c;
    ans->n = r.n;
    ans->realPart -= t.real();
    ans->imaginaryPart -= t.imag();
    return *ans;
};

ComplexVariable &solver::operator*(const solver::ComplexVariable &r, std::complex<double> t)
{
    ComplexVariable *ans = new ComplexVariable();
    ans->a = r.a;
    ans->b = r.b;
    ans->c = r.c;
    ans->n = r.n;
    ans->realPart = t.real();
    ans->imaginaryPart += t.imag();
    return *ans;
};

ComplexVariable &solver::operator/(const solver::ComplexVariable &r, std::complex<double> t)
{
    //need to implement
    ComplexVariable *ans = new ComplexVariable();
    return *ans;
};

ComplexVariable &solver::operator^(const solver::ComplexVariable &r, std::complex<double> t)
{
    //need to implement
    ComplexVariable *ans = new ComplexVariable();
    return *ans;
};

/* 
        overload operator int + complex<double>         
        from right and left side
    */
ComplexVariable &solver::operator+(int x, std::complex<double> &t)
{
    //need to implement
    ComplexVariable *temp = new ComplexVariable();
    return *temp;
};
ComplexVariable &solver::operator+(std::complex<double> &t, int x)
{
    //need to implement
    ComplexVariable *temp = new ComplexVariable();
    return *temp;
};
/* 
        overload operator complex<double> + ComplexVariable        
        from right  side
    */
ComplexVariable &solver::operator+(std::complex<double> &t, const solver::ComplexVariable &r)
{
    //need to implement
    ComplexVariable *ans = new ComplexVariable();
    return *ans;
};
/* 
        overload operator RealVariable+ ComplexVariable        
        from right and left side
    */
ComplexVariable &solver::operator+(solver::RealVariable r, solver::ComplexVariable c)
{
    //need to implement
    ComplexVariable *temp = new ComplexVariable();
    return *temp;
};
ComplexVariable &solver::operator+(solver::ComplexVariable c, solver::RealVariable r)
{
    //need to implement
    ComplexVariable *temp = new ComplexVariable();
    return *temp;
};
/* 
        overload operator double + complex<double>        
        from right and left side

    */

ComplexVariable &solver::operator+(double x, std::complex<double> &t)
{
    //need to implement
    ComplexVariable *temp = new ComplexVariable();
    return *temp;
};
ComplexVariable &solver::operator+(std::complex<double> &t, double x)
{
    //need to implement
    ComplexVariable *temp = new ComplexVariable();
    return *temp;
};
ComplexVariable &solver::operator*(std::complex<double> &t, double x)
{
    //need to implement
    ComplexVariable *temp = new ComplexVariable();
    return *temp;
}
/* 
    solve function
    return : double root1 or root2
    source : https://www.programiz.com/cpp-programming/examples/quadratic-roots
    */
double solver::solve(RealVariable &r)
{
    double a = r.a;
    double b = r.b;
    double c = r.c;
    double discriminant = b * b - 4 * a * c;
    if (r.n > 1)
    {
        if (discriminant > 0)
        {
            r.root1 = (-b + sqrt(discriminant)) / (2 * a);
            r.root2 = (-b - sqrt(discriminant)) / (2 * a);
        }

        else if (discriminant == 0)
        {
            r.root2 = (-b + sqrt(discriminant)) / (2 * a);
        }
        // roo1 and roo2 are complex
        else
        {
            throw std::invalid_argument("no solution");
        }
        return r.root1;
    }
    // if given equation is linear
    else if (r.n == 1)
    {
        return LinearEquation(r);
    }
    /*if given equation is not square or linear 
    throw exception
    */
    else
    {
        throw std::invalid_argument("unlegal");
    }
};
/* 
    solve function
    return : complex<double>  (realPart , imaginaryPart)
    source : https://www.programiz.com/cpp-programming/examples/quadratic-roots
           : https://fahad-cprogramming.blogspot.com/2017/07/complex-numbers-class-cpp-example.html 
    */
std::complex<double> solver::solve(ComplexVariable &c)
{
    double a1 = c.a;
    double b1 = c.b;
    double c1 = c.c;
    double discriminant = b1 * b1 - 4 * a1 * c1;
    /*
    if given equation is :  a + yi = b + xi
    then it return (c.realPart, c.imaginaryPart);
    source : https://www.mathsisfun.com/algebra/complex-number-multiply.html
    */
    if (c.imaginaryPart != 0 || c.realPart != 0)
    {
        return std::complex<double>(c.realPart, c.imaginaryPart);
    }
    if (c.n == 2)
    {
        if (discriminant > 0)
        {
            return std::complex<double>((-b1 + sqrt(discriminant)) / (2 * a1), 0.0);
        }
        else if (discriminant == 0)
        {
            return std::complex<double>((-b1 + sqrt(discriminant)) / (2 * a1), 0.0);
        }

        else
        {
            return std::complex<double>(-b1 / (2 * a1), sqrt(-discriminant) / (2 * a1));
        }
    }
    //linear equation
    else
    {
        double a = -(c1 / b1);
        if (b1 == 0)
            throw std::invalid_argument("no solution");
        return std::complex<double>(a, 0);
    }
    /*if given equation power is gretter then 2 or less then 1
    will throw "unlegal eqution" exception
    */
    throw std::invalid_argument("unlegal equation");
};
std::ostream &solver::operator<<(std::ostream &o, const solver::RealVariable &s)
{
    return o;
};
std::ostream &solver::operator<<(std::ostream &o, const solver::ComplexVariable &s)
{
    return o;
};