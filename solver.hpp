#include <iostream>
#include <complex>
#include <stdbool.h>
#include <vector>

using namespace std;
namespace solver
{
  /*
    class RealVariable
    root1 ,root2 result of square equation
    n = power of current x value
    ax^2 +bx + c = 0
 */
  class RealVariable
  {
  public:
    double a, b, c, root1, root2, discriminant;
    int n = 1;
    //constructor
    RealVariable();
    /* 
        overload operator  == , <= ,>=
                        + , - , * , / , ^ 
        from right and left side
    */
    friend RealVariable &operator==(const RealVariable &d, double x);
    friend RealVariable &operator==(double x, const RealVariable &d);
    friend RealVariable &operator==(const RealVariable &d, const int x);
    friend RealVariable &operator==(int x, const RealVariable &d);
    friend RealVariable &operator==(const RealVariable &r, const RealVariable &d);
    /*  overload Inequalities <=
             throw exception
            */
    friend RealVariable &operator<=(const RealVariable &d, double x);
    friend RealVariable &operator<=(double x, const RealVariable &d);
    friend RealVariable &operator<=(const RealVariable &d, int x);
    friend RealVariable &operator<=(int x, const RealVariable &d);
    friend RealVariable &operator<=(const RealVariable &r, const RealVariable &d);
    /*  overload Inequalities >=
             throw exception
            */
    friend RealVariable &operator>=(const RealVariable &d, double x);
    friend RealVariable &operator>=(double x, const RealVariable &d);
    friend RealVariable &operator>=(const RealVariable &d, int x);
    friend RealVariable &operator>=(int x, const RealVariable &d);
    friend RealVariable &operator>=(const RealVariable &r, const RealVariable &d);
    /*  overload  
                RealVariable  <+,-,*,/,^> int
              */
    friend RealVariable &operator+(const RealVariable &r, int x);
    friend RealVariable &operator-(const RealVariable &r, int x);
    friend RealVariable &operator*(const RealVariable &r, int x);
    friend RealVariable &operator/(const RealVariable &r, int x);
    friend RealVariable &operator^(const RealVariable &r, int x);
    /*  overload  
                int  <+,-,*,/,^> RealVariable
              */
    friend RealVariable &operator+(int x, const RealVariable &r);
    friend RealVariable &operator-(int x, const RealVariable &r);
    friend RealVariable &operator*(int x, const RealVariable &r);
    friend RealVariable &operator/(int x, const RealVariable &r);
    friend RealVariable &operator^(int x, const RealVariable &r);

    /*  overload  
                RealVariable <+,-,*,/,^> RealVariable
              */
    friend RealVariable &operator+(const RealVariable &r, const RealVariable &t);
    friend RealVariable &operator-(const RealVariable &r, const RealVariable &t);
    friend RealVariable &operator*(const RealVariable &r, const RealVariable &t);
    friend RealVariable &operator/(const RealVariable &r, const RealVariable &t);

    /*  overload  
                RealVariable <+,-,*,/,^>complex<double>
              */
    friend RealVariable &operator+(const RealVariable &r, std::complex<double> t);
    friend RealVariable &operator+(std::complex<double> &t, const RealVariable r);
    friend ::ostream &operator<<(std::ostream &o, RealVariable const &s);
  };
  /*
    class ComplexVariable
    ans result of square equation(complex).
    ans : (real , imag).
    realPart:  real part of the complex number.
    imaginaryPart:  imaginary part of the complex number.
    discriminant = b*b - 4*a*c;
 */
  class ComplexVariable
  {
  public:
    std::complex<double> ans;
    double a, b, c, discriminant,realPart, imaginaryPart;;
    int n = 1;
    //constructor
    ComplexVariable();
    /* 
        overload operator  == , <= ,>=
                        + , - , * , / , ^ 
        from right and left side
    */
    friend ComplexVariable &operator==(const ComplexVariable &d, double x);
    friend ComplexVariable &operator==(double x, const ComplexVariable &d);
    friend ComplexVariable &operator==(const ComplexVariable &d, int x);
    friend ComplexVariable &operator==(int x, const ComplexVariable &d);
    friend ComplexVariable &operator==(std::complex<double> &t, int x);
    friend ComplexVariable &operator==(const ComplexVariable &d, const ComplexVariable &s);
    /*  overload Inequalities <=
             throw exception
            */
    friend ComplexVariable &operator<=(const ComplexVariable &d, double x);
    friend ComplexVariable &operator<=(double x, const ComplexVariable &d);
    friend ComplexVariable &operator<=(const ComplexVariable &d, int x);
    friend ComplexVariable &operator<=(int x, const ComplexVariable &d);
    friend ComplexVariable &operator<=(std::complex<double> &t, int x);
    friend ComplexVariable &operator<=(const ComplexVariable d, const ComplexVariable s);
    /*  overload Inequalities >=
             throw exception
            */
    friend ComplexVariable &operator>=(const ComplexVariable &d, double x);
    friend ComplexVariable &operator>=(double x, const ComplexVariable &d);
    friend ComplexVariable &operator>=(const ComplexVariable &d, int x);
    friend ComplexVariable &operator>=(int x, const ComplexVariable &d);
    friend ComplexVariable &operator>=(std::complex<double> &t, int x);
    friend ComplexVariable &operator>=(const ComplexVariable d, const ComplexVariable s);
    /*  overload  
                ComplexVariable  <+,-,*,/,^> int
              */
    friend ComplexVariable &operator+(const ComplexVariable &r, int x);
    friend ComplexVariable &operator-(const ComplexVariable &r, int x);
    friend ComplexVariable &operator*(const ComplexVariable &r, int x);
    friend ComplexVariable &operator/(const ComplexVariable &r, int x);
    friend ComplexVariable &operator^(const ComplexVariable &r, int x);
    /*  overload  
                int  <+,-,*,/,^> ComplexVariable
              */
    friend ComplexVariable &operator+(int x, const ComplexVariable &r);
    friend ComplexVariable &operator-(int x, const ComplexVariable &r);
    friend ComplexVariable &operator*(int x, const ComplexVariable &r);
    friend ComplexVariable &operator/(int x, const ComplexVariable &r);
    friend ComplexVariable &operator^(int x, const ComplexVariable &r);

    /*  overload  
                ComplexVariable <+,-,*,/,^> ComplexVariable
              */
    friend ComplexVariable &operator+(const ComplexVariable &r, const ComplexVariable &t);
    friend ComplexVariable &operator-(const ComplexVariable &r, const ComplexVariable &t);
    friend ComplexVariable &operator*(const ComplexVariable &r, const ComplexVariable &t);
    friend ComplexVariable &operator/(const ComplexVariable &r, const ComplexVariable &t);
    friend ComplexVariable &operator^(const ComplexVariable &r, const ComplexVariable &t);

    /*  overload  
                ComplexVariable <+,-,*,/,^> complex<double>
              */
    friend ComplexVariable &operator+(const ComplexVariable &r, std::complex<double> t);
    friend ComplexVariable &operator-(const ComplexVariable &r, std::complex<double> t);
    friend ComplexVariable &operator*(const ComplexVariable &r, std::complex<double> t);
    friend ComplexVariable &operator/(const ComplexVariable &r, std::complex<double> t);
    friend ComplexVariable &operator^(const ComplexVariable &r, std::complex<double> t);

    /*  overload  
                complex<double> + ComplexVariable
              */
    friend ComplexVariable &operator+(std::complex<double> &t, const ComplexVariable &r);
    /*  overload  
                RealVariable + ComplexVariable
              */
    friend ComplexVariable &operator+(const RealVariable r, const ComplexVariable c);
    friend ComplexVariable &operator+(const ComplexVariable c, const RealVariable r);
    /*  overload  
                int + complex<double>
              */
    friend ComplexVariable &operator+(int x, std::complex<double> &t);
    friend ComplexVariable &operator+(std::complex<double> &t, int x);
    /*  overload  
                double + complex<double>
              */
    friend ComplexVariable &operator+(double x, std::complex<double> &t);
    friend ComplexVariable &operator+(std::complex<double> &t, double x);
    friend ComplexVariable &operator*(std::complex<double> &t, double x);


    friend ::ostream &operator<<(std::ostream &o, ComplexVariable const &s);

  };

  std::complex<double> solve(ComplexVariable &c);
  double solve(RealVariable &r);
  double LinearEquation(RealVariable &equation);
  /*  overload  operator <<
              */

}; // namespace solver
