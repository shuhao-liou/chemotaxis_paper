#ifndef _ytlu_class_simpson_
#define _ytlu_class_simpson_
#include <iostream>
using namespace std;
class simpson
{
    double a, b;
    int n;
    
    function intg;
    
    public:
        
        simpson(){;}
        simpson(double, double, int, function);
        simpson(double, double, function);
        simpson(function);
        
        double result();
        double result(int);
        double result(double, double);
        double result(double, double, int);
        
 inline double lower() const    { return this->a;}
 inline void   lower(double ra) { this->a = ra; }
 inline double upper() const    { return this->b ;}
 inline void   upper(double rb) { this->b = rb; }
 inline void   limit(double ra, double rb)
                                { this->a = ra; this->b=rb;}
 inline int    mesh() const     {return this->n; }
 inline void   mesh(int nn)     { this->n = nn; }
        void   show();
};

simpson::simpson(double ra, double rb, int nn, function f)
{ 
  this->a = ra;
  this->b = rb;         
  this->n = ((nn%2)==0)? nn : nn + 1;
  this->intg  = f;
}
simpson::simpson(double ra, double rb, function f)
{ 
  this->a = ra;
  this->b = rb;
  this->n = 1000;   //set default value for mesh number.
  this->intg  = f;
}
simpson::simpson(function f)
{  
  this->a = 0.0;    //set the default lower limit
  this->b = 1.0;    //set the default upper limit.
  this->n = 1000;   //set default value for mesh number.
  this->intg  = f;
}
double simpson::result()
{
    double h, sum, x4, x2;
    int i;
    
    h = (this->b - this->a) / (double)this->n;
    sum = intg.f(a) + 4.0* intg.f(a+h)+ intg.f(b);
    for (i=2; i < (this->n - 1); i+=2)
        { x2 = a + (double)i * h;
          x4 = a + (double)(i+1) * h;
          sum = sum + 2.0 * intg.f(x2) + 4.0 * intg.f(x4);
        }
    return (sum * h / 3.0);
}    
double simpson::result(int nn)
{
     this->n = ((nn%2)==0)? nn : nn + 1;
     return this->result();
} 
double simpson::result(double ra, double rb)
{
     this->a = ra;
     this->b = rb;
     return this->result();
}
double simpson::result(double ra, double rb, int nn)
{
     this->n = ((nn%2)==0)? nn : nn + 1;
     this->a = ra;
     this->b = rb;
     return this->result();
}
void simpson::show()
{  cout << endl;
   cout << "    Upper limit = " << this->a << endl;
   cout << "    Lower limit = " << this->b << endl;
   cout << "    Grid number = " << this->n << endl;
   cout << endl;
}   
#endif
