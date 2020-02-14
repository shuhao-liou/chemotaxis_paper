#ifndef _ytlu_class_function_
#define _ytlu_class_function_
#include <math.h>

class function
{
    double omega, phi;
    
    public:
        
        function(){;};
        function(double r1, double r2)
          { omega = r1;  phi = r2; }
        
        double f(double x)
          { return x*omega+phi; }
        
        double df(double x)
          { return omega; }

/*        double f(double x)
          { return exp(omega*x+phi); }
        
        double df(double x)
          { return exp(omega*x+phi)/omega; }
  */     
};    
#endif
