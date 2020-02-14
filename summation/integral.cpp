//2010.9.22 
//Using Simpson Method to integrate
//
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#define Pi acos(-1.0000)
#define phi Pi/3.0
#define aa 0.001
#define nn 80000
using namespace std;
double f(double x,double p)
{ 
  double steep = -p *cos(x - phi) / 2.00;
  return ( aa*exp(steep) /( 1.0+aa*exp(steep)) );
  }
//       return x*x;}
double simpson(double a, double b, int n, double p)
{
 
    double h, sum, x4, x2;
    int i;
    h = (b - a) / (double)n;
    sum = f(a,p) + 4.0* f(a+h,p)+ f(b,p);
    for (i=2; i < (n - 1); i+=2)
        { x2 = a + (double)i * h;
          x4 = a + (double)(i+1) * h;
          sum = sum + 2.0 * f(x2,p) + 4.0 * f(x4,p);
        }
    return (sum * h / 3.0);
   
}
double f2(double x,double p)
{ 
 double steep = -p *cos(x - phi) / 2.00;
 return ( aa*exp(steep) / pow(1.0+ aa*exp(steep),2.0) );
}
//       return 1/(1+x*x);}
double simpson2(double a, double b, int n, double p)
{
 
    double h, sum, x4, x2;
    int i;
    h = (b - a) / (double)n;
    sum = f2(a,p) + 4.0* f2(a+h,p)+ f2(b,p);
    for (i=2; i < (n - 1); i+=2)
        { x2 = a + (double)i * h;
          x4 = a + (double)(i+1) * h;
          sum = sum + 2.0 * f2(x2,p) + 4.0 * f2(x4,p);
        }
    return (sum * h / 3.0);
}  
main()
{
  double ans=0.0, ans2= 0.0;
  int iter = 100;
  double upper = 10000.0;
  double lower = 0.0;
  double h = (upper - lower) / (double)iter;
  double concen;
   
  FILE *fp;
  fp = fopen("sensing_appr.txt","w");
  for(int i=0;i<iter;i++)
  {
   ans = simpson(0.0,2.0*Pi,1000,lower+i*h);
   ans2 = simpson2(0.0,2.0*Pi,1000,lower+i*h); 
   fprintf(fp," %lf %.15lf %.15lf\n",(lower+i*h),nn - nn * ans / Pi,2.0* nn *ans2/ Pi);   
  }
  fclose(fp);
  cout << ans<< endl;
  system("PAUSE");
  return 0;
      }
