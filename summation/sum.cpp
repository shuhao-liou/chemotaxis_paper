//Summation Method
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#define Pi acos(-1.0000)
#define a 0.001
#define phi Pi/6.0
#define n 80000
using namespace std;
double sn(double p,int nn)
{
 double steep = -p *cos(2.0*Pi*nn/n - phi) / 2.00;
 return (a*exp(steep)/(1.0 +a *exp(steep)));
// return p;
       }
double desn(double p,int nn)
{
 double steep = -p *cos(2.0*Pi*nn/n - phi) / 2.00;
 return (a*exp(steep)/pow((1.0 +a *exp(steep)),2.0));
       }
main()
{
 double sum_sn = 0.0, sum_desn = 0.0;
 double p = 1.0; 
 int i,j;
 
 int iter = 100;
 double upper = 10.0;
 double lower = 0.0;
 double h = (upper - lower) / (double)iter;
 
 FILE *fp;
 fp = fopen("summation.txt","w");
 for(j=0;j<iter;j++)
 {
  sum_sn =0.0;
  sum_desn = 0.0;
  for(i=0;i<n;i++) 
  {
   sum_sn+= sn(lower+j*h,i);
 //  cout << sum_sn << endl;
   sum_desn+= desn(lower+j*h,i);
   }
   fprintf(fp,"%lf %.15lf %.15lf\n",lower+j*h,n-2.0*sum_sn,4.0*sum_desn);
 }
 fclose(fp);
 system("PAUSE");
 return 0;
}
