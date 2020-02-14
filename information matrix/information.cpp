//2011.1.29 To correct the information matrix
//Two purposes:
//1.<Sncos>,<Snsin>,<SnCos>^2,<SnSin>^2
//2.deviation      
//
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include "jacobrot.h"
#include "cvector.h"
#include "rvector.h"
#include "cmatrix.h"
#include "symMatrix.h"
#define Pi acos(-1.0)
#define aa 0.01
#define nn 80000
using namespace std;
int GG(int n, int m){return (2*n+m);}
double sncos(double x,double p1,double p2,int q)
{ 
  if(q==0) return cos(x)/(1.0+aa *exp(-p1*cos(x)/2.0 -p2*sin(x)/2.0  ));
  if(q==1) return sin(x)/(1.0+aa *exp(-p1*cos(x)/2.0 -p2*sin(x)/2.0  ));
  }

double simpson_sncos(double a, double b, int n, double p1, double p2,int q)
{ 
    double h, sum, x4, x2;
    int i;
    h = (b - a) / (double)n;
    sum = sncos(a,p1,p2,q) + 4.0* sncos(a+h,p1,p2,q)+ sncos(b,p1,p2,q);
    for (i=2; i < (n - 1); i+=2)
        { x2 = a + (double)i * h;
          x4 = a + (double)(i+1) * h;
          sum = sum + 2.0 * sncos(x2,p1,p2,q) + 4.0 * sncos(x4,p1,p2,q);
        }
    return (sum * h / 3.0);
    
}

double first(double x,double p1, double p2,int g)
{
 if(g==0) return pow(cos(x),2.0)/(1.0+aa *exp(-p1*cos(x)/2.0 -p2*sin(x)/2.0  ));
 if(g==1) return pow(sin(x),2.0)/(1.0+aa *exp(-p1*cos(x)/2.0 -p2*sin(x)/2.0  ));      }
double second(double x,double p1, double p2,int g)
{
 if(g==0) return pow(cos(x)/(1.0+aa *exp(-p1*cos(x)/2.0 -p2*sin(x)/2.0)),2.0);
 if(g==1) return pow(sin(x)/(1.0+aa *exp(-p1*cos(x)/2.0 -p2*sin(x)/2.0)),2.0);     
 }

double simpson_sncos2(double a, double b, int n, double p1, double p2,int q)
{ 
    double h, x4, x2;
    double sum_first =0.0, sum_second = 0.0;
    int i;
    h = (b - a) / (double)n;
    //first term
    sum_first = first(a,p1,p2,q) + 4.0* first(a+h,p1,p2,q)+ first(b,p1,p2,q);
    for (i=2; i < (n - 1); i+=2)
        { x2 = a + (double)i * h;
          x4 = a + (double)(i+1) * h;
          sum_first = sum_first + 2.0 * first(x2,p1,p2,q) + 4.0 * first(x4,p1,p2,q);
        }
    //second term
    sum_second = second(a,p1,p2,q) + 4.0* second(a+h,p1,p2,q)+ second(b,p1,p2,q);
    for (i=2; i < (n - 1); i+=2)
        { x2 = a + (double)i * h;
          x4 = a + (double)(i+1) * h;
          sum_second = sum_second + 2.0 * second(x2,p1,p2,q) + 4.0 * second(x4,p1,p2,q);
        }
    //return (sum_first * h / 3.0-sum_second * h / 3.0);
    return (sum_first * h / 3.0-sum_second * h / 3.0);
}

main()
{
  double ans=0.0, ans_phi= 0.0;
  int iter = 50;
  int i=0, j=0;
  double upper = 1.0;
  double lower = 0.02;
  double h = (upper - lower) / (double)iter;
  
  double phi = 0.0;
  double h_phi = 2.0* Pi / iter;
  
  double *sncos, *dsncos, *dsncos_phi;
  sncos = new double [2000];
  dsncos = new double [iter];
  dsncos_phi = new double [iter];
  double *sncos2, *dsncos2, *dsncos2_phi;
  sncos2 = new double [2000];
  dsncos2 = new double [iter];
  dsncos2_phi = new double [iter]; 
  
  double *snsin, *dsnsin, *dsnsin_phi;
  snsin = new double [2000];
  dsnsin = new double [iter];
  dsnsin_phi = new double [iter];
  double *snsin2, *dsnsin2, *dsnsin2_phi;
  snsin2 = new double [2000];
  dsnsin2 = new double [iter];
  dsnsin2_phi = new double [iter]; 
 
  //information process:
  //1. assign the matrix
  //2. To calculate the eigenvalue of the matrix
 
  FILE *fp;
  fp = fopen("sigmaphi a001.txt","w");
  for(i=0;i<iter;i++)
  {    
  //1. obtain the expectation value and fluctuation
    for(j=0;j<2;j++)
    {
      sncos[GG(i,j)] = nn*simpson_sncos(0.0,2.0*Pi,1000,(lower+i*h)*cos(phi+j*h_phi),(lower+i*h)*sin(phi+j*h_phi),0)/(2.0*Pi); 
      sncos2[GG(i,j)]= nn*simpson_sncos2(0.0,2.0*Pi,1000,(lower+i*h)*cos(phi+j*h_phi),(lower+i*h)*sin(phi+j*h_phi),0)/(2.0*Pi); 
      snsin[GG(i,j)] = nn*simpson_sncos(0.0,2.0*Pi,1000,(lower+i*h)*cos(phi+j*h_phi),(lower+i*h)*sin(phi+j*h_phi),1)/(2.0*Pi); 
      snsin2[GG(i,j)]= nn*simpson_sncos2(0.0,2.0*Pi,1000,(lower+i*h)*cos(phi+j*h_phi),(lower+i*h)*sin(phi+j*h_phi),1)/(2.0*Pi); 
     }
  }
  
  //2. to differetiate the expectation value and fluctuation
  for(int j=2;j<iter-1;j++) 
   {
     dsncos[j] = (sncos[GG(j+1,0)]-sncos[GG(j,0)])/h;
     dsncos2[j] = (sncos2[GG(j+1,0)]-sncos2[GG(j,0)])/h;
     dsncos_phi[j] = (sncos[GG(j,1)]-sncos[GG(j,0)])/h_phi;
     dsncos2_phi[j] = (sncos2[GG(j,1)]-sncos2[GG(j,0)])/h_phi;
     
     dsnsin[j] = (snsin[GG(j+1,0)]-snsin[GG(j,0)])/h;
     dsnsin2[j] = (snsin2[GG(j+1,0)]-snsin2[GG(j,0)])/h;
     dsnsin_phi[j] = (snsin[GG(j,1)]-snsin[GG(j,0)])/h_phi;
     dsnsin2_phi[j] = (snsin2[GG(j,1)]-snsin2[GG(j,0)])/h_phi;
  
     ans = pow(dsncos[j],2.0)/sncos2[GG(j,0)]+pow(dsncos2[j]/sncos2[GG(j,0)],2.0)/2.0+pow(dsnsin[j],2.0)/snsin2[GG(j,0)]+pow(dsnsin2[j]/snsin2[GG(j,0)],2.0)/2.0;
     ans_phi = pow(dsncos_phi[j],2.0)/sncos2[GG(j,1)]+pow(dsncos2_phi[j]/sncos2[GG(j,1)],2.0)/2.0+pow(dsnsin_phi[j],2.0)/snsin2[GG(j,1)]+pow(dsnsin2_phi[j]/snsin2[GG(j,1)],2.0)/2.0;
     fprintf(fp,"%.15lf %.15lf\n",lower+j*h,1.0/ans_phi);
   }  
   
  fclose(fp);
//  cout << ans<< endl;
  system("PAUSE");
  return 0;
      }
