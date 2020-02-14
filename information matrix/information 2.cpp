// 2011.05.28 The two competitively ligand simulation
// Two different ligands inside the concentration pool
// 1. derive the <sncos >, and its variation
// 2. Derive the minimum estimator
// The further model is contain in other document.
// In here, all the calculation is obtained by directly summation.
// 2011.06.24 The information matrix for single ligand
// 2011.06.24 check the results by using the information matrix
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#define Pi acos(-1.0)
#define nn 80000
#define bb1  1.0  //d1
#define phi_1 Pi/5.0
using namespace std;
int GG(int n, int m){return (2*n+m);}
double sncos(int n,double p1, double phi1, int q)
{  
  double phi_n = (2*Pi*n/(nn));
  //if(q==0) return 0.5*cos(phi_n)*(-1.0+bb1*exp( p1/2.0*cos(phi_n-phi1) ) )/(1.0+bb1*exp( p1/2.0*cos(phi_n-phi1) ) );
  if(q==0) return -cos(phi_n)/(1.0+ bb1*exp( p1/2.0*cos(phi_n-phi1)) );
  if(q==1) return 0.5*sin(phi_n)*(-1.0+bb1*exp( p1/2.0*cos(phi_n-phi1) ) )/(1.0+bb1*exp( p1/2.0*cos(phi_n-phi1) ) );
  }

double sum_sncos(double p1, double phi1, int q)
{
 double sum=0.0;
 int n=0;
 for(n=0;n<nn;n++)sum += sncos(n,p1,phi1,q);
 return sum;
       }


double sncos2(int n, double p1, double phi1, int q)
{
  double phi_n = (2*Pi*n/(nn));
  if(q==0) return 0.25*cos(phi_n)*cos(phi_n)*(1.0- pow((-1.0+bb1*exp( p1/2.0*cos(phi_n-phi1) ) )/(1.0+bb1*exp( p1/2.0*cos(phi_n-phi1) ) ),2.0)  );
  if(q==1) return 0.25*sin(phi_n)*sin(phi_n)*(1.0- pow((-1.0+bb1*exp( p1/2.0*cos(phi_n-phi1) ) )/(1.0+bb1*exp( p1/2.0*cos(phi_n-phi1) ) ),2.0)  );
         }

double sum_sncos2(double p1, double phi1, int q)
{
 double sum=0.0;
 int n=0;
 for(n=0;n<nn;n++)sum += sncos2(n,p1,phi1,q);
 return sum;
       }

main()
{
  double ans=0.0, ans_phi= 0.0;
  double matrix11, matrix12, matrix21, matrix22;
  int iter = 50;
  int i=0, j=0;
  double upper = 1.0;
  double lower = 0.02;
  double h = (upper - lower) / (double)iter;
  double h_p= 0.0001;
  double h_phi = 0.001;
  
  double *dsncos_phi, *dsncos;
  dsncos_phi = new double [iter];
  dsncos = new double [iter];
  double *dsncos2_phi, *dsncos2;
  dsncos2_phi = new double [iter]; 
  dsncos2 = new double [iter];
  
  double *dsnsin_phi, *dsnsin;
  dsnsin_phi = new double [iter];
  dsnsin = new double [iter];
  double *dsnsin2_phi, *dsnsin2;
  dsnsin2_phi = new double [iter]; 
  dsnsin2 = new double [iter];
  double sncos2, snsin2;
 
  double p1= lower;
  FILE *fp;
  fp = fopen("d1 1 Pi5.txt","w");
  
  for(j=1;j<iter-1;j++) 
   {                             
     sncos2 = sum_sncos2(p1,phi_1,0);
     snsin2 = sum_sncos2(p1,phi_1,1);
     
     dsncos[j] = ( sum_sncos(p1+h_p,phi_1,0)-sum_sncos(p1-h_p,phi_1,0) )/(2.0*h_p);
     dsncos2[j] = ( sum_sncos2(p1+h_p,phi_1,0)-sum_sncos2(p1-h_p,phi_1,0))/(2.0*h_p);
     dsncos_phi[j] = ( sum_sncos(p1,phi_1+h_phi,0)-sum_sncos(p1,phi_1,0) )  /(h_phi);
     dsncos2_phi[j] = (sum_sncos2(p1,phi_1+h_phi,0)-sum_sncos2(p1,phi_1-h_phi,0))/(2.0*h_phi);
    
     dsnsin[j] = ( sum_sncos(p1+h_p,phi_1,1)-sum_sncos(p1-h_p,phi_1,1) )/(2.0*h_p);
     dsnsin2[j] = ( sum_sncos2(p1+h_p,phi_1,1)-sum_sncos2(p1-h_p,phi_1,1))/(2.0*h_p);
     dsnsin_phi[j] = ( sum_sncos(p1,phi_1+h_phi,1)-sum_sncos(p1,phi_1-h_phi,1) )  /(2.0*h_phi);
     dsnsin2_phi[j] = (sum_sncos2(p1,phi_1+h_phi,1)-sum_sncos2(p1,phi_1-h_phi,1))/(2.0*h_phi);
  
     matrix11 = pow(dsncos[j],2.0)/sncos2+pow(dsncos2[j]/sncos2,2.0)/2.0+pow(dsnsin[j],2.0)/snsin2+pow(dsnsin2[j]/snsin2,2.0)/2.0;
     matrix22 = pow(dsncos_phi[j],2.0)/sncos2+pow(dsncos2_phi[j]/sncos2,2.0)/2.0+pow(dsnsin_phi[j],2.0)/snsin2+pow(dsnsin2_phi[j]/snsin2,2.0)/2.0;
     matrix12 = dsncos[j]*dsncos_phi[j]/sncos2+ dsnsin[j]*dsnsin_phi[j]/snsin2 + dsncos2[j]*dsncos2_phi[j]/2.0/sncos2+ dsnsin2[j]*dsnsin2_phi[j]/2.0/snsin2;
     matrix21 = matrix12;
     
   //  cout << sum_sncos(p1+h_p,phi_1,0) << " " << sum_sncos(p1,phi_1,0) << endl;
//     fprintf(fp,"c%.15lf %.15lf\n",dsncos_phi[j],dsnsin_phi[j]);
     fprintf(fp,"coef %.15lf %.15lf %.15lf\n",matrix11,matrix12,matrix22);
     ans= matrix22/(matrix11*matrix22-matrix12*matrix21);
     ans_phi= matrix11/(matrix11*matrix22-matrix12*matrix21);
     
     fprintf(fp,"%.15lf %.15lf %.15lf\n",p1,ans,(ans_phi));

     p1+=h;
     
   }  
   
  fclose(fp);

  system("PAUSE");
  return 0;
}
