//===========================================================
// fundtions for real symmtric matrix
// Lanczos Tridiagonal : LanczosTrid( inMatrix, Q-Matrix);
// jacobrot : jacobian iterative rotations
//----------------------------------------------------------
#ifndef _real_symMatrix_functions__
#define _real_symMatrix_functions__
#include <cmath>
using namespace std;
#include "sqrMatrix.h"
#include "rvector.h"
//==========================================================
// Jacobian iterative rotation for a real sysmmetric Matrix
//==========================================================
#define MAX_DIM 20                          /* Maximum allowed dimension */
#define MAX_ITER_NO 2000                    /* set maximum No of iteration */
int jacobrot(const Matrix &mtx0, colVector &val, Matrix &egv, const double crit)
{ double vt1, vt2, ut1, ut2, tmpx[MAX_DIM], tmpy[MAX_DIM];
  double ra, rb, rc, ax;
  double e1, e2, rnn, ramb, rapb, rax;
  unsigned int i, j, nx, ny, ni, ndim;
  Matrix mtx(mtx0);
  ndim = mtx.ncol();
  if (ndim > MAX_DIM) return(-1);            /* Dimension overflow */
// ============================================================
// Initail matrix egv to be a unit matrix for eigen vector
//-------------------------------------------------------------
  for (i=0; i<ndim; i++)
     { for (j=0; j<ndim; j++)
	  { egv.cell(i, j) = (i==j)? 1.0: 0.0; }
     }
//==============================================================
// Starting the loop for itrative rotation
//--------------------------------------------------------------
  for (ni=0; ni<MAX_ITER_NO; ni++)
     {
// Find the non-diagonal element with maximum absolute value.
       rax = 0.0;
       ny=1;
       nx=0;
       for (i=0; i<ndim; i++)
	  { for (j=0; j<i; j++)
	       {  ax = mtx.cell(i, j);
	          rnn = abs(ax);
		     if (rax < rnn)
                { rax = rnn;
				  ny = i;
				  nx = j;
                }
           }
	  }
// Check if the convergent condition met?
       if (rax < crit) break;            /* job completed */
// 2x2 rotation to eliminate the non-diagonal element at (nx, ny)
       else { ra = mtx.cell(nx, nx);
	          rc = mtx.cell(nx,ny);
	          rb = mtx.cell(ny,ny);
	         ramb = 0.5*(ra-rb);
	         rapb = 0.5*(ra+rb);
	         rnn = sqrt( ramb*ramb + rc*rc );
	         e1 = rapb - rnn;  // eigen value of the 2 x 2 matrix
	         e2 = rapb + rnn;  // eigen value of the 2 x 2 matrix
	         vt1 = rc;
             vt2 = e1 - ra;
	         rnn = sqrt(vt1*vt1 + vt2*vt2);
	         vt1 = vt1 / rnn;
	         vt2 = vt2 / rnn;  // eigen vector of 2x2 matrix
	         ut1 = -vt2;
	         ut2 = vt1;  // eigen vector of 2x2 matrix
// ==================================================================
// Two vectors to update the next matrix rows and columnes of nx, ny
// ------------------------------------------------------------------
              for (i=0; i<ndim; i++)
                 { tmpx[i] = vt1*mtx.get(i,nx) + vt2* mtx.get(i,ny); }
              for (i=0; i<ndim; i++)
                 { tmpy[i] = ut1*mtx.get(i,nx) + ut2*mtx.get(i,ny); }
              tmpx[nx] = e1;
              tmpx[ny] = 0.0;
              tmpy[nx] = 0.0;
              tmpy[ny] = e2;

	      for (i=0; i<ndim; i++) { mtx.cell(nx,i) = tmpx[i]; }
	      for (i=0; i<ndim; i++) { mtx.cell(i,nx) = tmpx[i]; }
	      for (i=0; i<ndim; i++) { mtx.cell(ny,i) = tmpy[i]; }
	      for (i=0; i<ndim; i++) { mtx.cell(i,ny) = tmpy[i]; }
/*=======================================================================
 *   Update the eigen vector matrix.
 *----------------------------------------------------------------------*/
	      for (i=0; i<ndim; i++)
		      { tmpx[i] = vt1*egv.get(i,nx) + vt2*egv.get(i,ny); }
          for (i=0; i<ndim; i++)
		      { tmpy[i] = ut1*egv.get(i,nx) + ut2*egv.get(i,ny); }
	      for (i=0; i<ndim; i++) { egv.cell(i,nx) = tmpx[i]; }
          for (i=0; i<ndim; i++) { egv.cell(i,ny) = tmpy[i]; }
	    } // end of else
//      cout << "iteration = " << ni <<"\n";
//      cout << egv;
//      system("Pause");
      } // end of for(ni)
  if (ni >= MAX_ITER_NO) return (-2);  // reached maximum iteration limit
  val = diagonal(mtx);
  return(ni);
}
//=================================================================
// Tridiagonal a real symmetric Matrix using Lanczos vector space
//----------------------------------------------------------------
void LanczosTrid(const Matrix &amtx, Matrix &qmtx)
{   int i, ndim;
 //   double crt = 1.0e-5;
    double alpha, beta;
//
    Matrix mtx(amtx);
    ndim = mtx.ncol();
    colVector r1(ndim), r2(ndim), r3(ndim), q(ndim);
    q = col(mtx,0);
    beta = q.abs();
    r1 = colVector(ndim, 0.0);
    r2 = q / beta;
    putcol(qmtx, r2, 0);
    for (i=1; i<ndim; i++)
        {  q = mtx * r2;
           alpha = transport(r2) * q;
           q -= ((alpha*r2) + (beta*r1));
           beta = q.abs();
           r3 = q / beta;
           putcol(qmtx, r3, i);
           r1 = r2;
           r2 = r3;
        }
   return;
}
/**********************************************************
 * Using Household transformation  to tridiagonal a real
 * symmetric matrix.
 * "Numerical Analysis, 6th edition" p.569, by R.L. Burden
 * & J.D. Faires, (Brook/Cole, 1997).
 **********************************************************/

Matrix HouseholdTrid(const Matrix &amtx)
 { int i, j, ndim, nk;
   double tt, qsum, akk1, alpha, rsq;
   double *vctv, *vctu, *vctz;

   Matrix symx(amtx);
   ndim = symx.ncol();
   vctv = new double[ndim];
   vctu = new double[ndim];
   vctz = new double[ndim];

   for (nk=0; nk < (ndim-2); nk++)
     { qsum = (double)0.0;
       for (j=nk+1; j<ndim; j++)
         { tt = symx.cell(nk, j);
           qsum += tt*tt;
         }
       akk1 = symx.cell(nk, nk+1);
       alpha = (akk1==(double)0)? -(sqrt(qsum)) : -(sqrt(qsum)*akk1/fabs(akk1));
       rsq = alpha * ( alpha - akk1);

       vctv[nk] = 0.0;
       vctv[nk+1] = akk1 - alpha;

       for (j=nk+2; j<ndim; j++) vctv[j] =  symx.cell(nk, j);

       for (i=nk; i<ndim; i++)
         { qsum = 0.0;
           for (j=nk+1; j<ndim; j++) qsum += (symx.cell(i,j)*vctv[j]);
           vctu[i] = qsum / rsq;
         }

       qsum = (double)0.0;
       for (i=nk+1; i<ndim; i++) qsum += (vctu[i]*vctv[i]);

       for (i=nk; i<ndim; i++) vctz[i] = vctu[i] - qsum*vctv[i]/(2.0*rsq);

       for (i=nk+1; i<ndim; i++)
         { symx.cell(i,i) = symx.cell(i,i)- (double)2 * vctv[i] * vctz[i];
           for (j=i+1; j<ndim; j++)
             { tt = symx.cell(i,j) - vctv[i]*vctz[j]-vctv[j]*vctz[i];
               symx.cell(i,j) = tt;
               symx.cell(j,i) = tt;
             }
          }
       for (j=nk+2; j<ndim; j++)
         { symx.cell(nk,j) = 0.0;
           symx.cell(j,nk) = 0.0;
         }

       tt = akk1 - vctv[nk+1] * vctz[nk];
       symx.cell(nk, nk+1) =  tt;
       symx.cell(nk+1, nk) =  tt;
    }
//   for (i=0; i<ndim; i++) dpt[i] = symx.get(i,i);
//    for (i=0; i<(ndim-1); i++) opt[i] = symx.get(i, i+1);
//    symx.~Matrix();
    delete [] vctv;
    delete [] vctu;
    delete [] vctz;
    return symx;
}
#endif
