#ifndef _Jacobian_Iterative_Rotation__
#define _Jacobian_Iterative_Rotation__
#include <fstream>
#include "cmatrix.h"
using namespace std;
//==========================================================
// Jacobian iterative rotation for HERMITIAM Matrix
//==========================================================
#define MAX_DIM 20                          /* Maximum allowed dimension */
#define MAX_ITER_NO 2000                    /* set maximum No of iteration */
int jacobrot(C_Matrix & mtx, C_Matrix & egv, double crit)
{ COMPLX vt1, vt2, ut1, ut2, tmpx[MAX_DIM], tmpy[MAX_DIM];
  COMPLX ra, rb, rc, ax;
  double e1, e2, rnn, ramb, rapb, rax;
  unsigned int i, j, nx, ny, ni, ndim;

  ndim = mtx.ncol();
  if (ndim > MAX_DIM) return(-1);            /* Dimension overflow */
// ============================================================
// Initail matrix egv to be a unit matrix for eigen vector
//-------------------------------------------------------------
  for (i=0; i<ndim; i++)
     { for (j=0; j<ndim; j++)
	  { egv.cell(i, j) = (i==j)? ONE: ZERO; }
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
       if (rax < crit) return(ni);            /* job completed */
// 2x2 rotation to eliminate the non-diagonal element at (nx, ny)
       else { ra = mtx.cell(nx, nx);
	          rc = mtx.cell(nx,ny);
	          rb = mtx.cell(ny,ny);
	         ramb = 0.5*(ra-rb).real();
	         rapb = 0.5*(ra+rb).real();
	         rnn = sqrt( ramb*ramb + norm(rc) );
	         e1 = rapb - rnn;  // eigen value of the 2 x 2 matrix
	         e2 = rapb + rnn;  // eigen value of the 2 x 2 matrix
	         vt1 = rc;
             vt2 = e1 - ra;
	         rnn = sqrt(norm(vt1) + norm(vt2));
	         vt1 = vt1 / rnn;
	         vt2 = vt2 / rnn;  // eigen vector of 2x2 matrix
	         ut1 = conj(-vt2);
	         ut2 = conj(vt1);  // eigen vector of 2x2 matrix
// ==================================================================
// Two vectors to update the next matrix rows and columnes of nx, ny
// ------------------------------------------------------------------
              for (i=0; i<ndim; i++)
                 { tmpx[i] = vt1*mtx.get(i,nx) + vt2* mtx.get(i,ny); }
              for (i=0; i<ndim; i++)
                 { tmpy[i] = ut1*mtx.get(i,nx) + ut2*mtx.get(i,ny); }
              tmpx[nx] = COMPLX(e1, 0.0);
              tmpx[ny] = ZERO;
              tmpy[nx] = ZERO;
              tmpy[ny] = COMPLX(e2, 0.0);

	      for (i=0; i<ndim; i++) { mtx.cell(nx,i) = conj(tmpx[i]); }
	      for (i=0; i<ndim; i++) { mtx.cell(i,nx) = tmpx[i]; }
	      for (i=0; i<ndim; i++) { mtx.cell(ny,i) = conj(tmpy[i]); }
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
  return(ni);
}
#endif
