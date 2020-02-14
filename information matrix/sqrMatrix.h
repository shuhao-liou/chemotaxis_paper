#ifndef _real_squareMatrix_class__
#define _real_squareMatrix_class__
#include "rmatrix.h"
#include "rvector.h"
#include <cmath>
using namespace std;

double determinant(const Matrix& mtx)
{ double fct, tt, tx;
  int i, j, idx, imax, ndim;
  Matrix aa(mtx);        // Copy THIS matrix to aa for manuipul.

  fct = (double)1.0;        // initial determinat to 1
  ndim = mtx.ncol();
  for (idx=0; idx<ndim; idx++)
    {
//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//   Switch the larger element to diagonal.
//   Uppon switching, determinant chages sign.
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
      tx = fabs(aa.get(idx, idx));
      imax = idx;
      for (i=idx+1; i<ndim; i++)
        { tt = abs(aa.get(i, idx));
          if (tt > tx) { tx = tt;
                         imax = i; }
        }
      if (tx < (double)1.0e-6) return (double)0.0;
      if (imax != idx)
       { fct = -(fct);                // determinant change sign upon
                                      // switching two rows.
         for (j=idx; j<ndim; j++)
           { tt = aa.get(idx, j);
             aa.cell(idx, j) = aa.get(imax,j);
             aa.cell(imax, j) = tt;
           }
       }

//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    Normalize the diagonal element and change the row as well.
//    Determinant is multiplied by the value of diagonal element.
//  ___________________________________________________________ 

       tt = aa.get(idx, idx);
       fct *= tt;                      // det multiplied by
                                       // diagonal element.
       for (j=idx; j<ndim; j++)
         { aa.cell(idx, j) /= tt; }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    Linearly combine the idx row and other i row to nullify
//    the (i, idx) elements.
// _________________________________________________________ 

       for (i=idx+1; i<ndim; i++)
         { tt = aa.get(i, idx);
           for (j=idx; j<ndim; j++) aa.cell(i, j) -= (aa.get(idx,j) * tt);
         }
    }
  aa.~Matrix();  //deconstruct the working space.
  return fct;
}

///***************************************************
//** Find Inverse using Guass-Jordan elimination
//***************************************************/

Matrix invMatrix(const Matrix &mtx)
{ double tt, tx;
  int i, j, idx, imax, ndim;
  Matrix aa(mtx);                  // Copy THIS matrix to aa for manuipul.
  ndim = aa.ncol();
  Matrix bb(ndim, ndim, 1);       // creat a unit matrix to store
                                   // the inverse matrix.

  for (idx=0; idx<ndim; idx++)       // the main loop.
    {

//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//   Switch the larger element to diagonal.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
      tx = abs(aa.get(idx, idx));
      imax = idx;
      for (i=idx+1; i<ndim; i++)
        { tt = abs(aa.get(i, idx));
          if (tt > tx) { tx = tt;
                         imax = i; }
        }
      if (imax != idx)                       // Switch row(imax) with
                                             // row(idx)
       {
         for (j=idx; j<ndim; j++)
           { tt = aa.get(idx, j);
             aa.cell(idx,  j) = aa.get(imax,j);
             aa.cell(imax, j) = tt;
           }
         for (j=0; j<ndim; j++)
           { tt = bb.get(idx, j);
             bb.cell(idx,  j) = bb.get(imax,j);
             bb.cell(imax, j) = tt;
           }
        }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//   Normalize the diagonal element and change the row as well.
// ___________________________________________________________ 

       tt = aa.get(idx, idx);
       for (j=idx; j<ndim; j++) aa.cell(idx, j) /= tt;
       for (j=0;   j<ndim; j++) bb.cell(idx, j) /= tt;

  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Linearly combine the idx row and other i row to zero
    the (i, idx) elements.
    _________________________________________________________  */
       for (i=0; i<idx; i++)
         { tt = aa.get(i, idx);
           for (j=idx; j<ndim; j++) aa.cell(i, j) -= aa.get(idx,j)*tt;
           for (j=0;   j<ndim; j++) bb.cell(i, j) -= bb.get(idx,j)*tt;
         }

       for (i=idx+1; i<ndim; i++)
         { tt = aa.get(i, idx);
           for (j=idx; j<ndim; j++) aa.cell(i, j) -= aa.get(idx,j)*tt;
           for (j=0; j<ndim; j++) bb.cell(i, j) -= bb.get(idx,j)*tt;
         }
    }
  aa.~Matrix(); 
  return bb;
}
colVector linearSolve(const Matrix&mtx, const colVector& bvct)
{ double tt, tx;
  int i, j, idx, imax, ndim;
  Matrix aa(mtx);                  // Copy THIS matrix to aa for manuipul.
  ndim = aa.ncol();
  colVector bb(bvct);              // copy the right-hand-side vector
                                   // the inverse matrix.

  for (idx=0; idx<ndim; idx++)       // the main loop.
    {

//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//   Switch the larger element to diagonal.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
      tx = abs(aa.get(idx, idx));
      imax = idx;
      for (i=idx+1; i<ndim; i++)
        { tt = abs(aa.get(i, idx));
          if (tt > tx) { tx = tt;
                         imax = i; }
        }
      if (imax != idx)                       // Switch row(imax) with
                                             // row(idx)
       {
         for (j=idx; j<ndim; j++)
           { tt = aa.get(idx, j);
             aa.cell(idx,  j) = aa.get(imax,j);
             aa.cell(imax, j) = tt;
           }
         tt = bb.get(idx);
         bb[idx]  = bb[imax];
         bb[imax] = tt;
        }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//   Normalize the diagonal element and change the row as well.
// ___________________________________________________________ 

       tt = aa.get(idx, idx);
       for (j=idx; j<ndim; j++) aa.cell(idx, j) /= tt;
       bb[idx] /= tt;

  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Linearly combine the idx row and other i row to zero
    the (i, idx) elements.
    _________________________________________________________  */
       for (i=0; i<idx; i++)
         { tt = aa.get(i, idx);
           for (j=idx; j<ndim; j++) aa.cell(i, j) -= aa.get(idx,j)*tt;
           bb[i] -= bb[idx] * tt;
         }

       for (i=idx+1; i<ndim; i++)
         { tt = aa.get(i, idx);
           for (j=idx; j<ndim; j++) aa.cell(i, j) -= aa.get(idx,j)*tt;
           bb[i] -= bb[idx] * tt;
         }
    }
  aa.~Matrix(); 
  return bb;
}
#endif
