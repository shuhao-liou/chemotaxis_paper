#ifndef _Complex_Vector_Class__
#define _Complex_Vector_Class__
#ifndef _Complex_General_Matrix__
#define _Complex_General_Matrix__
#include "cmatrix.h"
#endif
#ifndef _ytlu_randomgenerator_1_
#define _ytlu_randomgenerator_1_
#define RANDOM_RHO  0x7a35ef05UL  // rho = 8t-3, 8-byte unsigned long
#define RANDOM_MAX  ((double)0x40000000)  // maximan of congruence sequence
#define RANDOM_HALF ((double)0x20000000)
#define NRANDOM(nn) ((nn *= RANDOM_RHO) >> 2)
#define RRANDOM(nn) ( (((nn *= RANDOM_RHO) >> 2)-RANDOM_HALF) / RANDOM_HALF)
#endif
#include <complex>
using namespace std;
class C_rowVector : public C_Matrix
{    
    public:
    // constructors ****/ 
        C_rowVector (): C_Matrix(){;}   
        C_rowVector (const int n): C_Matrix(1,n){};
        C_rowVector (const int n, COMPLX *dp): C_Matrix(1, n, dp){};
        C_rowVector (const C_Matrix & mtx);
        C_rowVector (const int n, const COMPLX x);
        C_rowVector (const C_rowVector & vct):C_Matrix(1, vct.n_col, vct.cptr){};
// member functions
       COMPLX operator[](const int i) const {return this->cptr[i];}
       void    random(unsigned long int&);
       double  vnorm() const 
                   { int i;
                     COMPLX tt;
                     double sum = 0.0;
                     for (i=0; i<n_dim; i++)
                         { tt = cptr[i];
                           sum += norm(tt);
                         }
                     return sum;
                   }  
       double  abs() const { return sqrt(this->vnorm()); }
       C_rowVector unit(){ return(*this/this->abs());}
//
//  inline rowVector ortho(const rowVector v1) const
//          { return rowVector((*this) - (v1*(*this).transport()) * v1);};
//
// friends

   friend class C_Matrix;
   friend class C_colVector;
};   // end of C_rowVector class definition  


class C_colVector : public C_Matrix
{    
    public:
// constructors ****/
        C_colVector (): C_Matrix(){;}   
        C_colVector (const int n): C_Matrix(n,1){};
        C_colVector (const int n, COMPLX *dp): C_Matrix(n, 1, dp){};
        C_colVector (const C_Matrix & mtx);
        C_colVector (const int n, const COMPLX x);
        C_colVector (const C_colVector & vct):C_Matrix(vct.n_row, 1, vct.cptr){};
// member functions
       COMPLX operator[](const int i)const {return this->cptr[i];}
       void    random(unsigned long int&);
       double  vnorm() const 
                   { int i;
                     COMPLX tt;
                     double sum = 0.0;
                     for (i=0; i<n_dim; i++)
                         { tt = cptr[i];
                           sum += norm(tt);
                         }
                     return sum;
                   }  
       double  abs() const { return sqrt(this->vnorm()); }
       C_colVector unit(){ return(*this/this->abs());}
//
//  inline colVector ortho (const colVector v1) const
//          { return colVector((*this) - v1 * ((*this).transport()*v1));};
//
// Friends

   friend class C_Matrix;
   friend class C_rowVector;
//   friend colVector operator- (colVector v1, colVector v2) 
//                { return colVector(v1-v2); }
}; 
// end of colVector class definition  
//
//  
C_rowVector::C_rowVector (const C_Matrix & mtx):C_Matrix(mtx)
{ if (n_row != 1) 
   { cerr << "Equate matrix to row_Vector: Dimesion error.\n";
     this->n_col = this->n_row = 0;
     delete [] this->cptr;
   }
}     
C_rowVector::C_rowVector (const int n, const COMPLX x):C_Matrix(1,n)
{ int i;
  for (i=0; i<n_dim; i++) this->cptr[i] = x;
}//end of C_rowVector::C_rowVector(int, COMPLX)  
C_colVector::C_colVector (const C_Matrix & mtx):C_Matrix(mtx)
{ if (n_col != 1) 
   { fprintf(stderr,"Equate matrix to col_Vector: Dimension error.\n");
     this->n_col = this->n_row = 0;
     delete [] this->cptr;
   }
} // end of  C_colVector::C_colVector(C_Matrix)    
C_colVector::C_colVector (const int n, const COMPLX x):C_Matrix(n,1)
{ int i;
  for (i=0; i<n_dim; i++) this->cptr[i] = x;
}//end of C_colVector::C_colVector(int, COMPLX)

void C_colVector::random(unsigned long int& seed)
{  int i;
   for (i=0; i<this->n_dim; i++) this->cptr[i] = COMPLX(RRANDOM(seed), RRANDOM(seed));
   double sum;
   sum = 0.0;
   for (i=0; i<this->n_dim; i++) sum += norm(this->cptr[i]);
   sum = sqrt(sum);
   for (i=0; i<this->n_dim; i++) this->cptr[i] /= sum;
   return;   
}
void C_rowVector::random(unsigned long int& seed)
{  int i;
   for (i=0; i<this->n_dim; i++) this->cptr[i] = COMPLX(RRANDOM(seed), RRANDOM(seed));
   double sum;
   sum = 0.0;
   for (i=0; i<this->n_dim; i++) sum += norm(this->cptr[i]);
   sum = sqrt(sum);
   for (i=0; i<this->n_dim; i++) this->cptr[i] /= sum;
   return;   
}
C_rowVector row(const C_Matrix & mx, const int row)
{ int i, nc;
  nc = mx.ncol();
  C_rowVector vrow(nc);
  for (i=0; i<nc; i++) vrow[i] = mx.get(row, i); 
  return vrow; 
}
C_colVector col(const C_Matrix & mx, const int col)
{ int i, nc;
  nc = mx.nrow();
  C_colVector vcol(nc);
  for (i=0; i<nc; i++) vcol[i] = mx.get(i, col); 
  return vcol; 
}
C_rowVector diagonal(const C_Matrix& mx)
{ int i, nd;
  nd = (mx.nrow() > mx.ncol())? mx.ncol() : mx.nrow();
  C_rowVector vct(nd);
  for (i=0; i<nd; i++) vct[i] = mx.get(i,i);
  return vct;
}
COMPLX operator*(C_rowVector& r, C_colVector& c)
{ COMPLX sum;
  int i, n;
  n = c.ndim();
  if (n == r.ndim())
     {     sum = complex<double>(0.0, 0.0);
           for (i=0; i<n; i++) sum += (r[i] * c[i]); 
     }
  else cerr << "Dimension error in cvector -- operator<>\n";
  return sum;
}         
#endif

