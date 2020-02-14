#ifndef _Real_Vector_Class__
#define _Real_Vector_Class__
#include "rmatrix.h"
#include <iostream>
#include <iomanip>
using namespace std;
//==============================================================
// Class row Vector
//
//--------------------------------------------------------------
class rowVector : public Matrix
{    
    public:
    // constructors ****/ 
        rowVector (): Matrix(){;}   
        rowVector (const int n): Matrix(1,n){};
        rowVector (const int n, double *dp): Matrix(1, n, dp){};
        rowVector (const Matrix&);
        rowVector (const int n, const double x);
        rowVector (const rowVector &vct):Matrix(1, vct.n_col, vct.rptr){};
        rowVector (ifstream&);
 // member functions
//==============================================================
// The following three operators cast results into row vectors
//--------------------------------------------------------------
      rowVector& operator=(const rowVector &v)
           { this->Matrix::operator=( v ); 
             return *this;}
      rowVector& operator+=(const rowVector&v)                                      
           { this->Matrix::operator+=( v );
             return *this;}                                      
      rowVector& operator-=(const rowVector&v)                                      
           { this->Matrix::operator-=( v );
             return *this;}                                      
      rowVector& operator*=(const double r)                                      
           {  int i;
              for (i=0; i<this->n_dim; i++) this->rptr[i] *= r;
              return *this;}                                      
      rowVector& operator/=(const double r)                                      
           {  int i;
              for (i=0; i<this->n_dim; i++) this->rptr[i] /= r;
              return *this;}                                      
//
// friends

   friend class Matrix;
   friend class colVector;
};   // end of C_rowVector class definition  

rowVector::rowVector(ifstream&f)
{ int i;
  this->n_row = 1;
  f >> this->n_col;
  this->n_dim = this->n_col;
  this->rptr = new double [this->n_dim];
  for (i=0; i<this->n_dim; i++) f >> this->rptr[i];
  return;
}
rowVector::rowVector (const Matrix &mtx):Matrix(mtx)
{ if (n_row != 1) 
   { cerr << "Equate matrix to row_Vector: Dimesion error.\n";
     this->n_col = this->n_row = 0;
     delete [] this->rptr;
   }
}     
rowVector::rowVector (const int n, const double x):Matrix(1,n)
{ int i;
  for (i=0; i<n_dim; i++) this->rptr[i] = x;
}//end of C_rowVector::C_rowVector(int, COMPLX)  
//==============================================================
// The following two operators cast results into row vectors
//--------------------------------------------------------------
rowVector operator+(const rowVector&v1, const rowVector&v2)
       {
          rowVector v3(v1);
          v3 += v2;
          return v3;  }
rowVector operator-(const rowVector&v1, const rowVector&v2)     
       {
          rowVector v3(v1);
          v3 -= v2;
          return v3;  }
//==============================================================
// Class column Vector
//
//--------------------------------------------------------------
class colVector : public Matrix
{    
    public:
// constructors ****/
        colVector (): Matrix(){;}   
        colVector (const int n): Matrix(n,1){};
        colVector (const int n, double *dp): Matrix(n, 1, dp){};
        colVector (const Matrix &mtx);
        colVector (const int n, const double x);
        colVector (const colVector &vct):Matrix(vct.n_row, 1, vct.rptr){};
        colVector (ifstream&fpt);
// member functions
//==============================================================
// The following three operators cast results into column vectors
//--------------------------------------------------------------
      colVector& operator=(const colVector &v)
           { this->Matrix::operator=( v ); 
             return *this;}
      colVector& operator+=(const colVector&v)
           { this->Matrix::operator+=( v );
             return *this;}                                      
      colVector& operator-=(const colVector&v)
           { this->Matrix::operator-=( v );
             return *this;}                                      
      colVector& operator*=(const double r)                                      
           {  int i;
              for (i=0; i<this->n_dim; i++) this->rptr[i] *= r;
              return *this;}                                      
      colVector& operator/=(const double r)                                      
           {  int i;
              for (i=0; i<this->n_dim; i++) this->rptr[i] /= r;
              return *this;}                                      
//
//  inline colVector ortho (const colVector v1) const
//          { return colVector((*this) - v1 * ((*this).transport()*v1));};
//
// Friends

   friend class Matrix;
   friend class rowVector;
   friend class Tridiagonal;
//   friend colVector operator- (colVector v1, colVector v2) 
//                { return colVector(v1-v2); }
}; 
// end of colVector class definition  
//
//  
colVector::colVector(ifstream&f)
{ int i;
  this->n_col = 1;
  f >> this->n_row;
  this->n_dim = this->n_row;
  this->rptr = new double [this->n_dim];
  for (i=0; i<this->n_dim; i++) f >> this->rptr[i];
  return;
}
colVector::colVector (const Matrix & mtx):Matrix(mtx)
{ if (n_col != 1) 
   { fprintf(stderr,"Equate matrix to col_Vector: Dimension error.\n");
     this->n_col = this->n_row = 0;
     delete [] this->rptr;
   }
} // end of  C_colVector::C_colVector(C_Matrix)    
colVector::colVector (const int n, const double x):Matrix(n,1)
{ int i;
  for (i=0; i<n_dim; i++) this->rptr[i] = x;
}//end of C_colVector::C_colVector(int, COMPLX)
//==============================================================
// The following two operators redirect + & - into column vector
//--------------------------------------------------------------
colVector operator+(const colVector&v1, const colVector&v2)
       {
          colVector v3(v1);
          v3 += v2;
          return v3;  }
colVector operator-(const colVector&v1, const colVector&v2)     
       {
          colVector v3(v1);
          v3 -= v2;
          return v3;  }
//======================================================
// extract a row or a column from a matrix
// and put int a row or a col
//------------------------------------------------------
rowVector row(const Matrix & mx, const int nrow)
{ int i, nc;
  nc = mx.ncol();
  rowVector vrow(nc);
  for (i=0; i<nc; i++) vrow[i] = mx.get(nrow, i); 
  return vrow; 
}
colVector col(const Matrix & mx, const int ncol)
{ int i, nc;
  nc = mx.nrow();
  colVector vcol(nc);
  for (i=0; i<nc; i++) vcol[i] = mx.get(i, ncol); 
  return vcol; 
}
void putrow(Matrix & mx, const rowVector & r, const int nrow)
{ int i, nc;
  nc = mx.ncol(); 
  for (i=0; i<nc; i++) mx.cell(nrow, i) = r.get(i); 
  return; 
}
void putcol(Matrix & mx, const colVector & c, const int ncol)
{ int i, nr;
  nr = mx.nrow(); 
  for (i=0; i<nr; i++) mx.cell(i, ncol) = c.get(i); 
  return; 
}
//======================================================
// extract the diagonal matrix element from a matrix
//------------------------------------------------------
colVector diagonal(const Matrix& mx)
{ int i, nd;
  nd = (mx.nrow() > mx.ncol())? mx.ncol() : mx.nrow();
  colVector vct(nd);
  for (i=0; i<nd; i++) vct[i] = mx.get(i,i);
  return vct;
}
//
//
double operator*(const rowVector& r, const colVector& c)
{ double sum;
  int i, n;
  n = c.ndim();
  if (n == r.ndim())
     {     sum = 0.0;
           for (i=0; i<n; i++) sum += (r.get(i) * c.get(i)); 
     }
  else cerr << "Dimension error in cvector -- operator<>\n";
  return sum;
}
//=========================================================
// cast operators into vectors
//---------------------------------------------------------
colVector operator*(const Matrix &mx, const colVector &v)
           { 
             Matrix m2(mx);
             m2 *= v;
             return colVector(m2);
           }    
rowVector operator*(const rowVector&v, const Matrix &mx)
       {
          rowVector v3(v);
          v3.Matrix::operator*=( mx);
          return v3;  }
rowVector operator*(const double a, const rowVector &v)
           { int i;
             rowVector vp(v);
             for (i=0; i<v.ndim(); i++) vp[i] *= a; 
             return vp; }        
rowVector operator*(const rowVector &v, const double a)
           { return (a * v); }                
rowVector operator/(const rowVector &v, const double a)
           { int i;
             colVector vp(v);
             for (i=0; i<v.ndim(); i++) vp[i] /= a; 
             return vp; }
colVector operator*(const double a, const colVector &v)
           { int i;
             colVector vp(v);
             for (i=0; i<v.ndim(); i++) vp[i] *= a; 
             return vp; }               
colVector operator/(const colVector &v, const double a)
           { int i;
             colVector vp(v);
             for (i=0; i<v.ndim(); i++) vp[i] /= a; 
             return vp; }
colVector transport(rowVector &v)
           { return colVector(v.transport());}        
rowVector transport(colVector &v)
           { return rowVector(v.transport());}        
#endif

