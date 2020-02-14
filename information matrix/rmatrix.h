// ==============================================
//  Class for n x m matrix of real elements
//                   1st Edition  2010, May 20
//  written by Yan-Ten Lu
//             Department of Physics,
//             National Cheng Kung University,
//             Tainan, Taiwan 70101
//             ytlu@mail.ncku.edu.tw
// ==============================================
#ifndef _Real_General_Matrix__
#define _Real_General_Matrix__
#define CV(i, j) (i * this->n_col + j)
#define DefaultPrintMode  0
#define DefaultPrecision  16
#define DefaultWidth      7
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;
#include "random.h"
class Matrix
{
      protected:
                int n_row, n_col; // number of rows, and columns of the matrix
                int n_dim;        // dimension = row * col;
                double *rptr;     // complex pointer for the storing array
                                  // array: dynamically allocated in realization
                int PrintMode;
      public:
      // constructors
                Matrix(){PrintMode = DefaultPrintMode;}     // empty constructor
                Matrix(const int, const int);
                Matrix(const Matrix&); // copy constructor
                // copy a non-object matrix, given the address of the element [0][0]
                Matrix(const int, const int, double*); // copy a non-object matrix
                Matrix(const int, const int, const int);// special matrixes
                        // n,n,1 unit matrix
                        // n,m,0 zero matrix
                Matrix(const int n) // square matrix.
                        { *this = Matrix(n, n); }
                Matrix(ifstream&);
                ~Matrix()
                        { if (n_dim>0)delete [] this->rptr;
                          n_col = n_row = n_dim = 0;
                        }
      // operators of member function
 inline void    printmode(int p) { PrintMode = p;}
      Matrix& operator=(const Matrix&);
      Matrix& operator+=(const Matrix&);
      Matrix& operator-=(const Matrix&);
      Matrix& operator*=(const Matrix&);
      Matrix& operator*=(const double);
      Matrix& operator/=(const double);
      // utility function
 inline double& operator[](const int i) {return this->rptr[i];}
 inline double  get(const int i, const int j) const
                    { return this->rptr[i*n_col+j]; }
 inline double  get(const int i) const
                    { return this->rptr[i]; }
 inline double& cell(const int i, const int j)
                    { return this->rptr[ i*n_col + j ]; }
 inline double& cell(const int i)
                    { return this->rptr[i]; }
 inline     int nrow() const {return this->n_row;}
 inline     int ncol() const {return this->n_col;}
 inline     int ndim() const {return this->n_dim;}
           void random(unsigned long int&);
           double norm()const // sum of square of all elements
                 { int i;
                   double sum, tt;
                   sum =0.0;
                   for (i=0; i<this->n_dim; i++)
                       { tt = this->rptr[i];
                         sum += tt*tt;
                       }
                   return sum;
                 }
 inline    double abs() const { return sqrt( this->norm() );}  // square roor of norm();
           void unit()
                     { int i;
                       double vabs;
                       vabs = this->abs();
                       for (i=0; i<this->n_dim; i++) this->rptr[i] /= vabs;
                       return;
                     }
       Matrix transport()const;
//       Matrix adjoint()const; // for complex matrix only
           void readfile(ifstream&); // read matrix from a ifstream object
           void writefile(char*)const;  // write matrix to file (name)
      // friend functions
      friend ostream& operator<<(ostream&, const Matrix&);
      friend Matrix operator+(const Matrix&, const Matrix&);
      friend Matrix operator-(const Matrix&, const Matrix&);
      friend Matrix operator*(const Matrix&, const Matrix&);
      friend Matrix operator*(const double, const Matrix&);
      friend Matrix operator*(const Matrix&, const double);
      friend Matrix operator/(const Matrix&, const double);
      friend class colVector;
      friend class rowVector;
}; // end definition of class Matrix
// ==================================================
// construct an nxm empty matrix
// --------------------------------------------------
Matrix::Matrix(const int n, const int m)
{
  this->PrintMode = DefaultPrintMode;
  this->n_row = n;
  this->n_col = m;
  this->n_dim = n * m;
  if (n_dim > 0)
     {  if ( (this->rptr = new double[n_dim]) == NULL )
           { cerr << "Failed in memory allocation in Matrix [E002]~!\n";
             n_dim = n_row = n_col = 0;
             return;
           }
     }
  else
     { cerr << "Matirx dimesion error~[E001]!\n";
       n_dim = n_row = n_col = 0;
       return;
     }
} // end of an n x m emmpty matrix
//=======================================================
// constructor by copying another Matrix
//-------------------------------------------------------
Matrix::Matrix(const Matrix &mx1) // copy constructor
{ int i;
  this->PrintMode = mx1.PrintMode;
  this->n_row = mx1.n_row;
  this->n_col = mx1.n_col;
  this->n_dim = this->n_row * this->n_col;
  if (this->n_dim > 0)
     {  if ( (this->rptr = new double[n_dim]) == NULL )
           { cerr << "Failed in memory allocation in Matrix [E004]~!\n";
             n_dim = n_row = n_col = 0;
             return;
           }
        else // copy the elements of matrix to new array.
           { for (i=0; i<n_dim; i++) this->rptr[i] = mx1.rptr[i];
             return;
           }
     }
  else
     { cerr << "Matirx dimesion error~[E003]!\n";
       n_dim = n_row = n_col = 0;
       return;
     }
} // end of copy constructor
//======================================================
// initail a unit matrix (id = 1) or a zero matrix (id = 0)
//=======================================================
Matrix::Matrix(const int n, const int m, const int id)   // special matrixes
{ int i, j, jstart;
  this->PrintMode = DefaultPrintMode;
  this->n_row = n;
  this->n_col = m;
  this->n_dim = this->n_row * this->n_col;
  if (this->n_dim > 0)
     {  if ( (this->rptr = new double[n_dim]) == NULL )
           { cerr << "Failed in memory allocation in Matrix [E008]~!\n";
             n_dim = n_row = n_col = 0;
             return;
           }
     }
  else
     { cerr << "Matirx dimesion error~[E007]!\n";
       n_dim = n_row = n_col = 0;
       return;
     }
    if (id == 0)
       { for (i=0; i<n_dim; i++) this->rptr[i] = 0.0;
         return;
       }
    else
       { for (i=0; i<n_row; i++)
             {  jstart = i * this->n_col;
                for (j=0; j<n_col; j++) rptr[jstart+j] = (i==j)? 1.0 : 0.0;
             }
       }
   return;
} // end of unit matrix and zero matrix
//======================================================
// initial a matri by n x m 2-dim array -- giving the address of first element
//------------------------------------------------------
Matrix::Matrix(const int n, const int m, double *pt) // copy a non-object matrix
{ int i;
  this->PrintMode = DefaultPrintMode;
  this->n_row = n;
  this->n_col = m;
  this->n_dim = n * m;
  if (n_dim > 0)
     {   this->rptr = new double [n_dim];
         if ( this->rptr == NULL )
           { cout << "Failed in memory allocation in Matrix [E006]~!\n";
             n_dim = n_row = n_col = 0;
             return;
           }
        else // copy matix
           {
                for (i=0; i<n_dim; i++)
                    {  this->rptr[i] = pt[i];
                    }
           }
     }
  else
     { cout << "Matirx dimesion error~[E005]!\n";
       n_dim = n_row = n_col = 0;
       return;
     }
}// end of input from n x m array.
Matrix Matrix::transport()const
{ int i, j;
  Matrix mx(n_col, n_row);
  for (i=0; i<mx.n_row; i++)
      { for (j=0; j<mx.n_col; j++) mx.cell(i,j)=this->get(j,i);
      }
  return mx;
}
// Matrix Matrix::adjoint()const
// { int i, j;
//   Matrix mx(n_col, n_row);
//   for (i=0; i<mx.n_row; i++)
//       { for (j=0; j<mx.n_col; j++) mx.cell(i,j)=conj(this->get(j,i));
//       }
//   return mx;
//}
// ======================================================================
// Overload ostream to print out the Matrix object
// PrintMode 0: print each element a line.
//           1: pinrt in a matrix block (test viewing).
// =======================================================================
ostream& operator<<(ostream& odev, const Matrix& mx)
{ int i, j;
  if (mx.PrintMode == 0)
  {
    odev << " row = " << mx.n_row << ", col = " << mx.n_col << "\n";
    for (i=0; i<mx.n_row; i++)
        {  for (j=0; j<mx.n_col; j++)
            odev << "Cell (" << i <<", " << j <<") = " << mx.get(i,j) <<"\n";
        }
  }
  else
  { cout << setprecision(DefaultPrecision) << fixed;
       cout << " row = " << mx.n_row << ", col = " << mx.n_col << "\n";
     for (i=0; i<mx.n_row; i++)
         {  for (j=0; j<mx.n_col; j++) cout << setw(7) << mx.get(i,j);
            cout <<"\n";
         }
  }
  return odev;
}// end of overlaod cout.
istream& operator>>(istream& idev, Matrix& mx)
{ int i;
  for (i=0; i<mx.ndim(); i++) idev >> mx.cell(i);
  return idev;
}// end of overlaod cout.
//==================================================
// member function operator=
//==================================================
Matrix& Matrix::operator=(const Matrix &mx) //overlaod equal sign
{ int i;
// deconstruct the original array, if it existed.
  if (this->n_dim >0) delete [] this->rptr;
// the following is same as the copy constructor.
  this->PrintMode = mx.PrintMode;
  this->n_row = mx.n_row;
  this->n_col = mx.n_col;
  this->n_dim = this->n_row * this->n_col;
  if (this->n_dim > 0)
     {  if ( (this->rptr = new double[n_dim]) == NULL )
           { cerr << "Failed in memory allocation in Matrix [E010]~!\n";
             n_dim = n_row = n_col = 0;
             return *this;
           }
        else // copy the elements of matrix to new array.
           { for (i=0; i<n_dim; i++) this->rptr[i] = mx.rptr[i];
             return *this;  // this is different from the copy constructor.
           }
     }
  else
     { cerr << "Matirx dimesion error~[E009]!\n";
       n_dim = n_row = n_col = 0;
       return *this;
     }
} // end of operator=
//==================================================
// member function operator+=
//==================================================
Matrix& Matrix::operator+=(const Matrix& mx) //overlaod equal sign
{ int i;
// check validity of + operation.
  if ((this->n_row!=mx.n_row) || (this->n_col!=mx.n_col))
     { cerr << "Operator += invalid~! [E011]\n";
       return *this;
     }
// adding up two matrixes.
  for (i=0; i<n_dim; i++) this->rptr[i] += mx.rptr[i];
  return *this;
} // end of operator+=
//==================================================
// member function operator-=
//==================================================
Matrix& Matrix::operator-=(const Matrix& mx) //overlaod equal sign
{ int i;
// check validity of - operation.
  if ((this->n_row!=mx.n_row) || (this->n_col!=mx.n_col))
     { cerr << "Operator += invalid~! [E012]\n";
       return *this;
     }
// remove mx from this.
  for (i=0; i<n_dim; i++) this->rptr[i] -= mx.rptr[i];
  return *this;
} // end of operator-=
Matrix& Matrix::operator*=(const Matrix &mx)
{ int i, j, k, n, r, c; double sum;
// check validity of * operation.
  n = mx.n_row;
  if (this->n_col != n)
     { cerr << "Operator *= invalid~! [E019]\n";
       delete [] this->rptr;
       n_dim = n_row = n_col = 0;
       return *this;
     }
  r = this->n_row;
  c = mx.n_col;
  Matrix mp(r, c);
  for (i=0; i<r; i++)
      {  for (j=0; j<c; j++)
              {  sum = 0.0;
                 for (k=0; k<n; k++) sum += this->get(i,k)*mx.get(k,j);
                 mp.cell(i,j) = sum;
              }
      }
  *this = mp;
  return *this;
} // end of operator*=
Matrix& Matrix::operator*=(const double r)
{ int i;
  for (i=0; i<this->n_dim; i++)
      { this->rptr[i] *= r;
      }
  return *this;
} // end of operator*=
Matrix& Matrix::operator/=(const double r)
{ int i;
  for (i=0; i<this->n_dim; i++)
      { this->rptr[i] /= r;
      }
  return *this;
} // end of operator/=

//==================================================
// friend function operator+
//==================================================
Matrix operator+(const Matrix& m1, const Matrix& m2)
{ int i;
  Matrix mt(m1);
// check validity of + operation.
  if ((mt.n_row!=m2.n_row) || (mt.n_col!=m2.n_col))
     { cerr << "Operator + invalid~! [E013]\n";
       return mt;
     }
// adding up two matrixes.
  for (i=0; i<mt.n_dim; i++) mt.rptr[i] += m2.rptr[i];
  return mt;
} // end of operator+
//==================================================
// friend function operator-
//==================================================
Matrix operator-(const Matrix& m1, const Matrix& m2)
{ int i;
  Matrix mt(m1);
// check validity of + operation.
  if ((mt.n_row!=m2.n_row) || (mt.n_col!=m2.n_col))
     { cerr << "Operator - invalid~! [E014]\n";
       return mt;
     }
// adding up two matrixes.
  for (i=0; i<mt.n_dim; i++) mt.rptr[i] -= m2.rptr[i];
  return mt;
} // end of operator-
//==================================================
// friend function operator*
//==================================================
Matrix operator*(const Matrix& m1, const Matrix& m2)
{ int i, j, k;
  double sum;
  Matrix mt(m1.n_row, m2.n_col);
// check validity of + operation.
  if ((m1.n_col!=m2.n_row))
     { cerr << "Operator * invalid~! [E015]\n";
       return mt;
     }
// hook up two matrixes.
  for (i=0; i<mt.n_row; i++)
      { for (j=0; j<mt.n_col; j++)
            { sum = 0.0;
              for (k=0; k<m1.n_col; k++) sum += ( m1.get(i,k)*m2.get(k,j) );
              mt.cell(i,j) = sum;
            }
      }
  return mt;
} // end of operator* a c-number
Matrix operator*(const double c2, const Matrix& m1)
{ int k, n;
  Matrix mt(m1);
  n = mt.ndim();
  for (k=0; k<n; k++) mt[k] *= c2;
  return mt;
} // end of operator*
Matrix operator*(const Matrix& m1, const double c2)
{ int k, n;
  Matrix mt(m1);
  n = mt.ndim();
  for (k=0; k<n; k++) mt[k] *= c2;
  return mt;
} // end of operator*
Matrix operator/(const Matrix& m1, const double c2)
{ int k, n;
  Matrix mt(m1);
  n = mt.ndim();
  for (k=0; k<n; k++) mt[k] /= c2;
  return mt;
} // end of operator*
//========================================================
// member function; read a matrix from ifstream
//--------------------------------------------------------
void Matrix::readfile(ifstream& fpt)
{  int i;
//   if (this->n_dim > 0) delete [] this->rptr;
//   fpt >> this->n_row >> this->n_col;  // extract n_row and n_col
//   this->n_dim = this->n_row * this->n_col;
//   if ((this->rptr = new double[this->n_dim]) != NULL)
      for (i=0; i<this->n_dim; i++) fpt >> this->rptr[i];
        return;
}
// end of void Matrix::readfile(ifstream&)
//========================================================
//constructor: read a matrix from ifstream
//--------------------------------------------------------
Matrix::Matrix(ifstream& fpt)
{  int i;
   this->PrintMode = DefaultPrintMode;
   fpt >> this->n_row >> this->n_col;
   this->n_dim = this->n_row * this->n_col;
   if ((this->rptr = new double[this->n_dim]) != NULL)
      {   for (i=0; i<this->n_dim; i++)
              { fpt >> this->rptr[i];
              }
          return;
      }
   else { cerr << "Memory fault in Matrix~! [E017]\n";
          this->n_dim = this->n_col = this->n_row = 0;
          return;
        }
}
// constructor: read a matrix from ifstream
//========================================================
// member function; write a matrix to a file (name)
//--------------------------------------------------------
void Matrix::writefile(char*fname) const
{  int i, j;
   ofstream fpt(fname);
   fpt << this->n_row << "  " << this->n_col <<"\n";
   for(i=0; i<n_row; i++)
      { for (j=0; j<n_col;j++)
            { fpt <<  "  " << this->get(i,j);
            }
        fpt << "\n";
      }
      fpt.close();
}
// end of Matrix::writefile(char*) const
//========================================================
// member function: generate a matrix of random matrix
//-------------------------------------------------------
void Matrix::random(unsigned long int&seed)
{ int i;
  double vabs;
  for (i=0; i<this->n_dim; i++) { this->rptr[i] = RRANDOM(seed); }
  vabs = this->abs();
  for (i=0; i<this->n_dim; i++) { this->rptr[i] /= vabs; }

}
// end of Matrix::random(unsigned long int&seed)
#endif
