// ==============================================
//  Class for n x m matrix of complex elements
//                   1st Edition  2010, May 20
//  written by Yan-Ten Lu
//             Department of Physics,
//             National Cheng Kung University,
//             Tainan, Taiwan 70101
//             ytlu@mail.ncku.edu.tw
// ==============================================
#ifndef _Complex_General_Matrix__
#define _Complex_General_Matrix__
#define COMPLX complex<double>
#define CV(i, j) (i * this->n_col + j)
#define DefaultPrintMode  0
#define ZERO (complex<double>(0.0, 0.0))
#define ONE (complex<double>(1.0, 0.0))
#include <iostream>
#include <complex>
#include <iomanip>
#include <fstream>
using namespace std;
class C_Matrix
{
      protected:
                int n_row, n_col; // number of rows, and columns of the matrix
                int n_dim;        // dimension = row * col;
                COMPLX *cptr;     // complex pointer for the storing array
                                  // array: dynamically allocated in realization
                int PrintMode;
      public:
      // constructors
                C_Matrix(){PrintMode = DefaultPrintMode;}     // empty constructor
                C_Matrix(const int, const int);
                C_Matrix(const C_Matrix&); // copy constructor
                C_Matrix(const int n) // square matrix.
                        { C_Matrix(n, n); }
                // copy a non-object matrix, given the address of the element [0][0]
                C_Matrix(const int, const int, COMPLX*); // copy a non-object matrix
                C_Matrix(const int, const int, const int);// special matrixes
                        // n,n,1 unit matrix
                        // n,m,0 zero matrix
                C_Matrix(ifstream&);
                ~C_Matrix()
                        { if (n_dim>0)delete [] this->cptr;
                          n_col = n_row = n_dim = 0;
                        }
      // operators of member function
      C_Matrix& operator=(const C_Matrix&);
      C_Matrix& operator+=(const C_Matrix&);
      C_Matrix& operator-=(const C_Matrix&);
      // utility function
         COMPLX get(const int i, const int j) const
                    { return this->cptr[i*n_col+j]; }
         COMPLX get(const int i) const
                    { return this->cptr[i]; }
        COMPLX& cell(const int i, const int j)
                    { return this->cptr[i*n_col+j]; }
        COMPLX& cell(const int i)
                    { return this->cptr[i]; }
  inline    int nrow() const {return this->n_row;}
  inline    int ncol() const {return this->n_col;}
  inline    int ndim() const {return this->n_dim;}
       C_Matrix transport()const;
       C_Matrix adjoint()const;
           void readfile(ifstream&);
           void writefile(char*);
      // friend functions
      friend ostream& operator<<(ostream&, const C_Matrix&);
      friend C_Matrix operator+(const C_Matrix&, const C_Matrix&);
      friend C_Matrix operator-(const C_Matrix&, const C_Matrix&);
      friend C_Matrix operator*(const C_Matrix&, const C_Matrix&);
      friend class C_colVector;
      friend class C_rowVector;
}; // end definition of class C_Matrix
// ==================================================
// construct an nxm empty matrix
// --------------------------------------------------
C_Matrix::C_Matrix(const int n, const int m)
{
  this->PrintMode = DefaultPrintMode;
  this->n_row = n;
  this->n_col = m;
  this->n_dim = n * m;
  if (n_dim > 0)
     {  if ( (this->cptr = new COMPLX[n_dim]) == NULL )
           { cerr << "Failed in memory allocation in C_Matrix [E002]~!\n";
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
C_Matrix::C_Matrix(const C_Matrix &mx1) // copy constructor
{ int i;
  this->PrintMode = mx1.PrintMode;
  this->n_row = mx1.n_row;
  this->n_col = mx1.n_col;
  this->n_dim = this->n_row * this->n_col;
  if (this->n_dim > 0)
     {  if ( (this->cptr = new COMPLX[n_dim]) == NULL )
           { cerr << "Failed in memory allocation in C_Matrix [E004]~!\n";
             n_dim = n_row = n_col = 0;
             return;
           }
        else // copy the elements of matrix to new array.
           { for (i=0; i<n_dim; i++) this->cptr[i] = mx1.cptr[i];
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
C_Matrix::C_Matrix(const int n, const int m, const int id)   // special matrixes
{ int i, j, jstart;
  this->PrintMode = DefaultPrintMode;
  this->n_row = n;
  this->n_col = m;
  this->n_dim = this->n_row * this->n_col;
  if (this->n_dim > 0)
     {  if ( (this->cptr = new COMPLX[n_dim]) == NULL )
           { cerr << "Failed in memory allocation in C_Matrix [E008]~!\n";
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
       { for (i=0; i<n_dim; i++) this->cptr[i] = ZERO;
         return;
       }
    else
       { for (i=0; i<n_row; i++)
             {  jstart = i * this->n_col;
                for (j=0; j<n_col; j++) cptr[jstart+j] = (i==j)? ONE : ZERO;
             }
       }
   return;
} // end of unit matrix and zero matrix
//======================================================
// initial a matri by n x m 2-dim array
//------------------------------------------------------
C_Matrix::C_Matrix(const int n, const int m, COMPLX *pt) // copy a non-object matrix
{ int i;
  this->PrintMode = DefaultPrintMode;
  this->n_row = n;
  this->n_col = m;
  this->n_dim = n * m;
  if (n_dim > 0)
     {   this->cptr = new complex<double>[n_dim];
         if ( this->cptr == NULL )
           { cout << "Failed in memory allocation in C_Matrix [E006]~!\n";
             n_dim = n_row = n_col = 0;
             return;
           }
        else // copy matix
           {
                for (i=0; i<n_dim; i++)
                    {  this->cptr[i] = pt[i];
                    }
           }
     }
  else
     { cout << "Matirx dimesion error~[E005]!\n";
       n_dim = n_row = n_col = 0;
       return;
     }
}// end of input from n x m array.
C_Matrix C_Matrix::transport()const
{ int i, j;
  C_Matrix mx(n_col, n_row);
  for (i=0; i<mx.n_row; i++)
      { for (j=0; j<mx.n_col; j++) mx.cell(i,j)=this->get(j,i);
      }
  return mx;
}
C_Matrix C_Matrix::adjoint()const
{ int i, j;
  C_Matrix mx(n_col, n_row);
  for (i=0; i<mx.n_row; i++)
      { for (j=0; j<mx.n_col; j++) mx.cell(i,j)=conj(this->get(j,i));
      }
  return mx;
}
// ======================================================================
// Overload ostream to print out the Matrix object
// PrintMode 0: print each element a line.
//           1: pinrt in a matrix block (test viewing).
// =======================================================================
ostream& operator<<(ostream& odev, const C_Matrix& mx)
{ int i, j;
  if (mx.PrintMode == 0)
  {
    odev << " row = " << mx.n_row << ", col = " << mx.n_col << "\n";
    for (i=0; i<mx.n_row; i++)
        {  for (j=0; j<mx.n_col; j++)
            odev << "(" << i <<", " << j <<") = " << mx.get(i,j) <<"\n";
        }
  }
  else
  {  cout << " row = " << mx.n_row << ", col = " << mx.n_col << "\n";
     for (i=0; i<mx.n_row; i++)
         {  for (j=0; j<mx.n_col; j++) cout << setw(10) << mx.get(i,j);
            cout <<"\n";
         }
  }
  return odev;
}// end of overlaod cout.
//==================================================
// member function operator=
//==================================================
C_Matrix& C_Matrix::operator=(const C_Matrix& mx) //overlaod equal sign
{ int i;
// deconstruct the original array, if it existed.
  if (this->n_dim >0) delete [] this->cptr;
// the following is same as the copy constructor.
  this->PrintMode = mx.PrintMode;
  this->n_row = mx.n_row;
  this->n_col = mx.n_col;
  this->n_dim = this->n_row * this->n_col;
  if (this->n_dim > 0)
     {  if ( (this->cptr = new COMPLX[n_dim]) == NULL )
           { cerr << "Failed in memory allocation in C_Matrix [E010]~!\n";
             n_dim = n_row = n_col = 0;
             return *this;
           }
        else // copy the elements of matrix to new array.
           { for (i=0; i<n_dim; i++) this->cptr[i] = mx.cptr[i];
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
C_Matrix& C_Matrix::operator+=(const C_Matrix& mx) //overlaod equal sign
{ int i;
// check validity of + operation.
  if ((this->n_row!=mx.n_row) || (this->n_col!=mx.n_col))
     { cerr << "Operator += invalid~! [E011]\n";
       return *this;
     }
// adding up two matrixes.
  for (i=0; i<n_dim; i++) this->cptr[i] += mx.cptr[i];
  return *this;
} // end of operator+=
//==================================================
// member function operator-=
//==================================================
C_Matrix& C_Matrix::operator-=(const C_Matrix& mx) //overlaod equal sign
{ int i;
// check validity of - operation.
  if ((this->n_row!=mx.n_row) || (this->n_col!=mx.n_col))
     { cerr << "Operator += invalid~! [E012]\n";
       return *this;
     }
// remove mx from this.
  for (i=0; i<n_dim; i++) this->cptr[i] -= mx.cptr[i];
  return *this;
} // end of operator-=
//==================================================
// friend function operator+
//==================================================
C_Matrix operator+(const C_Matrix& m1, const C_Matrix& m2)
{ int i;
  C_Matrix mt(m1);
// check validity of + operation.
  if ((mt.n_row!=m2.n_row) || (mt.n_col!=m2.n_col))
     { cerr << "Operator + invalid~! [E013]\n";
       return mt;
     }
// adding up two matrixes.
  for (i=0; i<mt.n_dim; i++) mt.cptr[i] += m2.cptr[i];
  return mt;
} // end of operator+
//==================================================
// friend function operator-
//==================================================
C_Matrix operator-(const C_Matrix& m1, const C_Matrix& m2)
{ int i;
  C_Matrix mt(m1);
// check validity of + operation.
  if ((mt.n_row!=m2.n_row) || (mt.n_col!=m2.n_col))
     { cerr << "Operator - invalid~! [E014]\n";
       return mt;
     }
// adding up two matrixes.
  for (i=0; i<mt.n_dim; i++) mt.cptr[i] -= m2.cptr[i];
  return mt;
} // end of operator-
//==================================================
// friend function operator*
//==================================================
C_Matrix operator*(const C_Matrix& m1, const C_Matrix& m2)
{ int i, j, k;
  COMPLX sum;
  C_Matrix mt(m1.n_row, m2.n_col);
// check validity of + operation.
  if ((m1.n_col!=m2.n_row))
     { cerr << "Operator * invalid~! [E015]\n";
       return mt;
     }
// hook up two matrixes.
  for (i=0; i<mt.n_row; i++)
      { for (j=0; j<mt.n_col; j++)
            { sum = ZERO;
              for (k=0; k<m1.n_col; k++) sum += ( m1.get(i,k)*m2.get(k,j) );
              mt.cell(i,j) = sum;
            }
      }
  return mt;
} // end of operator* a c-number
C_Matrix operator*(const COMPLX c2, const C_Matrix& m1)
{ int k, n;
  C_Matrix mt(m1);
  n = mt.ndim();
  for (k=0; k<n; k++) mt.cell(k) *= c2;
  return mt;
} // end of operator*
C_Matrix operator*(const C_Matrix& m1, const COMPLX c2)
{ int k, n;
  C_Matrix mt(m1);
  n = mt.ndim();
  for (k=0; k<n; k++) mt.cell(k) *= c2;
  return mt;
} // end of operator*
C_Matrix operator/(const C_Matrix& m1, const COMPLX c2)
{ int k, n;
  C_Matrix mt(m1);
  n = mt.ndim();
  for (k=0; k<n; k++) mt.cell(k) /= c2;
  return mt;
} // end of operator*
void C_Matrix::readfile(ifstream& fpt)
{  int i, r, c;
   if (this->n_dim > 0) delete [] this->cptr;
   fpt >> this->n_row >> this->n_col;
   this->n_dim = this->n_row * this->n_col;
   if ((this->cptr = new COMPLX[this->n_dim]) != NULL)
      {   for (i=0; i<this->n_dim; i++)
              { fpt >> r >> c >> this->cptr[r*n_col + c];
              }
          return;
      }
   else { cerr << "Memory fault in readfile~! [E016]\n";
          this->n_dim = this->n_col = this->n_row = 0;
          return;
        }
}
C_Matrix::C_Matrix(ifstream& fpt)
{  int i, r, c;
   this->PrintMode = DefaultPrintMode;
   fpt >> this->n_row >> this->n_col;
//   cout << "Chk1 " << n_row << "  " << n_col << endl;
   this->n_dim = this->n_row * this->n_col;
//   cout << "NDIM =" << n_dim << endl;
   this->cptr = new complex<double> [this->n_dim];
//   cout << "Chk memorail allocation\n";
 //  if (( cptr != NULL)
         for (i=0; i<this->n_dim; i++)
              { fpt >> r >> c;
                fpt >> cptr[r*n_col + c];
              }
          return;

 //  else { cerr << "Memory fault in C_Matrix~! [E017]\n";
 //         this->n_dim = this->n_col = this->n_row = 0;
 //         return;
 //       }
}
void C_Matrix::writefile(char*fname)
{  int i, j, c;
   ofstream fpt(fname);
   fpt << this->n_row << "  " << this->n_col <<"\n";
   for(i=0; i<n_row; i++)
      { for (j=0; j<n_col;j++)
            { fpt << i << "  " << j << "  " << this->get(i,j)<<"\n";
            }
      }
      fpt.close();
}
//
#endif
