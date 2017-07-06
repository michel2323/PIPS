/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#include "EmtlDenGenMatrix.h"
#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdio>
#include "SimpleVector.h"
#include "EmtlVector.h"
#include "SparseGenMatrix.h"
#include "SparseSymMatrix.h"

#include "DoubleMatrixTypes.h"

EmtlDenGenMatrix::EmtlDenGenMatrix( int m, int n, const EmtlContext &ctx_) :
  ctx(ctx_)
{
  this->m = m;
  this->n = n;

  if (!ctx.noop()) {
    A = new DistMatrix<double,MC,MR>(m, n, ctx.grid());

    nr = A->LocalHeight();
    nc = A->LocalWidth();
    data = A->Matrix().Buffer();
    assert(A->Matrix().LDim() == nr);
  } else {
    data = 0;
    nr = nc = 0;
  }

#ifdef TIMING
  if (ctx.noop()) {
    printf("EMTL NOOP PROC %d\n",ctx.mype());
  } else {
    printf("EMTL LOCSIZE PROC %d NR %d NC %d TOTAL %d MYROW %d MYCOL %d\n",
      ctx.mype(),nr,nc,nr*nc,ctx.myrow(),ctx.mycol());
  }
#endif
}

EmtlDenGenMatrix::~EmtlDenGenMatrix()
{
  if (!ctx.noop()) {
    delete A;
  }
}

int EmtlDenGenMatrix::isKindOf( int type )
{
  return (type & kEmtlDenGenMatrix);
}

void EmtlDenGenMatrix::getDiagonal( OoqpVector& vec )
{
  this->fromGetDiagonal(0, vec);
}

void EmtlDenGenMatrix::setToDiagonal( OoqpVector& vec )
{
  memset(data, 0, nr*nc*sizeof(double));
  this->atPutDiagonal(0, vec);
}

void EmtlDenGenMatrix::atPutDense( int  row,     int  col,
				   double *  A , int  lda,
				   int  rowExtent ,
				   int  colExtent  )
{
  assert( "Not implemented" && 0 );
}

void EmtlDenGenMatrix::fromGetDense( int  row,    int  col,
				     double *  A, int  lda,
				     int  rowExtent,
				     int  colExtent  )
{
  assert( "Not implemented" && 0 );
}

void EmtlDenGenMatrix::putSparseTriple( int irow[], int len,
					    int jcol[], double A[], 
					    int& info )
{
  assert( "Not implemented" && 0 );
}

void EmtlDenGenMatrix::fromGetSpRow( int row, int col,
					 double A[], int lenA,
					 int jcolA[], int& nnz,
					 int colExtent, int& info )
{
  assert( "Not implemented" && 0 );
}

void EmtlDenGenMatrix::atPutSpRow( int col, double A[],
				       int lenA, int jcolA[], int& info )
{
  assert( "Not implemented" && 0 );
}

void EmtlDenGenMatrix::getSize( int& m, int& n )
{
  m = this->m;
  n = this->n;
}

void EmtlDenGenMatrix::getSize( long long& m, long long& n )
{
  m = this->m;
  n = this->n;
}

void EmtlDenGenMatrix::atPutSubmatrix( int destRow, int destCol,
					   DoubleMatrix& M,
					   int srcRow, int srcCol,
					   int rowExtent, int colExtent )
{
  this->atPutSubmatrix(destRow, destCol, M, srcRow, srcCol, rowExtent, colExtent, false);
}

void EmtlDenGenMatrix::atPutSubmatrix( int destRow, int destCol,
    DoubleMatrix& mat_,
    int srcRow, int srcCol,
    int rowExtent, int colExtent, bool transpose )
{
  //if (noop) return;

  int *krowM, *jcolM;
  int nrows, ncols, row, col;
  double *M;
  mat_.getSize(nrows, ncols);
  assert((srcRow == 0) && (srcCol == 0));
  assert((!transpose && (nrows == rowExtent) && (ncols == colExtent)) ||
      (transpose && (nrows == colExtent) && (ncols == rowExtent)));
  assert(destRow >= 0 && destRow+rowExtent <= m);
  assert(destCol >= 0 && destCol+colExtent <= n);
  if (mat_.isKindOf(kSparseGenMatrix)) {
    SparseGenMatrix & mat = dynamic_cast<SparseGenMatrix&> (mat_);
    krowM = mat.krowM();
    jcolM = mat.jcolM();
    M = mat.M();
    for (int i = 0; i < nrows; i++) {
      for (int k = krowM[i]; k < krowM[i+1]; k++) {
        if (!transpose) {
          row = i; col = jcolM[k];
        } else {
          row = jcolM[k]; col = i;
        }
        putVal(destRow+row, destCol+col, M[k]);

      }
    }

  } else if (mat_.isKindOf(kSparseSymMatrix)) {
    SparseSymMatrix & mat = dynamic_cast<SparseSymMatrix&> (mat_);
    krowM = mat.krowM();
    jcolM = mat.jcolM();
    M = mat.M();
    for (int i = 0; i < nrows; i++) {
      for (int k = krowM[i]; k < krowM[i+1]; k++) {
        row = i; col = jcolM[k];
        putVal(destRow+row,destCol+col, M[k]);
        row = jcolM[k]; col = i;
        putVal(destRow+row,destCol+col, M[k]);

      }

    }

  } else {
    assert( "Not implemented" && 0 );
  }
}


void EmtlDenGenMatrix::mult ( double beta,  OoqpVector& y_in,
				  double alpha, OoqpVector& x_in )
{
assert( "Not implemented" && 0 );

}


void EmtlDenGenMatrix::transMult ( double beta,  OoqpVector& y_in,
				       double alpha, OoqpVector& x_in )
{

assert( "Not implemented" && 0 );
}




double EmtlDenGenMatrix::abmaxnorm()
{
#ifdef DEBUG
  printf("[EmtlDenGenMatrix::abmaxnorm] Wrong abmaxnorm\n");
#endif
  return El::FrobeniusNorm(*A);
  // return El::MaxNorm(*A);
}



void EmtlDenGenMatrix::writeToStream(ostream& out) const
{
  assert( "Not implemented" && 0 );
}


void EmtlDenGenMatrix::randomize(double /* alpha */, double /* beta */,
					double * /* seed */)
{
  assert( "Not implemented" && 0 );
}


void EmtlDenGenMatrix::atPutDiagonal( int idiag, OoqpVector& vec )
{
  assert(vec.length() <= m - idiag);
  
  //if (noop) return;
  
  if (vec.isKindOf(kSimpleVector)) {
    SimpleVector &v = dynamic_cast<SimpleVector&>(vec);
    for (int i = idiag; i - idiag < v.length(); i++) {
      putVal(i,i, v[i-idiag]);
    }
    
  } else {
    assert( "Not implemented" && 0 );
  }
  
}

void EmtlDenGenMatrix::fromGetDiagonal( int idiag, OoqpVector& vec )
{
  assert( "Not implemented" && 0 );
  
}

void EmtlDenGenMatrix::atPutZeros( int row, int col,
  			   int rowExtent, int colExtent )
{
  assert(row >= 0 && row+rowExtent <= m);
  assert(col >= 0 && col+colExtent <= n);
  
  for (int j = col; j < col+colExtent; j++) {
    for (int i = row; i < row+rowExtent; i++) {
      putVal(i,j, 0.);
    }
  }
  
}


void EmtlDenGenMatrix::SymmetricScale( OoqpVector& vec )
{
  assert( "Not implemented" && 0 );
}

void EmtlDenGenMatrix::ColumnScale( OoqpVector& vec )
{
  assert( "Not implemented" && 0 );
}

void EmtlDenGenMatrix::RowScale( OoqpVector& vec )
{
  assert( "Not implemented" && 0 );
}

void EmtlDenGenMatrix::scalarMult( double num )
{
  if (num != 0.) {
    for (int i = 0; i < nr*nc; i++) {
      data[i] *= num;
    }
  } else {
    memset(data, 0, nr*nc*sizeof(double));
  }
}

void EmtlDenGenMatrix::matTransDMultMat(OoqpVector&, SymMatrix**)
{
  assert( "Not implemented" && 0 );
}

void EmtlDenGenMatrix::matTransDinvMultMat(OoqpVector&, SymMatrix**)
{
  assert( "Not implemented" && 0 );
}
