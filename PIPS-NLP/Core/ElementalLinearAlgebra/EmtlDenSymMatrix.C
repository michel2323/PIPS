/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#include "EmtlDenSymMatrix.h"
#include "EmtlDenGenMatrix.h"
#include <cassert>
#include <cstdio>
#include <cmath>
#include "EmtlVector.h"
#include "DoubleMatrixTypes.h"

int EmtlDenSymMatrix::isKindOf( int type )
{
   return (type & kEmtlDenSymMatrix);

}

EmtlDenSymMatrix::EmtlDenSymMatrix( int size, const EmtlContext &ctx_ ) :
  ctx(ctx_)
{
  m = size;
  mat = EmtlDenGenMatrixHandle(new EmtlDenGenMatrix(size,size,ctx));
}
/*
EmtlDenSymMatrix::EmtlDenSymMatrix( int lm, int ln, int size, int nnz )
{

}
*/



void EmtlDenSymMatrix::getDiagonal( OoqpVector& vec )
{
  mat->getDiagonal(vec);
}

void EmtlDenSymMatrix::setToDiagonal( OoqpVector& vec )
{
  mat->setToDiagonal(vec);
}

void EmtlDenSymMatrix::atPutDense( int /* row */, int /* col */,
				   double * /* A */, int /* lda */,
				   int /* rowExtent */,
				   int /* colExtent */ )
{
  assert( "Not implemented" && 0 );
}

void EmtlDenSymMatrix::fromGetDense( int row, int col,
				     double * A, int lda,
				     int rowExtent,
				     int colExtent )
{
  assert( "Not implemented" && 0 );
}

void EmtlDenSymMatrix::randomizePSD(double * /* seed */)
{
  assert( "Not implemented" && 0 );
}
  
void EmtlDenSymMatrix::symAtPutSpRow( int row, double A[],
					  int lenA, int jcolA[],
					  int& info )
{
  int ierr;
 
}

void EmtlDenSymMatrix::fsymAtPutSpRow( int row, double A[],
					   int lenA, int jcolA[],
					   int& info )
{

}

void EmtlDenSymMatrix::fromGetSpRow( int row, int col,
					 double A[], int lenA,
					 int jcolA[], int& nnz,
					 int colExtent, int& info )
{

}


void EmtlDenSymMatrix::symAtPutSubmatrix( int destRow, int destCol,
					      DoubleMatrix& M,
					      int srcRow, int srcCol,
					      int rowExtent, int colExtent )
{
  //int a,b;
  //M.getSize(a,b);
  //printf("%d %d\n", a, b); 
  assert(srcRow == 0 && srcCol == 0);
  if (destRow == destCol) {
    assert(M.isKindOf(kSymMatrix));
    mat->atPutSubmatrix(destRow,destCol,M,srcRow,srcCol,rowExtent,colExtent);
  } else {
    //printf("%d %d %d %d %d\n", m, destRow, destCol, rowExtent, colExtent);

    //printf("%d %d\n", a, b); 
    assert(!(( (destCol < destRow) && (destCol+colExtent > destRow) ) ||
	   ( (destRow < destCol) && (destRow+rowExtent > destCol) )));
    mat->atPutSubmatrix(destRow,destCol,M,srcRow,srcCol,rowExtent,colExtent,false);
    mat->atPutSubmatrix(destCol,destRow,M,srcCol,srcRow,colExtent,rowExtent,true);
  }
}

void EmtlDenSymMatrix::putSparseTriple( int irow[], int len,
					    int jcol[], double A[], 
					    int& info )
{

}

// Pass these to storage
void EmtlDenSymMatrix::getSize( int& m, int& n )
{
  mat->getSize(m,n);
}

int EmtlDenSymMatrix::size()
{
  return m;
}

void EmtlDenSymMatrix::atPutZeros( int row, int col,
				       int rowExtent, int colExtent )
{
    assert( "Not implemented" && 0 );
}
  
void EmtlDenSymMatrix::transMult ( double beta,  OoqpVector& y_in,
				       double alpha, OoqpVector& x_in )
{
  mat->transMult(beta, y_in, alpha, x_in);
}

void EmtlDenSymMatrix::mult ( double beta,  OoqpVector& y_in,
				  double alpha, OoqpVector& x_in )
{
  mat->mult(beta, y_in, alpha, x_in);
}

double EmtlDenSymMatrix::abmaxnorm()
{ 
  return mat->abmaxnorm();

}

// print lower triangle, naively, in row-major order to match sparse output
// debugging only
void EmtlDenSymMatrix::writeToStream(ostream& out) const
{
  MPI_Request req;
  double buf;
  
  for (int i = 0; i < m; i++) {
    for (int j = 0; j <= i; j++) {
      int proc = ctx.index2proc(i,j);
      if (ctx.mype() == proc) {
        int lrow, lcol;
        ctx.global2local(i,j,lrow,lcol);
        MPI_Isend(&mat->data[lrow+mat->getNR()*lcol],1,MPI_DOUBLE,0,0,ctx.comm(),&req);
        MPI_Request_free(&req);
      }
      if (ctx.mype() == 0) {
        MPI_Recv(&buf,1,MPI_DOUBLE,proc,0,ctx.comm(),MPI_STATUS_IGNORE);
        if (buf != 0.) {
          out << i << "\t" << j << "\t" << buf << endl;
        }
      }
    }
  }
  // do this instead of keeping track of requests
  MPI_Barrier(ctx.comm());
}


void EmtlDenSymMatrix::atPutDiagonal( int idiag, OoqpVector& v )
{
  mat->atPutDiagonal(idiag, v);
}

void EmtlDenSymMatrix::fromGetDiagonal( int idiag, OoqpVector& v )
{
  mat->atPutDiagonal(idiag, v);
}


void EmtlDenSymMatrix::SymmetricScale( OoqpVector& vec )
{
  assert( "Not implemented" && 0 );
}

void EmtlDenSymMatrix::ColumnScale( OoqpVector& vec )
{
  assert( "Not implemented" && 0 );
}

void EmtlDenSymMatrix::RowScale( OoqpVector& vec )
{
  assert( "Not implemented" && 0 );
}

void EmtlDenSymMatrix::scalarMult( double num )
{
  mat->scalarMult(num);
}

void EmtlDenSymMatrix::symAtPutZeros( int row, int col,
  			   int rowExtent, int colExtent )
{
  if (row == col) {
    // just check the block is symmetric
    assert(rowExtent == colExtent);
    mat->atPutZeros(row, col, rowExtent, colExtent);
  } else {
    // zero out the transpose also
    // if block crosses the diagonal, this probably isn't intended
    assert(!(( (col < row) && (col+colExtent > row) ) ||
	   ( (row < col) && (row+rowExtent > col) )));
    mat->atPutZeros(row, col, rowExtent, colExtent);
    mat->atPutZeros(col, row, colExtent, rowExtent);
  }
  
}


bool EmtlDenSymMatrix::isNoop() const 
{ 
  return mat->isNoop(); 
}

/*int EmtlDenSymMatrix::getNR() const 
{ 
  return mat->getNR(); 
}*/

// debugging only
double EmtlDenSymMatrix::getVal(int row, int col) const
{
  double v1, v2;
  int proc, lrow, lcol;
  proc = ctx.index2proc(row,col);
  if (ctx.mype() == proc) {
    ctx.global2local(row, col, lrow, lcol);
    v1 = mat->data[lrow+getNR()*lcol];
  }
  MPI_Bcast(&v1, 1, MPI_DOUBLE, proc, ctx.comm());
  proc = ctx.index2proc(col,row);
  if (ctx.mype() == proc) {
    ctx.global2local(col, row, lrow, lcol);
    v2 = mat->data[lrow+getNR()*lcol];
  }
  MPI_Bcast(&v2, 1, MPI_DOUBLE, proc, ctx.comm());
  assert( fabs(v1 - v2) < 10E-6 );
  return v1;
}
