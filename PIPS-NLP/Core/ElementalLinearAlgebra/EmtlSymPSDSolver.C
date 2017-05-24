/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#include "EmtlSymPSDSolver.h"
#include "EmtlVector.h"
#include <cassert>
#include <sys/resource.h>
#include <iostream>


EmtlSymPSDSolver::EmtlSymPSDSolver( EmtlDenSymMatrix *mat_ )
{
  // three different variables named mat here
  EmtlDenGenMatrix *mat = mat_->mat;
  SpReferTo(this->mat, mat);
  mat->getSize(m,n);
  assert(m==n);
  
}

void EmtlSymPSDSolver::diagonalChanged( int idiag, int extent )
{
  this->matrixChanged();
}


void EmtlSymPSDSolver::matrixChanged()
{
  if (mat->isNoop()) return;
  
  DistMatrix<double,MC,MR> &A = *mat->A;
  lapack::Chol(Lower, A);
}

void EmtlSymPSDSolver::solve( OoqpVector& v )
{
  if (mat->isNoop()) return;
  EmtlVector &vec = dynamic_cast<EmtlVector&>(v);
  DistMatrix<double,MC,MR> &A = *mat->A;
  DistMatrix<double,MC,MR> &x = *vec.A;

  /*
  rusage before_solve;
  getrusage( RUSAGE_SELF, &before_solve );
  */

  blas::Trsv(Lower, Normal, NonUnit, A, x);
  blas::Trsv(Lower, Transpose, NonUnit, A, x);
  
  /*
  rusage  after_solve;
  getrusage( RUSAGE_SELF, &after_solve );
  
  double solve_time =
	  (after_solve.ru_utime.tv_sec - before_solve.ru_utime.tv_sec)
	  + (after_solve.ru_utime.tv_usec - before_solve.ru_utime.tv_usec)
	  / 1000000.0;
    
  if (mat->cinfo.mype == 0) {
    cout << "solve time: " << solve_time << endl;
  }*/
}

EmtlSymPSDSolver::~EmtlSymPSDSolver()
{
}
