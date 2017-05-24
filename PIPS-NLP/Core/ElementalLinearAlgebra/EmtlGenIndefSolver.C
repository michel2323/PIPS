/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#include "EmtlGenIndefSolver.h"
#include "EmtlVector.h"
#include <cassert>
#include <sys/resource.h>
#include <iostream>
#include "elemental/lapack_internal.hpp"

EmtlGenIndefSolver::EmtlGenIndefSolver( EmtlDenGenMatrix *mat)
{
  SpReferTo(this->mat, mat);
  
  mat->getSize(m,n);
  assert(m==n);
  if (!mat->isNoop()) { 
    p = new DistMatrix<int,VC,Star>(m,1,mat->ctx.grid());
  } 
}

EmtlGenIndefSolver::EmtlGenIndefSolver( EmtlDenSymMatrix *mat_ )
{
  EmtlDenGenMatrix *mat = mat_->mat;
  SpReferTo(this->mat, mat);
  
  mat->getSize(m,n);
  assert(m==n);
  if (!mat->isNoop()) { 
    p = new DistMatrix<int,VC,Star>(m,1,mat->ctx.grid());
  } 
}

void EmtlGenIndefSolver::diagonalChanged( int idiag, int extent )
{
  this->matrixChanged();
}


void EmtlGenIndefSolver::matrixChanged()
{
  if (mat->isNoop()) return;
  DistMatrix<double,MC,MR> &A = *mat->A;

  lapack::LU(A,*p);
  needtocompose = true;
  // don't time compose pivots
  /*
  DistMatrix<int,Star,Star> p_Star_Star(mat->ctx.grid());
  p_Star_Star = *p;
  lapack::internal::ComposePivots( p_Star_Star, image, preimage, 0 );
  */
}

void EmtlGenIndefSolver::solve( OoqpVector& v )
{
  if (mat->isNoop()) return;
  if (needtocompose) {
    DistMatrix<int,Star,Star> p_Star_Star(mat->ctx.grid());
    p_Star_Star = *p;
    lapack::internal::ComposePivots( p_Star_Star, image, preimage, 0 );
    needtocompose = false;
  }
  
  EmtlVector &vec = dynamic_cast<EmtlVector&>(v);
  DistMatrix<double,MC,MR> &A = *mat->A;
  DistMatrix<double,MC,MR> &x = *vec.A;

  /*
  rusage before_solve;
  getrusage( RUSAGE_SELF, &before_solve );
  */
  
  lapack::internal::ApplyRowPivots( x, image, preimage, 0 );
  blas::Trsv( Lower, Normal, Unit, A, x );
  blas::Trsv( Upper, Normal, NonUnit, A, x );
  
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

EmtlGenIndefSolver::~EmtlGenIndefSolver()
{
  if (!mat->isNoop()) {
    delete p;
  }
}
