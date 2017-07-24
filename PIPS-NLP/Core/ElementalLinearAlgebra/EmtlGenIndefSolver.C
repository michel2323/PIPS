/* PIPS
   Authors: Miles Lubin, Michel Schanen
   See license and copyright information in the documentation */

#include "EmtlGenIndefSolver.h"
#include "EmtlVector.h"
#include <cassert>
#include <sys/resource.h>
#include <iostream>

EmtlGenIndefSolver::EmtlGenIndefSolver( EmtlDenGenMatrix *mat)
{
  SpReferTo(this->mat, mat);
  
  mat->getSize(m,n);
  assert(m==n);
// The diagonal of the LDL fact. has to be on the same grid
  DistMatrix<double,MC,MR> &A = *mat->A;
  dSub.SetGrid(A.Grid());
}

EmtlGenIndefSolver::EmtlGenIndefSolver( EmtlDenSymMatrix *mat_ )
{
  EmtlDenGenMatrix *mat = mat_->mat;
  SpReferTo(this->mat, mat);
  
  mat->getSize(m,n);
  assert(m==n);
// The diagonal of the LDL fact. has to be on the same grid
  DistMatrix<double,MC,MR> &A = *mat->A;
  dSub.SetGrid(A.Grid());
}

void EmtlGenIndefSolver::diagonalChanged( int idiag, int extent )
{
  this->matrixChanged();
}


int EmtlGenIndefSolver::matrixChanged()
{
  if (mat->isNoop()) return 0;
  DistMatrix<double,MC,MR> &A = *mat->A;
#ifdef DEBUG
  // printf("[EmtlGenIndefSolver]\n");
  // El::Display(A);
#endif
  double before=mat->abmaxnorm();
  El::LDL( A, dSub, p, false, El::BUNCH_KAUFMAN_A);
  double after=mat->abmaxnorm();
  GetDiagonal(A,d);
  InertiaType inertia=El::ldl::Inertia(d,dSub);
  negEigVal=inertia.numNegative;
#ifdef DEBUG
  printf("Inertia: %d %lf %lf %lf\n", negEigVal, before, after, El::FrobeniusNorm(A));
#endif
  return negEigVal;
}

void EmtlGenIndefSolver::solve( OoqpVector& v )
{
  if (mat->isNoop()) return;
  
  EmtlVector &vec = dynamic_cast<EmtlVector&>(v);
  DistMatrix<double,MC,MR> &A = *mat->A;
  DistMatrix<double,MC,MR> &x = *vec.A;
  #ifdef DEBUG
  // Display(x);
  printf("[DeSymIndefSolver::solve 1] norm: %f %lf\n", El::FrobeniusNorm(x), El::FrobeniusNorm(A));
  #endif
  El::ldl::SolveAfter(A, dSub, p, x, false);
  #ifdef DEBUG
  printf("[DeSymIndefSolver::solve 2] norm: %f %f\n", mat->abmaxnorm(), vec.onenorm());
  #endif
}

EmtlGenIndefSolver::~EmtlGenIndefSolver()
{
}
