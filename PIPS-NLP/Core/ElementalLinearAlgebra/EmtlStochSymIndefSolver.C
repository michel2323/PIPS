/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#include "EmtlStochSymIndefSolver.h"
#include "EmtlVector.h"
#include <cassert>
#include <sys/resource.h>
#include <iostream>

// defined in MehrotraStochSolver.C
extern double gmu;

EmtlStochSymIndefSolver::EmtlStochSymIndefSolver( EmtlDenSymMatrix *mat_, int nx_ ) :
  nx(nx_), Acopy(mat_->mat->A->Grid()), rhscopy(mat_->mat->A->Grid())
{
  EmtlDenGenMatrix *mat = mat_->mat;
  SpReferTo(this->mat, mat);
  mat->getSize(m,n);
  assert(m==n);
}

void EmtlStochSymIndefSolver::diagonalChanged( int idiag, int extent )
{
  this->matrixChanged();
}

/*
  [ Q * ] becomes [ L     *  ] 
  [ A 0 ]         [ AL^-T ~L ]
  Q = LL^T
  AQ^(-1)A^T = (~L)(~L)^T
*/

void EmtlStochSymIndefSolver::matrixChanged()
{
  if (mat->isNoop()) return;
  
  DistMatrix<double,MC,MR> &A = *mat->A;
 
  perturb = (gmu < 1e-5);
  Acopy = A;
  if (perturb) {
    // perturb the diagonal
    for (int i = 0; i < nx; i++) {
      mat->addAt(i,i,1e-7);
    }
  }
  
  const Grid &g = A.Grid();
  DistMatrix<double,MC,MR> ATL(g), ATR(g), ABL(g), ABR(g); 
  // get references to blocks
  PartitionDownDiagonal(A, 
                        ATL, ATR, 
                        ABL, ABR,
                        nx);
  // replace Q with L
  lapack::Chol(Lower, ATL);
  // replace A with AL^-T
  blas::Trsm(Right, Lower, Transpose, NonUnit, 1., ATL, ABL);
  // replace bottom right with (AL^-T)(AL^-T)^T = AQ^(-1)A^T
  blas::Syrk(Lower, Normal, 1., ABL, 0., ABR);
  // factorize bottom right into ~L
  lapack::Chol(Lower, ABR);
 
}

// separate solve routine so member solve() can manage iterative refinement
static void solveLDLT(DistMatrix<double,MC,MR> &A, DistMatrix<double,MC,MR> &x, int nx)
{
  const Grid &g = A.Grid(); 
 
  // L solve 
  blas::Trsv(Lower, Normal, NonUnit, A, x);

  // D solve
  // D = [ I  0 ]
  //     [ 0 -I ]
  // so we just negate elements of the RHS after nx
  DistMatrix<double,MC,MR> xT(g), xB(g);
  PartitionDown(x, 
                xT, 
                xB, 
                nx);
  blas::Scal(-1., xB);

  // L^t solve
  blas::Trsv(Lower, Transpose, NonUnit, A, x);
}

// x solution to Ax=b
// resid is returned in r, nothing else is modified
static void calcResid(DistMatrix<double,MC,MR> &A,
  DistMatrix<double,MC,MR> &x, DistMatrix<double,MC,MR> &b, 
  DistMatrix<double,MC,MR> &r, bool print)
{
  r = b;
  double n1 = lapack::OneNorm(b), n2 = lapack::InfinityNorm(b);

  blas::Symv(Lower, -1., A, x, 1., r);

  double n1_ = lapack::OneNorm(r), n2_ = lapack::InfinityNorm(r);
  if (print) printf("SC solve rel. ||Resid||_1 = %e, ||Resid||_inf = %e\n", n1_/n1,n2_/n2);

}

// we do have separate Lsolve, Dsolve, Ltsolve parts, but put them
// here so we don't need to modify Linsys which expects everything to happen
// at Dsolve (manages distributing right hand side)
void EmtlStochSymIndefSolver::solve( OoqpVector& v )
{
  if (mat->isNoop()) return;
  EmtlVector &vec = dynamic_cast<EmtlVector&>(v);
  DistMatrix<double,MC,MR> &A = *mat->A;
  DistMatrix<double,MC,MR> &x = *vec.A;
  const Grid &g = A.Grid(); 

  rhscopy = x;

  solveLDLT(A,x,nx);

  //if (!perturb) return;


  bool pe0 = mat->ctx.mype()==0;

  DistMatrix<double,MC,MR> r(g);
  calcResid(Acopy,x,rhscopy,r,pe0);

  // perform a step of iterative refinement
  if (perturb) {
    solveLDLT(A,r,nx);
    blas::Axpy(1.,r,x);
    if (pe0) printf("After refinement: ");
    calcResid(Acopy,x,rhscopy,r,pe0);
  }

  

}




EmtlStochSymIndefSolver::~EmtlStochSymIndefSolver()
{
 
}
