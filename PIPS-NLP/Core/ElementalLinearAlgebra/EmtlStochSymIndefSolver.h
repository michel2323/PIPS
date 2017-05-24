/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#ifndef EMTLSTOCHINDEFSOLVER_H
#define EMTLSTOCHINDEFSOLVER_H

#include "DoubleLinearSolver.h"
#include "EmtlDenGenMatrix.h"
#include "EmtlDenSymMatrix.h"

class EmtlStochSymIndefSolver : public DoubleLinearSolver {
  protected:
  int m,n;
  int nr;
  int nx;
  bool perturb;
  DistMatrix<double,MC,MR> Acopy,rhscopy;
  public:
  EmtlDenGenMatrixHandle mat;
    
  EmtlStochSymIndefSolver( EmtlDenSymMatrix *mat, int nx);
  
  virtual void diagonalChanged( int idiag, int extent );
  virtual void matrixChanged();
  virtual void solve ( OoqpVector& vec );
  virtual ~EmtlStochSymIndefSolver();


};


#endif
