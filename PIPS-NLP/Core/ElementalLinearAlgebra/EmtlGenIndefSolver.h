/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#ifndef EMTLGENINDEFSOLVER_H
#define EMTLGENINDEFSOLVER_H

#include "DoubleLinearSolver.h"
#include "EmtlDenGenMatrix.h"
#include "EmtlDenSymMatrix.h"

class EmtlGenIndefSolver : public DoubleLinearSolver {
  protected:
  // DistMatrix<int,VC,STAR> *p;
  DistPermutation p;
  DistMatrix<double,MD,STAR> dSub;
  DistMatrix<double,MC,MR> d;

  vector<int> image, preimage;
  int m,n;
  int nr;
  public:
  EmtlDenGenMatrixHandle mat;
    
  EmtlGenIndefSolver( EmtlDenGenMatrix *mat);
  EmtlGenIndefSolver( EmtlDenSymMatrix *mat);
  
  virtual void diagonalChanged( int idiag, int extent );
  virtual int matrixChanged();
  virtual void solve ( OoqpVector& vec );
  virtual ~EmtlGenIndefSolver();


};


#endif
