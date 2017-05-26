

#ifndef EMTLSYMPSDSOLVER_H
#define EMTLSYMPSDSOLVER_H

#include "DoubleLinearSolver.h"
#include "EmtlDenSymMatrix.h"

class EmtlSymPSDSolver : public DoubleLinearSolver {
  protected:
  int m,n;
  public:
  EmtlDenGenMatrixHandle mat;
    
  EmtlSymPSDSolver( EmtlDenSymMatrix *mat);
  
  virtual void diagonalChanged( int idiag, int extent );
  virtual int matrixChanged();
  virtual void solve ( OoqpVector& vec );
  virtual ~EmtlSymPSDSolver();


};


#endif
