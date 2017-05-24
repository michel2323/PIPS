/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#ifndef EMTLLINEARALGEBRA
#define EMTLLINEARALGEBRA

#include "LinearAlgebraPackage.h"
#include "EmtlContext.h"

class EmtlLinearAlgebraPackage : public LinearAlgebraPackage {
protected:
  const EmtlContext &ctx;
public:
  EmtlLinearAlgebraPackage(const EmtlContext &ctx_);
  virtual ~EmtlLinearAlgebraPackage() {};
  virtual SymMatrix * newSymMatrix( int size, int nnz );
  virtual GenMatrix * newGenMatrix( int m, int n, int nnz );
  virtual OoqpVector * newVector( int n );
  virtual void whatami( char type[32] );
};

#endif
