/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#include "EmtlLinearAlgebraPackage.h"
#include "EmtlDenSymMatrix.h"
#include "EmtlDenGenMatrix.h"

#include "EmtlVector.h"

EmtlLinearAlgebraPackage::EmtlLinearAlgebraPackage(const EmtlContext &ctx_) :
ctx(ctx_)
{
}

SymMatrix * EmtlLinearAlgebraPackage::newSymMatrix( int size, int nnz )
{
  return new EmtlDenSymMatrix(size, ctx);
}

GenMatrix * EmtlLinearAlgebraPackage::newGenMatrix( int m, int n,
							 int nnz )
{
  return new EmtlDenGenMatrix( m, n, ctx );
}

OoqpVector * EmtlLinearAlgebraPackage::newVector( int n )
{
  return new EmtlVector( n, ctx );
}

void EmtlLinearAlgebraPackage::whatami( char type[] )
{
  char type_[] = "EmtlLinearAlgebraPackage";

  int i = 0;
  
  type[0] = type_[0];
  while( type[i] != 0 && i < 31 ) {
    ++i;
    type[i] = type_[i];
  }
}

