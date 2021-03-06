/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/   

#include "sFactoryAug.h"

#include "sData.h"

#include "StochTree.h"
#include "StochInputTree.h"

#include "sLinsysRootAug.h"



sFactoryAug::sFactoryAug( StochInputTree* inputTree)
  : sFactory(inputTree)
{ };

sFactoryAug::sFactoryAug( stochasticInput& in, MPI_Comm comm)
  : sFactory(in,comm)
{ }


sFactoryAug::sFactoryAug( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : sFactory(nx, my, mz, nnzQ, nnzA, nnzC)
{ };

sFactoryAug::sFactoryAug()
{ };

sFactoryAug::~sFactoryAug()
{ };


sLinsysRoot* sFactoryAug::newLinsysRoot()
{
#ifdef DEBUG
  printf("sFactoryAug::newLinsysRoot 1\n");
#endif
  return new sLinsysRootAug(this, data);
}

sLinsysRoot* 
sFactoryAug::newLinsysRoot(sData* prob,
			   OoqpVector* dd,OoqpVector* dq,
			   OoqpVector* nomegaInv, OoqpVector* rhs, OoqpVector* additiveDiag)
{
#ifdef DEBUG
  printf("sFactoryAug::newLinsysRoot 2\n");
#endif
  return new sLinsysRootAug(this, prob, dd, dq, nomegaInv, rhs, additiveDiag);
}


LinearSystem* sFactoryAug::makeLinsys( Data * prob_in )
{  
#ifdef DEBUG
  printf("sFactoryAug::newLinsysRoot 3\n");
#endif
  linsys = NULL;
  linsys = newLinsysRoot();

  assert(linsys);
  return linsys; 
}


