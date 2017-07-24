/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */


#include "sFactoryAugEmtl.h"

#include "sData.h"

#include "StochTree.h"
#include "StochInputTree.h"

#include "sLinsysRootAugEmtl.h"
#include "StochResourcePlanner.h"

sFactoryAugEmtl::sFactoryAugEmtl( StochInputTree* inputTree)
  : sFactory(inputTree)
{ 
  MPI_Comm world = tree->commWrkrs;
  int emtlprocs = StochResourcePlanner::noScaProcesses;
  ctx = new EmtlContext(world, emtlprocs);

  if (ctx->mype() == 0) {
    int nb = Blocksize();
#ifdef TIMING
    printf("EMTL NPROCS %d GRID %d %d BLOCKSIZE %d\n",
      ctx->nprow()*ctx->npcol(), ctx->nprow(), ctx->npcol(), nb);
#else
    printf("Using %d processes for Elemental, %d by %d grid\n, blocksize %d", 
      ctx->nprow()*ctx->npcol(), ctx->nprow(), ctx->npcol(), nb);
#endif
  }
  
  
};

sFactoryAugEmtl::sFactoryAugEmtl( stochasticInput& in, MPI_Comm comm)
  : sFactory(in,comm)
{ 
  // MPI_Comm world = tree->commWrkrs
  int emtlprocs;
  MPI_Comm_size(comm, &emtlprocs);
  // int emtlprocs = StochResourcePlanner::noScaProcesses;
#ifdef DEBUG
  printf("[sFactoryAugEmtl::sFactoryAugEmtl] comm: %d\n", emtlprocs);
#endif
  ctx = new EmtlContext(comm, emtlprocs);
  

  if (ctx->mype() == 0) {
    int nb = Blocksize();
#ifdef TIMING
    printf("EMTL NPROCS %d GRID %d %d BLOCKSIZE %d\n",
      ctx->nprow()*ctx->npcol(), ctx->nprow(), ctx->npcol(), nb);
#else
    printf("Using %d processes for Elemental, %d by %d grid\n, blocksize %d", 
      ctx->nprow()*ctx->npcol(), ctx->nprow(), ctx->npcol(), nb);
#endif
  }
};

sFactoryAugEmtl::sFactoryAugEmtl( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : sFactory(nx, my, mz, nnzQ, nnzA, nnzC)
{ };

sFactoryAugEmtl::sFactoryAugEmtl()
{ 
  assert("Don't call this" && 0);
};

sFactoryAugEmtl::~sFactoryAugEmtl()
{ 
  delete ctx;
};


sLinsysRoot* sFactoryAugEmtl::newLinsysRoot()
{
#ifdef DEBUG
  printf("sFactoryAugEmtl::newLinsysRoot 1\n");
#endif
  return new sLinsysRootAugEmtl(this, data, *ctx);
}

sLinsysRoot* 
sFactoryAugEmtl::newLinsysRoot(sData* prob,
			   OoqpVector* dd,OoqpVector* dq,
			   OoqpVector* nomegaInv, OoqpVector* rhs, OoqpVector* additiveDiag)
{
#ifdef DEBUG
  printf("sFactoryAugEmtl::newLinsysRoot 2\n");
#endif
  return new sLinsysRootAugEmtl(this, prob, dd, dq, nomegaInv, rhs, additiveDiag, *ctx);
}
