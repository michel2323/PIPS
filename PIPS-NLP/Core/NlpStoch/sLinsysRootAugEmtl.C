/* PIPS
   Authors: Miles Lubin and Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysRootAugEmtl.h"
#include "sData.h"
#include "EmtlDenSymMatrix.h"
#include "EmtlGenIndefSolver.h"
#include "EmtlSymPSDSolver.h"
#include "EmtlVector.h"

#ifdef TIMING
#include "../../global_var.h"
#include "../PIPS-NLP/Core/Utilities/PerfMetrics.h"
#endif

#ifdef STOCH_TESTING
extern double g_iterNumber;
extern double g_scenNum;
#endif

extern int gInnerSCsolve;
extern int gOuterSolve;

sLinsysRootAugEmtl::sLinsysRootAugEmtl(sFactory * factory_, sData * prob_, const EmtlContext &ctx_)
  : sLinsysRoot(factory_, prob_), ctx(ctx_), CtDC(NULL)
{ 
  prob_->getLocalSizes(locnx, locmy, locmz);
  kkt = createKKT(prob_);
  solver = this->createSolver(prob_, kkt);
  redRhs = new SimpleVector(locnx+locmz+locmy+locmz);
  assert(gOuterSolve>=3);
#ifdef DEBUG
  printf("[sLinsysRootAugEmtl 1] sizes: %d %d %d\n",locnx,locmy,locmz);
#endif
  emtlRhs = new EmtlVector(locnx+locmz+locmy+locmz, ctx);
  iAmDistrib = 1;
};

sLinsysRootAugEmtl::sLinsysRootAugEmtl(sFactory* factory_,
			       sData* prob_,
			       OoqpVector* dd_, 
			       OoqpVector* dq_,
			       OoqpVector* nomegaInv_,
			       OoqpVector* rhs_,
             OoqpVector* additiveDiag_,
			       const EmtlContext &ctx_)
  : sLinsysRoot(factory_, prob_, dd_, dq_, nomegaInv_, additiveDiag_, rhs_), 
    ctx(ctx_), CtDC(NULL)
{ 
  prob_->getLocalSizes(locnx, locmy, locmz);
  kkt = createKKT(prob_);
  solver = this->createSolver(prob_, kkt);
#ifdef DEBUG
  printf("[sLinsysRootAugEmtl 2] sizes: %d %d %d\n",locnx,locmy,locmz);
#endif
  redRhs = new SimpleVector(locnx+locmz+locmy+locmz);
  emtlRhs = new EmtlVector(locnx+locmz+locmy+locmz, ctx);
  iAmDistrib = 1;
};

sLinsysRootAugEmtl::~sLinsysRootAugEmtl()
{
  if(CtDC) delete CtDC;
  delete redRhs;
  delete emtlRhs;
}


SymMatrix* 
sLinsysRootAugEmtl::createKKT(sData* prob)
{
  int n;

  if(gOuterSolve < 3){
    n = locnx+locmy;
    assert(locmz==0);
  }else{
    n = locnx+locmy+locmz+locmz;
  }
  return new EmtlDenSymMatrix(n, ctx);
}


DoubleLinearSolver*
sLinsysRootAugEmtl::createSolver(sData* prob, SymMatrix* kktmat_)
{

  EmtlDenSymMatrix* kktmat = dynamic_cast<EmtlDenSymMatrix*>(kktmat_);
  return new EmtlGenIndefSolver(kktmat);
}

void sLinsysRootAugEmtl::initializeKKT(sData* prob, Variables* vars)
{
  kkt->scalarMult(0.);
  
}



void sLinsysRootAugEmtl::solveReduced( sData *prob, SimpleVector& b)
{
  assert(locnx+locmz+locmy+locmz==b.length());
  SimpleVector& r = (*redRhs);
  assert(r.length() == b.length());
  SparseGenMatrix& C = prob->getLocalD();
#ifdef DEBUG
  printf("[sLinsysRootAugEmtl::solveReduced 1] b: %1.10e\n",b.onenorm());
#endif

  stochNode->resMon.recDsolveTmLocal_start();


  ///////////////////////////////////////////////////////////////////////
  // LOCAL SOLVE
  ///////////////////////////////////////////////////////////////////////
 
  ///////////////////////////////////////////////////////////////////////
  // b=[b1;b2;b3] is a locnx+locmy+locmz vector 
  // the new rhs should be 
  //           r = [b1-C^T*(zDiag)^{-1}*b3; b2]
  ///////////////////////////////////////////////////////////////////////

  r.copyFromArray(b.elements()); //will copy only as many elems as r has

  // aliases to parts (no mem allocations)
  SimpleVector r4(&r[locnx+locmz+locmy], locmz);
  SimpleVector r3(&r[locnx+locmz], locmy); //r3 is used as a temp buffer for b3
  SimpleVector r2(&r[locnx],       locmz);
  SimpleVector r1(&r[0],           locnx);

  ///////////////////////////////////////////////////////////////////////
  // compute r1 = b1 - C^T*(zDiag)^{-1}*b3
  ///////////////////////////////////////////////////////////////////////

  // Plug in elemental here
  EmtlVector& realRhs = dynamic_cast<EmtlVector&>(*emtlRhs);
  realRhs.copyFromArray(&r[0]);

  // end of elemental
  if(gInnerSCsolve==0) {
    // Option 1. - solve with the factors
    solver->Dsolve(realRhs);
  } else if(gInnerSCsolve==1) {
    printf("solveWithIterRef not implemented.\n");
    assert(0);
    // Option 2 - solve with the factors and perform iter. ref.
    // solveWithIterRef(prob, realRhs);
  } else {
    assert(gInnerSCsolve==2);
    printf("solveWithBiCGStab not implemented.\n");
    assert(0);
    // Option 3 - use the factors as preconditioner and apply BiCGStab
    // solveWithBiCGStab(prob, realRhs);
  }
  realRhs.copyIntoArray(&r[0]);

  ///////////////////////////////////////////////////////////////////////
  // r is the sln to the reduced system
  // the sln to the aug system should be 
  //      x = [r1; r2;  (zDiag)^{-1} * (b3-C*r1);
  ///////////////////////////////////////////////////////////////////////
  SimpleVector b1(&b[0],           locnx);
  SimpleVector b2(&b[locnx],       locmz);
  SimpleVector b3(&b[locnx+locmz], locmy);
  SimpleVector b4(&b[locnx+locmz+locmy], locmz);  
  b1.copyFrom(r1);
  b2.copyFrom(r2);
  b3.copyFrom(r3);
  b4.copyFrom(r4);
#ifdef DEBUG  
  printf("[sLinsysRootAugEmtl::solveReduced 2] b: %1.10e %d %lld\n",b.onenorm(), realRhs.getLocalSize(), r.length());
#endif
  stochNode->resMon.recDsolveTmLocal_stop();

  
}

void sLinsysRootAugEmtl::finalizeKKT(sData* prob, Variables* vars)
{
  int j, p, pend; double val;

  stochNode->resMon.recSchurMultLocal_start();

  EmtlDenSymMatrix * kktd = dynamic_cast<EmtlDenSymMatrix*>(kkt);
  int m__, n__ ; 
  kkt->getSize(m__,n__);
#ifdef DEBUG
  printf("[sLinsysRoot::finalizeKKT 1] kktd: %1.10e\n",kktd->abmaxnorm());
#endif
 

  //////////////////////////////////////////////////////
  // compute Q+diag(xdiag) - C' * diag(zDiag) * C 
  // and update the KKT
  //////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////
  // update the KKT with Q (DO NOT PUT DIAG)
  /////////////////////////////////////////////////////////////
  SparseSymMatrix& Q = prob->getLocalQ();
  int* krowQ=Q.krowM(); int* jcolQ=Q.jcolM(); double* dQ=Q.M();
  for(int i=0; i<locnx; i++) {
    pend = krowQ[i+1];
    for(p=krowQ[i]; p<pend; p++) {     
      j = jcolQ[p]; 
      if(i==j) continue;
      val = dQ[p];
      kktd->symAddAt(i,j,val);
    }
  }
#ifdef DEBUG
  printf("[sLinsysRoot::finalizeKKT 2] kktd: %1.10e\n",kktd->abmaxnorm());
#endif

  
  /////////////////////////////////////////////////////////////
  // update the KKT with the diagonals
  // xDiag is in fact diag(Q)+X^{-1}S
  /////////////////////////////////////////////////////////////
  //kktd->atPutDiagonal( 0, *xDiag );
  SimpleVector& sxDiag = dynamic_cast<SimpleVector&>(*xDiag);
  SimpleVector& syDiag = dynamic_cast<SimpleVector&>(*yDiag);
  
  for(int i=0; i<locnx; i++) kktd->symAddAt(i,i,sxDiag[i]);
#ifdef DEBUG
  printf("[sLinsysRoot::finalizeKKT 3] kktd: %1.10e\n",kktd->abmaxnorm());
#endif

  SimpleVector& ssDiag = dynamic_cast<SimpleVector&>(*sDiag);
  SimpleVector& szDiag = dynamic_cast<SimpleVector&>(*zDiag);

  /////////////////////////////////////////////////////////////
  // update the KKT with  S part
  /////////////////////////////////////////////////////////////
  if(locmz>0) {
    for(int i=locnx; i<locnx+locmz; i++) {
        kktd->symAddAt(i,i,ssDiag[i-locnx]);
    }
  } //~end if locmz>0
  /////////////////////////////////////////////////////////////
  // update the KKT with A (symmetric update forced)
  /////////////////////////////////////////////////////////////
  if(locmy>0){
    kktd->symAtPutSubmatrix( locnx+locmz, 0, prob->getLocalB(), 0, 0, locmy, locnx);
    for(int i=locnx+locmz; i<locnx+locmz+locmy; i++) {
      kktd->symAddAt(i,i,syDiag[i-locnx-locmz]);
    }
  }
  
  /////////////////////////////////////////////////////////////
  // update the KKT with C (symmetric update forced) ,  -I and dual reg
  /////////////////////////////////////////////////////////////  
  if(locmz>0){
    kktd->symAtPutSubmatrix( locnx+locmz+locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx);
    for(int i=0; i<locmz; i++){
      kktd->symAddAt(i+locnx+locmz+locmy,i+locnx,-1.0);
      kktd->symAddAt(i+locnx+locmz+locmy,i+locnx+locmz+locmy,szDiag[i]);
    }
  }
#ifdef DEBUG
  printf("[sLinsysRoot::finalizeKKT 4] kktd: %1.10e\n",kktd->abmaxnorm());
#endif
  stochNode->resMon.recSchurMultLocal_stop();
}

static inline int firstcol(const int mycol,const int startcol,const int npcol)
{
  return (mycol - startcol + (startcol/npcol+1)*npcol)%npcol;
}

const double MAX_MB_FOR_COL_BUFFERS = 100;

int sLinsysRootAugEmtl::factor2(sData *prob, Variables *vars)
{
  int negEVal=0, tempNegEVal=0;
  int return_NegEval=-1;
  int matIsSingular=0,matIsSingularAllReduce;
#ifdef TIMING
  double stime=MPI_Wtime();
  double stime1=MPI_Wtime();
  gprof.n_factor2++;
#endif

  EmtlDenSymMatrix& kktd = dynamic_cast<EmtlDenSymMatrix&>(*kkt);
  int nxP = locnx;
  if(locmz>0) nxP=nxP+2*locmz;
  if(locmy>0) nxP=nxP+locmy;
  
  const int BLOCKSIZE = MIN((1048576*MAX_MB_FOR_COL_BUFFERS/
                            (2*sizeof(double)*nxP)),nxP);

  initializeKKT(prob, vars);
  // we're only sending upper Q block,
  // count how many elements each processor has of this block per column
  // indexed by processor row
  int *nr_counts = new int[ctx.nprow()];
  for (int i = 0; i < ctx.nprow(); i++) {
    nr_counts[i] = Length(nxP, i, ctx.nprow());
  }
  int max_nr = MaxLength(nxP, ctx.nprow());

  DenseGenMatrix colbuffer(2*BLOCKSIZE, nxP);
  double *recvbuffer = new double[max_nr*BLOCKSIZE];
  double *sendbuffer = new double[nxP*BLOCKSIZE];
  int *recvcounts = new int[ctx.nprocs()];
  memset(recvcounts, 0, ctx.nprocs()*sizeof(int));

  //printf("got to factorize\n");
  
  // First tell children to factorize.
  for(unsigned int c=0; c<children.size(); c++) {
#ifdef STOCH_TESTING
    g_scenNum=c;
#endif
    tempNegEVal = children[c]->factor2(prob->children[c], vars);
    if(tempNegEVal<0) {
      matIsSingular = 1;
    }
    else {
      negEVal += tempNegEVal;
    }
  }
#ifdef TIMING
  gprof.t_initializeKKT+=MPI_Wtime()-stime;
#endif
  
  for (int startcol = 0; startcol < nxP; startcol += BLOCKSIZE) {
#ifdef TIMING
    stime=MPI_Wtime();
#endif
    int endcol = MIN(startcol+BLOCKSIZE, nxP); // exclusive
    int numcols = endcol-startcol;
    /*if (ctx.mype() == 0) {
      printf("startcol: %d endcol: %d\n", startcol, endcol);
    }*/
    memset(&colbuffer[0][0], 0, 2*BLOCKSIZE*nxP*sizeof(double));
#ifdef DEBUG  
  printf("[sLinsysRoot::factor2 -1] kktd: %1.2e\n",kktd.abmaxnorm());
#endif
    for (unsigned int c=0; c<children.size(); c++) {
      if(children[c]->mpiComm == MPI_COMM_NULL)
      	continue;
    
      children[c]->stochNode->resMon.recFactTmChildren_start();    
#ifdef DEBUG  
     printf("[sLinsysRoot::factor2 argsaddCols] %d %d\n",startcol, endcol);
#endif
      //---------------------------------------------
      children[c]->addColsToDenseSchurCompl(prob->children[c], colbuffer, startcol, endcol);
      //---------------------------------------------
      children[c]->stochNode->resMon.recFactTmChildren_stop();    
    }
#ifdef TIMING
    gprof.t_initializeKKT+=MPI_Wtime()-stime;
    stime=MPI_Wtime();
#endif
#ifdef DEBUG  
  printf("[sLinsysRoot::factor2 0] kktd: %1.10e\n",kktd.abmaxnorm());
#endif

    // only to improve timing of reduce 
    //MPI_Barrier(ctx.comm());

    //printf("pe %d got columns\n", ctx.mype());
    stochNode->resMon.recReduceTmLocal_start(); 
    // now the fun part, first rearrange the elements so that
    // we can call reducescatter
    // that is, all elements belonging to the first proc go first, etc
    // TODO: this won't work when using the torus
    assert(!ctx.usingTorus());
    int desti = 0;
    for (int pcol = 0; pcol < ctx.npcol(); pcol++)
    for (int prow = 0; prow < ctx.nprow(); prow++) {
      for(int j = firstcol(pcol,startcol,ctx.npcol()); 
            j < numcols; j+= ctx.npcol()) {
        //printf("pe %d for %d %d at col %d\n",ctx.mype(),prow,pcol,j);
        for (int i = prow; i < nxP; i += ctx.nprow()) {
          sendbuffer[desti++] = colbuffer[j][i];
        }
      }
      //printf("pe %d loaded buffer for %d %d\n",ctx.mype(),prow,pcol); 
    }
    assert(desti == nxP*numcols);
    //printf("pe %d loaded send buffer\n", ctx.mype());
    
    int destproc = 0;
    for (int pcol = 0; pcol < ctx.npcol(); pcol++) {
      for (int prow = 0; prow < ctx.nprow(); prow++) {
        // how many columns are we sending to this processor column
        int ncols = Length(
                      numcols-firstcol(pcol,startcol,ctx.npcol()), 
                      0, ctx.npcol());
        recvcounts[destproc++] = nr_counts[prow]*ncols;
      }
    }

    stochNode->resMon.recReduceScatterTmLocal_start();
    
    MPI_Reduce_scatter(sendbuffer, recvbuffer, recvcounts, MPI_DOUBLE, 
      MPI_SUM, ctx.comm());
#ifdef TIMING
    gprof.t_reduceKKTonnode+=MPI_Wtime()-stime;
    stime=MPI_Wtime();
#endif
    MPI_Allreduce(MPI_IN_PLACE, recvbuffer, recvcounts[ctx.mype()], MPI_DOUBLE, MPI_SUM, ctx.intercomm());
    stochNode->resMon.recReduceScatterTmLocal_stop();


    // now unpack on the local processor
    // each column is already in continuous memory
    // if zero equality constraints, actually the whole block is continous
    if ( recvcounts[ctx.mype()] > 0 ) {
      desti = 0;
      int lcol = (startcol+firstcol(ctx.mycol(),startcol,ctx.npcol()))
                    /ctx.npcol();
      int local_nr = nr_counts[ctx.myrow()];
      int ncols = recvcounts[ctx.mype()]/local_nr;
      //printf("%d: %d %d %d %d\n", ctx.mype(), lcol, local_nr, numcols,ncols); 
      for (int j = 0; j < ncols; j++) {
        //printf("%d %d\n", ctx.mype(), j);
        memcpy(kktd.mat->data+(j+lcol)*kktd.getNR(),recvbuffer+j*local_nr,
          local_nr*sizeof(double));
      }
    }

    stochNode->resMon.recReduceTmLocal_stop(); 
#ifdef TIMING
    gprof.t_reduceKKTinternode+=MPI_Wtime()-stime;
#endif

    
  }
  //printf("done factorizing\n");

  delete [] recvcounts;
  delete [] recvbuffer;
  delete [] sendbuffer;
  delete [] nr_counts;
#ifdef TIMING
  stime=MPI_Wtime();
#endif
  
#ifdef DEBUG  
  printf("[sLinsysRoot::factor2 1] kktd: %1.10e\n",kktd.abmaxnorm());
#endif
  finalizeKKT(prob, vars);
#ifdef DEBUG  
  printf("[sLinsysRoot::factor2 2] kktd: %1.10e\n",kktd.abmaxnorm());
#endif
#ifdef TIMING
  gprof.t_finalizeKKT+=MPI_Wtime()-stime;
  stime=MPI_Wtime();
#endif
  
  //double val = kktd.getVal(PROW,PCOL);
  //if (cinfo.mype == 0) {
  //  printf("(%d,%d) --- %f\n", PROW, PCOL, val);
  //}
  MPI_Allreduce(&matIsSingular, &matIsSingularAllReduce, 1, MPI_INT, MPI_SUM, mpiComm);

  if(0==matIsSingularAllReduce){
  	// all the diag mat is nonsingular
  	MPI_Allreduce(&negEVal, &return_NegEval, 1, MPI_INT, MPI_SUM, mpiComm);
#ifdef DEBUG  
    printf("[sLinsysRoot::factor2 3] kktd: %1.10e\n",kktd.abmaxnorm());
#endif
    negEVal = factorizeKKT();
#ifdef DEBUG  
    printf("[sLinsysRoot::factor2 4] kktd: %1.10e\n",kktd.abmaxnorm());
#endif
    if(negEVal<0) {
      return_NegEval = -1;
    } else {
      return_NegEval += negEVal;
    }
  }


#ifdef TIMING
  afterFactor();
  gprof.t_factorizeKKT+=MPI_Wtime()-stime;
#endif

  return return_NegEval;
}

