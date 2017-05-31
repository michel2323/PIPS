/* PIPS
   Authors: Miles Lubin and Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYSEMTLSYM
#define SAUGLINSYSEMTLSYM

#include "sLinsysRoot.h"
#include "sLinsysRootAugEmtl.h"
#include "EmtlContext.h"


class sData;
/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAugEmtlSym : public sLinsysRootAugEmtl {
 protected:
  //sLinsysRootAugEmtlSym() {};

  virtual DoubleLinearSolver* 
                       createSolver  (sData* prob, 
				      SymMatrix* kktmat);

  //virtual void         createChildren(sData* prob) 
  //{sLinsysRoot::createChildren(prob);};
 public:

  sLinsysRootAugEmtlSym(sFactory * factory_, sData * prob_, const EmtlContext &ctx_);
  sLinsysRootAugEmtlSym(sFactory* factory,
			     sData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_,
           OoqpVector* additiveDiag_,
			     const EmtlContext &ctx_);
  virtual ~sLinsysRootAugEmtlSym();

 public:
  virtual int factor2(sData *prob, Variables *vars);
 protected:

};


#endif

