/* PIPS
   Authors: Miles Lubin and Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYSEMTL
#define SAUGLINSYSEMTL

#include "sLinsysRoot.h"
#include "EmtlContext.h"

class sData;
/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAugEmtl : public sLinsysRoot {
 protected:
  //sLinsysRootAugEmtl() {};

  virtual SymMatrix*   createKKT     (sData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (sData* prob, 
				      SymMatrix* kktmat);

  //virtual void         createChildren(QpGenStochData* prob) 
  //{sLinsysRoot::createChildren(prob);};
 public:

  sLinsysRootAugEmtl(sFactory * factory_, sData * prob_, const EmtlContext &ctx_);
  sLinsysRootAugEmtl(sFactory* factory,
			     sData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_,
           OoqpVector* additiveDiag_,
			     const EmtlContext &ctx_);
  virtual ~sLinsysRootAugEmtl();

 public:
  virtual void finalizeKKT(sData* prob, Variables* vars);
  virtual void initializeKKT(sData* prob, Variables* vars);
  virtual int factor2(sData *prob, Variables *vars);
 protected:
  virtual void solveReduced( sData *prob, SimpleVector& b);

  const EmtlContext &ctx;
  SymMatrix* CtDC;
  SimpleVector* redRhs;
  OoqpVector* emtlRhs;
};


#endif

