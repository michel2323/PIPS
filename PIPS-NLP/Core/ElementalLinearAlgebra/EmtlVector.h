/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */
/* based on PetscVector.h */


#ifndef EMTLVECTOR_H
#define EMTLVECTOR_H

#include <cstring>
#include <iostream>
#include "EmtlVectorHandle.h"
#include "assert.h"
#include "OoqpVector.h"
#include "EmtlContext.h"
#include "SimpleVector.h"

class EmtlVector : public OoqpVector {
protected:
  int m;
  int nr,nc;
public:
  DistMatrix<double,MC,MR> *A;
  const EmtlContext &ctx;
  double *data;
  static int instances;

  EmtlVector( int m_, const EmtlContext &ctx_ );

  virtual ~EmtlVector();

  virtual int isKindOf( int kind );
  virtual void setToZero();
  virtual void setToConstant( double c );
  virtual void randomize( double alpha, double beta, double *ix );
  virtual void copyFrom( OoqpVector& v );
  virtual double infnorm();
  virtual double onenorm();
  virtual double twonorm();
  virtual void min( double& m, int& index );
  virtual double stepbound(OoqpVector & v, double maxStep );
  virtual double findBlocking(OoqpVector & wstep_vec, 
			      OoqpVector & u_vec, 
			      OoqpVector & ustep_vec, 
			      double maxStep,
			      double *w_elt, 
			      double *wstep_elt,
			      double *u_elt, 
			      double *ustep_elt,
			      int& first_or_second);

  virtual void componentMult( OoqpVector& v );
  virtual void componentDiv ( OoqpVector& v );
  virtual void scalarMult( double num);
  virtual void writeToStream(ostream& out) const;
  virtual void writefToStream( ostream& out,
			       const char format[] ) const;

  virtual void scale( double alpha );

  virtual void axpy  ( double alpha, OoqpVector& x );
  virtual void axzpy ( double alpha, OoqpVector& x, OoqpVector& z );
  virtual void axdzpy( double alpha, OoqpVector& x, OoqpVector& z );

  virtual void addConstant( double c );
  virtual void gondzioProjection( double rmin, double rmax );
  virtual double dotProductWith( OoqpVector& v );
  virtual double shiftedDotProductWith( double alpha, OoqpVector& mystep,
					OoqpVector& yvec,
					double beta,  OoqpVector& ystep );
  virtual void negate();
  virtual void invert();
  virtual int allPositive();

  virtual int matchesNonZeroPattern( OoqpVector& select );
  virtual void selectNonZeros( OoqpVector& select );
  virtual long long numberOfNonzeros();
  virtual void addSomeConstants( double c, OoqpVector& select );
  virtual void writefSomeToStream( ostream& out,
				   const char format[],
				   OoqpVector& select ) const;
  virtual void axdzpy( double alpha, OoqpVector& x,
		       OoqpVector& z, OoqpVector& select );

  virtual int somePositive( OoqpVector& select );
  virtual void divideSome( OoqpVector& div, OoqpVector& select );
  // output is the same on all processors
  virtual void copyIntoArray( double v[] ) const;
  // same array input on all processors, no communication
  virtual void copyFromArray( double v[] );
  virtual void copyFromArray( char v[] );
  

  int getLocalSize() { return nr; };
  //void printout(int tag);
  
  /* TO BE IMPLEMENTED */
  virtual void copyIntoArray( double v[] ) { assert( "Not implemented" && 0 ); };
  
  virtual void findBlockingPD(OoqpVector & wstep_vec, 
            OoqpVector & u_vec, 
            OoqpVector & ustep_vec, 
            double *w_elt, 
            double *wstep_elt,
            double *u_elt, 
            double *ustep_elt,
            int& first_or_second, double * alphaPri, double * alphaDual) { assert( "Not implemented" && 0 ); };  
  
  virtual void print() { assert( "Not implemented" && 0 ); };
  
  virtual double sumLog(OoqpVector* select) { assert( "Not implemented" && 0 ); };
  
  virtual double sumPowElt(const double pow_in) { assert( "Not implemented" && 0 ); };
  
  virtual double sumElt() { assert( "Not implemented" && 0 ); };
  
  virtual void MinComponentOrConstant( OoqpVector* vec_in, double minVal ) { assert( "Not implemented" && 0 ); };
  
  virtual void MinComponentBetween( OoqpVector* vec_in, OoqpVector *select=NULL) { assert( "Not implemented" && 0 ); }; 

  virtual void correctLargeVal( const double testVal, const double corVal, const int absFlag) { assert( "Not implemented" && 0 ); };

  virtual void setToConstantFromTo( double c, int start, int end )  { assert( "Not implemented" && 0 ); }; 

  virtual void SetComponentFromMaxXorY( OoqpVector* x_in, OoqpVector *y_in, OoqpVector *select=NULL) { assert( "Not implemented" && 0 ); };

  virtual void SetComponentFromMaxXorY( OoqpVector* x_in, OoqpVector *y_in, 
			int Start, int End,int xStart, int xEnd, int yStart, int yEnd) { assert( "Not implemented" && 0 ); };

  virtual void SetComponentFromMaxXorConstant( OoqpVector* x_in, const double y_in) { assert( "Not implemented" && 0 ); };

  virtual void SetComponentFromMaxXorConstant( OoqpVector* x_in, const double y_in, 
		  int Start, int End,int xStart, int xEnd) { assert( "Not implemented" && 0 ); };
      
  virtual void copyFromFromTo( OoqpVector* vec_in, int VStart, int VEnd, int VinStart, int VinEnd) { assert( "Not implemented" && 0 ); };

  virtual void copyFromArrayFromTo( double *w, int VStart, int VEnd, int Wstart, int Wend) { assert( "Not implemented" && 0 ); };


  /** Return the absolute value */
  virtual void absVal(OoqpVector *vec_in)  { assert( "Not implemented" && 0 ); };
  


};

#endif

