/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef EMTLSPSYMMATRIX_H
#define EMTLSPSYMMATRIX_H

#include "DoubleMatrix.h"
#include "EmtlDenGenMatrixHandle.h"
#include "EmtlDenGenMatrix.h"

// DEBUG
//#define PROW 10
//#define PCOL 0


class EmtlDenSymMatrix : public SymMatrix {
protected:
  long long m;
public:
  EmtlDenGenMatrixHandle mat;
  const EmtlContext &ctx;

  EmtlDenSymMatrix( long long size, const EmtlContext &ctx_);
  //EmtlDenSymMatrix( int lm, int ln, int size, int nnz );



  virtual int isKindOf( int type );
  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent );
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent );

  virtual void symAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			      int& info );

  virtual void fsymAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			       int& info );

  virtual void getSize( long long& m, long long& n );
  virtual void getSize( int& m, int& n );

  virtual long long size();

  virtual void symAtPutSubmatrix( int destRow, int destCol,
				  DoubleMatrix& M,
				  int srcRow, int srcCol,
				  int rowExtent, int colExtent );

  virtual void fromGetSpRow( int row, int col,
                             double A[], int lenA, int irowA[], int& nnz,
                             int rowExtent, int& info );

  virtual void atPutZeros( int row, int col,
			   int rowExtent, int colExtent );
  virtual void mult ( double beta,  OoqpVector& y,
		      double alpha, OoqpVector& x );
  virtual void transMult ( double beta,  OoqpVector& y,
			   double alpha, OoqpVector& x );
  
  virtual double abmaxnorm();
  
  virtual void writeToStream(ostream& out) const;

  virtual void randomizePSD(double * seed);
  
  virtual void getDiagonal( OoqpVector& vec );
  virtual void setToDiagonal( OoqpVector& vec );
  virtual void atPutDiagonal( int idiag, OoqpVector& v );
  virtual void fromGetDiagonal( int idiag, OoqpVector& x );

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info );

  virtual void SymmetricScale ( OoqpVector& vec );
  virtual void ColumnScale ( OoqpVector& vec );
  virtual void RowScale ( OoqpVector& vec );
  virtual void scalarMult( double num);
  
  virtual void symAtPutZeros( int row, int col,
  			   int rowExtent, int colExtent );
           
  virtual void symAtSetSubmatrix( int destRow, int destCol, DoubleMatrix& M,
				 int srcRow, int srcCol,
				 int rowExtent, int colExtent,bool firstCall, std::map<int,int> &ValIdxMap ); 

  inline int getNR() const { return mat->getNR(); }
  bool isNoop() const;
  // DEBUG USE ONLY:
  double getVal(int row, int col) const;
  // ---

  inline void symPutVal(int row, int col, double val) {
    mat->putVal(row, col, val);
    if (col != row) mat->putVal(col, row, val);
  }
  inline void symAddAt(int row, int col, double val, bool nosym = false) {
    mat->addAt(row, col, val);
    if (col != row && !nosym) mat->addAt(col, row, val);
  }
  
  virtual ~EmtlDenSymMatrix() {};
};

typedef SmartPointer<EmtlDenSymMatrix> EmtlDenSymMatrixHandle;

#endif
