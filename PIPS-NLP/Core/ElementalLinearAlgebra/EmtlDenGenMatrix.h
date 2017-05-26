/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#ifndef EMTLSPGENMATRIXBODY 
#define EMTLSPGENMATRIXBODY

#include "DoubleMatrix.h"
#include "EmtlDenGenMatrixHandle.h"
#include "EmtlContext.h"


class EmtlDenGenMatrix : public GenMatrix {
protected:
  int m,n;
  int nr,nc;
public:
  DistMatrix<double,MC,MR> *A;
  const EmtlContext &ctx;
  // we need fast access to the data
  double *data;
  EmtlDenGenMatrix(int m, int n, const EmtlContext &ctx_);


  virtual int isKindOf( int type );

  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent );
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent );
  virtual void fromGetSpRow( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int colExtent, int& info );
  virtual void atPutSpRow( int row, double A[], int lenA, int jcolA[],
			   int& info );

  virtual void getSize( int& m, int& n );
  virtual void getSize( long long& m, long long& n );

  virtual void atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent );
			       
  virtual void atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent, bool transpose);

  virtual void mult ( double beta,  OoqpVector& y,
		      double alpha, OoqpVector& x );

  virtual void transMult ( double beta,  OoqpVector& y,
			   double alpha, OoqpVector& x );

  virtual void getDiagonal( OoqpVector& vec );
  virtual void setToDiagonal( OoqpVector& vec );

  virtual void atPutDiagonal( int idiag, OoqpVector& v );
  virtual void fromGetDiagonal( int idiag, OoqpVector& v );

  virtual double abmaxnorm();

  virtual void writeToStream(ostream& out) const;

  virtual void randomize(double alpha, double beta, double * seed);

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info );
				

  virtual void SymmetricScale ( OoqpVector& vec );
  virtual void ColumnScale ( OoqpVector& vec );
  virtual void RowScale ( OoqpVector& vec );
  virtual void scalarMult( double num);
  
  virtual void atPutZeros( int row, int col,
  			   int rowExtent, int colExtent );

  virtual void matTransDMultMat(OoqpVector&, SymMatrix**);
  virtual void matTransDinvMultMat(OoqpVector&, SymMatrix**);
  
  // needed for Solver
  int getNR() const { return nr; }

  bool isNoop() const { return ctx.noop(); }
  
  // called with the same value from all procs
  // store the value if it belongs to this proc
  inline void putVal(int row, int col, double val) {
    int proc = ctx.index2proc(row,col);
    if (ctx.mype() == proc) {
      int lrow, lcol;
      ctx.global2local(row,col,lrow,lcol);
      data[lrow + nr*lcol] = val;
    }
  }
  inline void addAt(int row, int col, double val) {
    int proc = ctx.index2proc(row,col);
    if (ctx.mype() == proc) {
      int lrow, lcol;
      ctx.global2local(row,col,lrow,lcol);
      data[lrow + nr*lcol] += val;
    }
  }
  
  virtual ~EmtlDenGenMatrix();
};


#endif

