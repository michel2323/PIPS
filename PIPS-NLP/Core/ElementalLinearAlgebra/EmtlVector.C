/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */


#include "EmtlVector.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "VectorUtilities.h"
#include <cassert>
#include <cfloat>
#include <cstring>
#include <cmath>


struct ooqp_mpi_double_int{ double val; int rank; };

int EmtlVector::instances = 0;
EmtlVector::EmtlVector( int m_, const EmtlContext &ctx_ ) : OoqpVector( m_ ), ctx(ctx_)
{
  this->m = m_;

  // note the loops over i < nr*nc won't do anything if nc == 0
  
  if (!ctx.noop()) {
    A = new DistMatrix<double,MC,MR>(m, 1, ctx.grid());

    nr = A->LocalHeight();
    nc = A->LocalWidth();
    data = A->LocalMatrix().Buffer();
    assert(A->LocalMatrix().LDim() == nr);
  } else {
    data = 0;
    nr = nc = 0;
  }


 EmtlVector::instances++;
}

EmtlVector::~EmtlVector()
{
  if (!ctx.noop()) {
    delete A;
  }
  EmtlVector::instances--;
}

int EmtlVector::numberOfNonzeros()
{
  int mycount = 0;
  for (int i = 0; i < nr*nc; i++)
  {
    if (data[i] != 0) mycount++;
  }
    
  int allcount;
  MPI_Allreduce(&mycount, &allcount, 1, MPI_INT, MPI_SUM, ctx.comm());
  
  return allcount;
}
int EmtlVector::isKindOf( int kind )
{
  return (kEmtlVector == kind );
}

void EmtlVector::setToZero()
{
  for (int i = 0; i < nr*nc; i++) data[i] = 0.0;
}

void EmtlVector::setToConstant( double c )
{
  for (int i = 0; i < nr*nc; i++) data[i] = c;
}

void EmtlVector::randomize( double alpha, double beta, double * /* ix */ )
{
  assert( "Not implemented" && 0 );
}

void EmtlVector::copyIntoArray( double v[] ) const
{
  int maxsize = utilities::MaxLocalLength(m, ctx.nprow());
  double *buf = new double[maxsize];
  
  
  for (int row = 0; row < ctx.nprow(); row++) {
    int nr_row = utilities::LocalLength(m, row, ctx.nprow());
    if (nr_row < 1) continue;
    int proc = ctx.get_pnum(row, 0);
    if (ctx.mype() == proc) {
      assert(nr_row == nr && nc == 1);
      memcpy(buf, data, nr_row*sizeof(double));
    }
    MPI_Bcast(buf, nr_row, MPI_DOUBLE, proc, ctx.comm());
    for (int i = 0; i < nr_row; i++) {
      int gindex = row + i*ctx.nprow();
      v[gindex] = buf[i];
    }
  }
  delete [] buf;
}

void EmtlVector::copyFromArray( double v[] )
{
  for (int i = 0; i < nr*nc; i++) {
    int gindex = ctx.myrow() + i*ctx.nprow();
    data[i] = v[gindex];
  }
}

void EmtlVector::copyFromArray( char v[] )
{
  assert( "Not implemented" && 0 );
}

void EmtlVector::copyFrom( OoqpVector& v_in )
{
  if (v_in.isKindOf(kEmtlVector)) {
    EmtlVector &vec = (EmtlVector&) v_in;
    assert(this->m == vec.m && this->nr == vec.nr);
    for (int i = 0; i < nr*nc; i++)
    {
      data[i] = vec.data[i];
    }
  } else {
    assert( "Not implemented" && 0 );
  }
}

double EmtlVector::infnorm()
{
  double m = 0, mall;
  for (int i = 0; i < nr*nc; i++) {
    m = MAX(abs(data[i]),m);
  }
  MPI_Allreduce(&m, &mall, 1, MPI_DOUBLE, MPI_MAX, ctx.comm());
  
  return mall;
  
  //int one = 1;
  //return pdlange_("Max", &m, &one, data, &one, &one, desc, NULL);
}


void EmtlVector::componentMult( OoqpVector& v_in )
{
  if (v_in.isKindOf(kEmtlVector)) {
    EmtlVector &vec = (EmtlVector&) v_in;
    assert(this->m == vec.m && this->nr == vec.nr);
    for (int i = 0; i < nr*nc; i++)
    {
      data[i] *= vec.data[i];
    }
  } else {
    assert( "Not implemented" && 0 );
  }
}

void EmtlVector::componentDiv ( OoqpVector& v_in )
{
  if (v_in.isKindOf(kEmtlVector)) {
    EmtlVector &vec = (EmtlVector&) v_in;
    assert(this->m == vec.m && this->nr == vec.nr);
    for (int i = 0; i < nr*nc; i++)
    {
      data[i] /= vec.data[i];
    }
  } else {
    assert( "Not implemented" && 0 );
  }
}

void EmtlVector::writeToStream(ostream& out) const
{
  assert( "Not implemented" && 0 );
}

void EmtlVector::writefToStream( ostream& out,
				      const char format[] ) const
{
  assert( "Not implemented" && 0 );
}


void EmtlVector::scale( double alpha )
{
  for (int i = 0; i < nr*nc; i++) data[i] *= alpha;
}


void EmtlVector::axpy  ( double alpha, OoqpVector& x )
{
  if (x.isKindOf(kEmtlVector)) {
    EmtlVector &vec = (EmtlVector&) x;
    assert(this->m == vec.m && this->nr == vec.nr);
    for (int i = 0; i < nr*nc; i++)
    {
      data[i] += alpha*vec.data[i];
    }
  } else {
    assert( "Not implemented" && 0 );
  }
}

void EmtlVector::axzpy ( double alpha, OoqpVector& vx_in, OoqpVector& vz_in )
{
  if (vx_in.isKindOf(kEmtlVector) && vz_in.isKindOf(kEmtlVector)) {
    EmtlVector &x = (EmtlVector&) vx_in;
    EmtlVector &z = (EmtlVector&) vz_in;
    assert(this->m == x.m && this->nr == x.nr && this->m == z.m && this->nr == z.nr);
    for (int i = 0; i < nr*nc; i++)
    {
      data[i] += alpha * x.data[i] * z.data[i];
    }
  } else {
    assert( "Not implemented" && 0 );
  }
}

void EmtlVector::axdzpy( double alpha, OoqpVector& vx_in,
			      OoqpVector& vz_in )
{
  if (vx_in.isKindOf(kEmtlVector) && vz_in.isKindOf(kEmtlVector)) {
    EmtlVector &x = (EmtlVector&) vx_in;
    EmtlVector &z = (EmtlVector&) vz_in;
    assert(this->m == x.m && this->nr == x.nr && this->m == z.m && this->nr == z.nr);
    for (int i = 0; i < nr*nc; i++)
    {
      data[i] += alpha * x.data[i] / z.data[i];
    }
  } else {
    assert( "Not implemented" && 0 );
  }
}

//#undef __FUNC__
//#define __FUNC__ EmtlVector::ADDCONSTANT
void EmtlVector::addConstant( double c )
{
  for (int i = 0; i < nr*nc; i++) data[i] += c;
}

void EmtlVector::gondzioProjection( double rmin, double rmax )
{
  assert( "Not implemented" && 0 );
}

double EmtlVector::dotProductWith( OoqpVector& v )
{
  EmtlVector &vec = dynamic_cast<EmtlVector&>(v);
  double dot = 0., alldot;
  assert(this->m == vec.m && this->nr == vec.nr);
  
  for (int i = 0; i < nr*nc; i++) dot += data[i]*vec.data[i];
  
  MPI_Allreduce(&dot, &alldot, 1, MPI_DOUBLE, MPI_SUM, ctx.comm());
  
  return alldot;
}

double EmtlVector::shiftedDotProductWith( double alpha, OoqpVector& pvec_in,
					       OoqpVector& yvec_in,
					       double beta,  OoqpVector& qvec_in )
{

  double dot = 0., alldot;
  
  EmtlVector &mystep = dynamic_cast<EmtlVector&>(pvec_in);
  EmtlVector &yvec = dynamic_cast<EmtlVector&>(yvec_in);
  EmtlVector &ystep = dynamic_cast<EmtlVector&>(qvec_in);
  assert(this->m == mystep.m && this->nr == mystep.nr);
  assert(this->m == yvec.m && this->nr == yvec.nr);
  assert(this->m == ystep.m && this->nr == ystep.nr);
  
  for (int i = 0; i < nr*nc; i++)
    dot += (data[i] + alpha*mystep.data[i])*(yvec.data[i]+beta*ystep.data[i]);
  
  MPI_Allreduce(&dot, &alldot, 1, MPI_DOUBLE, MPI_SUM, ctx.comm());
  
  return alldot;
}

void EmtlVector::negate()
{
  for (int i = 0; i < nr*nc; i++) data[i] *= -1.0;
}

void EmtlVector::invert()
{
  for (int i = 0; i < nr*nc; i++) data[i] = 1.0/data[i];
}

int EmtlVector::allPositive()
{
  int good = 1;
  for (int i = 0; i < nr*nc; i++)
  {
    if (data[i] <= 0.0) {
      good = 0;
      break;
    }
  }
  int allgood;
    
  MPI_Allreduce(&good, &allgood, 1, MPI_INT, MPI_LAND, ctx.comm());

  return allgood;

}


int EmtlVector::matchesNonZeroPattern( OoqpVector& select )
{
  EmtlVector &sel = dynamic_cast<EmtlVector&>(select);
  int good = 1;
  assert(this->m == sel.m && this->nr == sel.nr);
  
  for (int i = 0; i < nr*nc ; i++)
  {
    // only care about nonzeros in this that aren't in select?
    if (data[i] != 0.0 && sel.data[i] == 0.0) 
    {
      good = 0;
      break;
    }
  }

  int allgood;
  MPI_Allreduce(&good, &allgood, 1, MPI_INT, MPI_LAND, ctx.comm());

  return allgood;
}

void EmtlVector::selectNonZeros( OoqpVector& select )
{
  EmtlVector &sel = dynamic_cast<EmtlVector&>(select);
  assert(this->m == sel.m && this->nr == sel.nr);
  for (int i = 0; i < nr*nc; i++) 
  {
    if (sel.data[i] == 0.0) data[i] = 0.0;
  }
  
}

void EmtlVector::addSomeConstants( double c, OoqpVector& select )
{
  EmtlVector &sel = dynamic_cast<EmtlVector&>(select);
  assert(this->m == sel.m && this->nr == sel.nr);
  for (int i = 0; i < nr*nc; i++)
  {
    if (sel.data[i] != 0.0) data[i] += c;
  }

}

void EmtlVector::writefSomeToStream( ostream& out,
					  const char format[],
					  OoqpVector& select ) const
{
  assert( "Not implemented" && 0 );
}

void EmtlVector::axdzpy( double alpha, OoqpVector& vx_in,
			      OoqpVector& vz_in, OoqpVector& select )
{
  if (vx_in.isKindOf(kEmtlVector) && vz_in.isKindOf(kEmtlVector)) {
    EmtlVector &x = (EmtlVector&) vx_in;
    EmtlVector &z = (EmtlVector&) vz_in;
    EmtlVector &sel = dynamic_cast<EmtlVector&>(select);
    assert(this->m == sel.m && this->nr == sel.nr);
    assert(this->m == x.m && this->nr == x.nr && this->m == z.m && this->nr == z.nr);
    for (int i = 0; i < nr*nc; i++)
    {
      if (sel.data[i] != 0.0) data[i] += alpha * x.data[i] / z.data[i];
    }
  } else {
    assert( "Not implemented" && 0 );
  }

}

int EmtlVector::somePositive( OoqpVector& select ) {
  assert( "Not implemented" && 0 );
  return 0;
}

void EmtlVector::divideSome( OoqpVector& vecb, OoqpVector& select )
{
    if (vecb.isKindOf(kEmtlVector)) {
    EmtlVector &x = (EmtlVector&) vecb;
    EmtlVector &sel = dynamic_cast<EmtlVector&>(select);
    assert(this->m == sel.m && this->nr == sel.nr);
    assert(this->m == x.m && this->nr == x.nr);
    for (int i = 0; i < nr*nc; i++)
    {
      if (sel.data[i] != 0.0) data[i] /= x.data[i];
    }
  } else {
    assert( "Not implemented" && 0 );
  }
}


double EmtlVector::stepbound(OoqpVector & svec_in,
				  double bound )
{
  EmtlVector &vec = dynamic_cast<EmtlVector&>(svec_in);
  assert(this->m == vec.m && this->nr == vec.nr);
  
  bound = ::stepbound(data, nr*nc, 1, vec.data, 1, bound);
  
  double allbound;
  MPI_Allreduce(&bound, &allbound, 1, MPI_DOUBLE, MPI_MIN,
		ctx.comm());

  return allbound;
}


double EmtlVector::onenorm()
{
  assert( "Not implemented" && 0 );
}

double EmtlVector::twonorm()
{
  double n2 = 0, nall;
  for (int i = 0; i < nr*nc; i++) {
    n2 += data[i]*data[i];
  }
  MPI_Allreduce(&n2, &nall, 1, MPI_DOUBLE, MPI_SUM, ctx.comm());
  
  return sqrt(nall);
}

double EmtlVector::findBlocking(OoqpVector & wstep_vec, 
				     OoqpVector & u_vec, 
				     OoqpVector & ustep_vec, 
				     double maxStep,
				     double *w_elt, 
				     double *wstep_elt,
				     double *u_elt, 
				     double *ustep_elt,
				     int& first_or_second)
{
  int minrank;
  double bound;
  
  EmtlVector &wstep = dynamic_cast<EmtlVector&>(wstep_vec);
  EmtlVector &u = dynamic_cast<EmtlVector&>(u_vec);
  EmtlVector &ustep = dynamic_cast<EmtlVector&>(ustep_vec);
  
  bound = ::find_blocking(data, nr*nc, 1, wstep.data, 1, u.data, 1, ustep.data, 
			  1,maxStep, w_elt, wstep_elt, u_elt, ustep_elt, first_or_second);
		
  ooqp_mpi_double_int mini, allmini;
  
  mini.val = bound; mini.rank = ctx.mype();
  allmini.val = 0; allmini.rank = 0;
  MPI_Allreduce(&mini, &allmini, 1, MPI_DOUBLE_INT, MPI_MINLOC, ctx.comm());
  bound = allmini.val; minrank = allmini.rank;

  double steps[4] = { *w_elt, *wstep_elt, *u_elt, *ustep_elt };
  MPI_Bcast(steps, 4, MPI_DOUBLE, minrank, ctx.comm());
  
  *w_elt = steps[0];     *wstep_elt = steps[1];
  *u_elt = steps[2];     *ustep_elt = steps[3];
  
  MPI_Bcast(&first_or_second, 1, MPI_INT, minrank, ctx.comm());
  
//    cout << " fos2 " << first_or_second << " rank " << selfrank << endl;
  return bound;
}

void EmtlVector::min( double& min, int& index )
{
  ooqp_mpi_double_int di, alldi;
  di.val = DBL_MAX;
  di.rank = -1;
  for (int i = 0; i < nr*nc; i++) {
    if (data[i] < di.val) {
      di.val = data[i];
      di.rank = i;
    }
  }
  // convert to global index
  di.rank = ctx.myrow() + di.rank*ctx.nprow();
  MPI_Allreduce(&di, &alldi, 1, MPI_DOUBLE_INT, MPI_MINLOC, ctx.comm());
  min = alldi.val;
  index = alldi.rank;

}

void EmtlVector::scalarMult( double num)
{
  this->scale(num);
}





