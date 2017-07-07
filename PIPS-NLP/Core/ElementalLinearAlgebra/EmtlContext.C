/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#include "../../global_var.h"
#include "EmtlContext.h"
#include <cstdio>
#include <vector>

EmtlContext::EmtlContext(MPI_Comm comm, int emtlprocs) : mpicomm(comm), _usingTorus(false)
{
  _noop=false;
  MPI_Comm_split_type(comm,MPI_COMM_TYPE_SHARED,gmyid,MPI_INFO_NULL,&emtlcomm);
  MPI_Comm_rank(emtlcomm, &_mype);
  MPI_Comm_size(emtlcomm, &emtlprocs);
  if(gmyid==0) {
    printf("%d emtlprocs\n", emtlprocs);
  }
  MPI_Comm_size(emtlcomm, &_nprocs);
  int dims[2] = {0,0};
  MPI_Dims_create(emtlprocs, 2, dims);
  _nprow = dims[0];
  _npcol = dims[1];
  _grid = new Grid(emtlcomm, _nprow);
  assert(_nprow == _grid->Height());
  assert(_npcol == _grid->Width());
  _myrow = _grid->MCRank();
  _mycol = _grid->MRRank();
  SetBlocksize(224);
  return;
}

EmtlContext::~EmtlContext() 
{
  if (!_noop) {
    delete _grid;
  }
  if (!_usingTorus) {
    //MPI_Comm_free(&emtlcomm);
  } else {
    delete [] procmap;
  }
  // elemental::Finalize();
}
