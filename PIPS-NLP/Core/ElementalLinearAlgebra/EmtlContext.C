/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#include "../../global_var.h"
#include "EmtlContext.h"
#include <cstdio>
#include <vector>

EmtlContext::EmtlContext(MPI_Comm comm, int emtlprocs) : _usingTorus(false)
{
  _noop=false;
  // Build intranode communicator
  MPI_Comm_split_type(comm,MPI_COMM_TYPE_SHARED,gmyid,MPI_INFO_NULL,&intracomm);
  // MPI_Comm_split(comm,gmyid%2,gmyid,&intracomm);
  // intracomm=MPI_COMM_WORLD;
  MPI_Comm_rank(intracomm, &_mype);
  MPI_Comm_size(intracomm, &emtlprocs);
  
  // Build internode communicator
  int color=0;
  color=_mype; 
  MPI_Comm_split(comm,color,gmyid,&_intercomm);
  
  if(gmyid==0) {
    printf("%d emtlprocs\n", emtlprocs);
  }
  MPI_Comm_size(intracomm, &_nprocs);
  int dims[2] = {0,0};
  MPI_Dims_create(emtlprocs, 2, dims);
  _nprow = dims[0];
  _npcol = dims[1];
  _grid = new Grid(intracomm, _nprow);
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
