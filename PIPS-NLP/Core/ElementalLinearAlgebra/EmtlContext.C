/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#include "EmtlContext.h"
#include <cstdio>
#include <vector>

EmtlContext::EmtlContext(MPI_Comm comm, int emtlprocs) : mpicomm(comm), _usingTorus(false)
{
  MPI_Comm_size(mpicomm, &_nprocs);
  MPI_Comm_rank(mpicomm, &_mype);
  assert(emtlprocs <= _nprocs);
  if (emtlprocs < 1) emtlprocs = _nprocs;

  // may have to modify this if using elemental at multiple levels
  // i.e. check if already initialized
  // elemental::Init(0,0);

  _noop = (_mype >= emtlprocs);

  #ifdef __bg__
  if (mpicomm != MPI_COMM_WORLD || emtlprocs != _nprocs) {
    if (_mype == 0) {
      printf("Not using torus; need MPI_COMM_WORLD and emtlprocs == nprocs\n");
    }
  } else {
    // make sure we can use the torus with this topology
    MPI_Comm c;
    int ret = MPIX_Cart_comm_create(&c);
    if (ret == MPI_ERR_TOPOLOGY) {
      if (_mype == 0) {
        printf("Not using torus; job needs to run on full partition\n");
      }
    } else {
      assert(ret == MPI_SUCCESS);
      MPI_Comm_free(&c);
      if (_mype == 0) {
        printf("Using torus\n");
      }
      _grid = new Grid(true, true, false, false);
      _usingTorus = true;
      _nprow = _grid->Height();
      _npcol = _grid->Width();
    }
  }
  #endif

  if (!_usingTorus) {

    // MPI_Group allgroup, emtlgroup;
    // MPI_Comm_group(mpicomm, &allgroup);
    // 
    // int *ranks = new int[emtlprocs];
    // for (int i = 0; i < emtlprocs; i++) {
    //   ranks[i] = i;
    // }
    // MPI_Group_incl(allgroup, emtlprocs, ranks, &emtlgroup);
    // MPI_Comm_create(mpicomm, emtlgroup, &emtlcomm);
    // 
    // MPI_Group_free(&emtlgroup);
    // MPI_Group_free(&allgroup);
    // delete [] ranks;
   
    int dims[2] = {0,0};
    MPI_Dims_create(emtlprocs, 2, dims);
    _nprow = dims[0];
    _npcol = dims[1];
  #ifdef DEBUG
    printf("[EmtlContext::EmtlContext] _nprow, _npcol: %d %d\n", _nprow, _npcol);
  #endif

    if (!_noop) {
      // _grid = new Grid(emtlcomm, _nprow, _npcol);
      _grid = new Grid(mpicomm, _nprow);
      assert(_nprow == _grid->Height());
      assert(_npcol == _grid->Width());

    }
  }
  
  if (!_noop) {
    
    _myrow = _grid->MCRank();
    _mycol = _grid->MRRank();
    
    if (!_usingTorus) {
      assert(_myrow == _mype % _nprow);
      assert(_mycol == _mype / _nprow);
    }
    // This is the algorithmic blocksize
    // tune it to the local system, multiples of 32
    SetBlocksize(224);
  }

  if (_usingTorus) {
    // set up map from process grid to global process id
    int my[2] = {_myrow, _mycol};
    std::vector<int> recv(_nprocs*2);
    MPI_Allgather(my, 2, MPI_INT, &recv[0], 2, MPI_INT, mpicomm);
    procmap = new int[_nprocs];
    for (int i = 0; i < _nprocs; i++) {
      procmap[recv[2*i+1] + recv[2*i]*_npcol] = i;
    }
  }

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
