/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#ifndef EMTLCONTEXT_H
#define EMTLCONTEXT_H

#define WITHOUT_COMPLEX // weird g++ bug
#include "El.hpp"
#include "mpi.h"
#include <cassert>

using namespace El;

#ifndef MAX
#define MAX(a,b) ( (a>b) ? a : b )
#define MIN(a,b) ( (a>b) ? b : a )
#endif

class EmtlContext {
  public:
  EmtlContext(MPI_Comm comm, int emtlprocs);
  ~EmtlContext();

  // local processors
  inline int mype() const { return _mype; }
  // total processors in communicator
  inline int nprocs() const { return _nprocs; }
  // row in processor grid
  inline int myrow() const { assert(!_noop); return _myrow; }
  // col in processor grid
  inline int mycol() const { assert(!_noop); return _mycol; }
  // number of cols in processor grid
  inline int npcol() const { return _npcol; }
  // number of rows in processor grid
  inline int nprow() const { return _nprow; }
  inline MPI_Comm comm() const { return mpicomm; }
  inline Grid& grid() const { assert(!_noop); return *_grid; }
  // No-op, true if this processor isn't used for elemental
  inline bool noop() const { return _noop; }

  inline int get_pnum(const int prow, const int pcol) const { 
    if (!_usingTorus) {
      // column-major
      return prow + pcol*_nprow;
    } else {
      return procmap[pcol + _npcol*prow];
    }
  }
  // returns row/column index
  inline int index2proc(const int row, const int col) const {
    int prow = row % _nprow;
    int pcol = col % _npcol;
    return get_pnum(prow, pcol);
  }
  inline void global2local(const int row, const int col, int &lrow, int &lcol) const {
    lrow = row / _nprow;
    lcol = col / _npcol;
  }

  inline bool usingTorus() const { return _usingTorus; }

  protected:
  MPI_Comm mpicomm, emtlcomm;
  int _mype;
  int _nprocs;
  int _myrow;
  int _mycol;
  int _npcol;
  int _nprow;
  bool _noop;
  Grid *_grid;
  // row-major nprow x npcol matrix
  int *procmap;
  bool _usingTorus;

};




#endif
