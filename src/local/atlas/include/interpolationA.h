
#pragma once

#include "atlas/grid.h"
#include "atlas/functionspace.h"
#include "atlas/util/ObjectHandle.h"

#include <vector>

class interpolationA
{
public:

  interpolationA () = default;
  interpolationA (const atlas::grid::Distribution &, const atlas::functionspace::StructuredColumns &,
                  const atlas::grid::Distribution &, const atlas::functionspace::StructuredColumns &);

  atlas::FieldSet shuffle (const atlas::FieldSet &) const;

  int getCnt (int jloc2) const
  {
    return desc[jloc2].icnt;
  }

  int getOff (int jloc2) const
  {
    return desc[jloc2].ioff;
  }

private:
  size_t isize_recv;
  size_t isize_send;

  struct recv_t
  {
    // Here we describe how to use a message from a remote processor
    int iproc;                 // Rank of processor we receive data from
    size_t isize;              // Number of points we receive
    struct loccntoff_t
    {
      atlas::idx_t iloc;     // Local indice
      int          irem_cnt; // Remote count for current local indice
      int          irem_off; // Offset in remote buffer current local index
    };  
    std::vector<loccntoff_t> desc; // Allocated to number of points we have on this task for grid #2
  };

  struct send_t
  {
    int iproc;
    size_t isize;
    std::vector<atlas::idx_t> iloc; // Local points to send
  };

  std::vector<recv_t> yl_recv;
  std::vector<send_t> yl_send;

  std::vector<int> isort; // Sort received points (size = all points from grid #1 we receive on this task)

  // For each point on grid #2 : length/offset of points from grid #1
  
  struct loc2cntoff_t
  {
    int ioff; 
    int icnt;
  };

  std::vector<loc2cntoff_t> desc; // Size = number of points of grid #2 on this task
};

