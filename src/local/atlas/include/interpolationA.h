
#pragma once

#include "atlas/grid.h"
#include "atlas/functionspace.h"
#include "atlas/util/ObjectHandle.h"

#include <vector>

class interpolationAimpl
{
public:

  interpolationAimpl () = default;
  interpolationAimpl (const atlas::grid::Distribution &, const atlas::functionspace::StructuredColumns &,
                      const atlas::grid::Distribution &, const atlas::functionspace::StructuredColumns &);

  ~interpolationAimpl ();

  template <typename T>
  atlas::FieldSet shuffle (const atlas::FieldSet &) const;

  template <typename T> atlas::FieldSet
  interpolate (const atlas::FieldSet &) const;

  int getCnt (int jloc2) const
  {
    return desc[jloc2].icnt;
  }

  int getOff (int jloc2) const
  {
    return desc[jloc2].ioff;
  }

  int getMissing () const
  {
    return isize_miss;
  }

private:
  size_t isize_recv = 0;
  size_t isize_send = 0;
  size_t isize_miss = 0;
  bool verbose = false;

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

  const atlas::StructuredGrid grid1;
  const atlas::StructuredGrid grid2;
  const atlas::grid::Distribution dist1;
  const atlas::grid::Distribution dist2;

  const atlas::functionspace::StructuredColumns fs1;
  const atlas::functionspace::StructuredColumns fs2;

};


class interpolationA : public atlas::util::ObjectHandle<interpolationAimpl>
{
public:
  interpolationA () = default;

  interpolationA 
  (const atlas::grid::Distribution & dist1, const atlas::functionspace::StructuredColumns & fs1,
   const atlas::grid::Distribution & dist2, const atlas::functionspace::StructuredColumns & fs2)
   : Handle::Handle (new interpolationAimpl (dist1, fs1, dist2, fs2))
  {
  }

  template <typename T> atlas::FieldSet
  shuffle (const atlas::FieldSet & pgp1) const
  {
    return get ()->shuffle<T> (pgp1);
  }

  int getCnt (int jloc2) const
  {
    return get ()->getCnt (jloc2);
  }

  int getOff (int jloc2) const
  {
    return get ()->getOff (jloc2);
  }

  template <typename T> atlas::FieldSet
  interpolate (const atlas::FieldSet & pgp1) const
  {
    return get ()->interpolate<T> (pgp1);
  }
};

extern "C"
{
interpolationAimpl * interpolationA__new 
  (const atlas::grid::DistributionImpl *, const atlas::functionspace::detail::StructuredColumns *,
   const atlas::grid::DistributionImpl *, const atlas::functionspace::detail::StructuredColumns *);
atlas::field::FieldSetImpl * interpolationA__interpolate (interpolationAimpl *, atlas::field::FieldSetImpl *);
};

