
#pragma once

#include "atlas/grid.h"
#include "atlas/functionspace.h"
#include "atlas/util/ObjectHandle.h"

#include <vector>

class interpolationAimpl : public atlas::util::Object
{
public:

  enum class opt_t
  {
    OPT_AVG=0,
    OPT_MIN=1,
    OPT_MAX=2,
    OPT_SUM=3
  };

  interpolationAimpl () = default;
  interpolationAimpl (const atlas::grid::Distribution &, const atlas::functionspace::StructuredColumns &,
                      const atlas::grid::Distribution &, const atlas::functionspace::StructuredColumns &,
                      const bool ldopenmp);

  template <typename T>
  atlas::FieldSet shuffle (const atlas::FieldSet &) const;

  template <typename T> atlas::FieldSet
  interpolate (const atlas::FieldSet &, const opt_t = opt_t::OPT_AVG) const;

  int getCnt (int jloc2) const
  {
    return desc[jloc2].icnt;
  }

  int getOff (int jloc2) const
  {
    return desc[jloc2].ioff;
  }

  int getLen () const
  {
    return desc.size ();
  }

  int getMissing () const
  {
    return isize_miss;
  }

private:

  template <typename T, typename O, typename E>
  void reduce (atlas::array::ArrayView<T,1> & v2, atlas::array::ArrayView<T,1> & v2e, 
               const size_t size2, T t0, T zundef, E eval, O op) const;

  size_t isize_recv = 0;
  size_t isize_send = 0;
  size_t isize_miss = 0;
  bool verbose = false;
  bool llopenmp = true;  // Use OpenMP when possible

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
  using opt_t = interpolationAimpl::opt_t;

  interpolationA () = default;

  interpolationA 
  (const atlas::grid::Distribution & dist1, const atlas::functionspace::StructuredColumns & fs1,
   const atlas::grid::Distribution & dist2, const atlas::functionspace::StructuredColumns & fs2,
   const bool ldopenmp)
   : Handle::Handle (new interpolationAimpl (dist1, fs1, dist2, fs2, ldopenmp))
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

  int getLen () const
  {
    return get ()->getLen ();
  }

  template <typename T> atlas::FieldSet
  interpolate (const atlas::FieldSet & pgp1, const opt_t = opt_t::OPT_AVG) const
  {
    return get ()->interpolate<T> (pgp1, opt_t);
  }
};

extern "C"
{
interpolationAimpl * interpolationA__new 
  (const atlas::grid::DistributionImpl *, const atlas::functionspace::detail::StructuredColumns *,
   const atlas::grid::DistributionImpl *, const atlas::functionspace::detail::StructuredColumns *,
   const int);
atlas::field::FieldSetImpl * interpolationA__interpolate (interpolationAimpl *, atlas::field::FieldSetImpl *, const int);
atlas::field::FieldSetImpl * interpolationA__shuffle (interpolationAimpl *, atlas::field::FieldSetImpl *);
int interpolationA__getlen (const interpolationAimpl *);
void interpolationA__getcnt (const interpolationAimpl *, int cnt[]);
void interpolationA__getoff (const interpolationAimpl *, int off[]);
};

