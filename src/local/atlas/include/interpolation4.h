#pragma once

#include "atlas/grid.h"
#include "atlas/functionspace.h"
#include "atlas/util/ObjectHandle.h"
#include "atlas/field.h"

#include <vector>

class interpolation4impl
{
public:
  
  interpolation4impl () = default;

  interpolation4impl (const atlas::grid::Distribution &, const atlas::functionspace::StructuredColumns &,
                      const atlas::grid::Distribution &, const atlas::functionspace::StructuredColumns &);

  template <typename T> atlas::FieldSet
  shuffle (const atlas::FieldSet &) const;

  template <typename T> atlas::FieldSet
  interpolate (const atlas::FieldSet &) const;

  enum
  {
    ISW = -3, 
    ISE = -2, 
    INW = -1, 
    INE = -0 
  };


private:

  class weights4_t
  {
  public:
    weights4_t () = default;
    std::vector<double> values;
  };

  void create_weights ();

  class recv_t
  {
  public:
    int iprc, icnt, ioff;
  };

  class send_t
  {
  public:
    int iprc, icnt, ioff;
    std::vector<atlas::gidx_t> iglo;
    std::vector<atlas::idx_t> iloc;
  };

  size_t isize_recv;
  size_t isize_send;

  std::vector<int> isort;
  std::vector<recv_t> yl_recv;
  std::vector<send_t> yl_send;
  
  // Should be const
  const atlas::StructuredGrid grid1;
  const atlas::StructuredGrid grid2;
  const atlas::grid::Distribution dist1;
  const atlas::grid::Distribution dist2;

  size_t size1;
  size_t size2;

  const atlas::functionspace::StructuredColumns fs1;
  const atlas::functionspace::StructuredColumns fs2;


  class weights4_t weights4;
};


class interpolation4 : public atlas::util::ObjectHandle<interpolation4impl> 
{
public:
  interpolation4 () = default;

  interpolation4 
  (const atlas::grid::Distribution & dist1, const atlas::functionspace::StructuredColumns & fs1,
   const atlas::grid::Distribution & dist2, const atlas::functionspace::StructuredColumns & fs2)
   : Handle::Handle (new interpolation4impl (dist1, fs1, dist2, fs2))
  {
  }

  template <typename T> atlas::FieldSet
  shuffle (const atlas::FieldSet & pgp1) const
  {
    return get ()->shuffle<T> (pgp1);
  }

  template <typename T> atlas::FieldSet
  interpolate (const atlas::FieldSet & pgp1) const
  {
    return get ()->interpolate<T> (pgp1);
  }

};


extern "C"
{
interpolation4impl * interpolation4__new 
  (const atlas::grid::DistributionImpl *, const atlas::functionspace::detail::StructuredColumns *,
   const atlas::grid::DistributionImpl *, const atlas::functionspace::detail::StructuredColumns *);
atlas::field::FieldSetImpl * interpolation4__interpolate (interpolation4impl *, atlas::field::FieldSetImpl *);
};

