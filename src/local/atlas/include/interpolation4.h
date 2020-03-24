#pragma once

#include "atlas/grid.h"
#include "atlas/functionspace.h"
#include "atlas/field.h"

#include <vector>

class interpolate4
{
public:
  
  interpolate4 () = default;

  interpolate4 (const atlas::grid::Distribution &, const atlas::functionspace::StructuredColumns &,
              const atlas::grid::Distribution &, const atlas::functionspace::StructuredColumns &);

  template <typename T> atlas::FieldSet
  shuffle (const atlas::FieldSet &);

  template <typename T> atlas::FieldSet
  interpolate (const atlas::FieldSet &);

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

  atlas::functionspace::StructuredColumns fs1;
  atlas::functionspace::StructuredColumns fs2;

  static const int ISW, ISE, INW, INE;

  class weights4_t weights4;
};





