#include "atlas/library/Library.h"
#include "atlas/util/Config.h"
#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"

#include <algorithm>
#include <vector>

#include <sys/time.h>
#include <stdio.h>
#include <math.h>

namespace
{

const double deg2rad = M_PI / 180.0;
const double rad2deg = 180.0 / M_PI;

};

atlas::FieldSet 
getXYZ (const atlas::functionspace::StructuredColumns & fs)
{
  atlas::FieldSet xyz;

  auto t = atlas::array::DataType::kind<double> ();
  auto s = atlas::array::make_shape (fs.sizeOwned ());

  xyz.add (atlas::Field (std::string ("x"), t, s));
  xyz.add (atlas::Field (std::string ("y"), t, s));
  xyz.add (atlas::Field (std::string ("z"), t, s));
  
  auto x = atlas::array::make_view<double,1> (xyz[0]);
  auto y = atlas::array::make_view<double,1> (xyz[1]);
  auto z = atlas::array::make_view<double,1> (xyz[2]);

  auto i = atlas::array::make_view<int,1> (fs.index_i ());
  auto j = atlas::array::make_view<int,1> (fs.index_j ());

  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {
      atlas::PointLonLat ll = fs.grid ().StructuredGrid::lonlat (i (jloc)-1, j (jloc)-1);
      double coslon = cos (deg2rad * ll.lon ()), sinlon = sin (deg2rad * ll.lon ());
      double coslat = cos (deg2rad * ll.lat ()), sinlat = sin (deg2rad * ll.lat ());
      x (jloc) = coslon * coslat; 
      y (jloc) = sinlon * coslat; 
      z (jloc) = sinlat;
    }

  return xyz;
}


atlas::FieldSet
gradient (const atlas::functionspace::StructuredColumns & fs, const atlas::FieldSet & pgp)
{
  atlas::FieldSet pgrad;




  return pgrad;
}


int main (int argc, char * argv[]) 
{
  atlas::Library::instance ().initialise (argc, argv);

  auto & comm = atlas::mpi::comm ();

  int irank = comm.rank (), nproc = comm.size ();

  atlas::StructuredGrid grid ("N16");
  atlas::grid::Distribution dist (grid, atlas::util::Config ("type", "equal_regions"));


  printf ("---- START ----\n");

  atlas::functionspace::StructuredColumns fs (grid, dist, atlas::util::Config ("halo", 1));

  printf (" fs.size () = %d,  fs.sizeOwned () = %d\n",  fs.size (),  fs.sizeOwned ());

  atlas::FieldSet xyz = getXYZ (fs);

  gradient (fs, xyz);

  atlas::Library::instance ().finalise ();

  std::cout << "----- STOP -----" << std::endl;

  return 0;
}

