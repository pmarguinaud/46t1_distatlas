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

template <typename T>
T modulo (T a, T b)
{
  a = a - int (a / b) * b;
  return a >= 0 ? a : a + b;
}

};

atlas::FieldSet 
getXYZ (const atlas::functionspace::StructuredColumns & fs)
{
  atlas::FieldSet xyz;

  auto t = atlas::array::DataType::kind<double> ();
  auto s = atlas::array::make_shape (fs.size ());

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
  auto & comm = atlas::mpi::comm ();
  int irank = comm.rank (), nproc = comm.size ();

  atlas::FieldSet pgrad;

  const atlas::StructuredGrid & grid = fs.grid ();

  fs.haloExchange (pgp);

  auto iv = atlas::array::make_view<int,1> (fs.index_i ());
  auto jv = atlas::array::make_view<int,1> (fs.index_j ());

  bool glob  = grid.domain ().global ();

  const auto & xspc = grid.xspace ();

  printf (" grid.ny () = %8d\n", grid.ny ());

  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {
      int i = iv (jloc) - 1, j = jv (jloc) - 1;

      atlas::PointXY xy = grid.xy (i, j);

      if ((j == 0) || (j == grid.ny () - 1))
        continue;

      int jm = j - 1, j0 = j + 0, jp = j + 1;

      auto xj_to_imip = [&] (double x, int j, int & im, int & ip)
      {
        
        int inx = xspc.nx ()[j];
        double dx = xspc.dx ()[j];
        double xmin = xspc.xmin ()[j];
        im = floor ((x - xmin) / dx); 
        ip = im + 1;
        if (! glob)
          {
            if ((im < 0) || (im >= inx)) im = -1;
            if ((ip < 0) || (ip >= inx)) ip = -1;
          }
      };

      int inw, jnw = jm, ine, jne = jm,  
          isw, jsw = jp, ise, jse = jp;

      xj_to_imip (xy.x (), jnw, inw, ine);
      xj_to_imip (xy.x (), jsw, isw, ise);

      auto check_ij = [&] (int i, int j)
      {
        if ((i < 0) || (j < 0))
          return;
        if (j <  fs.j_begin_halo ( )) { printf ("%s:%d\n", __FILE__, __LINE__); abort (); }
        if (j >= fs.j_end_halo   ( )) { printf ("%s:%d\n", __FILE__, __LINE__); abort (); }
        if (i <  fs.i_begin_halo (j)) { printf ("%s:%d\n", __FILE__, __LINE__); abort (); }
        if (i >= fs.i_end_halo   (j)) { printf ("%s:%d\n", __FILE__, __LINE__); abort (); }
      };

      auto ij_to_index_and_xy = [&] (bool w, bool e, int i, int j, atlas::PointXY & xy)
      {
        if ((i < 0) || (j < 0))
          return -1;
        xy = grid.StructuredGrid::xy (i, j);

        if (glob)
          {
            if (w && (i == xspc.nx ()[j]-1)) 
              xy = atlas::PointXY (-xspc.dx ()[j], xy.y ());
            if (e && (i == 0))  // TODO: check xmax
              xy = atlas::PointXY (xspc.xmax ()[j] + xspc.dx ()[j], xy.y ());
          }
        return fs.index (i, j);
      };

      atlas::PointXY xynw, xyne, xysw, xyse;

      atlas::idx_t 
            jlocnw = ij_to_index_and_xy (true, false, jnw, inw, xynw), 
            jlocne = ij_to_index_and_xy (false, true, jne, ine, xyne),
            jlocsw = ij_to_index_and_xy (true, false, jsw, isw, xysw), 
            jlocse = ij_to_index_and_xy (false, true, jse, ise, xyse);


      printf ("| %8d, %8d | %8d, %8d | %8d, %8d | %8d, %8d |\n", 
              i, j, inw, jnw, ine, jne, isw, jsw, ise, jse);


    }


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

  printf (" fs.j_begin () = %8d,  fs.j_end () = %8d\n", fs.j_begin (), fs.j_end ());
  for (int j = fs.j_begin (); j < fs.j_end (); j++)
    printf (" %8d > %8d .. %8d [%8d]\n", j, fs.i_begin (j), fs.i_end (j), grid.nx (j));

  printf (" fs.j_begin_halo () = %8d,  fs.j_end_halo () = %8d\n", fs.j_begin_halo (), fs.j_end_halo ());
  for (int j = fs.j_begin_halo (); j < fs.j_end_halo (); j++)
    printf (" %8d > %8d .. %8d [%8d]\n", j, fs.i_begin_halo (j), fs.i_end_halo (j), grid.nx (j));


  atlas::FieldSet xyz = getXYZ (fs);

  gradient (fs, xyz);

  atlas::Library::instance ().finalise ();

  std::cout << "----- STOP -----" << std::endl;

  return 0;
}

