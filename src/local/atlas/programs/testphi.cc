#include "interpolation4.h"

#include "atlas/library/Library.h"
#include "atlas/util/Config.h"
#include "atlas/runtime/Log.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/array.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Trace.h"

#include "atlas/util/CoordinateEnums.h" 
#include "atlas/interpolation.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <string>  
#include <type_traits>
		

#include <sys/time.h>
#include <stdio.h>

const double deg2rad = M_PI / 180.0;
const double rad2deg = 180.0 / M_PI;

void pp (const std::string & file, int line)
{
  auto & comm = atlas::mpi::comm ();
  comm.barrier ();
  std::cout << file << ":" << line << std::endl;
}
#define PP() pp (__FILE__, __LINE__)


using namespace atlas;
using namespace atlas::grid;
using namespace atlas::mesh;
using namespace atlas::util;
using namespace atlas::array;
using namespace atlas::meshgenerator;
using namespace atlas::output;

FieldSet
getLonLat (const functionspace::StructuredColumns & fs)
{
  FieldSet lonlat;

  auto t = atlas::array::DataType::kind<double> ();
  auto s = atlas::array::make_shape (fs.sizeOwned ());

  lonlat.add (Field (std::string ("lon"), t, s));
  lonlat.add (Field (std::string ("lat"), t, s));
  
  auto lon = array::make_view<double,1> (lonlat[0]);
  auto lat = array::make_view<double,1> (lonlat[1]);

  auto i = array::make_view<int,1> (fs.index_i ());
  auto j = array::make_view<int,1> (fs.index_j ());

  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {
      PointLonLat ll = fs.grid ().StructuredGrid::lonlat (i (jloc)-1, j (jloc)-1);
      lon (jloc) = ll.lon ();
      lat (jloc) = ll.lat ();
    }

  return lonlat;
}

FieldSet 
getXYZ (const functionspace::StructuredColumns & fs)
{
  FieldSet xyz;

  auto t = atlas::array::DataType::kind<double> ();
  auto s = atlas::array::make_shape (fs.sizeOwned ());

  xyz.add (Field (std::string ("x"), t, s));
  xyz.add (Field (std::string ("y"), t, s));
  xyz.add (Field (std::string ("z"), t, s));
  
  auto x = array::make_view<double,1> (xyz[0]);
  auto y = array::make_view<double,1> (xyz[1]);
  auto z = array::make_view<double,1> (xyz[2]);

  auto i = array::make_view<int,1> (fs.index_i ());
  auto j = array::make_view<int,1> (fs.index_j ());

  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {
      PointLonLat ll = fs.grid ().StructuredGrid::lonlat (i (jloc)-1, j (jloc)-1);
      double coslon = cos (deg2rad * ll.lon ()), sinlon = sin (deg2rad * ll.lon ());
      double coslat = cos (deg2rad * ll.lat ()), sinlat = sin (deg2rad * ll.lat ());
      x (jloc) = coslon * coslat; 
      y (jloc) = sinlon * coslat; 
      z (jloc) = sinlat;
    }

  return xyz;
}



int main (int argc, char * argv[]) 
{
  atlas::Library::instance ().initialise (argc, argv);

  auto & comm = atlas::mpi::comm ();

  int myproc = comm.rank (), nproc = comm.size ();

  Config config;
  
  StructuredGrid grid1 ("L80x40");
  grid::Distribution dist1 = grid::Distribution (grid1, Config ("type", "checkerboard") | Config ("nbands", nproc));
  ReducedGaussianGrid grid2 ("N16");
  grid::Distribution dist2 (grid2, Config ("type", "equal_regions"));


if(0){
//const auto & part1 = dist1.partition ();
//
//printf ("-- part1 --\n");
//for (int i = 0; i < part1.size (); i++)
//  printf (" %8d > %8d\n", i, part1[i]);
//
//printf ("-- grid1 --\n");

  for (int iy = 0; iy < grid1.ny (); iy++)
  for (int ix = 0; ix < grid1.nx (iy); ix++)
    {
      double lonlat[2];
      grid1.lonlat (ix, iy, lonlat);
      printf (" %8d > %12.4f %12.4f\n", iy, lonlat[0], lonlat[1]);
    }
}

if(0)
{
  int iglo = 0;
  for (int iy = 0; iy < grid1.ny (); iy++)
  for (int ix = 0; ix < grid1.nx (iy); ix++)
    {
      idx_t ij[2];
      grid1.gidx2ij (iglo, ij);
      printf (" ix, iy, iglo = %8d, %8d, %8d, ij = %8d, %8d\n", ix, iy, iglo, ij[0], ij[1]);
      iglo++;
    }
}


  comm.barrier ();

  printf ("---- START ----\n");

  functionspace::StructuredColumns fs1 {grid1, dist1};
  functionspace::StructuredColumns fs2 {grid2, dist2};

  std::string f = std::string ("out2.") + std::to_string (myproc) + ".txt";

  FILE * fp = fopen (f.c_str (), "w");

  interpolation4 int4 (dist1, fs1, dist2, fs2);

  auto xyz1 = getXYZ (fs1);
  auto xyz2 = getXYZ (fs2);

  auto xyz2i = int4.interpolate<double> (xyz1);

  auto x2  = array::make_view<double,1> (xyz2 [0]);
  auto y2  = array::make_view<double,1> (xyz2 [1]);
  auto z2  = array::make_view<double,1> (xyz2 [2]);

  auto x2i = array::make_view<double,1> (xyz2i[0]);
  auto y2i = array::make_view<double,1> (xyz2i[1]);
  auto z2i = array::make_view<double,1> (xyz2i[2]);

  for (int jloc2 = 0; jloc2 < fs2.sizeOwned (); jloc2++)
    {
      double dx = x2  (jloc2) - x2i (jloc2);
      double dy = y2  (jloc2) - y2i (jloc2);
      double dz = z2  (jloc2) - z2i (jloc2);

      double n = sqrt (x2i (jloc2) * x2i (jloc2) + y2i (jloc2) * y2i (jloc2) + z2i (jloc2) * z2i (jloc2));
      x2i (jloc2) /= n;
      y2i (jloc2) /= n;
      z2i (jloc2) /= n;

      double a = acos (std::min (+1.0, std::max (-1.0, x2 (jloc2) * x2i (jloc2) + y2 (jloc2) * y2i (jloc2) + z2 (jloc2) * z2i (jloc2))));

      fprintf (fp, " (%12.4f, %12.4f, %12.4f) (%12.4f, %12.4f, %12.4f) %12.4f %12.4f\n",
              x2  (jloc2), y2  (jloc2), z2  (jloc2),
              x2i (jloc2), y2i (jloc2), z2i (jloc2),
              sqrt (dx * dx + dy * dy + dz * dz) * rad2deg, a * rad2deg);
    }

  fclose (fp);
    
  atlas::Library::instance ().finalise ();

  std::cout << "----- STOP -----" << std::endl;

  return 0;
}

