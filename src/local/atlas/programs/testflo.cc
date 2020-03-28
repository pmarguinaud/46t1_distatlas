#include "interpolationA.h"

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
#include <math.h>

using namespace atlas;
using namespace atlas::grid;
using namespace atlas::mesh;
using namespace atlas::util;
using namespace atlas::array;
using namespace atlas::meshgenerator;
using namespace atlas::output;


namespace
{

const double deg2rad = M_PI / 180.0;
const double rad2deg = 180.0 / M_PI;

};


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

static 
double distlonlat (double lon1, double lat1, double lon2, double lat2)
{
  double coslon1  = cos (lon1 * deg2rad), sinlon1  = sin (lon1 * deg2rad);
  double coslat1  = cos (lat1 * deg2rad), sinlat1  = sin (lat1 * deg2rad);
  double x1 = coslon1 * coslat1, y1 = sinlon1 * coslat1, z1 = sinlat1; 

  double coslon2  = cos (lon2 * deg2rad), sinlon2  = sin (lon2 * deg2rad);
  double coslat2  = cos (lat2 * deg2rad), sinlat2  = sin (lat2 * deg2rad);
  double x2 = coslon2 * coslat2, y2 = sinlon2 * coslat2, z2 = sinlat2; 

  return acos (x1 * x2 + y1 * y2 + z1 * z2) * rad2deg;
}

FieldSet
getJGlo (const functionspace::StructuredColumns & fs)
{
  FieldSet jglo;

  auto t = atlas::array::DataType::kind<double> ();
  auto s = atlas::array::make_shape (fs.sizeOwned ());

  jglo.add (Field (std::string ("jglo"), t, s));
  
  auto v = array::make_view<double,1> (jglo[0]);

  auto i = atlas::array::make_view<int,1> (fs.index_i ());
  auto j = atlas::array::make_view<int,1> (fs.index_j ());

 
  const auto & grid = fs.grid ();

  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {
      atlas::gidx_t iglo = grid.ij2gidx (i (jloc)-1, j (jloc)-1);
      v (jloc) = double (iglo);
    }

  return jglo;
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
  StructuredGrid grid2 ("N16");
  grid::Distribution dist2 (grid2, Config ("type", "equal_regions"));


if(0){
  
  printf ("-- grid1 --\n");
  for (int iy = 0, iglo = 0; iy < grid1.ny (); iy++)
  for (int ix = 0; ix < grid1.nx (iy); ix++)
    {
      double lonlat[2];
      grid1.lonlat (ix, iy, lonlat);
      printf (" %8d > %12.4f %12.4f\n", iglo, lonlat[0], lonlat[1]);
      iglo++;
    }
}

if(0){
  
  printf ("-- grid2 --\n");
  for (int iy = 0, iglo = 0; iy < grid2.ny (); iy++)
  for (int ix = 0; ix < grid2.nx (iy); ix++)
    {
      double lonlat[2];
      grid2.lonlat (ix, iy, lonlat);
      printf (" %8d > %12.4f %12.4f\n", iglo, lonlat[0], lonlat[1]);
      iglo++;
    }
}


  comm.barrier ();

  printf ("---- START ----\n");

  functionspace::StructuredColumns fs1 {grid1, dist1};
  functionspace::StructuredColumns fs2 {grid2, dist2};


  interpolationA intA (dist1, fs1, dist2, fs2);

  FieldSet lonlat1 = getLonLat (fs1);
  FieldSet lonlat2 = getLonLat (fs2);

  FieldSet lonlat2e = intA.shuffle<double> (lonlat1);

  {
    auto lon2  = array::make_view<double,1> (lonlat2 [0]);
    auto lat2  = array::make_view<double,1> (lonlat2 [1]);
    auto lon2e = array::make_view<double,1> (lonlat2e[0]);
    auto lat2e = array::make_view<double,1> (lonlat2e[1]);

    for (int jloc2 = 0; jloc2 < fs2.sizeOwned (); jloc2++)
      {
        int icnt = intA.getCnt (jloc2); 
        printf (" %8d > (%12.4f, %12.4f) < %8d\n", jloc2, lon2[jloc2], lat2[jloc2], icnt);

        for (int jj = 0; jj < icnt; jj++)
          {
            int ioff = intA.getOff (jloc2); 
            printf ("         -> (%12.4f, %12.4f) %12.4f \n", 
                    lon2e[ioff+jj], lat2e[ioff+jj], distlonlat (lon2[jloc2], lat2[jloc2],
                    lon2e[ioff+jj], lat2e[ioff+jj]));
          }
      } 
  }

  FieldSet jglo1 = getJGlo (fs1);
  FieldSet jglo2 = getJGlo (fs2);

  FieldSet jglo2e = intA.shuffle<double> (jglo1);

  {
    auto v2  = array::make_view<double,1> (jglo2 [0]);
    auto v2e = array::make_view<double,1> (jglo2e[0]);

if(0)
    for (int jloc2 = 0; jloc2 < fs2.sizeOwned (); jloc2++)
      {
        int icnt = intA.getCnt (jloc2); 
        printf (" %8d > %8d\n", jloc2, int (v2[jloc2]));

        for (int jj = 0; jj < icnt; jj++)
          {
            int ioff = intA.getOff (jloc2); 
            printf ("         -> %8d\n", int (v2e[ioff+jj]));
          }
      } 
  }


  atlas::Library::instance ().finalise ();

  std::cout << "----- STOP -----" << std::endl;

  return 0;
}

