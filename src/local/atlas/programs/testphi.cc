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

#include "atlas/util/CoordinateEnums.h" //to have LON LAT
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

using namespace atlas;
using namespace atlas::grid;
using namespace atlas::mesh;
using namespace atlas::util;
using namespace atlas::array;
using namespace atlas::meshgenerator;
using namespace atlas::output;


int main (int argc, char * argv[]) 
{
  atlas::Library::instance ().initialise (argc, argv);
  auto & comm = atlas::mpi::comm ();
  int myproc = comm.rank ();
  int nproc = comm.size ();

  const int Nx = 64, Ny = 64, Nux = 53, Nuy = 53;

  const double DxInMetres = 50000.;
  const double DyInMetres = 50000.;
  const double LaDInDegrees = 46.2;
  const double Latin1InDegrees = 46.2;
  const double Latin2InDegrees = 46.2;
  const double LoVInDegrees = 2.0;
  const double XMinInMetres = -Nux / 2 * DxInMetres;
  const double YMinInMetres = -Nuy / 2 * DyInMetres;

  std::vector<Spacing> spacings (Ny);

  for (int i = 0; i < Ny; i++)
    spacings[i] = Spacing (Config ("type", "linear") | Config ("N", Nx) 
                         | Config ("start", XMinInMetres) | Config ("end", XMinInMetres + (Nx - 1) * DxInMetres));

  StructuredGrid::XSpace xspace (spacings);
  StructuredGrid::YSpace yspace (Config ("type", "linear") | Config ("N", Ny) | Config ("start", YMinInMetres) | Config ("end", YMinInMetres + (Ny - 1 ) * DyInMetres));
  Projection proj (Config ("type", "lambert_conformal_conic") | Config ("longitude0", LoVInDegrees)
                 | Config ("latitude0", LaDInDegrees) | Config ("latitude1", Latin1InDegrees) 
                 | Config ("latitude2", Latin2InDegrees));

  atlas::StructuredGrid grid (xspace, yspace, proj, Domain ());

  for (int j = 0, jglo = 0; j < grid.ny (); j++)
  for (int i = 0; i < grid.nx (j); i++, jglo++)
    {
      double lonlat[2];
      grid.lonlat (i, j, lonlat);
      printf (" i, j = %8d, %8d > %12.4f, %12.4f\n", i, j, lonlat[0], lonlat[1]);
    }


  return 0;
}

