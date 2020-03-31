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

  Config config;
  
  config.set ("nx", 80);
  config.set ("ny", 40);
  config.set ("type", "regular_lonlat");
  config.set ("domain", Config ("type", "rectangular") | Config ("xmin", 0.) | Config ("xmax", 360.) | Config ("ymin", -90.) | Config ("ymax", +90.) | Config ("units", "degrees"));

  StructuredGrid grid1 (config);
  grid::Distribution dist1 = grid::Distribution (grid1, Config ("type", "checkerboard") | Config ("nbands", nproc));

  util::Config periodic_halo;
  periodic_halo.set( "periodic_points", true );
  periodic_halo.set("halo",1);
  functionspace::StructuredColumns fs1 (grid1, dist1, periodic_halo);
  Field field = fs1.createField<double> (option::name ("test"));
  auto view = array::make_view<double,1> (field);
  for (int j = fs1.j_begin_halo (); j < fs1.j_end_halo (); j++){
    for (int i = fs1.i_begin_halo (j); i < fs1.i_end_halo (j); i++){
       gidx_t idx=fs1.index(i,j);
       view(idx)=double(myproc);
    }
  }
  fs1.haloExchange(field);
  printf("FIN DU CAS TEST N>4\n");

  atlas::Library::instance ().finalise ();

  std::cout << "----- STOP -----" << std::endl;



  return 0;
}

