#include "atlas/library/Library.h"
#include "atlas/util/Config.h"
#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"

#include <stdlib.h>
#include <iostream>
#include <sstream>  


int main (int argc, char * argv[]) 
{
  const int Nx = 129, Ny = 105;
  const double xmax = 17, xmin = -15, ymax = 61, ymin = 35;

  std::vector<atlas::grid::Spacing> spacings (Ny);
    
  for (int i = 0; i < Ny; i++)
    spacings[i] = atlas::grid::Spacing (atlas::util::Config ("type", "linear") | atlas::util::Config ("N", Nx)
                                       | atlas::util::Config ("start", xmin) | atlas::util::Config ("end", xmax));
  
  atlas::StructuredGrid::XSpace xspace (spacings);
  atlas::StructuredGrid::YSpace yspace (atlas::util::Config ("type", "linear") | atlas::util::Config ("N", Ny) 
                                      | atlas::util::Config ("start", ymin) | atlas::util::Config ("end", ymax));
  atlas::Projection proj (atlas::util::Config ("type", "lonlat"));

  atlas::StructuredGrid grid (xspace, yspace, proj, atlas::Domain ());

  printf (" grid.ny () = %8d\n", grid.ny ());  




  return 0;
}

