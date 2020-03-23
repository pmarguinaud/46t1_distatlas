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

const double deg2rad = M_PI / 180.0;

int partition2[] = 
{
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 
};

using namespace atlas;
using namespace atlas::grid;
using namespace atlas::mesh;
using namespace atlas::util;
using namespace atlas::array;
using namespace atlas::meshgenerator;
using namespace atlas::output;


namespace
{

template <typename T>
T modulo (T a, T b)
{
  a = a - int (a / b) * b;
  return a >= 0 ? a : a + b;
}

template <typename T, typename I>
T reorder (const T & vec, const I & ord)
{
  T v (vec.size ());
  for (typename I::value_type i = 0; i < vec.size (); i++)
    v[i] = vec[ord[i]];
  return v;
}

template <typename I>
I reverse (const I & ord)
{
  I rev (ord.size ());
  for (typename I::value_type i = 0; i < ord.size (); i++)
    rev[ord[i]] = i;
  return rev;
}


}

struct shuffle4_t
{

  struct shuffle4_recv_t
  {
    int iprc, icnt, ioff;
  };

  struct shuffle4_send_t
  {
    int iprc, icnt, ioff;
    std::vector<gidx_t> iglo;
    std::vector<idx_t> iloc;
  };

  size_t isize_recv;
  size_t isize_send;

  std::vector<int> isort;
  std::vector<shuffle4_recv_t> yl_recv;
  std::vector<shuffle4_send_t> yl_send;
  
  // Should be const
  StructuredGrid grid1;
  StructuredGrid grid2;
  grid::Distribution dist1;
  grid::Distribution dist2;
  size_t size1;
  size_t size2;
  functionspace::StructuredColumns fs1;
  functionspace::StructuredColumns fs2;
  

  static const int ISW, ISE, INW, INE;
};

const int shuffle4_t::ISW = -3, shuffle4_t::ISE = -2, shuffle4_t::INW = -1, shuffle4_t::INE = -0;

struct weights4_t
{
  std::vector<double> values;
};



shuffle4_t
create_shuffle4 (const StructuredGrid & grid1, 
                 const StructuredGrid & grid2,
                 const grid::Distribution & dist1, 
                 const grid::Distribution & dist2,
                 const functionspace::StructuredColumns & fs1,
                 const functionspace::StructuredColumns & fs2)
{
  shuffle4_t shuffle4;

  shuffle4.grid1 = grid1;
  shuffle4.grid2 = grid2;
  shuffle4.dist1 = dist1;
  shuffle4.dist2 = dist2;
  shuffle4.size1 = fs1.sizeOwned ();
  shuffle4.size2 = fs2.sizeOwned ();
  shuffle4.fs1   = fs1;
  shuffle4.fs2   = fs2;

  


  auto & comm = atlas::mpi::comm ();

  int myproc = comm.rank (), nproc = comm.size ();

  std::vector<gidx_t> iglo1all (4 * fs2.sizeOwned ());

  const auto & proj1 = grid1.projection ();
  const auto & xspc1 = grid1.xspace ();
  const auto & yspc1 = grid1.yspace ();

  auto i2 = array::make_view<int,1> (fs2.index_i ());
  auto j2 = array::make_view<int,1> (fs2.index_j ());

// TODO: Use OpenMP on jloc2
  for (int jloc2 = 0; jloc2 < fs2.sizeOwned (); jloc2++)
    {
      PointLonLat lonlat2 = grid2.StructuredGrid::lonlat (i2 (jloc2)-1, j2 (jloc2)-1);
      PointXY xy1 = proj1.xy (lonlat2);

      int iny1 = grid1.ny (), iy1a, iy1b;

      // Search along Y axis
// TODO : increasing Y coordinate
      if (xy1.y () > yspc1.front ())
        iy1a = -1;
      else if (xy1.y () < yspc1.back ())
        iy1a = iny1-1;
      else
        {
          iy1a = 0;
          iy1b = iny1-1;
          while (1)
            {
              // Dichotomy
              int iy1m = (iy1a + iy1b) / 2;
              if ((yspc1[iy1a] >= xy1.y ()) && (xy1.y () >= yspc1[iy1m]))
                iy1b = iy1m;
              else
                iy1a = iy1m;
              if (abs (iy1b - iy1a) <= 1)
                break;
            }
        }
       
      iy1b = iy1a + 1;

      // Search along X axis

      int ix1a, ix1b;

      if (iy1a > -1)
        {
// TODO : handle non global domains (grid1.domain.global () == false)
// TODO : handle shifted longitudes ??
          int inx1 = xspc1.nx ()[iy1a];
          double dx = inx1 * xspc1.dx ()[iy1a];
          int ix1a = modulo (int (inx1 * xy1.x () / dx), inx1);
          int ix1b = modulo (ix1a + 1, inx1);
          iglo1all[4 * (jloc2 + 1) + shuffle4_t::ISW - 1] = grid1.ij2gidx (ix1a, iy1a);
          iglo1all[4 * (jloc2 + 1) + shuffle4_t::ISE - 1] = grid1.ij2gidx (ix1b, iy1a);
        }
      else
        {
          // No points were found
          iglo1all[4 * (jloc2 + 1) + shuffle4_t::ISW - 1] = -1;
          iglo1all[4 * (jloc2 + 1) + shuffle4_t::ISE - 1] = -1;
        }
      
      if (iy1b < iny1)
        {
          int inx1 = xspc1.nx ()[iy1b];
          double dx = inx1 * xspc1.dx ()[iy1b];
          int ix1a = modulo (int (inx1 * xy1.x () / dx), inx1);
          int ix1b = modulo (ix1a + 1, inx1);
          iglo1all[4 * (jloc2 + 1) + shuffle4_t::INW - 1] = grid1.ij2gidx (ix1a, iy1b);
          iglo1all[4 * (jloc2 + 1) + shuffle4_t::INE - 1] = grid1.ij2gidx (ix1b, iy1b);
        }
      else
        {
          // No points were found
          iglo1all[4 * (jloc2 + 1) + shuffle4_t::INW - 1] = 0;
          iglo1all[4 * (jloc2 + 1) + shuffle4_t::INE - 1] = 0;
        }

    }

  std::vector<int> 
       iord_by_glo1 (iglo1all.size ()), // Sort by glo1
       irev_by_glo1 (iglo1all.size ()), // Reverse sort
       ired_by_glo1 (iglo1all.size ()); // Reduction of sorted array

  std::iota (std::begin (iord_by_glo1), std::end (iord_by_glo1), 0);

  std::sort (std::begin (iord_by_glo1), std::end (iord_by_glo1), 
             [&iglo1all] (int a, int b) { return iglo1all[a] < iglo1all[b]; });

  // Reverse indices array

  irev_by_glo1 = reverse (iord_by_glo1);
  
  // Reorder indices of grid #1

  iglo1all = reorder (iglo1all, iord_by_glo1);

  // Reduce list of indices on grid #1 : ired_by_glo1 contains the 
  // indice of the first occurrence of a given value

  ired_by_glo1[0] = 0;
  int inpt1 = 1; // Number of distinct points required by this task on grid #1
  for (int i = 1; i < iglo1all.size (); i++)
    if (iglo1all[i-1] == iglo1all[i])
      ired_by_glo1[i] = ired_by_glo1[i-1];
    else
      ired_by_glo1[i] = inpt1++;
 
  // Build the list of indices of grid #1 (unique indices)

  struct prcglo_t
  {
    int iprc;
    gidx_t iglo;
  };

  std::vector<prcglo_t> prcglo1 (inpt1); // Processor, jglo distinct list

  for (int i = 0; i < iglo1all.size (); i++)
    prcglo1[ired_by_glo1[i]].iglo = iglo1all[i];

  // Find MPI tasks of points on grid #1

  for (int i = 0; i < inpt1; i++)
    if (prcglo1[i].iglo == -1)
      prcglo1[i].iprc = -1;
    else
      prcglo1[i].iprc = dist1.partition (prcglo1[i].iglo);

  // Sort by (MPI task, global index)

  std::vector<int> 
     iord_by_prc1glo1 (inpt1),
     irev_by_prc1glo1 (inpt1);

  std::iota (std::begin (iord_by_prc1glo1), std::end (iord_by_prc1glo1), 0);

  std::sort (std::begin (iord_by_prc1glo1), std::end (iord_by_prc1glo1), 
             [&prcglo1] (int a, int b) { 
                if (prcglo1[a].iprc == prcglo1[b].iprc)
                  return prcglo1[a].iglo < prcglo1[b].iglo; 
                return prcglo1[a].iprc < prcglo1[b].iprc; 
            });

  irev_by_prc1glo1 = reverse (iord_by_prc1glo1);

  // Reorder

  prcglo1 = reorder (prcglo1, iord_by_prc1glo1);

  // Final sort array

  shuffle4.isort = reorder (irev_by_prc1glo1, reorder (ired_by_glo1, irev_by_glo1));

  // Send/recv counts

  std::vector<int> isendcnt (nproc), irecvcnt (nproc);

  int iskip = 0; // Number of points not found; they are supposed to be at the begining of the list

  std::fill (std::begin (irecvcnt), std::end (irecvcnt), 0);

  for (int i = 0; i < prcglo1.size (); i++)
    {
      int iproc = prcglo1[i].iprc;
      if (iproc > -1) 
        irecvcnt[iproc] = irecvcnt[iproc] + 1;
      else
        iskip = iskip + 1;
    }

  // Copy indices to contiguous array for sending

  std::vector<gidx_t> iglobal1 (prcglo1.size ());
  for (int i = 0; i < prcglo1.size (); i++)
    iglobal1[i] = prcglo1[i].iglo;

  prcglo1.clear ();

  // Exchange send/recv counts

  comm.allToAll (irecvcnt, isendcnt);

  int insend = std::count_if (std::begin (isendcnt), std::end (isendcnt), [] (int k) { return k > 0; });
  int inrecv = std::count_if (std::begin (irecvcnt), std::end (irecvcnt), [] (int k) { return k > 0; });

  // Create send/recv descriptors

  shuffle4.yl_send.resize (insend);
  shuffle4.yl_recv.resize (inrecv);
  
  insend = 0;

  for (int iproc = 0; iproc < nproc; iproc++) 
    if (isendcnt[iproc] > 0)
      {
        shuffle4.yl_send[insend].iprc = iproc;
        shuffle4.yl_send[insend].icnt = isendcnt[iproc];
        shuffle4.yl_send[insend].iglo.resize (isendcnt[iproc]);
        insend++;
      }
   
  // Send offsets

  shuffle4.yl_send[0].ioff = 0;
  for (int ii = 1; ii < insend; ii++)
    shuffle4.yl_send[ii].ioff = shuffle4.yl_send[ii-1].ioff + shuffle4.yl_send[ii-1].icnt;


  inrecv = 0;

  for (int iproc = 0; iproc < nproc; iproc++)
    if (irecvcnt[iproc] > 0)
      {
        shuffle4.yl_recv[inrecv].iprc = iproc;
        shuffle4.yl_recv[inrecv].icnt = irecvcnt[iproc];
        inrecv++;
      }

  // Recv offsets

  shuffle4.yl_recv[0].ioff = iskip;
  for (int ii = 1; ii < inrecv; ii++)
    shuffle4.yl_recv[ii].ioff = shuffle4.yl_recv[ii-1].ioff + shuffle4.yl_recv[ii-1].icnt;

// Exchange global indices of grid #1

  std::vector<eckit::mpi::Request> reqsend (insend), reqrecv (inrecv);

// Send indices we need to send values for

  for (int i = 0; i < insend; i++)
// Receive global indices in ILOCAL; they will be translated to local indices later
    reqsend[i] = comm.iReceive (&shuffle4.yl_send[i].iglo[0],
                                shuffle4.yl_send[i].iglo.size (), 
                                shuffle4.yl_send[i].iprc, 101);

  comm.barrier ();

// Receive indices we shall send values for

  for (int i = 0; i < inrecv; i++)
    reqrecv[i] = comm.iSend (&iglobal1[shuffle4.yl_recv[i].ioff], 
                             shuffle4.yl_recv[i].icnt, 
                             shuffle4.yl_recv[i].iprc, 101);

  for (int i = 0; i < insend; i++)
    comm.wait (reqsend[i]);

// We have received the global indices we should transmit; decode them into
// lat/lon indices, then into local indices

  for (int i = 0; i < insend; i++)
    {
      size_t sz = shuffle4.yl_send[i].iglo.size ();
      shuffle4.yl_send[i].iloc.resize (sz);
      for (int j = 0; j < sz; j++)
        {
          gidx_t jglo = shuffle4.yl_send[i].iglo[j];
          idx_t jloc;


          if (jglo < 0)
            {
              jloc = -1;
            } 
          else
            {
              idx_t ij[2];

              jloc = fs1.index (ij[0], ij[1]);

              grid1.gidx2ij (jglo, ij);

              // Check this (i, j) is held by current MPI task
              ATLAS_ASSERT ((fs1.j_begin () <= ij[1]) && (ij[1] < fs1.j_end ()));
              ATLAS_ASSERT ((fs1.i_begin (ij[1]) <= ij[0]) && (ij[0] < fs1.i_end (ij[1])));

              jloc = fs1.index (ij[0], ij[1]);
            }

          shuffle4.yl_send[i].iloc[j] = jloc;
        }
      shuffle4.yl_send[i].iglo.clear ();
    }

  // Wait for MPI requests to complete
  
  for (int i = 0; i < inrecv; i++)
    comm.wait (reqrecv[i]);

  shuffle4.isize_recv = shuffle4.yl_recv[0].ioff;
  for (const auto & r : shuffle4.yl_recv)
    shuffle4.isize_recv += r.icnt;

  shuffle4.isize_send = shuffle4.yl_send[0].ioff;
  for (const auto & r : shuffle4.yl_send)
    shuffle4.isize_send += r.icnt;

  return shuffle4;
}

template <typename T> FieldSet
do_shuffle4 (const shuffle4_t & shuffle4, const FieldSet & pgp1)
{
  FieldSet pgp2e;

  int infld = pgp1.size ();

  // Temporary buffers
  
  atlas::vector<T> 
     zbufr (shuffle4.isize_recv),
     zbufs (shuffle4.isize_send);

  auto & comm = atlas::mpi::comm ();

  int myproc = comm.rank (), nproc = comm.size ();

  // Requests for send/recv

  std::vector<eckit::mpi::Request> 
     reqrecv (shuffle4.yl_recv.size ()), 
     reqsend (shuffle4.yl_send.size ());

  // Check dimensions & type of pgp1 

  for (int jfld = 0; jfld < pgp1.size (); jfld++)
    {
      auto & f = pgp1[jfld];
      if (f.datatype () != atlas::array::DataType::kind<T> ())
        throw_Exception ("Datatype mismatch", Here ());
      if (f.size () < shuffle4.size1)
        throw_Exception ("Field too small", Here ());
    }
  
  std::fill (std::begin (zbufr), std::begin (zbufr) + shuffle4.yl_recv[0].ioff, 0);

  // Post receives

  for (int i = 0; i < shuffle4.yl_recv.size (); i++)
    reqrecv[i] = comm.iReceive (&zbufr[infld*shuffle4.yl_recv[i].ioff], 
                                shuffle4.yl_recv[i].icnt * infld,
                                shuffle4.yl_recv[i].iprc, 101);

  comm.barrier ();

  // Send data

  for (int i = 0; i < shuffle4.yl_send.size (); i++)
    {
      int ioff = shuffle4.yl_send[i].ioff;
      int icnt = shuffle4.yl_send[i].icnt;
      for (int jfld = 0; jfld < infld; jfld++)
        {
          auto v = array::make_view<T,1> (pgp1[jfld]);
          for (int k = 0; k < icnt; k++)
            zbufs[infld*ioff+jfld*icnt+k] = v (k);
        }
      reqsend[i] = comm.iSend (&zbufs[infld*ioff], icnt * infld,
                               shuffle4.yl_send[i].iprc, 101);
    }

  // Create fields in pgp2
  
  for (int jfld = 0; jfld < infld; jfld++)
    {
      pgp2e.add (Field (std::string ("#") + std::to_string (jfld), 
                        atlas::array::DataType::kind<T> (), 
                        atlas::array::make_shape (shuffle4.isize_recv)));
      auto v = array::make_view<T,1> (pgp2e[jfld]);
      for (int k = 0; k < shuffle4.yl_recv[0].ioff; k++)
        v (k) = 0;
    }

  for (auto & req : reqrecv)
    comm.wait (req);

  // Shuffle values in pgp2e
  for (int i = 0; i < shuffle4.yl_recv.size (); i++)
    {
      int ioff = shuffle4.yl_recv[i].ioff;
      int icnt = shuffle4.yl_recv[i].icnt;
      for (int jfld = 0; jfld < infld; jfld++)
        {
          auto v = array::make_view<T,1> (pgp2e[jfld]);
          for (int k = 0; k < icnt; k++)
            v (ioff + k) = zbufr[infld * ioff + jfld * icnt + k];
        }
    }

  for (auto & req : reqsend)
    comm.wait (req);

  return pgp2e;
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

weights4_t
create_weights4 (const shuffle4_t & shuffle4)
{
  weights4_t weights4;

  auto xyz1 = getXYZ (shuffle4.fs1);
  auto xyz2 = getXYZ (shuffle4.fs2);

  FieldSet xyz2e = do_shuffle4<double> (shuffle4, xyz1);

  weights4.values.resize (4 * shuffle4.size2);

  auto x2  = array::make_view<double,1> (xyz2 [0]);
  auto y2  = array::make_view<double,1> (xyz2 [1]);
  auto z2  = array::make_view<double,1> (xyz2 [2]);
  
  auto x2e = array::make_view<double,1> (xyz2e[0]);
  auto y2e = array::make_view<double,1> (xyz2e[1]);
  auto z2e = array::make_view<double,1> (xyz2e[2]);
  
  for (int jloc2 = 0; jloc2 < shuffle4.size2; jloc2++) 
    {
      // The weight is the inverse of the distance in radian between the target point 
      // and the points used for interpolation
      int c[4] = {shuffle4_t::ISW, shuffle4_t::ISE, shuffle4_t::INW, shuffle4_t::INE};
      for (int j = 0; j < 4; j++)
        {
          int jj = c[j];
          int jind1 = shuffle4.isort[4*(jloc2+1)+jj-1];
          // affect a zero value when there is no point
          if (jind1 < 0)
            weights4.values[4*(jloc2+1)+jj-1] = 0;
          else 
            // Prevent division by zero
            weights4.values[4*(jloc2+1)+jj-1] = 1.0 / std::max (1.0E-10, 
                       // Scalar product
                       x2[jloc2] * x2e[jind1] + 
                       y2[jloc2] * y2e[jind1] + 
                       z2[jloc2] * z2e[jind1]);
        }

      // Rebalance weights so that their sum be 1
      
      double s = 0;
      for (int j = 0; j < 4; j++)
        s += weights4.values[4*(jloc2+1)+c[j]-1];

      for (int j = 0; j < 4; j++)
        weights4.values[4*(jloc2+1)+c[j]-1] /= s;
  
    }

  return weights4;
}


int main (int argc, char * argv[]) 
{
  atlas::Library::instance ().initialise (argc, argv);

  auto & comm = atlas::mpi::comm ();

  int myproc = comm.rank (), nproc = comm.size ();

  Config config;
  
  config.set ("nx", 80);
  config.set ("ny", 40);
  config.set ("type", "regular_lonlat");
  config.set ("domain", Config ("type", "rectangular") | Config ("xmin", 0.) | Config ("xmax", 360.) | Config ("ymin", -90.) | Config ("ymax", +90.) | Config ("units", "degrees"));

  StructuredGrid grid1 (config);
  grid::Distribution dist1 = grid::Distribution (grid1, Config ("type", "checkerboard") | Config ("nbands", nproc));

if(0){
  const auto & part1 = dist1.partition ();

  printf ("-- part1 --\n");
  for (int i = 0; i < part1.size (); i++)
    printf (" %8d > %8d\n", i, part1[i]);

  printf ("-- grid1 --\n");

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

  printf (" --- start ---\n");

  const std::vector<int> pl2 = {20, 30, 36, 48, 54, 60, 64, 80, 80, 90, 90, 96, 
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 96, 90, 90, 80, 80, 
    64, 60, 54, 48, 36, 30, 20};

  ReducedGaussianGrid grid2 (pl2);
  grid::Distribution dist2 (nproc, grid2.size (), partition2, 1);

if(0){
  const auto & part2 = dist2.partition ();

  printf ("-- part2 --\n");
  for (int i = 0; i < part2.size (); i++)
    printf (" %8d > %8d\n", i, part2[i]);


  printf ("-- grid2 --\n");
  for (int iy = 0; iy < grid2.ny (); iy++)
  for (int ix = 0; ix < grid2.nx (iy); ix++)
    {
      double lonlat[2];
      grid2.StructuredGrid::lonlat (ix, iy, lonlat);
      printf (" %8d > %12.4f %12.4f\n", iy, lonlat[0], lonlat[1]);
    }


}

  functionspace::StructuredColumns fs1 {grid1, dist1};

  auto i1 = array::make_view<int,1> (fs1.index_i ());
  auto j1 = array::make_view<int,1> (fs1.index_j ());

  Field x1 = fs1.createField<double> (option::name ("x"));
  Field y1 = fs1.createField<double> (option::name ("y"));
  Field z1 = fs1.createField<double> (option::name ("z"));

  // Create xyz fields
  {
    auto vx1 = array::make_view<double,1> (x1);
    auto vy1 = array::make_view<double,1> (y1);
    auto vz1 = array::make_view<double,1> (z1);

    for (int k = 0; k < fs1.sizeOwned (); k++)
      {
        double lonlat[2];
        grid1.StructuredGrid::lonlat (i1 (k)-1, j1 (k)-1, lonlat);
        gidx_t gidx = grid1.ij2gidx (i1 (k)-1, j1 (k)-1);
        double coslon = cos (deg2rad * lonlat[0]), sinlon = sin (deg2rad * lonlat[0]);
        double coslat = cos (deg2rad * lonlat[0]), sinlat = sin (deg2rad * lonlat[0]);
        double x = coslon * coslat, y = sinlon * coslat, z = sinlat;
        vx1 (k) = x; vy1 (k) = y; vz1 (k) = z;
//      printf (" %8d, %8d > %8d > %12.4f, %12.4f\n", i1 (k)-1, j1 (k)-1, gidx, lonlat[1], lonlat[0]);
      }

  }
  
  functionspace::StructuredColumns fs2 {grid2, dist2};

if(0){
  printf (" j_begin, j_end = %8d, %8d\n", fs2.j_begin (), fs2.j_end ());

  for (int j = fs2.j_begin (); j < fs2.j_end (); j++)
    printf (" %8d > i_begin, i_end = %8d, %8d\n", j, fs2.i_begin (j), fs2.i_end (j));

  printf ("-- i,j > index --\n");


  for (int j = fs2.j_begin (); j < fs2.j_end (); j++)
  for (int i = fs2.i_begin (j); i < fs2.i_end (j); i++)
    printf (" %8d, %8d > %8d\n", i, j, fs2.index (i, j));
}
   
  create_shuffle4 (grid1, grid2, dist1, dist2, fs1, fs2);


#ifdef UNDEF

  bool verbose = atoi (argv[1]);
  const int NX = atoi (argv[2]), NY = atoi (argv[3]);
  bool light = atoi (argv[4]);

  grid::Distribution dist;

  if (light)
    dist = grid::Distribution (grid, Config ("light", light) | Config ("blocksize", NX));
  else

  pt ("dist");

if(verbose)
  for (int j = 0, g = 0; j < NY; j++)
  for (int i = 0; i < NX; i++, g++)
    {
      int p = dist.partition (g);
      printf (" %8d -> %8d\n", g, p);
    }
  
  functionspace::StructuredColumns fs {grid, dist, Config ("halo", 40)};

  pt ("fs");

  int k = 0;

if(verbose)
  for (const auto & xy : grid.xy ())
    {
      const PointLonLat ll = grid.projection ().lonlat (xy);
      printf (" %8d > %20.10f, %20.10f | %20.10f, %20.10f\n", k, ll.lon (), ll.lat (), xy.x (), xy.y ());
      k++;
    }

  auto lonlat       = array::make_view<double, 2> (fs.xy ());
  auto global_index = array::make_view<  long, 1> (fs.global_index ());
  auto partition    = array::make_view<   int, 1> (fs.partition ());

  Field field = fs.createField<double> (option::name ("myproc"));

  pt ("field");

if(verbose)
  printf (" fs.size () = %d\n", fs.size ());

if(verbose)
  for (int i = 0; i < fs.size (); i++)
    printf (" %8d > %8d | %c | %20.10f, %20.10f\n", 
            global_index (i)-1, partition (i),  
            partition (i) == myproc ? ' ' : 'G', 
            lonlat (i, LON), lonlat (i, LAT));

  auto view = array::make_view<double, 1> (field);

  for (int i = 0; i < fs.size (); i++)
    if (partition (i) == myproc)
      view (i) = double (myproc);
    else
      view (i) = -1.0;

  field.haloExchange ();
  pt ("halo");

if(verbose)
  printf (" -- field --\n");

if(verbose)
  for (int i = 0; i < fs.size (); i++)
    printf (" %8d > %8d | %c | %20.10f\n", 
            global_index (i)-1, partition (i),  
            partition (i) == myproc ? ' ' : 'G', 
            view (i));

  StructuredGrid grid2 ("N16");
  grid::Distribution dist2 (grid2, Config ("type", "equal_regions"));
//functionspace::StructuredColumns fs2 {grid2, dist2, Config ("halo", 1)};
  functionspace::StructuredColumns fs2 {grid2, dist2};
  Field field2 = fs2.createField<double> (option::name ("myproc1"));

  std::cout << " fs.size () = " << fs.size () << std::endl;
  std::cout << " grid.size () = " << grid.size () << std::endl;

  std::cout << " fs2.size () = " << fs2.size () << std::endl;
  std::cout << " grid2.size () = " << grid2.size () << std::endl;


  Interpolation interpolation_fwd (Config ("type", "structured-bilinear"), fs, fs2);
  interpolation_fwd.execute (field, field2);

  auto view2 = array::make_view<double, 1> (field2);

  for (int i = 0; i < fs2.size (); i++)
    printf (" %8d > %20.10f\n", i, view2 (i));


#endif

  atlas::Library::instance ().finalise ();

  std::cout << "----- STOP -----" << std::endl;



  return 0;
}

