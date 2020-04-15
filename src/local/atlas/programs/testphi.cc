#include "atlas/library/Library.h"
#include "atlas/util/Config.h"
#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"
#include "atlas/util/Constants.h"

#include <stdlib.h>
#include <iostream>
#include <sstream>  


namespace
{

static double D2R( const double x ) {                                                                              
    return atlas::util::Constants::degreesToRadians() * x;                                                         
}                                                                                                                  
static double R2D( const double x ) {                                                                              
    return atlas::util::Constants::radiansToDegrees() * x;                                                         
}                                                                                                                  


atlas::Projection::Jacobian
getJacobian (const atlas::Projection & proj, const atlas::PointLonLat & lonlat0)
{
  double dlon = 0.0001, dlat = 0.0001;

  double xy0[2] = {lonlat0.lon () + 0.00, lonlat0.lat () + 0.00};
  double xy1[2] = {lonlat0.lon () + dlon, lonlat0.lat () + 0.00};
  double xy2[2] = {lonlat0.lon () + 0.00, lonlat0.lat () + dlat};

  proj.lonlat2xy (xy0);
  proj.lonlat2xy (xy1);
  proj.lonlat2xy (xy2);

  atlas::Projection::Jacobian jac;

  jac[0] = {(xy1[0] - xy0[0]) / dlon, (xy2[0] - xy0[0]) / dlat};
  jac[1] = {(xy1[1] - xy0[1]) / dlon, (xy2[1] - xy0[1]) / dlat};

  return jac;
}

void printJacobian (const atlas::Projection::Jacobian & jac, const std::string & name)
{
  printf (" %s = \n", name.c_str ());
  printf ("  %12.4f, %12.4f | %12.4f\n", jac[0][0], jac[0][1], sqrt (jac[0][0] * jac[0][0] + jac[0][1] * jac[0][1]));
  printf ("  %12.4f, %12.4f | %12.4f\n", jac[1][0], jac[1][1], sqrt (jac[1][0] * jac[1][0] + jac[1][1] * jac[1][1]));
  printf ("  %12.4f, %12.4f\n", sqrt (jac[0][0] * jac[0][0] + jac[1][0] * jac[1][0]), 
                                sqrt (jac[0][1] * jac[0][1] + jac[1][1] * jac[1][1]));
}

void printJacobianRatio (const atlas::Projection::Jacobian & jac1, const atlas::Projection::Jacobian & jac2)
{
  printf ("  %12.4f, %12.4f\n", jac1[0][0] / jac2[0][0], jac1[0][1] / jac2[0][1]);
  printf ("  %12.4f, %12.4f\n", jac1[1][0] / jac2[1][0], jac1[1][1] / jac2[1][1]);
}

};

int main (int argc, char * argv[]) 
{
  std::vector<int> nx = 
  {
    20, 27, 32, 40, 45, 48, 60, 60, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 60, 60, 48, 45, 40, 32, 27, 20
  };

  double stretch = 2.4;
  double centre[2] = {2.0, 46.7};


  std::vector<atlas::grid::Spacing> spacings (nx.size ());

  for (int i = 0; i < nx.size (); i++)
    {   
      double lonmax = 360.0 * double (nx[i] - 1) / double (nx[i]);
      spacings[i] = atlas::grid::Spacing (atlas::util::Config ("type", "linear") | atlas::util::Config ("N", nx[i])
                                        | atlas::util::Config ("start", 0) | atlas::util::Config ("end", lonmax));
    }   

  atlas::StructuredGrid::XSpace xspace (spacings);
  atlas::StructuredGrid::YSpace yspace (atlas::util::Config ("type", "gaussian") | atlas::util::Config ("N", nx.size ()));

  atlas::Projection proj (atlas::util::Config ("type", "rotated_schmidt") 
                        | atlas::util::Config ("stretching_factor", stretch) 
                        | atlas::util::Config ("rotation_angle", 0.0)
                        | atlas::util::Config ("north_pole", std::vector<double>{centre[0], centre[1]})
                        | atlas::util::Config ("arpege", true));



  atlas::StructuredGrid grid (xspace, yspace, proj, atlas::Domain (atlas::util::Config ("type", "global")));


  int jglo = 0;
  for (int j = 0; j < grid.ny ( ); j++)
  for (int i = 0; i < grid.nx (j); i++)
    {
      if (jglo % 100 == 0)
        {
          double lonlat[2], xy[2];
          grid.lonlat (i, j, lonlat);
          auto jacA = proj.getJacobianAtLonLat (atlas::PointLonLat (lonlat[0], lonlat[1]));

          auto jacC = jacA.inverse ();

          auto jacD = jacA * jacC;

          auto jacB = getJacobian (proj, atlas::PointLonLat (lonlat[0], lonlat[1]));


          grid.xy (i, j, xy);
          printf ("------------------\n");
          printf (" %8d > %8d, %8d | %12.4f, %12.4f | %12.4f, %12.4f\n", jglo, i, j, xy[0], xy[2], lonlat[0], lonlat[1]);
          printf (" %12.4f, %12.4f | %12.4f, %12.4f\n", cos (D2R (lonlat[0])), sin (D2R (lonlat[0])), cos (D2R (lonlat[1])), sin (D2R (lonlat[1])));
          printf (" %12.4f, %12.4f | %12.4f, %12.4f\n", cos (D2R (xy    [0])), sin (D2R (xy    [0])), cos (D2R (xy    [1])), sin (D2R (xy    [1])));
          printf ("\n");
          printJacobian (jacA, "jacA");
          printf ("\n");
          printJacobian (jacB, "jacB");
          printf ("\n");
          printJacobian (jacC, "jacC");
          printf ("\n");
          printJacobian (jacD, "jacD");
          printf ("\n");
          printJacobianRatio (jacA, jacB);
          printf ("\n");
        }
      jglo++;
    }






  return 0;
}

