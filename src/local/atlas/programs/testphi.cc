#include "atlas/projection.h"
#include "atlas/grid.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"

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

void normalize (double coslat, atlas::Projection::Jacobian & jac)
{
  jac[0][0] /= coslat;
  jac[1][0] /= coslat;

  double n0 = sqrt (jac[0][0] * jac[0][0] + jac[0][1] * jac[0][1]);
  double n1 = sqrt (jac[1][0] * jac[1][0] + jac[1][1] * jac[1][1]);

  jac[0][0] /= n0; jac[0][1] /= n0;
  jac[1][0] /= n1; jac[1][1] /= n1;
}

const double deg2rad = M_PI / 180.0;
const double rad2deg = 180.0 / M_PI;

template <typename SKIP>
void doTaylorTest (const atlas::StructuredGrid & grid, double dmax, SKIP skip)
{
  const auto & proj = grid.projection ();

  for (int j = 0, jglo = 0; j < grid.ny ( ); j++)
  for (int i = 0; i < grid.nx (j); i++, jglo++)
  if (! skip (grid, i, j))
    {
      double lonlat[2];
      grid.lonlat (i, j, lonlat);

      auto jacA = proj.getJacobianAtLonLat (atlas::PointLonLat (lonlat[0], lonlat[1]));
      auto jacB = getJacobian (proj, atlas::PointLonLat (lonlat[0], lonlat[1]));

if (j == 0)
  {
    printf (" i, j = %8d, %8d\n", i, j);
    normalize (cos (deg2rad * lonlat[1]), jacB);
    printJacobian (jacB, "jacB");

    printf ("\n");
  }

      auto jacC = jacB * jacA.inverse () - atlas::Projection::Jacobian::Id ();

      double diff = jacC.norm ();

if (0)
      if (diff > dmax)
        {
          printf ("------------------\n");
          printf (" i, j = %8d, %8d\n", i, j);
          printJacobian (jacA, "jacA");
          printf ("\n");
          printJacobian (jacB, "jacB");
          printf ("\n");
          printf ("%12.4e\n", diff);
          printf ("\n");
        }

//    EXPECT (diff < dmax);
    }

}


void doTaylorTest (const atlas::StructuredGrid & grid, double dmax)
{
  auto noskip = [] (const atlas::StructuredGrid & grid, int i, int j) { return false; };
  doTaylorTest (grid, dmax, noskip);
}

int main (int argc, char * argv[])
{

  const std::vector<int> nx =
  { 
    20, 27, 32, 40, 45, 48, 60, 60, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 60, 60, 48, 45, 40, 32, 27, 20
  };
  
  const double stretch = 2.4, centre[2] = {2.0, 46.7};
  
  std::vector<atlas::grid::Spacing> spacings (nx.size ());
  
  for (int i = 0; i < nx.size (); i++)
    {   
      double lonmax = 360.0 * double (nx[i] - 1) / double (nx[i]);
      spacings[i] = atlas::grid::Spacing (atlas::util::Config ("type", "linear") | atlas::util::Config ("N", nx[i])
                                        | atlas::util::Config ("start", 0) | atlas::util::Config ("end", lonmax));
    }
  
  atlas::StructuredGrid::XSpace xspace (spacings);
  atlas::StructuredGrid::YSpace yspace (atlas::util::Config ("type", "gaussian") | atlas::util::Config ("N", nx.size ()));
  
  auto
  proj = atlas::Projection (atlas::util::Config ("type", "rotated_schmidt") 
                          | atlas::util::Config ("stretching_factor", stretch)
                          | atlas::util::Config ("rotation_angle", 0.0)
                          | atlas::util::Config ("north_pole", std::vector<double>{centre[0], centre[1]})
                          | atlas::util::Config ("arpege", true));
  

  doTaylorTest (atlas::StructuredGrid (xspace, yspace, proj, atlas::Domain (atlas::util::Config ("type", "global"))), 1e-3);

  return 0;

}

