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

static atlas::Domain domain( const atlas::Grid::Config& grid ) {                                                                      
    atlas::Grid::Config config;                                                                                                
    if ( grid.get( "domain", config ) ) {                                                                               
        return atlas::Domain( config );                                                                                        
    }                                                                                                                   
    return atlas::Domain();                                                                                                    
}                                                                                                                       


}

int main (int argc, char * argv[]) 
{

  atlas::StructuredGrid grid;

  atlas::util::Config config = atlas::util::Config ("nx", 40) | atlas::util::Config ("ny", 20) 
                              | atlas::util::Config ("type", "shifted_lonlat") 
                              | atlas::util::Config ("domain", atlas::util::Config ("type", "global") | atlas::util::Config ("west", -180.));
  

  grid = atlas::StructuredGrid (config);
  
  auto spec = grid.spec ();

  std::cout << spec << std::endl;


  auto ll0 = grid.lonlat (0, 0);
  auto xy0 = grid.xy (0, 0);

  std::cout << " ll0 = " << ll0.lon () << ", " << ll0.lat () << std::endl;
  std::cout << " xy0 = " << xy0.x   () << ", " << xy0.y   () << std::endl;




if (0)
{
    const int nlon = 40, nlat = 20;
    bool shifted_x = false, shifted_y = false;

    double start_x                   = -180 + ( shifted_x ? 0.5 : 0.0 ) * 360.0 / double( nlon );
    std::array<double, 2> interval_x = {start_x, start_x + 360.};
    bool no_endpoint                 = false;
    atlas::StructuredGrid::XSpace xspace( interval_x, std::vector<atlas::idx_t>( nlat, nlon ), no_endpoint );

    std::cout << xspace.xmin ()[0] << std::endl;
    std::cout << xspace.xmax ()[0] << std::endl;


    // spacing is uniform in y
    // If shifted_y, the whole interval is shifted by -dy/2, and last latitude
    // would be -90-dy/2 (below -90!!!), if endpoint=true.
    // Instead, we set endpoint=false so that last latitude is -90+dy/2 instead.
    atlas::StructuredGrid::YSpace yspace( [&] {
        atlas::Grid::Config config_spacing;
        config_spacing.set( "type", "linear" );
        config_spacing.set( "start", 90.0 - ( shifted_y ? 90.0 / double( nlat ) : 0.0 ) );
        config_spacing.set( "end", -90.0 - ( shifted_y ? 90.0 / double( nlat ) : 0.0 ) );
        config_spacing.set( "endpoint", shifted_y ? false : true );
        config_spacing.set( "N", nlat );
        return config_spacing;
    }() );

    atlas::Projection projection;
    atlas::Grid::Config config_projection;
    if ( config.get( "projection", config_projection ) ) {
        projection = atlas::Projection( config_projection );
    }

    std::string name;

    if ( shifted_x and shifted_y ) {
        name = "S";
    }
    else if ( shifted_x and not shifted_y ) {
        name = "Slon";
    }
    else if ( not shifted_x and shifted_y ) {
        name = "Slat";
    }
    else {
        name = "L";
    }

    name += std::to_string( nlon ) + "x" + std::to_string( nlat );

    atlas::StructuredGrid grid (xspace, yspace, projection, domain( config ));

  auto spec = grid.spec ();

  std::cout << spec << std::endl;


  auto ll0 = grid.lonlat (0, 0);
  auto xy0 = grid.xy (0, 0);

  std::cout << " ll0 = " << ll0.lon () << ", " << ll0.lat () << std::endl;
  std::cout << " xy0 = " << xy0.x   () << ", " << xy0.y   () << std::endl;

  std::cout << grid.xmin (0) << std::endl;

  auto xsp = grid.xspace ();

  std::cout << xsp.xmin () << std::endl;

  

}


  return 0;
}

