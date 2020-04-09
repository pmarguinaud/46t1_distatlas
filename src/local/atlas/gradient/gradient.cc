#include "atlas/util/Config.h"
#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"

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

void zero (atlas::Field & fld)
{
  auto v = atlas::array::make_view<double,1> (fld);
  for (int i = 0; i < fld.size (); i++)
    v (i) = 0.0;
}

};

atlas::FieldSet
gradient (const atlas::functionspace::StructuredColumns & fs, const atlas::FieldSet & pgp)
{
  auto & comm = atlas::mpi::comm ();
  int irank = comm.rank (), nproc = comm.size ();

  atlas::FieldSet pgpg;

  const atlas::StructuredGrid & grid = fs.grid ();

  fs.haloExchange (pgp);

  for (int jfld = 0; jfld < pgp.size (); jfld++)
    {
      auto f1 = pgp[jfld];

      auto f2x = atlas::Field (f1.name () + ".dx",
                               atlas::array::DataType::kind<double> (), 
                               atlas::array::make_shape (fs.size ()));
      auto f2y = atlas::Field (f1.name () + ".dy",
                               atlas::array::DataType::kind<double> (), 
                               atlas::array::make_shape (fs.size ()));

      f2x.metadata () = f1.metadata ();
      f2y.metadata () = f1.metadata ();

      std::string CLSUFF;

      if (f2x.metadata ().get ("CLSUFF", CLSUFF)) f2x.metadata ().set ("CLSUFF", CLSUFF + ".DX");
      if (f2y.metadata ().get ("CLSUFF", CLSUFF)) f2y.metadata ().set ("CLSUFF", CLSUFF + ".DY");

      pgpg.add (f2x);
      pgpg.add (f2y);

      zero (f2x);
      zero (f2y);

    }

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


      printf ("| %8d, %8d | %8d, %8d | %8d, %8d | %8d, %8d | %8d, %8d |\n", 
              i, j, inw, jnw, ine, jne, isw, jsw, ise, jse);


  


    }


  return pgpg;
}


atlas::field::FieldSetImpl * gradient__
(const atlas::functionspace::detail::StructuredColumns * fs, atlas::field::FieldSetImpl * pgp)
{
  atlas::FieldSet pgpg = gradient (atlas::functionspace::StructuredColumns (fs), atlas::FieldSet (pgp));
  atlas::field::FieldSetImpl * pgpg_ = pgpg.get ();
  pgpg_->attach ();
  return pgpg_;
}


