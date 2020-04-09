#include "gradient.h"

#include "atlas/util/Config.h"
#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"

#include <limits>
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
  atlas::Field fxy = fs.xy ();

  const atlas::StructuredGrid & grid = fs.grid ();



  fs.haloExchange (pgp);
  fs.haloExchange (fxy);

  printf (" pgp.size () = %8d\n", pgp.size ());


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

      printf (" fs.size (), fs.sizeOwned () = %8d, %8d\n", fs.size (), fs.sizeOwned ());

  auto iv = atlas::array::make_view<int,1> (fs.index_i ());
  auto jv = atlas::array::make_view<int,1> (fs.index_j ());

  bool glob  = grid.domain ().global ();

  const auto & xspc = grid.xspace ();

  printf (" glob = %d\n", glob);
  printf (" grid.ny () = %8d\n", grid.ny ());

  for (int j = 0; j < grid.ny (); j++)
    printf (" %8d < %8d\n", j, grid.nx (j));

  printf ("\n");

  printf ("i_begin, i_end\n");
  for (int j = fs.j_begin (); j < fs.j_end (); j++)
    printf ("%8d > %8d, %8d\n", j, fs.i_begin (j), fs.i_end (j));

  printf ("\n");

  printf ("i_begin_halo, i_end_halo\n");
  for (int j = fs.j_begin_halo (); j < fs.j_end_halo (); j++)
    printf ("%8d > %8d, %8d\n", j, fs.i_begin_halo (j), fs.i_end_halo (j));

  printf ("\n");

  printf (" %8s | %8s, %8s | %8s, %8s | %8s, %8s | %8s, %8s | %8s, %8s | %8s, %8s | %8s, %8s | %8s, %8s |\n", 
              "jloc", "i", "j", "inw", "jnw", "ine", "jne", "isw", "jsw", "ise", "jse", "i_w", "j_w", "i_e", "j_e", "jlocw", "jloce");

  int inv = std::numeric_limits<int>::max ();

  auto vxy = atlas::array::make_view<double,2> (fxy);

  for (int j = fs.j_begin_halo (); j < fs.j_end_halo (); j++)
  for (int i = fs.i_begin_halo (j); i < fs.i_end_halo (j); i++)
    {
      atlas::idx_t idx = fs.index (i, j);
      printf (" %8d, %8d > %12.4f, %12.4f\n", i, j, vxy (idx, 0), vxy (idx, 1));
    }

  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {

bool dbg = (irank == 0) && (jloc == 0);

      int i = iv (jloc) - 1, j = jv (jloc) - 1;

      atlas::PointXY xy = atlas::PointXY (vxy (jloc, 0), vxy (jloc, 1));

      int jm = j - 1, j0 = j + 0, jp = j + 1;
      int im = i - 1, i0 = i + 0, ip = i + 1;

      auto xj_to_imip = [&] (double x, int j, int & im, int & ip)
      {
        int js = glob ? std::min (grid.ny () - 1, std::max (0, j)) : j; // Pole
        double dx = xspc.dx ()[js];
        double xmin = xspc.xmin ()[js];
        im = floor ((x - xmin) / dx); 
        ip = im + 1;
      };

      int inw, jnw = jm, ine, jne = jm,  
          isw, jsw = jp, ise, jse = jp,
          i_w, j_w = j0, i_e, j_e = j0;

      xj_to_imip (xy.x (), jnw, inw, ine);
      xj_to_imip (xy.x (), jsw, isw, ise);

if (dbg)
printf (" inw, jnw = %8d, %8d\n", inw, jnw);

      i_w = i - 1;
      i_e = i + 1;

      auto ij_to_index_and_xy = [&] (int i, int j)
      {
        if ((j < fs.j_begin_halo ()) || (fs.j_end_halo () <= j))
          return inv;
        if ((i < fs.i_begin_halo (j)) || (fs.i_end_halo (j) <= i))
          return inv;
        atlas::idx_t idx = fs.index (i, j);
        return idx;
      };

      atlas::idx_t 
           jlocnw = ij_to_index_and_xy (inw, jnw), 
           jlocne = ij_to_index_and_xy (ine, jne),
           jlocsw = ij_to_index_and_xy (isw, jsw), 
           jlocse = ij_to_index_and_xy (ise, jse),
           jloc_w = ij_to_index_and_xy (i_w, j_w),
           jloc_e = ij_to_index_and_xy (i_e, j_e);

if (0)
if (j == 0)
{
  printf (" jloc = %8d\n", jloc);

  printf (" jlocnw, jlocne, jlocsw, jlocse, jloc_w, jloc_e = %8d, %8d, %8d, %8d, %8d, %8d\n", 
            jlocnw, jlocne, jlocsw, jlocse, jloc_w, jloc_e);

  continue;

}

      for (int jfld = 0; jfld < pgp.size (); jfld++)
        {
          auto v  = atlas::array::make_view<double,1> (pgp [  jfld  ]);
          auto vx = atlas::array::make_view<double,1> (pgpg[2*jfld+0]);
          auto vy = atlas::array::make_view<double,1> (pgpg[2*jfld+1]);

          double vw = v (jloc_w);
          double ve = v (jloc_e);
          double xw = vxy (jloc_w, 0);
          double xe = vxy (jloc_e, 0);

          vx (jloc) = (vw - ve) / (xw - xe);

          double an = (vxy (jloc, 0) - vxy (jlocne, 0)) / (vxy (jlocnw, 0) - vxy (jlocne, 0));
          double as = (vxy (jloc, 0) - vxy (jlocse, 0)) / (vxy (jlocsw, 0) - vxy (jlocse, 0));
   
          double vn = an * v (jlocnw) + (1.0 - an) * v (jlocne);
          double vs = as * v (jlocsw) + (1.0 - as) * v (jlocse);
          double yn = vxy (jlocnw, 1);
          double ys = vxy (jlocsw, 1);

          if (glob && (j == 0))
            yn = +180.0 - yn;
          if (glob && (j == grid.ny () - 1))
            ys = -180.0 + ys;

          vy (jloc) = (vn - vs) / (yn - ys);
        }


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


