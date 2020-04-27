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

template <typename T>
void 
rotate (const atlas::functionspace::StructuredColumns & fs, 
        const atlas::StructuredGrid & grid, atlas::FieldSet & pgp)
{
  atlas::FieldSet pgpr;

  const auto & grid1 = grid;
  const auto & grid2 = fs.grid ();
  const auto & proj1 = grid1.projection ();
  const auto & proj2 = grid2.projection ();
  
  auto vxy = atlas::array::make_view<double,2> (fs.xy ());

  auto normalize = [] (double coslat, atlas::Projection::Jacobian & jac)
  {
    jac[0][0] /= coslat;
    jac[1][0] /= coslat;

    double n0 = sqrt (jac[0][0] * jac[0][0] + jac[0][1] * jac[0][1]);
    double n1 = sqrt (jac[1][0] * jac[1][0] + jac[1][1] * jac[1][1]);
 
    jac[0][0] /= n0; jac[0][1] /= n0;
    jac[1][0] /= n1; jac[1][1] /= n1;
  };

  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {
      atlas::PointXY xy = atlas::PointXY (vxy (jloc, 0), vxy (jloc, 1));
      atlas::PointLonLat lonlat = proj2.lonlat (xy);

      auto dir1 = proj1.getJacobianAtLonLat (lonlat);
      auto dir2 = proj2.getJacobianAtLonLat (lonlat);

      double coslat = cos (deg2rad * lonlat.lat ());

      normalize (cos (deg2rad * lonlat.lat ()), dir1);
      normalize (cos (deg2rad * lonlat.lat ()), dir2);

      auto inv1 = dir1.inverse ();

      auto mat = dir2 * inv1;

      for (int jfld = 0; jfld < pgp.size () / 2; jfld++)
        {
          T zundef = std::numeric_limits<T>::max ();
          bool llundef = pgp[jfld].metadata ().get ("undef", zundef);
          auto u = atlas::array::make_view<T,1> (pgp[2*jfld+0]);
          auto v = atlas::array::make_view<T,1> (pgp[2*jfld+1]);

          double ux = u (jloc), uy = v (jloc);

          if (llundef)
            {
              if ((ux == zundef) || (uy == zundef))
                {
                  u (jloc) = zundef;
                  v (jloc) = zundef;
                }
              else
                {
                  u (jloc) = mat[0][0] * ux + mat[0][1] * uy;
                  v (jloc) = mat[1][0] * ux + mat[1][1] * uy;
                }
            }
          else
            {
               u (jloc) = mat[0][0] * ux + mat[0][1] * uy;
               v (jloc) = mat[1][0] * ux + mat[1][1] * uy;
            }
        }
      
    }


}

template <typename T>
atlas::FieldSet
gradient (const atlas::functionspace::StructuredColumns & fs, const atlas::FieldSet & pgp)
{
  auto & comm = atlas::mpi::comm ();
  int irank = comm.rank (), nproc = comm.size ();

  atlas::FieldSet pgpg;
  atlas::Field fxy = fs.xy ();

  const atlas::StructuredGrid & grid = fs.grid ();
  const atlas::Projection & proj = grid.projection ();

  fs.haloExchange (pgp);
  fs.haloExchange (fxy);

  for (int jfld = 0; jfld < pgp.size (); jfld++)
    {
      auto f1 = pgp[jfld];

      auto f2x = atlas::Field (f1.name () + ".DX",
                               atlas::array::DataType::kind<T> (), 
                               atlas::array::make_shape (fs.size ()));
      auto f2y = atlas::Field (f1.name () + ".DY",
                               atlas::array::DataType::kind<T> (), 
                               atlas::array::make_shape (fs.size ()));

      f2x.metadata () = f1.metadata ();
      f2y.metadata () = f1.metadata ();

      f2x.rename (f1.name () + ".DX");
      f2y.rename (f1.name () + ".DY");

      pgpg.add (f2x);
      pgpg.add (f2y);

      zero (f2x);
      zero (f2y);
    }

  auto iv = atlas::array::make_view<int,1> (fs.index_i ());
  auto jv = atlas::array::make_view<int,1> (fs.index_j ());

  bool glob  = grid.domain ().global ();

  const auto & xspc = grid.xspace ();

  int nai = std::numeric_limits<int>::max ();

  auto vxy = atlas::array::make_view<double,2> (fxy);

  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {
      int i = iv (jloc) - 1, j = jv (jloc) - 1;

      atlas::PointXY xy = atlas::PointXY (vxy (jloc, 0), vxy (jloc, 1));
      atlas::PointLonLat lonlat = proj.lonlat (xy);

      auto dir = proj.getJacobianAtLonLat (lonlat);
      auto inv = dir.inverse ();

      double xdx = cos (deg2rad * lonlat.lat ()) * inv.dlon_dx (), xdy = inv.dlat_dx ();
      double ydx = cos (deg2rad * lonlat.lat ()) * inv.dlon_dy (), ydy = inv.dlat_dy ();

      double bx = sqrt (xdx * xdx + xdy * xdy);
      double by = sqrt (ydx * ydx + ydy * ydy);

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

      i_w = i - 1;
      i_e = i + 1;

      auto ij_to_index_and_xy = [&] (int i, int j)
      {
        if ((j < fs.j_begin_halo ()) || (fs.j_end_halo () <= j))
          return nai;
        if ((i < fs.i_begin_halo (j)) || (fs.i_end_halo (j) <= i))
          return nai;
        atlas::idx_t idx = fs.index (i, j);
        return idx;
      };

      atlas::idx_t 
           jlocnw = ij_to_index_and_xy (inw, jnw), 
           jlocne = ij_to_index_and_xy (ine, jne),
           jlocsw = ij_to_index_and_xy (isw, jsw), 
           jlocse = ij_to_index_and_xy (ise, jse),
           jloc_w = ij_to_index_and_xy (i_w, j_w),
           jloc_e = ij_to_index_and_xy (i_e, j_e),
           jloc_0 = jloc;

      for (int jfld = 0; jfld < pgp.size (); jfld++)
        {
          T zundef = std::numeric_limits<T>::max ();
          bool llundef = pgp[jfld].metadata ().get ("undef", zundef);


          auto v  = atlas::array::make_view<T,1> (pgp [  jfld  ]);
          auto vx = atlas::array::make_view<T,1> (pgpg[2*jfld+0]);
          auto vy = atlas::array::make_view<T,1> (pgpg[2*jfld+1]);

          T v0 = v (jloc_0);
          T vw = v (jloc_w);
          T ve = v (jloc_e);
          double x0 = vxy (jloc_0, 0);
          double xw = vxy (jloc_w, 0);
          double xe = vxy (jloc_e, 0);

          if (glob && ((xe - xw) < -180.0))
            xw = xw - 360.0;
          if (glob && ((xe - x0) < -180.0))
            x0 = x0 - 360.0;

          double an = (vxy (jloc, 0) - vxy (jlocne, 0)) / (vxy (jlocnw, 0) - vxy (jlocne, 0));
          double as = (vxy (jloc, 0) - vxy (jlocse, 0)) / (vxy (jlocsw, 0) - vxy (jlocse, 0));
   
          double vn = an * v (jlocnw) + (1.0 - an) * v (jlocne);
          double vs = as * v (jlocsw) + (1.0 - as) * v (jlocse);
          T y0 = vxy (jloc_0, 1);
          T yn = vxy (jlocnw, 1);
          T ys = vxy (jlocsw, 1);

          if (glob && (j == 0))
            yn = +180.0 - yn;
          if (glob && (j == grid.ny () - 1))
            ys = -180.0 + ys;

          if (llundef)
            {
              int ifx = 0, ify = 0;

              if (v0 != zundef)
                {
                  if (ve == zundef) { ve = v0; xe = x0; ifx++; }
                  if (vw == zundef) { vw = v0; xw = x0; ifx++; }
                  if (vn == zundef) { vn = v0; yn = y0; ify++; }
                  if (vs == zundef) { vs = v0; ys = y0; ify++; }
                }

              if ((ve == zundef) || (vw == zundef)) ifx = 2;
              if ((vn == zundef) || (vs == zundef)) ify = 2;

              vx (jloc) = ifx == 2 ? zundef : (ve - vw) / (bx * (xe - xw));
              vy (jloc) = ify == 2 ? zundef : (vn - vs) / (by * (yn - ys));
            }
          else
            {
              vx (jloc) = (ve - vw) / (bx * (xe - xw));
              vy (jloc) = (vn - vs) / (by * (yn - ys));
            }

        }


    }

  return pgpg;
}


atlas::field::FieldSetImpl * gradient__
(const atlas::functionspace::detail::StructuredColumns * fs, atlas::field::FieldSetImpl * pgp)
{
  atlas::FieldSet pgpg = gradient<double> (atlas::functionspace::StructuredColumns (fs), atlas::FieldSet (pgp));
  atlas::field::FieldSetImpl * pgpg_ = pgpg.get ();
  pgpg_->attach ();
  return pgpg_;
}

void rotate__
(const atlas::functionspace::detail::StructuredColumns * fs, 
 const atlas::grid::detail::grid::Structured * grid, atlas::field::FieldSetImpl * pgp)
{
  atlas::FieldSet pgp_ (pgp);
  rotate<double> (atlas::functionspace::StructuredColumns (fs), atlas::StructuredGrid (grid), pgp_);
}

