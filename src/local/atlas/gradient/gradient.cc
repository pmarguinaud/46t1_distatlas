#include "gradient.h"

#include "atlas/runtime/Trace.h"
#include "atlas/util/Config.h"
#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"
#include "atlas/runtime/Exception.h"

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

      auto mat = dir1 * dir2.inverse ();

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
  atlas::FieldSet pgpg;

  ATLAS_TRACE_SCOPE ("gradient")
  {

  auto & comm = atlas::mpi::comm ();
  int irank = comm.rank (), nproc = comm.size ();

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


  // Prepare vectors of views

  std::vector<atlas::array::ArrayView<T,1>> pv;
  std::vector<atlas::array::ArrayView<T,1>> pvx;
  std::vector<atlas::array::ArrayView<T,1>> pvy;

  std::vector<bool> lundef;
  std::vector<T> zundef;

  pv .reserve (pgp.size ());
  pvx.reserve (pgp.size ());
  pvy.reserve (pgp.size ());

  lundef.reserve (pgp.size ());
  zundef.reserve (pgp.size ());

  for (int jfld = 0; jfld < pgp.size (); jfld++)
    {
      pv .emplace_back (atlas::array::make_view<T,1> (pgp [  jfld  ]));
      pvx.emplace_back (atlas::array::make_view<T,1> (pgpg[2*jfld+0]));
      pvy.emplace_back (atlas::array::make_view<T,1> (pgpg[2*jfld+1]));

      zundef.push_back (std::numeric_limits<T>::max ());
      lundef.push_back (pgp[jfld].metadata ().get ("undef", zundef.back ()));
    }

  auto iv = atlas::array::make_view<int,1> (fs.index_i ());
  auto jv = atlas::array::make_view<int,1> (fs.index_j ());

  bool glob  = grid.domain ().global ();

  const auto & xspc = grid.xspace ();

  int nai = std::numeric_limits<int>::max ();

  auto vxy = atlas::array::make_view<double,2> (fxy);

#pragma omp parallel for
  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {
      int i = iv (jloc) - 1, j = jv (jloc) - 1;

      atlas::PointXY xy = atlas::PointXY (vxy (jloc, 0), vxy (jloc, 1));
      atlas::PointLonLat lonlat = proj.lonlat (xy);

      auto dir = proj.getJacobianAtLonLat (lonlat);
      auto inv = dir.inverse ();

      double xdx = cos (deg2rad * lonlat.lat ()) * inv.dlon_dx (), xdy = inv.dlat_dx ();
      double ydx = cos (deg2rad * lonlat.lat ()) * inv.dlon_dy (), ydy = inv.dlat_dy ();

      double bx = deg2rad * sqrt (xdx * xdx + xdy * xdy);
      double by = deg2rad * sqrt (ydx * ydx + ydy * ydy);

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
          T zundf = zundef[jfld];
          bool lundf = lundef[jfld];

          auto & v  = pv [jfld];
          auto & vx = pvx[jfld];
          auto & vy = pvy[jfld];

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

          // For global domain, assume x, y in degrees
          if (glob && (j == 0))
            yn = +180.0 - yn;
          if (glob && (j == grid.ny () - 1))
            ys = -180.0 + ys;

          if (lundf)
            {
              int ifx = 0, ify = 0;

              if (v0 != zundf)
                {
                  if (ve == zundf) { ve = v0; xe = x0; ifx++; }
                  if (vw == zundf) { vw = v0; xw = x0; ifx++; }
                  if (vn == zundf) { vn = v0; yn = y0; ify++; }
                  if (vs == zundf) { vs = v0; ys = y0; ify++; }
                }

              if ((ve == zundf) || (vw == zundf)) ifx = 2;
              if ((vn == zundf) || (vs == zundf)) ify = 2;

              vx (jloc) = ifx == 2 ? zundf : (ve - vw) / (bx * (xe - xw));
              vy (jloc) = ify == 2 ? zundf : (vn - vs) / (by * (yn - ys));

              if (v0 != zundf)
                {
                  if (vx (jloc) == zundf) vx (jloc) = 0.0;
                  if (vy (jloc) == zundf) vy (jloc) = 0.0;
                }

            }
          else
            {
              vx (jloc) = (ve - vw) / (bx * (xe - xw));
              vy (jloc) = (vn - vs) / (by * (yn - ys));
            }

        }

    }

  }

  return pgpg;
}


template <typename T>
atlas::FieldSet
halfdiff (const atlas::functionspace::StructuredColumns & fs, const atlas::FieldSet & pgp)
{
  atlas::FieldSet pgpg;

  ATLAS_TRACE_SCOPE ("halfdiff")
  {

  auto & comm = atlas::mpi::comm ();
  int irank = comm.rank (), nproc = comm.size ();

  atlas::Field fxy = fs.xy ();

  const atlas::StructuredGrid & grid = fs.grid ();
  const atlas::Projection & proj = grid.projection ();

  fs.haloExchange (pgp);
  fs.haloExchange (fxy);

  auto t = atlas::array::DataType::kind<T> ();
  auto s = atlas::array::make_shape (fs.size ());
  
  auto addf = [&] (const atlas::Field & f, const std::string & name)
  {
    std::string n = f.name () + name;
    auto g = atlas::Field (n, t, s);
    g.metadata () = f.metadata (); 
    g.rename (n);
    pgpg.add (g);
  };

  for (int jfld = 0; jfld < pgp.size (); jfld++)
    {
      addf (pgp[jfld], ".XP"); addf (pgp[jfld], ".XM");
      addf (pgp[jfld], ".YP"); addf (pgp[jfld], ".YM");
    }


  auto d = atlas::array::DataType::kind<double> (); 

  auto xp = atlas::Field ("XP", d, s); auto vxp = atlas::array::make_view<double,1> (xp);
  auto xm = atlas::Field ("XM", d, s); auto vxm = atlas::array::make_view<double,1> (xm);
  auto yp = atlas::Field ("YP", d, s); auto vyp = atlas::array::make_view<double,1> (yp);
  auto ym = atlas::Field ("YM", d, s); auto vym = atlas::array::make_view<double,1> (ym);

  pgpg.add (xp);
  pgpg.add (xm);
  pgpg.add (yp);
  pgpg.add (ym);

  // Prepare vectors of views

  std::vector<atlas::array::ArrayView<T,1>> pv, pvxp, pvxm, pvyp, pvym;

  std::vector<bool> lundef;
  std::vector<T> zundef;

  pv  .reserve (pgp.size ());
  pvxp.reserve (pgp.size ());
  pvxm.reserve (pgp.size ());
  pvyp.reserve (pgp.size ());
  pvym.reserve (pgp.size ());

  lundef.reserve (pgp.size ());
  zundef.reserve (pgp.size ());

  for (int jfld = 0; jfld < pgp.size (); jfld++)
    {
      pv  .emplace_back (atlas::array::make_view<T,1> (pgp [  jfld  ]));
      pvxp.emplace_back (atlas::array::make_view<T,1> (pgpg[4*jfld+0]));
      pvxm.emplace_back (atlas::array::make_view<T,1> (pgpg[4*jfld+1]));
      pvyp.emplace_back (atlas::array::make_view<T,1> (pgpg[4*jfld+2]));
      pvym.emplace_back (atlas::array::make_view<T,1> (pgpg[4*jfld+3]));

      zundef.push_back (std::numeric_limits<T>::max ());
      lundef.push_back (pgp[jfld].metadata ().get ("undef", zundef.back ()));
    }

  auto iv = atlas::array::make_view<int,1> (fs.index_i ());
  auto jv = atlas::array::make_view<int,1> (fs.index_j ());

  bool glob  = grid.domain ().global ();

  const auto & xspc = grid.xspace ();

  int nai = std::numeric_limits<int>::max ();

  auto vxy = atlas::array::make_view<double,2> (fxy);

#pragma omp parallel for
  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {
      int i = iv (jloc) - 1, j = jv (jloc) - 1;

      atlas::PointXY xy = atlas::PointXY (vxy (jloc, 0), vxy (jloc, 1));
      atlas::PointLonLat lonlat = proj.lonlat (xy);

      auto dir = proj.getJacobianAtLonLat (lonlat);
      auto inv = dir.inverse ();

      double xdx = cos (deg2rad * lonlat.lat ()) * inv.dlon_dx (), xdy = inv.dlat_dx ();
      double ydx = cos (deg2rad * lonlat.lat ()) * inv.dlon_dy (), ydy = inv.dlat_dy ();

      double bx = deg2rad * sqrt (xdx * xdx + xdy * xdy);
      double by = deg2rad * sqrt (ydx * ydx + ydy * ydy);

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

      double x0 = vxy (jloc_0, 0);
      double xw = vxy (jloc_w, 0);
      double xe = vxy (jloc_e, 0);

      if (glob && ((xe - xw) < -180.0))
        xw = xw - 360.0;
      if (glob && ((xe - x0) < -180.0))
        x0 = x0 - 360.0;

      double an = (vxy (jloc, 0) - vxy (jlocne, 0)) / (vxy (jlocnw, 0) - vxy (jlocne, 0));
      double as = (vxy (jloc, 0) - vxy (jlocse, 0)) / (vxy (jlocsw, 0) - vxy (jlocse, 0));
   
      T y0 = vxy (jloc_0, 1);
      T yn = vxy (jlocnw, 1);
      T ys = vxy (jlocsw, 1);

      // For global domain, assume x, y in degrees
      if (glob && (j == 0))
        yn = +180.0 - yn;
      if (glob && (j == grid.ny () - 1))
        ys = -180.0 + ys;

      vxp (jloc) = (bx * (xe - x0));
      vxm (jloc) = (bx * (xw - x0));
      vyp (jloc) = (by * (yn - y0));
      vym (jloc) = (by * (ys - y0));

      for (int jfld = 0; jfld < pgp.size (); jfld++)
        {
          T zundf = zundef[jfld];
          bool lundf = lundef[jfld];

          auto & fv   = pv  [jfld];
          auto & fvxp = pvxp[jfld];
          auto & fvxm = pvxm[jfld];
          auto & fvyp = pvyp[jfld];
          auto & fvym = pvym[jfld];

          T v0 = fv (jloc_0);

          if (lundf && (v0 == zundf))
            {
              fvxp (jloc) = zundf;
              fvxm (jloc) = zundf;
              fvyp (jloc) = zundf;
              fvym (jloc) = zundf;
              continue;
            }

          fvxp (jloc) = 0.0;
          fvxm (jloc) = 0.0;
          fvyp (jloc) = 0.0;
          fvym (jloc) = 0.0;

          T vw = fv (jloc_w);
          T ve = fv (jloc_e);

          double vn = an * fv (jlocnw) + (1.0 - an) * fv (jlocne);
          double vs = as * fv (jlocsw) + (1.0 - as) * fv (jlocse);

          if (lundf)
            {
              if (ve != zundf) fvxp (jloc) = (ve - v0);
              if (vw != zundf) fvxm (jloc) = (vw - v0);
              if (vn != zundf) fvyp (jloc) = (vn - v0);
              if (vs != zundf) fvym (jloc) = (vs - v0);
            }
          else
            {
              fvxp (jloc) = (ve - v0);
              fvxm (jloc) = (vw - v0);
              fvyp (jloc) = (vn - v0);
              fvym (jloc) = (vs - v0);
            }

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

atlas::field::FieldSetImpl * halfdiff__ (atlas::functionspace::detail::StructuredColumns * fs, atlas::field::FieldSetImpl * pgp)
{
  atlas::FieldSet pgpg = halfdiff<double> (atlas::functionspace::StructuredColumns (fs), atlas::FieldSet (pgp));
  atlas::field::FieldSetImpl * pgpg_ = pgpg.get (); 
  pgpg_->attach (); 
  return pgpg_;
}

