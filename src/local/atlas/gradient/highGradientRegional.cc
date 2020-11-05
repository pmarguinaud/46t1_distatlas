#include "highGradientRegional.h"

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

template <int O>
struct coef
{
  static const double * c ()
  {
    static_assert (false, "Coefficients not available");
  }
};

template <> struct coef<1>
{
  static const double * c () { static const double _c[] = {1./2.}; return &_c[0]; }
};

template <> struct coef<2>
{
  static const double * c () { static const double _c[] = {2./3., -1./12.}; return &_c[0]; }
};

template <> struct coef<3>
{
  static const double * c () { static const double _c[] = {3./4., -3./20., 1./60.}; return &_c[0]; }
};

template <> struct coef<4>
{
  static const double * c () { static const double _c[] = {4./5., -1./5., 4./105., -1./280.}; return &_c[0]; }
};

template <typename T, int O, int D>
struct diff
{
  static T d (int i, int j, const atlas::array::ArrayView<T,1> & v, const atlas::functionspace::StructuredColumns & fs)
  {
    T t = 0;
    const double * c = coef<O>::c ();
    for (int k = 0; k < O; k++)
      {
        int l = k + 1;
        int ii0 = D == 0 ? i - l : i, ii1 = D == 0 ? i + l : i;
        int jj0 = D == 0 ? j : j - l, jj1 = D == 0 ? j : j + l;
        t += c[k] * (v (fs.index (ii1, jj1)) - v (fs.index (ii0, jj0)));
      }
    return t;
  }
};

};




template <typename T, int O>
atlas::FieldSet
highGradientRegional (const atlas::functionspace::StructuredColumns & fs, const atlas::FieldSet & pgp)
{
  atlas::FieldSet pgpg;

  ATLAS_TRACE_SCOPE ("highGradientRegional")
  {

  atlas::Field fxy = fs.xy ();

  for (int j = 0; j < fs.grid ().ny (); j++)
    ATLAS_ASSERT (fs.grid ().nx (j) == fs.grid ().nx (0));


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

  pv .reserve (pgp.size ());
  pvx.reserve (pgp.size ());
  pvy.reserve (pgp.size ());

  for (int jfld = 0; jfld < pgp.size (); jfld++)
    {
      pv .emplace_back (atlas::array::make_view<T,1> (pgp [  jfld  ]));
      pvx.emplace_back (atlas::array::make_view<T,1> (pgpg[2*jfld+0]));
      pvy.emplace_back (atlas::array::make_view<T,1> (pgpg[2*jfld+1]));
    }

  auto iv = atlas::array::make_view<int,1> (fs.index_i ());
  auto jv = atlas::array::make_view<int,1> (fs.index_j ());

#pragma omp parallel for
  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {
      int i = iv (jloc) - 1, j = jv (jloc) - 1;

      for (int jfld = 0; jfld < pgp.size (); jfld++)
        {
          auto & v  = pv [jfld];
          auto & vx = pvx[jfld];
          auto & vy = pvy[jfld];

          vx (jloc) = diff<T,O,0>::d (i, j, v, fs);
          vy (jloc) = diff<T,O,1>::d (i, j, v, fs);
        }

    }

  }

  return pgpg;
}

atlas::field::FieldSetImpl * highGradientRegional__
(const atlas::functionspace::detail::StructuredColumns * fs, atlas::field::FieldSetImpl * pgp, const int order)
{ 
  atlas::FieldSet pgpg;
  ATLAS_ASSERT (order < 5);
  switch (order)
    {
      case 1: pgpg = highGradientRegional<double,1> (atlas::functionspace::StructuredColumns (fs), atlas::FieldSet (pgp)); break;
      case 2: pgpg = highGradientRegional<double,2> (atlas::functionspace::StructuredColumns (fs), atlas::FieldSet (pgp)); break;
      case 3: pgpg = highGradientRegional<double,3> (atlas::functionspace::StructuredColumns (fs), atlas::FieldSet (pgp)); break;
      case 4: pgpg = highGradientRegional<double,4> (atlas::functionspace::StructuredColumns (fs), atlas::FieldSet (pgp)); break;
    }
  atlas::field::FieldSetImpl * pgpg_ = pgpg.get ();
  pgpg_->attach ();
  return pgpg_;
}

