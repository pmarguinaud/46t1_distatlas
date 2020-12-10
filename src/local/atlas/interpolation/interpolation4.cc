#include "interpolation4.h"
#include "ompsort.h"


#include "atlas/runtime/Exception.h"
#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"
#include "atlas/runtime/Trace.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <string>  
#include <type_traits>

#include <sys/time.h>
#include <stdio.h>

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

template <typename T, typename I>
T reorder (const T & vec, const I & ord)
{
  T v (ord.size ());
#pragma omp parallel for
  for (typename I::value_type i = 0; i < ord.size (); i++)
    v[i] = vec[ord[i]];
  return v;
}

template <typename I>
I reverse (const I & ord)
{
  I rev (ord.size ());
#pragma omp parallel for
  for (typename I::value_type i = 0; i < ord.size (); i++)
    rev[ord[i]] = i;
  return rev;
}

template <typename I, typename J>
void ompiota (I b, I e, J j)
{
#pragma omp parallel for
  for (int i = 0; i < e - b; i++)
    b[i] = i + j;
}

atlas::FieldSet 
getXYZ (const atlas::functionspace::StructuredColumns & fs)
{
  ATLAS_TRACE_SCOPE ("getXYZ")
  {

  atlas::FieldSet xyz;

  auto t = atlas::array::DataType::kind<double> ();
  auto s = atlas::array::make_shape (fs.sizeOwned ());

  xyz.add (atlas::Field (std::string ("x"), t, s));
  xyz.add (atlas::Field (std::string ("y"), t, s));
  xyz.add (atlas::Field (std::string ("z"), t, s));
  
  auto x = atlas::array::make_view<double,1> (xyz[0]);
  auto y = atlas::array::make_view<double,1> (xyz[1]);
  auto z = atlas::array::make_view<double,1> (xyz[2]);

  auto i = atlas::array::make_view<int,1> (fs.index_i ());
  auto j = atlas::array::make_view<int,1> (fs.index_j ());

#pragma omp parallel for
  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {
      atlas::PointLonLat ll = fs.grid ().StructuredGrid::lonlat (i (jloc)-1, j (jloc)-1);
      double coslon = cos (deg2rad * ll.lon ()), sinlon = sin (deg2rad * ll.lon ());
      double coslat = cos (deg2rad * ll.lat ()), sinlat = sin (deg2rad * ll.lat ());
      x (jloc) = coslon * coslat; 
      y (jloc) = sinlon * coslat; 
      z (jloc) = sinlat;
    }

  return xyz;
  }
}


}

interpolation4impl::interpolation4impl 
(const atlas::grid::Distribution & _dist1, const atlas::functionspace::StructuredColumns & _fs1,
 const atlas::grid::Distribution & _dist2, const atlas::functionspace::StructuredColumns & _fs2,
 const bool ldopenmp)
: dist1 (_dist1), dist2 (_dist2), grid1 (_fs1.grid ()), grid2 (_fs2.grid ()), fs1 (_fs1), fs2 (_fs2), llopenmp (ldopenmp)
{
  ATLAS_TRACE_SCOPE ("interpolation4impl::interpolation4impl")
  {

  size1 = fs1.sizeOwned ();
  size2 = fs2.sizeOwned ();

  auto & comm = atlas::mpi::comm ();
  int irank = comm.rank (), nproc = comm.size ();
  const bool llmpi = nproc > 1;

  std::vector<atlas::gidx_t> iglo1all (4 * fs2.sizeOwned ());

  const auto & proj1 = grid1.projection ();
  const auto & xspc1 = grid1.xspace ();
  const auto & yspc1 = grid1.yspace ();

  auto i2 = atlas::array::make_view<int,1> (fs2.index_i ());
  auto j2 = atlas::array::make_view<int,1> (fs2.index_j ());

  bool glob  = grid1.domain ().global ();

  bool yincr = yspc1.front () < yspc1.back ();
  auto lt = [yincr] (const double a, const double b) { return yincr ? a <  b : a >  b; };
  auto le = [yincr] (const double a, const double b) { return yincr ? a <= b : a >= b; };
  auto gt = [yincr] (const double a, const double b) { return yincr ? a >  b : a <  b; };

  ATLAS_TRACE_SCOPE ("Coordinates mapping")
  {
#pragma omp parallel for if (llopenmp)
  for (int jloc2 = 0; jloc2 < fs2.sizeOwned (); jloc2++)
    {

      atlas::PointLonLat lonlat2 = grid2.StructuredGrid::lonlat (i2 (jloc2)-1, j2 (jloc2)-1);
      atlas::PointXY xy1 = proj1.xy (lonlat2);

      int iny1 = grid1.ny (), iy1a, iy1b;

      // Search along Y axis
      if (lt (xy1.y (), yspc1.front ()))
        {
          iy1a = -1;
          iy1b =  0;
        }
      else if (gt (xy1.y (), yspc1.back ()))
        {
          iy1a = iny1-1;
          iy1b = -1;
        }
      else
        {
          iy1a = 0;
          iy1b = iny1-1;
          while (1)
            {
              // Dichotomy
              int iy1m = (iy1a + iy1b) / 2;
              if (le (yspc1[iy1a], xy1.y ()) && le (xy1.y (), yspc1[iy1m]))
                iy1b = iy1m;
              else
                iy1a = iy1m;
              if (abs (iy1b - iy1a) <= 1)
                break;
            }
          iy1b = iy1a + 1;
        }
       
      // Search along X axis

      auto x1iy1_to_iwie = [&] (double x1, int iy1, atlas::gidx_t & iw, atlas::gidx_t & ie)
      {
        int ix1a = -1, ix1b = -1;
        if (iy1 > -1)
          {
            int inx1 = xspc1.nx ()[iy1];
            double dx = xspc1.dx ()[iy1];
            double xmin = xspc1.xmin ()[iy1];
            ix1a = floor ((x1 - xmin) / dx); 
            ix1b = ix1a + 1;
            if (glob)
              {
                ix1a = modulo (ix1a, inx1);
                ix1b = modulo (ix1b, inx1);
              }
            else
              {
                if ((ix1a < 0) || (ix1a >= inx1)) ix1a = -1;
                if ((ix1b < 0) || (ix1b >= inx1)) ix1b = -1;
              }
          }

        iw = (iy1 >= 0) && (ix1a >= 0) ? grid1.index (ix1a, iy1) : -1;
        ie = (iy1 >= 0) && (ix1b >= 0) ? grid1.index (ix1b, iy1) : -1;
      };

      x1iy1_to_iwie (xy1.x (), iy1a, 
                     iglo1all[4 * (jloc2 + 1) + INW - 1], 
                     iglo1all[4 * (jloc2 + 1) + INE - 1]);

      x1iy1_to_iwie (xy1.x (), iy1b, 
                     iglo1all[4 * (jloc2 + 1) + ISW - 1], 
                     iglo1all[4 * (jloc2 + 1) + ISE - 1]);


    }
  }


  std::vector<int> 
       iord_by_glo1 (iglo1all.size ()), // Sort by glo1
       irev_by_glo1 (iglo1all.size ()), // Reverse sort
       ired_by_glo1 (iglo1all.size ()); // Reduction of sorted array

  ATLAS_TRACE_SCOPE ("Sort indices")
  {

  ompiota (std::begin (iord_by_glo1), std::end (iord_by_glo1), 0);

  ompsort (std::begin (iord_by_glo1), std::end (iord_by_glo1), 
           [&iglo1all] (int a, int b) { return iglo1all[a] < iglo1all[b]; }, llopenmp);
   

  // Reverse indices array

  irev_by_glo1 = reverse (iord_by_glo1);
  
  // Reorder indices of grid #1

  iglo1all = reorder (iglo1all, iord_by_glo1);

  };

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
    atlas::gidx_t iglo;
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

  ATLAS_TRACE_SCOPE ("Sort by (MPI task, global index)")
  {

  std::vector<int> 
     iord_by_prc1glo1 (inpt1),
     irev_by_prc1glo1 (inpt1);

  ompiota (std::begin (iord_by_prc1glo1), std::end (iord_by_prc1glo1), 0);

  ompsort (std::begin (iord_by_prc1glo1), std::end (iord_by_prc1glo1), 
           [&prcglo1] (int a, int b) { 
              if (prcglo1[a].iprc == prcglo1[b].iprc)
                return prcglo1[a].iglo < prcglo1[b].iglo; 
              return prcglo1[a].iprc < prcglo1[b].iprc; 
           }, llopenmp);

  irev_by_prc1glo1 = reverse (iord_by_prc1glo1);

  // Reorder

  prcglo1 = reorder (prcglo1, iord_by_prc1glo1);

  // Final sort array

  isort = reorder (irev_by_prc1glo1, reorder (ired_by_glo1, irev_by_glo1));

  };

  // Send/recv counts

  std::vector<int> isendcnt (nproc), irecvcnt (nproc);

  iskip = 0; // Number of points not found; they are supposed to be at the begining of the list

  std::fill (std::begin (irecvcnt), std::end (irecvcnt), 0);

  for (int i = 0; i < prcglo1.size (); i++)
    {
      int iproc = prcglo1[i].iprc;
      if (iproc > -1) 
        irecvcnt[iproc] = irecvcnt[iproc] + 1;
      else
        iskip++;
    }

  // Copy indices to contiguous array for sending

  std::vector<atlas::gidx_t> iglobal1 (prcglo1.size ());
  for (int i = 0; i < prcglo1.size (); i++)
    iglobal1[i] = prcglo1[i].iglo;

  prcglo1.clear ();

  // Exchange send/recv counts

  if (llmpi)
    comm.allToAll (irecvcnt, isendcnt);
  else
    isendcnt[0] = irecvcnt[0];

  int insend = std::count_if (std::begin (isendcnt), std::end (isendcnt), [] (int k) { return k > 0; });
  int inrecv = std::count_if (std::begin (irecvcnt), std::end (irecvcnt), [] (int k) { return k > 0; });

  // Create send/recv descriptors

  yl_send.resize (insend);
  yl_recv.resize (inrecv);
  
  insend = 0;

  for (int iproc = 0; iproc < nproc; iproc++) 
    if (isendcnt[iproc] > 0)
      {
        yl_send[insend].iprc = iproc;
        yl_send[insend].icnt = isendcnt[iproc];
        yl_send[insend].iglo.resize (isendcnt[iproc]);
        insend++;
      }
   
  // Send offsets

  if (insend > 0)
    yl_send[0].ioff = 0;
  for (int ii = 1; ii < insend; ii++)
    yl_send[ii].ioff = yl_send[ii-1].ioff + yl_send[ii-1].icnt;


  inrecv = 0;

  for (int iproc = 0; iproc < nproc; iproc++)
    if (irecvcnt[iproc] > 0)
      {
        yl_recv[inrecv].iprc = iproc;
        yl_recv[inrecv].icnt = irecvcnt[iproc];
        inrecv++;
      }

  // Recv offsets

  if (inrecv > 0)
    yl_recv[0].ioff = iskip;
  for (int ii = 1; ii < inrecv; ii++)
    yl_recv[ii].ioff = yl_recv[ii-1].ioff + yl_recv[ii-1].icnt;

  // Exchange global indices of grid #1


  std::vector<eckit::mpi::Request> reqsend (insend), reqrecv (inrecv);

  if (llmpi)
    {
      // Send indices we need to send values for
      for (int i = 0; i < insend; i++)
        reqsend[i] = comm.iReceive (&yl_send[i].iglo[0],
                                    yl_send[i].iglo.size (), 
                                    yl_send[i].iprc, 101);

      comm.barrier ();

      // Receive indices we shall send values for

      for (int i = 0; i < inrecv; i++)
        reqrecv[i] = comm.iSend (&iglobal1[yl_recv[i].ioff], 
                                 yl_recv[i].icnt, 
                                 yl_recv[i].iprc, 101);

      for (int i = 0; i < insend; i++)
        comm.wait (reqsend[i]);

      for (int i = 0; i < inrecv; i++)
        comm.wait (reqrecv[i]);

    }
  else
    {
      for (int j = 0; j < yl_send[0].iglo.size (); j++)
        yl_send[0].iglo[j] = iglobal1[yl_recv[0].ioff+j];
    }

  // We have received the global indices we should transmit; decode them into
  // lat/lon indices, then into local indices

  for (int i = 0; i < insend; i++)
    {
      size_t sz = yl_send[i].iglo.size ();
      yl_send[i].iloc.resize (sz);
      for (int j = 0; j < sz; j++)
        {
          atlas::gidx_t jglo = yl_send[i].iglo[j];
          atlas::idx_t jloc;

          if (jglo < 0)
            {
              jloc = -1;
            } 
          else
            {
              atlas::idx_t ij[2];

              grid1.index2ij (jglo, ij[0], ij[1]);
              jloc = fs1.index (ij[0], ij[1]);

              // Check this (i, j) is held by current MPI task
              ATLAS_ASSERT ((fs1.j_begin () <= ij[1]) && (ij[1] < fs1.j_end ()));
              ATLAS_ASSERT ((fs1.i_begin (ij[1]) <= ij[0]) && (ij[0] < fs1.i_end (ij[1])));

              jloc = fs1.index (ij[0], ij[1]);
            }
          yl_send[i].iloc[j] = jloc;
        }
      yl_send[i].iglo.clear ();
    }

  
  if (inrecv > 0)
    isize_recv = yl_recv[0].ioff;
  for (const auto & r : yl_recv)
    isize_recv += r.icnt;

  isize_send = insend > 0 ? yl_send[0].ioff : 0;
  for (const auto & r : yl_send)
    isize_send += r.icnt;

  create_weights ();
  }
}

template <typename T> atlas::FieldSet
interpolation4impl::shuffle (const atlas::FieldSet & pgp1) const
{
  ATLAS_TRACE_SCOPE ("interpolation4impl::shuffle")
  {

  atlas::FieldSet pgp2e;

  int infld = pgp1.size ();

  // Temporary buffers
  
  atlas::vector<T> 
     zbufr (isize_recv * infld),
     zbufs (isize_send * infld);

  auto & comm = atlas::mpi::comm ();
  int irank = comm.rank (), nproc = comm.size ();
  const bool llmpi = nproc > 1;

  // Requests for send/recv

  std::vector<eckit::mpi::Request> 
     reqrecv (yl_recv.size ()), 
     reqsend (yl_send.size ());


  // Check dimensions & type of pgp1 

  for (int jfld = 0; jfld < pgp1.size (); jfld++)
    {
      auto & f = pgp1[jfld];
      if (f.datatype () != atlas::array::DataType::kind<T> ())
        atlas::throw_Exception ("Datatype mismatch", Here ());
      if (f.size () < size1)
        atlas::throw_Exception ("Field too small", Here ());
    }

  
  std::fill (std::begin (zbufr), std::begin (zbufr) + iskip, 0);

  if (llmpi)
    {
      // Post receives

      for (int i = 0; i < yl_recv.size (); i++)
        reqrecv[i] = comm.iReceive (&zbufr[infld*yl_recv[i].ioff], 
                                    yl_recv[i].icnt * infld,
                                    yl_recv[i].iprc, 101);

      comm.barrier ();
    }

  // Send data

  for (int i = 0; i < yl_send.size (); i++)
    {
      int ioff = yl_send[i].ioff;
      int icnt = yl_send[i].icnt;
// TODO: Collapse loops
      for (int jfld = 0; jfld < infld; jfld++)
        {
          auto v = atlas::array::make_view<T,1> (pgp1[jfld]);
#pragma omp parallel for if (llopenmp)
          for (int k = 0; k < icnt; k++)
            {
if (k < 0) abort ();
if (k >= pgp1[jfld].size ()) abort ();
              zbufs[infld*ioff+jfld*icnt+k] = v (yl_send[i].iloc[k]);
            }
        }

      if (llmpi)
        reqsend[i] = comm.iSend (&zbufs[infld*ioff], icnt * infld,
                                 yl_send[i].iprc, 101);
      else
        {
          for (int j = 0; j < infld * icnt; j++)
            zbufr[infld*yl_recv[0].ioff+j] = zbufs[j];
        }
    }


  // Create fields in pgp2
  
// TODO: Collapse loops
  for (int jfld = 0; jfld < infld; jfld++)
    {
      auto f1  = pgp1[jfld];
      auto f2e = atlas::Field (f1.name (),
                               atlas::array::DataType::kind<T> (), 
                               atlas::array::make_shape (isize_recv));
      f2e.metadata () = f1.metadata ();
      pgp2e.add (f2e);
      auto v = atlas::array::make_view<T,1> (f2e);
#pragma omp parallel for if (llopenmp)
      for (int k = 0; k < yl_recv[0].ioff; k++)
        v (k) = 0;
    }

  if (llmpi)
    for (auto & req : reqrecv)
      comm.wait (req);

  // Shuffle values in pgp2e
  for (int i = 0; i < yl_recv.size (); i++)
    {
      int ioff = yl_recv[i].ioff;
      int icnt = yl_recv[i].icnt;
// TODO: Collapse loops
      for (int jfld = 0; jfld < infld; jfld++)
        {
          auto v = atlas::array::make_view<T,1> (pgp2e[jfld]);
#pragma omp parallel for if (llopenmp)
          for (int k = 0; k < icnt; k++)
            v (ioff + k) = zbufr[infld * ioff + jfld * icnt + k];
        }
    }

  if (llmpi)
    for (auto & req : reqsend)
      comm.wait (req);

  return pgp2e;
  }
}

void
interpolation4impl::create_weights ()
{
  ATLAS_TRACE_SCOPE ("interpolation4impl::create_weights")
  {

  auto xyz1 = getXYZ (fs1);
  auto xyz2 = getXYZ (fs2);

  atlas::FieldSet xyz2e = shuffle<double> (xyz1);

  weights4.values.resize (4 * size2);

  auto x2  = atlas::array::make_view<double,1> (xyz2 [0]);
  auto y2  = atlas::array::make_view<double,1> (xyz2 [1]);
  auto z2  = atlas::array::make_view<double,1> (xyz2 [2]);
  
  auto x2e = atlas::array::make_view<double,1> (xyz2e[0]);
  auto y2e = atlas::array::make_view<double,1> (xyz2e[1]);
  auto z2e = atlas::array::make_view<double,1> (xyz2e[2]);
  
  ATLAS_TRACE_SCOPE ("Compute weights");
#pragma omp parallel for if (llopenmp)
  for (int jloc2 = 0; jloc2 < size2; jloc2++) 
    {
      // The weight is the inverse of the distance in radian between the target point 
      // and the points used for interpolation


      int c[4] = {ISW, ISE, INW, INE};
      for (int j = 0; j < 4; j++)
        {
          int jj = c[j];
          int jind1 = isort[4*(jloc2+1)+jj-1];
          // affect a zero value when there is no point
          if (jind1 < iskip)
            weights4.values[4*(jloc2+1)+jj-1] = 0;
          else 
            // Prevent division by zero
            weights4.values[4*(jloc2+1)+jj-1] = 1.0 / std::max (1.0E-10, 
              acos (x2[jloc2] * x2e[jind1] + y2[jloc2] * y2e[jind1] + z2[jloc2] * z2e[jind1])); // Scalar product

        }

      // Rebalance weights so that their sum be 1
      
      double s = 0;
      for (int j = 0; j < 4; j++)
        s += weights4.values[4*(jloc2+1)+c[j]-1];

      for (int j = 0; j < 4; j++)
        weights4.values[4*(jloc2+1)+c[j]-1] /= s;
    }
  };
}


template <typename T> atlas::FieldSet
interpolation4impl::interpolate (const atlas::FieldSet & pgp1) const
{
  ATLAS_TRACE_SCOPE ("interpolation4impl::interpolate")
  {

  atlas::FieldSet pgp2;

  int infld = pgp1.size ();

  size_t size2 = fs2.sizeOwned ();

  auto pgp2e = shuffle<T> (pgp1);
  
  for (int jfld = 0; jfld < infld; jfld++)
    {
      auto f1 = pgp1[jfld];
      auto f2 = atlas::Field (f1.name (),
                              atlas::array::DataType::kind<T> (), 
                              atlas::array::make_shape (fs2.size ()));
      f2.metadata () = f1.metadata ();
      pgp2.add (f2);
    }

// TODO: Collapse loops
  for (int jfld = 0; jfld < infld; jfld++)
    {
      T zundef = std::numeric_limits<T>::max ();
      bool llundef = pgp1[jfld].metadata ().get ("undef", zundef);
      auto v2  = atlas::array::make_view<T,1> (pgp2 [jfld]);
      auto v2e = atlas::array::make_view<T,1> (pgp2e[jfld]);

#pragma omp parallel for if (llopenmp)
      for (int jloc2 = 0; jloc2 < size2; jloc2++)
        {
          int kisw = 4*(jloc2+1)+ISW-1, jisw = isort[kisw];
          int kise = 4*(jloc2+1)+ISE-1, jise = isort[kise];
          int kinw = 4*(jloc2+1)+INW-1, jinw = isort[kinw];
          int kine = 4*(jloc2+1)+INE-1, jine = isort[kine];

          if (llundef)
            {
              T vs = static_cast<T> (0);
              double ws = 0.0;
              auto add = [&] (T v, double w) { if (v != zundef) { vs = vs + w * v; ws = ws + w; } };
              add (v2e[jisw], weights4.values[kisw]);
              add (v2e[jise], weights4.values[kise]);
              add (v2e[jinw], weights4.values[kinw]);
              add (v2e[jine], weights4.values[kine]);
              v2 (jloc2) = ws > 0.0 ? vs / ws : zundef;
            }
          else
            {
              v2 (jloc2) = 
                weights4.values[kisw] * v2e[jisw] + weights4.values[kise] * v2e[jise] + 
                weights4.values[kinw] * v2e[jinw] + weights4.values[kine] * v2e[jine];
            }
        }
    }

  return pgp2;
  }
}


#define DEF(T) \
template atlas::FieldSet interpolation4impl::interpolate<T>  (const atlas::FieldSet &) const; 

DEF (double);



interpolation4impl * interpolation4__new 
  (const atlas::grid::DistributionImpl * dist1, const atlas::functionspace::detail::StructuredColumns * fs1,
   const atlas::grid::DistributionImpl * dist2, const atlas::functionspace::detail::StructuredColumns * fs2,
   const int ldopenmp)
{
  interpolation4impl * int4 = new interpolation4impl
  (atlas::grid::Distribution (dist1), atlas::functionspace::StructuredColumns (fs1), 
   atlas::grid::Distribution (dist2), atlas::functionspace::StructuredColumns (fs2),
   ldopenmp > 0);
  return int4;
}

atlas::field::FieldSetImpl * interpolation4__interpolate 
(interpolation4impl * This , atlas::field::FieldSetImpl * pgp1)
{
  atlas::FieldSet pgp2 = This->interpolate<double> (atlas::FieldSet (pgp1));
  atlas::field::FieldSetImpl * pgp2_ = pgp2.get ();
  pgp2_->attach ();
  return pgp2_;
}

