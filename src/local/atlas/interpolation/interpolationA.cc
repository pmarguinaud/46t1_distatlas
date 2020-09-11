#include "interpolationA.h"
#include "ompsort.h"


#include "eckit/mpi/Comm.h"
#include "atlas/runtime/Trace.h"
#include "atlas/runtime/Exception.h"
#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <vector>
#include <map>
#include <string>  
#include <type_traits>
		

#include <sys/time.h>
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
  for (int i = 0; i < e - b; i++)
    b[i] = i + j;
}

};

interpolationAimpl::interpolationAimpl 
(const atlas::grid::Distribution & _dist1, const atlas::functionspace::StructuredColumns & _fs1,
 const atlas::grid::Distribution & _dist2, const atlas::functionspace::StructuredColumns & _fs2)
: dist1 (_dist1), dist2 (_dist2), grid1 (_fs1.grid ()), grid2 (_fs2.grid ()), fs1 (_fs1), fs2 (_fs2)
{
  ATLAS_TRACE_SCOPE ("interpolationAimpl::interpolationAimpl")
  {

  size_t size1 = fs1.sizeOwned ();
  size_t size2 = fs2.sizeOwned ();

  auto & comm = atlas::mpi::comm ();

  int myproc = comm.rank (), nproc = comm.size ();
  const bool llmpi = nproc > 1;

  std::vector<atlas::gidx_t> iglo1all (4 * fs2.sizeOwned ());

  const auto & proj2 = grid2.projection ();
  const auto & xspc2 = grid2.xspace ();
  const auto & yspc2 = grid2.yspace ();

  auto i1 = atlas::array::make_view<int,1> (fs1.index_i ());
  auto j1 = atlas::array::make_view<int,1> (fs1.index_j ());

  struct prcglo_t
  {
    atlas::gidx_t iglo1;  // Global indice on grid #1
    int iprc2;            // Task to send to on grid #2 (for averaging)
    atlas::gidx_t iglo2;  // Global indice on grid #2 (for averaging)
    atlas::idx_t  iloc1;  // Local indice on grid #1
  };

  std::vector<prcglo_t> prcglo1 (size1);

  bool glob  = grid2.domain ().global ();

  bool yincr = yspc2.front () < yspc2.back ();
  auto lt = [yincr] (const double a, const double b) { return yincr ? a <  b : a >  b; };
  auto le = [yincr] (const double a, const double b) { return yincr ? a <= b : a >= b; };
  auto gt = [yincr] (const double a, const double b) { return yincr ? a >  b : a <  b; };

  int iny2 = grid2.ny ();

  double dy2f = (yspc2[     1] - yspc2[     0]) / 2.0f;
  double dy2b = (yspc2[iny2-1] - yspc2[iny2-2]) / 2.0f;


  ATLAS_TRACE_SCOPE ("Coordinates mapping") 
  {

#pragma omp parallel for
  for (int jloc1 = 0; jloc1 < size1; jloc1++)
    {
      //  For all points of grid #1 (only in our region), we find the nearest point of grid #2. Here,
      //  "nearest" means nearest latitude and nearest longitude. We also find the task
      //  which holds this point of grid #2

      atlas::PointLonLat lonlat1 = grid1.StructuredGrid::lonlat (i1 (jloc1)-1, j1 (jloc1)-1);
      atlas::gidx_t iglo1 = grid1.ij2gidx (i1 (jloc1)-1, j1 (jloc1)-1);
      atlas::PointXY xy2 = proj2.xy (lonlat1);

      int iy2 = -1;

      // Search along Y axis
      if (lt (xy2.y (), yspc2[0]))
        {
          if (gt (xy2.y (), yspc2[0] - dy2f) || glob)
            iy2 = 0;
        }
      else if (gt (xy2.y (), yspc2[iny2-1]))
        {
          if (lt (xy2.y (), yspc2[iny2-1] + dy2b) || glob)
            iy2 = iny2-1;
        }
      else
        {
          int iy2a = 0, iy2b = iny2-1;
          while (1)
            {
              // Dichotomy
              int iy2m = (iy2a + iy2b) / 2;
              if (le (yspc2[iy2a], xy2.y ()) && le (xy2.y (), yspc2[iy2m]))
                iy2b = iy2m;
              else
                iy2a = iy2m;
              if (abs (iy2b - iy2a) <= 1)
                break;
            }
          if (fabs (yspc2[iy2a] - xy2.y ()) < fabs (yspc2[iy2b] - xy2.y ()))
            iy2 = iy2a;
          else
            iy2 = iy2b;
        }
       

      // Search along X axis
      
      int ix2 = -1;

      if (iy2 >= 0)
        {
          int inx2 = grid2.nx (iy2);
          double dx = xspc2.dx ()[iy2];
          double xmin = xspc2.xmin ()[iy2];
          ix2 = round ((xy2.x () - xmin) / dx);

          if (glob)
            {
              ix2 = modulo (ix2, inx2);
            }
          else
            {
               if ((ix2 < 0) || (ix2 >= inx2)) ix2 = -1;
            }
        }

      atlas::gidx_t iglo2 = (iy2 >= 0) && (ix2 >= 0) ? grid2.ij2gidx (ix2, iy2) : -1;

      int iprc2 = iglo2 >= 0 ? dist2.partition (iglo2) : -1;

      prcglo1[jloc1].iglo1 = iglo1;
      prcglo1[jloc1].iglo2 = iglo2;
      prcglo1[jloc1].iprc2 = iprc2;
      prcglo1[jloc1].iloc1 = jloc1;

    }

  }


  // Compute & exchange send/recv counts
  std::vector<size_t> isendcnt (nproc, 0), irecvcnt (nproc), isendoff (nproc);
  
  int iskip = 0; // Number of points on grid #1 that we will not use
  for (int jloc1 = 0; jloc1 < size1; jloc1++)
    if (prcglo1[jloc1].iprc2 >= 0)
      isendcnt[prcglo1[jloc1].iprc2]++;
    else
      iskip++;

  isendoff[0] = 0;
  for (int iproc = 1; iproc < nproc; iproc++)
    isendoff[iproc] = isendoff[iproc-1] + isendcnt[iproc-1];

  if (llmpi)
    comm.allToAll (isendcnt, irecvcnt);
  else
    irecvcnt[0] = isendcnt[0];

  ATLAS_TRACE_SCOPE ("Sort points by task, global indice")
  {

  std::vector<int> iord_by_prc2glo2 (prcglo1.size ());

  ompiota (std::begin (iord_by_prc2glo2), std::end (iord_by_prc2glo2), 0);

  if (llmpi)
     ompsort (std::begin (iord_by_prc2glo2), std::end (iord_by_prc2glo2), 
              [&prcglo1] (int a, int b) 
              { 
                if (prcglo1[a].iprc2 == prcglo1[b].iprc2)
                  return prcglo1[a].iglo2 < prcglo1[b].iglo2;
                return prcglo1[a].iprc2 < prcglo1[b].iprc2;
           });
  else
     std::stable_sort (std::begin (iord_by_prc2glo2), std::end (iord_by_prc2glo2), 
              [&prcglo1] (int a, int b) 
              { 
                if (prcglo1[a].iprc2 == prcglo1[b].iprc2)
                  return prcglo1[a].iglo2 < prcglo1[b].iglo2;
                return prcglo1[a].iprc2 < prcglo1[b].iprc2;
           });

  prcglo1 = reorder (prcglo1, iord_by_prc2glo2);

  }

  int insend = std::count_if (std::begin (isendcnt), std::end (isendcnt), [] (size_t k) { return k > 0; });
  int inrecv = std::count_if (std::begin (irecvcnt), std::end (irecvcnt), [] (size_t k) { return k > 0; });

  std::vector<eckit::mpi::Request> reqsend (insend), reqrecv (inrecv);

  struct exch_t
  {
    size_t size;
    std::vector<atlas::gidx_t> iglo1iglo2; 
  };

  std::vector<exch_t> yl_exch_recv (inrecv), yl_exch_send (insend);

  // Post receives

  yl_recv.resize (inrecv);

  inrecv = 0;
  for (int iproc = 0; iproc < nproc; iproc++)
    if (irecvcnt[iproc])
      {
        yl_recv[inrecv].iproc = iproc;
        yl_exch_recv[inrecv].size = irecvcnt[iproc];
        yl_exch_recv[inrecv].iglo1iglo2.resize (2 * irecvcnt[iproc]);

        if (llmpi)
          reqrecv[inrecv] = comm.iReceive (&yl_exch_recv[inrecv].iglo1iglo2[0],
                                          yl_exch_recv[inrecv].iglo1iglo2.size (),
                                          iproc, 101);
          
        inrecv++;
      }

  if (llmpi)
    comm.barrier ();

  // Send

  yl_send.resize (insend);

  insend = 0;
  for (int iproc = 0; iproc < nproc; iproc++)
    if (isendcnt[iproc])
      {
        yl_exch_send[insend].size = isendcnt[iproc];
        yl_exch_send[insend].iglo1iglo2.resize (2 * isendcnt[iproc]);

        yl_send[insend].iproc = iproc;
        yl_send[insend].isize = isendcnt[iproc];
        yl_send[insend].iloc.resize (isendcnt[iproc]);
     


        // List of points of grid #1 to send to iproc
#pragma omp parallel for
        for (int i = 0; i < isendcnt[iproc]; i++)
          {
            int jind1 = iskip + isendoff[iproc] + i;
            yl_send[insend].iloc[i]                = prcglo1[jind1].iloc1;
            yl_exch_send[insend].iglo1iglo2[2*i+0] = prcglo1[jind1].iglo1;
            yl_exch_send[insend].iglo1iglo2[2*i+1] = prcglo1[jind1].iglo2;
          }

        if (llmpi)
          reqsend[insend] = comm.iSend (&yl_exch_send[insend].iglo1iglo2[0],
                                       yl_exch_send[insend].iglo1iglo2.size (),
                                       iproc, 101);
        else
          yl_exch_recv[0].iglo1iglo2 = yl_exch_send[0].iglo1iglo2;

        insend++;
      }

  prcglo1.clear ();

  isize_send = 
    std::accumulate (std::begin (yl_send), std::end (yl_send),
                     0, [] (int a, const interpolationAimpl::send_t & s)
                     {
                       return a + s.isize;
                     });
  if (llmpi)
    {
      for (int i = 0; i < inrecv; i++)
        comm.wait (reqrecv[i]);

      // TODO : move wait send further
      for (int i = 0; i < insend; i++)
        comm.wait (reqsend[i]);
    }

  for (auto & yl_exch : yl_exch_send)
    yl_exch.iglo1iglo2.clear ();

  ATLAS_TRACE_SCOPE ("Prepare receive descriptors")
  {
#pragma omp parallel for
  for (int ii = 0; ii < yl_recv.size (); ii++) 
    {
      // Total number of points we get from this task
      yl_recv[ii].isize = yl_exch_recv[ii].iglo1iglo2.size () / 2;

      int inloc = 0; // Number of local indices we will process using points from this remote proc

      atlas::gidx_t jglo2 = -1;
      for (int jj = 0; jj < yl_exch_recv[ii].size; jj++)
        {
          atlas::gidx_t iglo2 = yl_exch_recv[ii].iglo1iglo2[2*jj+1];
          if (iglo2 != jglo2)
            {
              inloc++;
              jglo2 = iglo2;
            }
        }
        
      // Count number of different remote points for each local point on grid #2
      
      yl_recv[ii].desc.resize (inloc);

      jglo2 = -1;
      inloc = 0;
      for (int jj = 0; jj < yl_exch_recv[ii].size; jj++)
        {
          atlas::gidx_t iglo2 = yl_exch_recv[ii].iglo1iglo2[2*jj+1];
          if (iglo2 != jglo2)
            {

              atlas::idx_t ij[2];

              grid2.gidx2ij (iglo2, ij);
              atlas::idx_t jloc2 = fs2.index (ij[0], ij[1]);

              // Check this (i, j) is held by current MPI task

              ATLAS_ASSERT ((fs2.j_begin () <= ij[1]) && (ij[1] < fs2.j_end ()));
              ATLAS_ASSERT ((fs2.i_begin (ij[1]) <= ij[0]) && (ij[0] < fs2.i_end (ij[1])));

              yl_recv[ii].desc[inloc].iloc = jloc2;
              yl_recv[ii].desc[inloc].irem_off = jj;

              inloc++;
              jglo2 = iglo2;
            }
        }

      int isize = yl_recv[ii].desc.size ();
      for (int jj = 0; jj < isize-1; jj++)
        yl_recv[ii].desc[jj].irem_cnt = yl_recv[ii].desc[jj+1].irem_off 
                                      - yl_recv[ii].desc[jj+0].irem_off;

      yl_recv[ii].desc[isize-1].irem_cnt 
                                      = yl_recv[ii].isize 
                                      - yl_recv[ii].desc[isize-1].irem_off;
      
    }
  }

  // Total number of points received by this proc
  isize_recv = std::accumulate (std::begin (yl_recv), std::end (yl_recv), 
                                0, [] (int a, const interpolationAimpl::recv_t & r) { return a + r.isize; });

  desc.resize (fs2.sizeOwned ());

  for (auto & c : desc)
    c.icnt = 0;

  for (int ii = 0; ii < yl_recv.size (); ii++)
  for (int jj = 0; jj < yl_recv[ii].desc.size (); jj++)
    {
      atlas::idx_t jloc2 = yl_recv[ii].desc[jj].iloc;
      // Number of points for jloc2
      desc[jloc2].icnt += yl_recv[ii].desc[jj].irem_cnt;
    }

  // Offset of points for jloc2
  desc[0].ioff = 0;
  for (int jloc2 = 1; jloc2 < fs2.sizeOwned (); jloc2++)
    desc[jloc2].ioff = desc[jloc2-1].ioff + desc[jloc2-1].icnt;


  // Compute sort vector
  
  struct jloc2iglo1_t
  {
    atlas::idx_t  jloc2;
    atlas::gidx_t iglo1;
  };

  std::vector<jloc2iglo1_t> isortk (isize_recv); 

  for (int ii = 0, jl = 0; ii < yl_recv.size (); ii++)
    {
#pragma omp parallel for
      for (int jj = 0; jj < yl_recv[ii].desc.size (); jj++)
        {
          int jloc2 = yl_recv[ii].desc[jj].iloc         ;
          int ioff  = yl_recv[ii].desc[jj].irem_off + jl;
          int icnt  = yl_recv[ii].desc[jj].irem_cnt     ;
          for (int k = 0; k < icnt; k++)
            isortk[ioff+k].jloc2 = jloc2;
        }
      jl = jl + yl_recv[ii].isize;
    }

  for (int ii = 0, jl = 0; ii < yl_exch_recv.size (); ii++)
    {
      size_t isize = yl_exch_recv[ii].iglo1iglo2.size () / 2;
#pragma omp parallel for
      for (int jj = 0; jj < isize; jj++)
        isortk[jl+jj].iglo1 = yl_exch_recv[ii].iglo1iglo2[2*jj+0];
      jl = jl + isize;
    }

  ATLAS_TRACE_SCOPE ("Sort by local index, remote index")
  {
    std::vector<int> iord (isize_recv);
    ompiota (std::begin (iord), std::end (iord), 0);

    ompsort (std::begin (iord), std::end (iord), [&isortk] (int a, int b)
             {
               if (isortk[a].jloc2 == isortk[b].jloc2) 
                 return isortk[a].iglo1 < isortk[b].iglo1;
               return isortk[a].jloc2 < isortk[b].jloc2;
             });

    isort = reverse (iord);
  }

  for (int jloc2 = 0; jloc2 < fs2.sizeOwned (); jloc2++)  
    if (desc[jloc2].icnt == 0)
      {
        isize_miss++;
        if ((isize_miss < 20) && verbose)
          std::cerr << " shuffle : no points were found for local " << jloc2 << std::endl;
      }

  if ((isize_miss > 0) && verbose)
    std::cerr << " shuffle : no points were for " << isize_miss << " points " << std::endl;


  }

}

template <typename T>
atlas::FieldSet interpolationAimpl::shuffle (const atlas::FieldSet & pgp1) const
{
  ATLAS_TRACE_SCOPE ("interpolationAimpl::shuffle")
  {

  atlas::FieldSet pgp2e;

  auto & comm = atlas::mpi::comm ();
  const int myproc = comm.rank ();
  const bool llmpi = comm.size () > 1;

  const int infld = pgp1.size ();
  const int insend = yl_send.size ();
  const int inrecv = yl_recv.size ();

  // Create fields in pgp2e

  for (int jfld = 0; jfld < infld; jfld++)
    {
      auto f1  = pgp1[jfld];
      auto f2e = atlas::Field (f1.name (),
                               atlas::array::DataType::kind<T> (), 
                               atlas::array::make_shape (isize_recv));
      f2e.metadata () = f1.metadata ();
      pgp2e.add (f2e);
    }

  std::vector<eckit::mpi::Request> reqsend (insend), reqrecv (inrecv);

  atlas::vector<T> zbufr (infld * isize_recv), zbufs (infld * isize_send);

  for (int ii = 0, ioffr = 0; ii < inrecv; ii++)
    {
      int iproc = yl_recv[ii].iproc;
      size_t isize = yl_recv[ii].isize;
      if (llmpi)
        reqrecv[ii] = comm.iReceive (&zbufr[ioffr], isize * infld, iproc, 101);
      ioffr = ioffr + infld * isize;
    }

  if (llmpi)
    comm.barrier ();

  std::vector<size_t> ioffs_all (insend);

  for (int ii = 0, ioffs = 0; ii < insend; ii++)
    {
      size_t isize = yl_send[ii].isize;
      ioffs_all[ii] = ioffs;
      ioffs += infld * isize;
    }


  std::vector<atlas::array::ArrayView<T,1>> fv1, fv2e;

  fv1 .reserve (pgp1 .size ());
  fv2e.reserve (pgp2e.size ());

  for (auto & f : pgp1)
    fv1 .emplace_back (atlas::array::make_view<T,1> (f));
  for (auto & f : pgp2e)
    fv2e.emplace_back (atlas::array::make_view<T,1> (f));


  for (int ii = 0; ii < insend; ii++)
    {
      size_t ioffs = ioffs_all[ii];
      int iproc = yl_send[ii].iproc;
      size_t isize = yl_send[ii].isize;
#pragma omp parallel for collapse (2)
      for (int jfld = 0; jfld < infld; jfld++)
        for (int jj = 0; jj < isize; jj++)
          {
            auto & v = fv1[jfld];
            zbufs[ioffs+jfld*isize+jj] = v (yl_send[ii].iloc[jj]);
          }
      if (llmpi)
        reqsend[ii] = comm.iSend (&zbufs[ioffs], infld * isize, iproc, 101);
      else
        zbufr.assign (&zbufs[0], &zbufs[0] + infld * isize);
    }

  if (llmpi)
    {
// TODO: Use MPI_WAITANY
      for (int i = 0; i < inrecv; i++)
        comm.wait (reqrecv[i]);
    }

  std::vector<size_t> ioffr_all (inrecv);

  for (int ii = 0, ioffr = 0; ii < inrecv; ii++)
    {
      ioffr_all[ii] = ioffr;
      size_t isize = yl_recv[ii].isize;
      ioffr += isize;
    }

  for (int ii = 0; ii < inrecv; ii++)
    {
      size_t ioffr = ioffr_all[ii];
      size_t isize = yl_recv[ii].isize;

#pragma omp parallel for collapse (2)
      for (int jfld = 0; jfld < infld; jfld++)
        for (int jj = 0; jj < isize; jj++)
          {
            auto & v = fv2e[jfld];
            v (isort[ioffr+jj]) = zbufr[ioffr*infld+jfld*isize+jj];
          }
    }

  if (llmpi)
    {
      for (int i = 0; i < insend; i++)
        comm.wait (reqsend[i]);
    }

  return pgp2e;
  }
}

template <typename T, typename O, typename E>
void interpolationAimpl::reduce (atlas::array::ArrayView<T,1> & v2, atlas::array::ArrayView<T,1> & v2e, 
                                 const size_t size2, T t0, T zundef, E eval, O op) const
{
#pragma omp parallel for
  for (int jloc2 = 0; jloc2 < size2; jloc2++)
    {
      int ioff = getOff (jloc2); 
      int icnt = getCnt (jloc2); 

      int jcnt = 0;
      T t = t0;
  
      for (int jj = 0; jj < icnt; jj++)
        if (v2e (ioff+jj) != zundef)
          {
            t = op (t, v2e (ioff+jj));
            jcnt++;
          }
  
      v2 (jloc2) = eval (jcnt, t, zundef);  
    } 
}

template <typename T> atlas::FieldSet
interpolationAimpl::interpolate (const atlas::FieldSet & pgp1, const opt_t opt) const
{
  ATLAS_TRACE_SCOPE ("interpolationAimpl::interpolate")
  {

  atlas::FieldSet pgp2;
  atlas::FieldSet pgp2e = shuffle<T> (pgp1);

  size_t size2 = fs2.sizeOwned ();

  int infld = pgp1.size ();

  for (int jfld = 0; jfld < infld; jfld++)
    {
      auto f1 = pgp1[jfld];
      auto f2 = atlas::Field (f1.name (),
                              atlas::array::DataType::kind<T> (), 
                              atlas::array::make_shape (fs2.size ()));
      f2.metadata () = f1.metadata ();
      T zundef = std::numeric_limits<T>::max ();
      if (! f2.metadata ().get ("undef", zundef))
        f2.metadata ().set ("undef", zundef);
      pgp2.add (f2);
    }

  for (int jfld = 0; jfld < infld; jfld++)
    {
      T zundef = std::numeric_limits<T>::max ();
      bool llundef = pgp1[jfld].metadata ().get ("undef", zundef);
      auto v2  = atlas::array::make_view<T,1> (pgp2 [jfld]);
      auto v2e = atlas::array::make_view<T,1> (pgp2e[jfld]);

      switch (opt)
        {
          case opt_t::OPT_AVG:
            reduce (v2, v2e, size2, static_cast<T> (0), zundef, 
                   [] (int jcnt, T t, T zundef) { return jcnt > 0 ? t / static_cast<T> (jcnt) : zundef; },
                   [] (T t1, T t2) { return t1 + t2; });
          break;
          case opt_t::OPT_MIN:
            reduce (v2, v2e, size2, std::numeric_limits<T>::max (), zundef, 
                   [] (int jcnt, T t, T zundef) { return jcnt > 0 ? t : zundef; },
                   [] (T t1, T t2) { return std::min (t1, t2); });
          break;
          case opt_t::OPT_MAX:
            reduce (v2, v2e, size2, std::numeric_limits<T>::min (), zundef, 
                   [] (int jcnt, T t, T zundef) { return jcnt > 0 ? t : zundef; },
                   [] (T t1, T t2) { return std::max (t1, t2); });
          break;
        }
    } 

  return pgp2;
  }
}

#define DEF(T) \
template atlas::FieldSet interpolationAimpl::interpolate<T>  (const atlas::FieldSet &, const opt_t) const; 

DEF (double);


interpolationAimpl * interpolationA__new 
  (const atlas::grid::DistributionImpl * dist1, const atlas::functionspace::detail::StructuredColumns * fs1,
   const atlas::grid::DistributionImpl * dist2, const atlas::functionspace::detail::StructuredColumns * fs2)
{
  interpolationAimpl * intA = new interpolationAimpl
  (atlas::grid::Distribution (dist1), atlas::functionspace::StructuredColumns (fs1), 
   atlas::grid::Distribution (dist2), atlas::functionspace::StructuredColumns (fs2));
  return intA;
}

atlas::field::FieldSetImpl * interpolationA__interpolate 
(interpolationAimpl * This , atlas::field::FieldSetImpl * pgp1, const int opt)
{
  atlas::FieldSet pgp2 = This->interpolate<double> (atlas::FieldSet (pgp1), static_cast<const interpolationAimpl::opt_t> (opt));
  atlas::field::FieldSetImpl * pgp2_ = pgp2.get ();
  pgp2_->attach ();
  return pgp2_;
}

atlas::field::FieldSetImpl * interpolationA__shuffle 
(interpolationAimpl * This , atlas::field::FieldSetImpl * pgp1)
{
  atlas::FieldSet pgp2 = This->shuffle<double> (atlas::FieldSet (pgp1));
  atlas::field::FieldSetImpl * pgp2_ = pgp2.get ();
  pgp2_->attach ();
  return pgp2_;
}

int interpolationA__getlen (const interpolationAimpl * This)
{
  return This->getLen ();
}

void interpolationA__getcnt (const interpolationAimpl * This, int cnt[])
{
#pragma omp parallel for
  for (size_t i = 0; i < This->getLen (); i++)
    cnt[i] = This->getCnt (i);
}

void interpolationA__getoff (const interpolationAimpl * This, int off[])
{
#pragma omp parallel for
  for (size_t i = 0; i < This->getLen (); i++)
    off[i] = This->getOff (i);
}
