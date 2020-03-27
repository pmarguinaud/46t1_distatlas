#include "atlas/library/Library.h"
#include "atlas/util/Config.h"
#include "atlas/runtime/Log.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/array.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Trace.h"

#include "atlas/util/CoordinateEnums.h" 
#include "atlas/interpolation.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <string>  
#include <type_traits>
		

#include <sys/time.h>
#include <stdio.h>
#include <math.h>

void pp (const std::string & file, int line)
{
  auto & comm = atlas::mpi::comm ();
  comm.barrier ();
  std::cout << file << ":" << line << std::endl;
}
#define PP() pp (__FILE__, __LINE__)


using namespace atlas;
using namespace atlas::grid;
using namespace atlas::mesh;
using namespace atlas::util;
using namespace atlas::array;
using namespace atlas::meshgenerator;
using namespace atlas::output;


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

template <typename T, typename IDX>
T reorder (const T & vec, const IDX & ord)
{
  T v (ord.size ());
  for (typename IDX::value_type i = 0; i < ord.size (); i++)
    v[i] = vec[ord[i]];
  return v;
}

template <typename IDX>
IDX reverse (const IDX & ord)
{
  IDX rev (ord.size ());
  for (typename IDX::value_type i = 0; i < ord.size (); i++)
    rev[ord[i]] = i;
  return rev;
}

};

FieldSet
getLonLat (const functionspace::StructuredColumns & fs)
{
  FieldSet lonlat;

  auto t = atlas::array::DataType::kind<double> ();
  auto s = atlas::array::make_shape (fs.sizeOwned ());

  lonlat.add (Field (std::string ("lon"), t, s));
  lonlat.add (Field (std::string ("lat"), t, s));
  
  auto lon = array::make_view<double,1> (lonlat[0]);
  auto lat = array::make_view<double,1> (lonlat[1]);

  auto i = array::make_view<int,1> (fs.index_i ());
  auto j = array::make_view<int,1> (fs.index_j ());

  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {
      PointLonLat ll = fs.grid ().StructuredGrid::lonlat (i (jloc)-1, j (jloc)-1);
      lon (jloc) = ll.lon ();
      lat (jloc) = ll.lat ();
    }

  return lonlat;
}

FieldSet 
getXYZ (const functionspace::StructuredColumns & fs)
{
  FieldSet xyz;

  auto t = atlas::array::DataType::kind<double> ();
  auto s = atlas::array::make_shape (fs.sizeOwned ());

  xyz.add (Field (std::string ("x"), t, s));
  xyz.add (Field (std::string ("y"), t, s));
  xyz.add (Field (std::string ("z"), t, s));
  
  auto x = array::make_view<double,1> (xyz[0]);
  auto y = array::make_view<double,1> (xyz[1]);
  auto z = array::make_view<double,1> (xyz[2]);

  auto i = array::make_view<int,1> (fs.index_i ());
  auto j = array::make_view<int,1> (fs.index_j ());

  for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
    {
      PointLonLat ll = fs.grid ().StructuredGrid::lonlat (i (jloc)-1, j (jloc)-1);
      double coslon = cos (deg2rad * ll.lon ()), sinlon = sin (deg2rad * ll.lon ());
      double coslat = cos (deg2rad * ll.lat ()), sinlat = sin (deg2rad * ll.lat ());
      x (jloc) = coslon * coslat; 
      y (jloc) = sinlon * coslat; 
      z (jloc) = sinlat;
    }

  return xyz;
}

struct shuffle_t
{
  size_t isize_recv;
  size_t isize_send;

  struct recv_t
  {
    // Here we describe how to use a message from a remote processor
    int iproc;                 // Rank of processor we receive data from
    size_t isize;              // Number of points we receive
    struct loccntoff_t
    {
      atlas::idx_t iloc;        // Local indice
      int          irem_cnt; // Remote count for current local indice
      int          irem_off; // Offset in remote buffer current local index
    };  
    std::vector<loccntoff_t> desc; // Allocated to number of points we have on this task for grid #2
  };

  struct send_t
  {
    int iproc;
    size_t isize;
    std::vector<atlas::idx_t> iloc; // Local points to send
  };

  std::vector<recv_t> yl_recv;
  std::vector<send_t> yl_send;

  std::vector<int> isort; // Sort received points (size = all points from grid #1 we receive on this task)

  // For each point on grid #2 : length/offset of points from grid #1
  
  struct loc2cntoff_t
  {
    int ioff; 
    int icnt;
  };

  std::vector<loc2cntoff_t> desc; // Size = number of points of grid #2 on this task
  

};


shuffle_t
create_shuffle (const grid::Distribution & dist1, const functionspace::StructuredColumns & fs1,
                const grid::Distribution & dist2, const functionspace::StructuredColumns & fs2)
{
  shuffle_t shuffle;

  size_t size1 = fs1.sizeOwned ();
  size_t size2 = fs2.sizeOwned ();
  const auto & grid1 = fs1.grid ();
  const auto & grid2 = fs2.grid ();

  auto & comm = atlas::mpi::comm ();

  int myproc = comm.rank (), nproc = comm.size ();

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
  };
  
  std::vector<prcglo_t> prcglo1 (size1);


// TODO : Use OpenMP on jloc1
  for (int jloc1 = 0; jloc1 < size1; jloc1++)
    {
      //  For all points of grid #1 (only in our region), we find the nearest point of grid #2. Here,
      //  "nearest" means nearest latitude and nearest longitude. We also find the task
      //  which holds this point of grid #2

      atlas::PointLonLat lonlat1 = grid1.StructuredGrid::lonlat (i1 (jloc1)-1, j1 (jloc1)-1);
      atlas::gidx_t iglo1 = grid1.ij2gidx (i1 (jloc1)-1, j1 (jloc1)-1);
      atlas::PointXY xy2 = proj2.xy (lonlat1);

      int iny2 = grid2.ny (), iy2;

      // Search along Y axis: find nearest latitude
// TODO : increasing Y coordinate
      if (xy2.y () > yspc2.front ())
        iy2 = 0;
      else if (xy2.y () < yspc2.back ())
        iy2 = iny2-1;
      else
        {
          int iy2a = 0, iy2b = iny2-1;
          while (1)
            {
              // Dichotomy
              int iy2m = (iy2a + iy2b) / 2;
              if ((yspc2[iy2a] >= xy2.y ()) && (xy2.y () >= yspc2[iy2m]))
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
       
      // Find nearest longitude

      int inx2 = grid2.nx (iy2), ix2;

// TODO : handle non global domains (grid1.domain.global () == false)
// TODO : handle shifted longitudes ??
      double dx = inx2 * xspc2.dx ()[iy2];
      ix2 = modulo (int (round (inx2 * xy2.x ()) / dx), inx2);

      gidx_t iglo2 = grid2.ij2gidx (ix2, iy2);
      int iprc2 = dist2.partition (iglo2);

      
      prcglo1[jloc1].iglo1 = iglo1;
      prcglo1[jloc1].iglo2 = iglo2;
      prcglo1[jloc1].iprc2 = iprc2;
    }


  // Compute & exchange send/recv counts
  std::vector<size_t> isendcnt (nproc, 0), irecvcnt (nproc), isendoff (nproc);
  
  for (int jloc1 = 0; jloc1 < size1; jloc1++)
    isendcnt[prcglo1[jloc1].iprc2]++;

  isendoff[0] = 0;
  for (int iproc = 1; iproc < nproc; iproc++)
    isendoff[iproc] = isendoff[iproc-1] + isendcnt[iproc-1];

  comm.allToAll (isendcnt, irecvcnt);

  // Sort points by task, global indice

  std::vector<int> iord_by_prc2glo2 (prcglo1.size ());

  std::iota (std::begin (iord_by_prc2glo2), std::end (iord_by_prc2glo2), 0);

  std::sort (std::begin (iord_by_prc2glo2), std::end (iord_by_prc2glo2), 
             [&prcglo1] (int a, int b) 
             { 
               if (prcglo1[a].iprc2 == prcglo1[b].iprc2)
                 return prcglo1[a].iglo2 < prcglo1[b].iglo2;
               return prcglo1[a].iprc2 < prcglo1[b].iprc2;
             });

  prcglo1 = reorder (prcglo1, iord_by_prc2glo2);

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

  shuffle.yl_recv.resize (inrecv);

  inrecv = 0;
  for (int iproc = 0; iproc < nproc; iproc++)
    if (irecvcnt[iproc])
      {
        shuffle.yl_recv[inrecv].iproc = iproc;
        yl_exch_recv[inrecv].size = irecvcnt[iproc];
        yl_exch_recv[inrecv].iglo1iglo2.resize (2 * irecvcnt[iproc]);
        reqrecv[inrecv] = comm.iReceive (&yl_exch_recv[inrecv].iglo1iglo2[0],
                                        yl_exch_recv[inrecv].iglo1iglo2.size (),
                                        iproc, 101);
        inrecv++;
      }

  comm.barrier ();

  // Send

// TODO : Use OpenMP

  shuffle.yl_send.resize (insend);

  insend = 0;
  for (int iproc = 0; iproc < nproc; iproc++)
    if (isendcnt[iproc])
      {
        yl_exch_send[insend].size = isendcnt[iproc];
        yl_exch_send[insend].iglo1iglo2.resize (2 * isendcnt[iproc]);

        shuffle.yl_send[insend].iproc = iproc;
        shuffle.yl_send[insend].isize = isendcnt[iproc];
        shuffle.yl_send[insend].iloc.resize (isendcnt[iproc]);
     
        // List of points of grid #1 to send to iproc
        for (int i = 0; i < isendcnt[iproc]; i++)
          {
            int jloc1 = isendoff[iproc] + i;
            shuffle.yl_send[insend].iloc[i] = prcglo1[jloc1].iglo1;
            yl_exch_send[insend].iglo1iglo2[2*i+0] = prcglo1[jloc1].iglo1;
            yl_exch_send[insend].iglo1iglo2[2*i+1] = prcglo1[jloc1].iglo2;
          }

        reqsend[insend] = comm.iSend (&yl_exch_send[insend].iglo1iglo2[0],
                                     yl_exch_send[insend].iglo1iglo2.size (),
                                     iproc, 101);
        insend++;
      }

  shuffle.isize_send = 
    std::accumulate (std::begin (shuffle.yl_send), std::end (shuffle.yl_send),
                     0, [] (int a, const shuffle_t::send_t & s)
                     {
                       return a + s.isize;
                     });
 

  for (int i = 0; i < inrecv; i++)
    comm.wait (reqrecv[i]);

  for (int i = 0; i < insend; i++)
    comm.wait (reqsend[i]);

  for (auto & yl_exch : yl_exch_send)
    yl_exch.iglo1iglo2.clear ();


  for (int ii = 0; ii < shuffle.yl_recv.size (); ii++) 
    {
      // Total number of points we get from this task
      shuffle.yl_recv[ii].isize = yl_exch_recv[ii].iglo1iglo2.size () / 2;

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
      
      shuffle.yl_recv[ii].desc.resize (inloc);
      
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

              shuffle.yl_recv[ii].desc[inloc].iloc = jloc2;
              shuffle.yl_recv[ii].desc[inloc].irem_off = jj;

              inloc++;
              jglo2 = iglo2;
            }
        }

      int isize = shuffle.yl_recv[ii].desc.size ();
      for (int jj = 0; jj < isize-1; jj++)
        shuffle.yl_recv[ii].desc[jj].irem_cnt = shuffle.yl_recv[ii].desc[jj+1].irem_off 
                                              - shuffle.yl_recv[ii].desc[jj+0].irem_off;

      shuffle.yl_recv[ii].desc[isize-1].irem_cnt 
                                              = shuffle.yl_recv[ii].isize 
                                              - shuffle.yl_recv[ii].desc[isize-1].irem_off;
      
    }

  // Total number of points received by this proc
  shuffle.isize_recv = std::accumulate (std::begin (shuffle.yl_recv), std::end (shuffle.yl_recv), 
                                        0, [] (int a, const shuffle_t::recv_t & r) { return a + r.isize; });

  shuffle.desc.resize (shuffle.isize_recv);

  for (auto & c : shuffle.desc)
    c.icnt = 0;

  for (int ii = 0; ii < shuffle.yl_recv.size (); ii++)
  for (int jj = 0; jj < shuffle.yl_recv[ii].desc.size (); jj++)
    {
      atlas::idx_t jloc2 = shuffle.yl_recv[ii].desc[jj].iloc;
      // Number of points for jloc2
      shuffle.desc[jloc2].icnt += shuffle.yl_recv[ii].desc[jj].irem_cnt;
    }

  // Offset of points for jloc2
  shuffle.desc[0].ioff = 0;
  for (int jloc2 = 1; jloc2 < fs2.sizeOwned (); jloc2++)
    shuffle.desc[jloc2].ioff = shuffle.desc[jloc2-1].ioff + shuffle.desc[jloc2-1].ioff;


  // Compute sort vector
  
  struct jloc2iglo1_t
  {
    atlas::idx_t  jloc2;
    atlas::gidx_t iglo1;
  };

  std::vector<jloc2iglo1_t> isort (shuffle.isize_recv); 
  std::vector<int> iord (shuffle.isize_recv);
  std:iota (std::begin (iord), std::end (iord), 0);



  for (int ii = 0, jl = 0; ii < shuffle.yl_recv.size (); ii++)
    {
// TODO: Use OpenMP
      for (int jj = 0; jj < shuffle.yl_recv[ii].desc.size (); jj++)
        {
          int jloc2 = shuffle.yl_recv[ii].desc[jj].iloc         ;
          int ioff  = shuffle.yl_recv[ii].desc[jj].irem_off + jl;
          int icnt  = shuffle.yl_recv[ii].desc[jj].irem_cnt     ;
          for (int k = 0; k < icnt; k++)
            isort[ioff+k].jloc2 = jloc2;
        }
      jl = jl + shuffle.yl_recv[ii].desc.size ();
    }

  for (int ii = 0, jl = 0; ii < yl_exch_recv.size (); ii++)
    {
      size_t isize = yl_exch_recv[ii].iglo1iglo2.size () / 2;
// TODO: Use OpenMP
      for (int jj = 0; jj < isize; jj++)
        isort[jl+jj].iglo1 = yl_exch_recv[ii].iglo1iglo2[2*jj+0];
      jl = jl + isize;
    }

  // Sort by local index, remote index

  std::sort (std::begin (iord), std::end (iord), [&isort] (int a, int b)
             {
               if (isort[a].jloc2 == isort[b].jloc2) 
                 return isort[a].iglo1 < isort[b].iglo1;
               return isort[a].jloc2 < isort[b].jloc2;
             });

  shuffle.isort = reverse (iord);

  int icount = 0;

  for (int jloc2 = 0; jloc2 < fs2.sizeOwned (); jloc2++)  
    {
      if (shuffle.desc[jloc2].icnt == 0)
        {
          icount++;
          if (icount < 20)
            std::cerr << " shuffle : no points were found for local " << jloc2 << std::endl;
        }
    }
  if (icount > 0)
    std::cerr << " shuffle : no points were for " << icount << " points " << std::endl;



  return shuffle;
}

FieldSet
do_shuffle (const shuffle_t & shuffle, const FieldSet & pgp1)
{
  FieldSet pgp2e;

  auto & comm = atlas::mpi::comm ();

  const int infld = pgp1.size ();
  const int insend = shuffle.yl_send.size ();
  const int inrecv = shuffle.yl_recv.size ();

  // Create fields in pgp2e

  for (int jfld = 0; jfld < infld; jfld++)
    pgp2e.add (atlas::Field (pgp1[jfld].name (),
                             atlas::array::DataType::kind<double> (), 
                             atlas::array::make_shape (shuffle.isize_recv)));

  std::vector<eckit::mpi::Request> 
                     reqsend (insend), 
                     reqrecv (inrecv);

  std::vector<double> zbufr (infld * shuffle.isize_recv);
  std::vector<double> zbufs (infld * shuffle.isize_send);

  for (int ii = 0, ioffr = 0; ii < inrecv; ii++)
    {
      int iproc = shuffle.yl_recv[ii].iproc;
      size_t isize = shuffle.yl_recv[ii].isize;
      reqrecv[ii] = comm.iReceive (&zbufr[ioffr], isize * infld, iproc, 101);

printf (" RECV FROM %8d : %8d ELEMENTS\n", iproc, isize * infld);

      ioffr = ioffr + infld * isize;
    }

  comm.barrier ();

  std::vector<size_t> ioffs_all (insend);

  for (int ii = 0, ioffs = 0; ii < insend; ii++)
    {
      size_t isize = shuffle.yl_send[ii].isize;
      ioffs_all[ii] = ioffs;
      ioffs += infld * isize;
    }

  for (int ii = 0; ii < insend; ii++)
    {
      size_t ioffs = ioffs_all[ii];
      int iproc = shuffle.yl_send[ii].iproc;
      size_t isize = shuffle.yl_send[ii].isize;

// TODO: Use OpenMP
      for (int jfld = 0; jfld < infld; jfld++)
        {
          auto v = atlas::array::make_view<double,1> (pgp1[jfld]);
          for (int jj = 0; jj < isize; jj++)
            zbufs[ioffs+jfld*isize+jj] = v (shuffle.yl_send[ii].iloc[jj]);
        }
      reqsend[ii] = comm.iSend (&zbufs[ioffs], infld * isize, iproc, 101);
printf (" SEND  TO  %8d : %8d ELEMENTS\n", iproc, isize * infld);
    }

// TODO: Use MPI_WAITANY
  for (int i = 0; i < inrecv; i++)
    comm.wait (reqrecv[i]);

  std::vector<size_t> ioffr_all (inrecv);

  for (int ii = 0, ioffr = 0; ii < inrecv; ii++)
    {
      ioffr_all[ii] = ioffr;
      size_t isize = shuffle.yl_recv[ii].isize;
      ioffr += isize;
    }

  for (int ii = 0; ii < inrecv; ii++)
    {
// TODO : Use OpenMP
      size_t ioffr = ioffr_all[ii];
      size_t isize = shuffle.yl_recv[ii].isize;

      for (int jfld = 0; jfld < infld; jfld++)
        {
          auto v = atlas::array::make_view<double,1> (pgp2e[jfld]);
          for (int jj = 0; jj < isize; jj++)
            v (shuffle.isort[ioffr+jj]) = zbufr[ioffr*infld+jfld*isize+jj];
        }
    }

  for (int i = 0; i < insend; i++)
    comm.wait (reqsend[i]);

  return pgp2e;
}

int main (int argc, char * argv[]) 
{
  atlas::Library::instance ().initialise (argc, argv);

  auto & comm = atlas::mpi::comm ();

  int myproc = comm.rank (), nproc = comm.size ();

  Config config;
  
  StructuredGrid grid1 ("L80x40");
  grid::Distribution dist1 = grid::Distribution (grid1, Config ("type", "checkerboard") | Config ("nbands", nproc));
  ReducedGaussianGrid grid2 ("N16");
  grid::Distribution dist2 (grid2, Config ("type", "equal_regions"));


if(0){
//const auto & part1 = dist1.partition ();
//
//printf ("-- part1 --\n");
//for (int i = 0; i < part1.size (); i++)
//  printf (" %8d > %8d\n", i, part1[i]);
//
//printf ("-- grid1 --\n");

  for (int iy = 0; iy < grid1.ny (); iy++)
  for (int ix = 0; ix < grid1.nx (iy); ix++)
    {
      double lonlat[2];
      grid1.lonlat (ix, iy, lonlat);
      printf (" %8d > %12.4f %12.4f\n", iy, lonlat[0], lonlat[1]);
    }
}

if(0)
{
  int iglo = 0;
  for (int iy = 0; iy < grid1.ny (); iy++)
  for (int ix = 0; ix < grid1.nx (iy); ix++)
    {
      idx_t ij[2];
      grid1.gidx2ij (iglo, ij);
      printf (" ix, iy, iglo = %8d, %8d, %8d, ij = %8d, %8d\n", ix, iy, iglo, ij[0], ij[1]);
      iglo++;
    }
}


  comm.barrier ();

  printf ("---- START ----\n");

  functionspace::StructuredColumns fs1 {grid1, dist1};
  functionspace::StructuredColumns fs2 {grid2, dist2};


  shuffle_t shuffle = create_shuffle (dist1, fs1, dist2, fs2);

  FieldSet lonlat1 = getLonLat (fs1);
  FieldSet lonlat2 = getLonLat (fs2);

  FieldSet lonlat2e = do_shuffle (shuffle, lonlat1);

  atlas::Library::instance ().finalise ();

  std::cout << "----- STOP -----" << std::endl;

  return 0;
}

