#include "halfdiff.h"

#include "atlas/runtime/Exception.h"
#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <string>  
#include <type_traits>

#include <sys/time.h>
#include <stdio.h>



atlas::FieldSet halfdiff (const atlas::functionspace::StructuredColumns & fs, const atlas::FieldSet & pgp)
{
  atlas::FieldSet pgpg;

  //check if functionspace has halo >1  todo
  //
  fs.haloExchange(pgp);

  //mpi stuff 
  auto & comm = atlas::mpi::comm ();
  int myproc = comm.rank ();
  int nproc = comm.size ();


  //lon lat of functionspace
  atlas::array::ArrayView<double, 2> xy =atlas::array::make_view<double, 2>( fs.xy() );
  
   
  const atlas::StructuredGrid & grid = fs.grid();
  struct neighbours
  {
     int w,e,nw,ne,se,sw;
  };
  struct grad_t
  {
     double w,e,n,s,dw,de,dn,ds;
  };

  const double zundef= std::numeric_limits<double>::max (); //read it on field metadata

  auto t = atlas::array::DataType::kind<double> ();
  auto s = atlas::array::make_shape (fs.sizeOwned ());
  
  auto addf = [&] (const atlas::Field & f, const std::string & name)
  {
    auto g = atlas::Field (f.name () + name, t, s);
    g.metadata () = f.metadata (); 
    pgpg.add (g);
  };

  for (int jfld = 0; jfld < pgp.size (); jfld++)
    {
      addf (pgp[jfld], "IP"); addf (pgp[jfld], "IM");
      addf (pgp[jfld], "JP"); addf (pgp[jfld], "JM");
    }

  for (int jfld = 0; jfld < pgp.size (); jfld++)
    {
      auto view   = atlas::array::make_view<double,1> (pgp[jfld]);
      auto viewIP = atlas::array::make_view<double,1> (pgpg[4*jfld+0]);
      auto viewIM = atlas::array::make_view<double,1> (pgpg[4*jfld+1]);
      auto viewJP = atlas::array::make_view<double,1> (pgpg[4*jfld+2]);
      auto viewJM = atlas::array::make_view<double,1> (pgpg[4*jfld+3]);

      for (int j = fs.j_begin (); j < fs.j_end (); j++)
      for (int i = fs.i_begin (j); i < fs.i_end (j); i++)
        {
          const auto jl = fs.index (i, j);
          grad_t pgrad;
          
          if (view(jl)==zundef)
            {
              viewIP(jl)=zundef; viewIM(jl)=zundef;
              viewJP(jl)=zundef; viewJM(jl)=zundef;
              continue;
            }

          pgrad.w=0;
          pgrad.e=0;
          pgrad.n=0;
          pgrad.s=0;
          pgrad.dw=1;
          pgrad.de=1;
          pgrad.dn=1;
          pgrad.ds=1;

          //value of the field
          double val=view(jl); 
          neighbours voiz;//w,e,n,s,nw,ne,se,sw;
          
          //correction because not the same lon by lat band number
          int icorN=0;
          double dcorN=0;
          int icorS=0;
          double dcorS=0;
          if (myproc==0 and j==0){
              icorN=(i+(grid.nx(j)/2))%(grid.nx(j));
              dcorN=0; //ca tombe juste
              voiz.nw=fs.index(icorN+1,j-1); //il faut inverser au pole
              voiz.ne=fs.index(icorN,j-1);
          }else{
              //int et pas round donc c'est le nw
              dcorN=grid.nx(j-1)*float(i)/grid.nx(j);
              icorN=int(dcorN);
              dcorN=dcorN-icorN;
              voiz.nw=fs.index(icorN,j-1);
              voiz.ne=fs.index(icorN+1,j-1);
          }
          if(myproc==nproc-1 and j==(fs.j_end ()-1)){
              icorS=(i+(grid.nx(j)/2))%(grid.nx(j));
              dcorS=0;
              voiz.se=fs.index(icorS,j+1); //inversion au pole
              voiz.sw=fs.index(icorS+1,j+1);
          }else{
              //sw
              dcorS=(grid.nx(j+1)*float(i)/grid.nx(j));
              icorS=int(dcorS);
              dcorS=dcorS-icorS;
              voiz.se=fs.index(icorS+1,j+1);
              voiz.sw=fs.index(icorS,j+1);
          }
          voiz.w =fs.index(i-1,j);
          voiz.e =fs.index(i+1,j);
          
          pgrad.w=xy(jl,atlas::LON)-xy(voiz.w,atlas::LON);
          if(view(voiz.w) != zundef){
             pgrad.dw=val-view(voiz.w); 
          }
          pgrad.e=xy(jl,atlas::LON)-xy(voiz.e,atlas::LON);
          if(view(voiz.e) != zundef){
             pgrad.de=val-view(voiz.e); 
          }   
          
          if((view(voiz.nw) != zundef)and(view(voiz.ne) != zundef)){
              //on prend le barycentre des 2
              pgrad.dn=val-(dcorN*view(voiz.nw) + (1-dcorN)*view(voiz.ne));
              pgrad.n=xy(jl,atlas::LAT)-(dcorN*xy(voiz.nw,atlas::LAT)+(1-dcorN)*xy(voiz.ne,atlas::LAT)); //besoin de faire le barycentre ici aussi?
          } 
          
          if((view(voiz.sw) != zundef)and(view(voiz.se) != zundef)){
              //on prend le barycentre des 2
              pgrad.ds=val-(dcorS*view(voiz.sw) + (1-dcorS)*view(voiz.se));
              pgrad.s=xy(jl,atlas::LAT)-(dcorS*xy(voiz.sw,atlas::LAT)+(1-dcorS)*xy(voiz.se,atlas::LAT)); //besoin de faire le barycentre ici aussi?
          }
          
          viewIM(jl)= pgrad.w ? pgrad.dw/pgrad.w : zundef;
          viewIP(jl)= pgrad.e ? pgrad.de/pgrad.e : zundef;
          viewJP(jl)= pgrad.n ? pgrad.dn/pgrad.n : zundef;
          viewJM(jl)= pgrad.s ? pgrad.ds/pgrad.s : zundef;

        }
    }

  return pgpg;
 
}

atlas::field::FieldSetImpl * halfdiff__ (atlas::functionspace::detail::StructuredColumns * fs, atlas::field::FieldSetImpl * pgp)
{
  atlas::FieldSet pgpg = halfdiff (atlas::functionspace::StructuredColumns (fs), atlas::FieldSet (pgp));
  atlas::field::FieldSetImpl * pgpg_ = pgpg.get (); 
  pgpg_->attach (); 
  return pgpg_;
}


