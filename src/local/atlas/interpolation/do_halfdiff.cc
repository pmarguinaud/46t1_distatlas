#include "do_halfdiff.h"


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



atlas::FieldSet do_halfdiff (const atlas::functionspace::StructuredColumns & fs, const atlas::Field & field)
{

  //check if functionspace has halo >1  todo
  //
  fs.haloExchange(field);

  //mpi stuff 
  auto & comm = atlas::mpi::comm ();
  int myproc = comm.rank ();
  int nproc = comm.size ();


  //view of the field
  auto view = atlas::array::make_view<double,1> (field);
  //lon lat of functionspace
  atlas::array::ArrayView<double, 2> lonlat =atlas::array::make_view<double, 2>( fs.xy() );
  
   
  atlas::StructuredGrid grid=fs.grid();
  struct neighbours
  {
     int w,e,nw,ne,se,sw;
  };
  struct valneighbours
  {
     double w,e,nw,ne,se,sw;
  };
  struct grad_t
  {
     double w,e,n,s,dw,de,dn,ds;
  };
  std::vector<grad_t> pgrad (fs.sizeOwned());
  double zundef=0; //read it on field metadata

  auto t= atlas::array::DataType::kind<double> ();
  auto s=  atlas::array::make_shape (fs.sizeOwned ());
  auto AOSIP = atlas::Field(std::string("AOSIP"),t,s);
  auto AOSIM = atlas::Field(std::string("AOSIM"),t,s);
  auto AOSJP = atlas::Field(std::string("AOSJP"),t,s);
  auto AOSJM = atlas::Field(std::string("AOSJM"),t,s);
  auto AODE = atlas::Field(std::string("AODE"),t,s);
  auto AODW = atlas::Field(std::string("AODW"),t,s);
  auto AODN = atlas::Field(std::string("AODN"),t,s);
  auto AODS = atlas::Field(std::string("AODS"),t,s);
  AOSIP.metadata () = field.metadata (); AOSIP.rename(std::string("AOSIP"));AOSIP.metadata().set("CLSUFF","AOSIP");
  AOSIM.metadata () = field.metadata ();AOSIM.rename(std::string("AOSIM"));AOSIM.metadata().set("CLSUFF","AOSIM");
  AOSJP.metadata () = field.metadata ();AOSJP.rename(std::string("AOSJP"));AOSJP.metadata().set("CLSUFF","AOSJP");
  AOSJM.metadata () = field.metadata ();AOSJM.rename(std::string("AOSJM"));AOSJM.metadata().set("CLSUFF","AOSJM");
  AODE.metadata () = field.metadata ();AODE.rename(std::string("AODE"));AODE.metadata().set("CLSUFF","AODE");
  AODW.metadata () = field.metadata ();AODW.rename(std::string("AODW"));AODW.metadata().set("CLSUFF","AODW");
  AODN.metadata () = field.metadata ();AODN.rename(std::string("AODN"));AODN.metadata().set("CLSUFF","AODN");
  AODS.metadata () = field.metadata ();AODS.rename(std::string("AODS"));AODS.metadata().set("CLSUFF","AODS");
  auto viewAOSIP = atlas::array::make_view<double,1> (AOSIP);
  auto viewAOSIM = atlas::array::make_view<double,1> (AOSIM);
  auto viewAOSJP = atlas::array::make_view<double,1> (AOSJP);
  auto viewAOSJM = atlas::array::make_view<double,1> (AOSJM);
  auto viewAODE = atlas::array::make_view<double,1> (AODE);
  auto viewAODW = atlas::array::make_view<double,1> (AODW);
  auto viewAODN = atlas::array::make_view<double,1> (AODN);
  auto viewAODS = atlas::array::make_view<double,1> (AODS);
  int jl=0;
  for (int j = fs.j_begin (); j < fs.j_end (); j++){
    for (int i = fs.i_begin (j); i < fs.i_end (j); i++){
      jl= fs.index(i,j);
      if(view(jl)==zundef){
          pgrad[jl].w=zundef;
          pgrad[jl].e=zundef;
          pgrad[jl].n=zundef;
          pgrad[jl].s=zundef;
          pgrad[jl].dw=zundef;
          pgrad[jl].de=zundef;
          pgrad[jl].dn=zundef;
          pgrad[jl].ds=zundef;
          viewAOSIP(jl)=zundef;
          viewAOSIM(jl)=zundef;
          viewAOSJP(jl)=zundef;
          viewAOSJM(jl)=zundef;
      }else{
          pgrad[jl].w=0;
          pgrad[jl].e=0;
          pgrad[jl].n=0;
          pgrad[jl].s=0;
          pgrad[jl].dw=1;
          pgrad[jl].de=1;
          pgrad[jl].dn=1;
          pgrad[jl].ds=1;
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
      
      pgrad[jl].w=lonlat(jl,atlas::LON)-lonlat(voiz.w,atlas::LON);
      if(view(voiz.w) != zundef){
         pgrad[jl].dw=val-view(voiz.w); 
      }
      pgrad[jl].e=lonlat(jl,atlas::LON)-lonlat(voiz.e,atlas::LON);
      if(view(voiz.e) != zundef){
         pgrad[jl].de=val-view(voiz.e); 
      }   

      if((view(voiz.nw) != zundef)and(view(voiz.ne) != zundef)){
          //on prend le barycentre des 2
          pgrad[jl].dn=val-(dcorN*view(voiz.nw) + (1-dcorN)*view(voiz.ne));
          pgrad[jl].n=lonlat(jl,atlas::LAT)-(dcorN*lonlat(voiz.nw,atlas::LAT)+(1-dcorN)*lonlat(voiz.ne,atlas::LAT)); //besoin de faire le barycentre ici aussi?
      } 
      
      if((view(voiz.sw) != zundef)and(view(voiz.se) != zundef)){
          //on prend le barycentre des 2
          pgrad[jl].ds=val-(dcorS*view(voiz.sw) + (1-dcorS)*view(voiz.se));
          pgrad[jl].s=lonlat(jl,atlas::LAT)-(dcorS*lonlat(voiz.sw,atlas::LAT)+(1-dcorS)*lonlat(voiz.se,atlas::LAT)); //besoin de faire le barycentre ici aussi?
      }

      ///
      if(pgrad[jl].w==0){
        viewAOSIM(jl)=zundef;
      }else{
        viewAOSIM(jl)=pgrad[jl].dw/pgrad[jl].w;
      }
      if(pgrad[jl].e==0){
        viewAOSIP(jl)=zundef;
      }else{
        viewAOSIP(jl)=pgrad[jl].de/pgrad[jl].e;
      }
      if(pgrad[jl].n==0){
        viewAOSJP(jl)=zundef;
      }else{
        viewAOSJP(jl)=pgrad[jl].dn/pgrad[jl].n;
      }
      if(pgrad[jl].s==0){
        viewAOSJM(jl)=zundef;
      }else{
        viewAOSJM(jl)=pgrad[jl].ds/pgrad[jl].s;
      }
 
    viewAODE(jl)=pgrad[jl].de;
    viewAODW(jl)=pgrad[jl].dw;
    viewAODN(jl)=pgrad[jl].dn;
    viewAODS(jl)=pgrad[jl].ds;
      }
    }
  }

  atlas::FieldSet AOS;
  AOS.add(AOSIP);
  AOS.add(AOSIM);
  AOS.add(AOSJP);
  AOS.add(AOSJM);
  AOS.add(AODE);
  AOS.add(AODW);
  AOS.add(AODN);
  AOS.add(AODS);
  AOS.add(field);
  return AOS;
 
}

atlas::field::FieldSetImpl * do_halfdiff__do_halfdiff (atlas::functionspace::detail::StructuredColumns * fs, atlas::field::FieldImpl * field)
{
  atlas::FieldSet grad = do_halfdiff (atlas::functionspace::StructuredColumns (fs), atlas::Field (field));
  atlas::field::FieldSetImpl * grad_ = grad.get (); 
  grad_->attach (); 
  return grad_;
}


