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

namespace
{

const double deg2rad = M_PI / 180.0;
const double rad2deg = 180.0 / M_PI;
const double ra = 6378137;

double distance(const double lat1, const double lon1,const double lat2, const double lon2){
  printf("Distance\n");
  double d=ra * acos ((sin(lat1*deg2rad)*sin(lat2*deg2rad))+(cos(lat1*deg2rad)*cos(lat2*deg2rad)*cos((lon2-lon1)*deg2rad)));
  return d;
}

double dLon(const double lon1,const double lon2,const double lat){
  //avoid import math
  double d;
  if (lon1>lon2){
    d=ra*cos(lat)*(lon1-lon2)*deg2rad;
  }else{
    d=ra*cos(lat)*(lon2-lon1)*deg2rad;
  }
  return d;
}
double dLat(const double lat1,const double lat2){
  //avoid import math
  double d;
  if (lat1>lat2){
    d=ra*(lat1-lat2)*deg2rad;
  }else{
    d=ra*(lat2-lat1)*deg2rad;
  }
  return d;
}

//double distanceBis(const double lat1, const double lon1,const double lat2, const double lon2){ 
// return (distance(lat1,lon1,lat1,lon2),distance(lat1,lon2,lat2,lon2));
//}

}


atlas::FieldSet do_halfdiff (const atlas::functionspace::StructuredColumns & fs, const atlas::Field & field)
{
  printf("TOTO\n");
  auto view = atlas::array::make_view<double,1> (field);
  auto & comm = atlas::mpi::comm ();
  int myproc = comm.rank ();
  int nproc = comm.size ();
  printf("A\n");
  atlas::array::ArrayView<double, 2> lonlat =atlas::array::make_view<double, 2>( fs.xy() );
  printf("B\n");
  fs.haloExchange(field);
  atlas::StructuredGrid grid=fs.grid();
  struct neighbours
  {
     int w,e,n,s,nw,ne,se,sw;
  };
  struct valneighbours
  {
     double w,e,n,s,nw,ne,se,sw;
  };
  struct grad_t
  {
     double w,e,n,s,dw,de,dn,ds;
  };
  std::vector<grad_t> pgrad (fs.sizeOwned());
  double zundef=0;

  printf("C\n");
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
      }else{
          pgrad[jl].w=0;
          pgrad[jl].e=0;
          pgrad[jl].n=0;
          pgrad[jl].s=0;
          pgrad[jl].dw=1;
          pgrad[jl].de=1;
          pgrad[jl].dn=1;
          pgrad[jl].ds=1;
      }
      neighbours voiz;//w,e,n,s,nw,ne,se,sw;
              //en haut C  gauche
      voiz.w =fs.index(i-1,j);
      voiz.e =fs.index(i+1,j);

      //pas le meme nombre de lon par bande de lat
      int icorN=0;
      int icorS=0;
      if (myproc==0 and j==0){
          icorN=(i+(grid.nx(j)/2))%(grid.nx(j));
      }else{
          icorN=std::round(grid.nx(j)*float(i)/grid.nx(j-1));
      }
      if(myproc==nproc-1 and j==(fs.j_end ()-1)){
          icorS=(i+(grid.nx(j)/2))%(grid.nx(j));
      }else{
          icorS=std::round(grid.nx(j)*float(i)/grid.nx(j+1));
      }
      printf("%5d,%5d,%5d,%5d\n",i,j,icorN,icorS);

      voiz.n =fs.index(icorN,j-1);
      voiz.s =fs.index(icorS,j+1);
      voiz.nw=fs.index(icorN-1,j-1);
      voiz.ne=fs.index(icorN+1,j-1);
      voiz.se=fs.index(icorS+1,j+1);
      voiz.sw=fs.index(icorS-1,j+1);
      printf("%5d,%5d,%5d,%5d,%5d,%5d,%5d,%5d\n",voiz.w,voiz.e,voiz.n,voiz.s,voiz.nw,voiz.ne,voiz.se,voiz.sw);
      valneighbours valvoiz;//w,e,n,s,nw,ne,se,sw;
      valvoiz.w=view(voiz.w);
      valvoiz.e=view(voiz.e);
      valvoiz.n=view(voiz.n);
      valvoiz.s=view(voiz.s);
      valvoiz.nw=view(voiz.nw);
      valvoiz.ne=view(voiz.ne);
      valvoiz.se=view(voiz.se);
      valvoiz.sw=view(voiz.sw);
      double val=view(jl);
      //pgrad[jl].e=distance(lonlat(voiz.e,atlas::LAT),lonlat(voiz.e,atlas::LON),lonlat(jl,atlas::LAT),lonlat(jl,atlas::LON));
      pgrad[jl].e=distance(lonlat(voiz.e,atlas::LAT),lonlat(voiz.e,atlas::LON),lonlat(jl,atlas::LAT),lonlat(jl,atlas::LON));
      if(valvoiz.e != zundef){
        pgrad[jl].de=valvoiz.e-val;
      }else{
        if((valvoiz.ne!=zundef) and (valvoiz.se!=zundef)){
            //moyenne des 2
            pgrad[jl].de=((valvoiz.ne+valvoiz.se))/2-val;

        }else if(valvoiz.ne!=zundef) {
           //on dit que la valeur du ne est la valeur a l'e
           pgrad[jl].de=valvoiz.ne-val;
        }else{
           //soit se soit zundeff
           pgrad[jl].de=valvoiz.se-val;
        }
      }
      pgrad[jl].w=distance(lonlat(voiz.w,atlas::LAT),lonlat(voiz.w,atlas::LON),lonlat(jl,atlas::LAT),lonlat(jl,atlas::LON));
      if(valvoiz.w != zundef){
        pgrad[jl].dw=val-valvoiz.w;
      }else{
        if((valvoiz.nw!=zundef) and (valvoiz.sw!=zundef)){

            pgrad[jl].dw=val-((valvoiz.nw+valvoiz.sw))/2;

        }else if(valvoiz.nw!=zundef) {
           pgrad[jl].dw=val-valvoiz.nw;
        }else{
           pgrad[jl].dw=val-valvoiz.sw;
        }
      }
      pgrad[jl].n=distance(lonlat(voiz.n,atlas::LAT),lonlat(voiz.n,atlas::LON),lonlat(jl,atlas::LAT),lonlat(jl,atlas::LON));
      if(valvoiz.n != zundef){
        pgrad[jl].dn=valvoiz.n-val;
      }else{
        if((valvoiz.ne !=zundef) and (valvoiz.nw !=zundef)){
           pgrad[jl].dn=((valvoiz.nw+valvoiz.ne))/2-val;
        }else if(valvoiz.nw!=zundef) {
           pgrad[jl].dn=valvoiz.nw-val;
        }else{
           pgrad[jl].dn=valvoiz.ne-val;
        }
      }
      pgrad[jl].s=distance(lonlat(voiz.s,atlas::LAT),lonlat(voiz.s,atlas::LON),lonlat(jl,atlas::LAT),lonlat(jl,atlas::LON));
      if(valvoiz.s != zundef){
        pgrad[jl].ds=val-valvoiz.s;
      }else{
        if((valvoiz.se !=zundef) and (valvoiz.sw !=zundef)){
            pgrad[jl].ds=val-((valvoiz.sw+valvoiz.se))/2;
        }else if(valvoiz.sw!=zundef) {
           pgrad[jl].ds=val-valvoiz.sw;
        }else{
           pgrad[jl].ds=val-valvoiz.se;
        }
      }
      printf("%5d\n",jl);

    }
  }

  printf("D\n");
  atlas::Field AOSIP = fs.createField (atlas::option::name ("AOSIP"),atlas::option::type("double"));
  atlas::Field AOSIM = fs.createField (atlas::option::name ("AOSIM"),atlas::option::type("double"));
  atlas::Field AOSJP = fs.createField (atlas::option::name ("AOSJP"),atlas::option::type("double"));
  atlas::Field AOSJM = fs.createField (atlas::option::name ("AOSJM"),atlas::option::type("double"));
  auto viewAOSIP = atlas::array::make_view<double,1> (AOSIP);
  auto viewAOSIM = atlas::array::make_view<double,1> (AOSIM);
  auto viewAOSJP = atlas::array::make_view<double,1> (AOSJP);
  auto viewAOSJM = atlas::array::make_view<double,1> (AOSJM);
  atlas::Field AODE = fs.createField (atlas::option::name ("AODE"),atlas::option::type("double"));
  atlas::Field AODW = fs.createField (atlas::option::name ("AODW"),atlas::option::type("double"));
  atlas::Field AODN = fs.createField (atlas::option::name ("AODN"),atlas::option::type("double"));
  atlas::Field AODS = fs.createField (atlas::option::name ("AODS"),atlas::option::type("double"));
  auto viewAODE = atlas::array::make_view<double,1> (AODE);
  auto viewAODW = atlas::array::make_view<double,1> (AODW);
  auto viewAODN = atlas::array::make_view<double,1> (AODN);
  auto viewAODS = atlas::array::make_view<double,1> (AODS);
  jl=0;
  for (int j = fs.j_begin (); j < fs.j_end (); j++){
    for (int i = fs.i_begin (j); i < fs.i_end (j); i++){
      jl= fs.index(i,j);

    if(view(jl)==zundef){
        viewAOSIP(jl)=zundef;
        viewAOSIM(jl)=zundef;
        viewAOSJP(jl)=zundef;
        viewAOSJM(jl)=zundef;
    }else{
      if(pgrad[jl].w==0){
        viewAOSIM(jl)=zundef;
      }else{
        viewAOSIM(jl)=-std::min(0.,pgrad[jl].dw/pgrad[jl].w);
      }
      if(pgrad[jl].e==0){
        viewAOSIP(jl)=zundef;
      }else{
        viewAOSIP(jl)=std::max(0.,pgrad[jl].de/pgrad[jl].e);
      }
      if(pgrad[jl].n==0){
        viewAOSJP(jl)=zundef;
      }else{
        viewAOSJP(jl)=std::max(0.,pgrad[jl].dn/pgrad[jl].n);
      }
      if(pgrad[jl].s==0){
        viewAOSJM(jl)=zundef;
      }else{
        viewAOSJM(jl)=-std::min(0.,pgrad[jl].ds/pgrad[jl].s);
      }


    }
    viewAODE(jl)=pgrad[jl].de;
    viewAODW(jl)=pgrad[jl].dw;
    viewAODN(jl)=pgrad[jl].dn;
    viewAODS(jl)=pgrad[jl].ds;
    }
  }
  printf("E\n");
  atlas::FieldSet AOS;
  AOS.add(AOSIP);
  AOS.add(AOSIM);
  AOS.add(AOSJP);
  AOS.add(AOSJM);
  AOS.add(AODE);
  AOS.add(AODW);
  AOS.add(AODN);
  AOS.add(AODS);
  return AOS;
 
}

atlas::field::FieldSetImpl * do_halfdiff__do_halfdiff (atlas::functionspace::detail::StructuredColumns * fs, atlas::field::FieldImpl * field)
{
  atlas::FieldSet grad = do_halfdiff (atlas::functionspace::StructuredColumns (fs), atlas::Field (field));
  atlas::field::FieldSetImpl * grad_ = grad.get (); 
  grad_->attach (); 
  return grad_;
}


