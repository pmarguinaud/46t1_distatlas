PROGRAM AAPGD

USE MPL_END_MOD,       ONLY : MPL_END
USE MPL_INIT_MOD,      ONLY : MPL_INIT
USE MPL_MODULE,        ONLY : MPL_MYRANK, MPL_GROUPS_CREATE, MPL_NUMPROC
USE MPL_BARRIER_MOD,   ONLY : MPL_BARRIER
USE MPL_ALLREDUCE_MOD, ONLY : MPL_ALLREDUCE
USE PARKIND1,          ONLY : JPIM, JPRB
USE YOMHOOK,           ONLY : DR_HOOK, LHOOK
USE MODD_IO_SURF_ARO,  ONLY : SURFEX_FIELD_BUF_CACHE,      &
                            & SURFEX_FIELD_BUF_ADD,        &
                            & SURFEX_FIELD_BUF_WRITE_MISC, &
                            & SURFEX_FIELD_BUF_DEALLOC
USE SFXFLDDESC_MOD,    ONLY : SFXFLDDESC, SFXFLDDESC_LOOKUP
USE HADES
USE ATLAS_ARPEGE_MODULE,  only :  wfaatlas => wfa




USE atlas_module  

!

IMPLICIT NONE

#include "setup_trans0.h"
#include "trans_inq.h"
#include "dist_grid.h"
#include "gath_grid.h"
#include "dir_trans.h"
#include "inv_trans.h"
#include "posnam.intfb.h"
#include "abor1.intfb.h"

INTEGER(KIND=JPIM) :: NPRGPNS, NPRGPEW, NPRTRW, NPRTRV, NPRINTLEV
INTEGER(KIND=JPIM) :: MYPROC, NPROC

INTEGER (KIND=JPIM) :: IA, IB

TYPE (GRID_OPTIONS_t), TARGET :: YLGROPT1, YLGROPT2
TYPE (GRID_OPTIONS_t), POINTER :: YLGROPT
TYPE (GRID_t) :: YLGRID1, YLGRID2
TYPE (DIST_t) :: YLDIST1, YLDIST2


type (atlas_structuredgrid)   :: YLGRID1A, YLGRID2A
type (atlas_griddistribution) :: YLDIST1A, YLDIST2A
type (SHUFFLE_t)  :: ylshflatlas

type(atlas_Config)                    :: YLGRID2AConf
type(atlas_Config)                    :: YLGRID1AConf

TYPE HDR_t
  REAL (KIND=JPRB) :: RUNDEF = 0
  REAL (KIND=JPRB) :: ZWEST = -999, ZEAST = -999, ZNORTH = -999, ZSOUTH = -999
  INTEGER (KIND=JPIM) :: ICOLS = -999, IROWS = -999, IWIDTH = -999
  CHARACTER (LEN=256) :: CLCOMMENT = ''
END TYPE HDR_t

TYPE FIELD_t
  CHARACTER (LEN=16) :: CPREF, CSUFF, CTYPE
  INTEGER (KIND=JPIM) :: IRANK, ILEVG = 0
END TYPE FIELD_t

REAL (KIND=JPRB), ALLOCATABLE, TARGET :: ZFLD1 (:,:), ZFLD2 (:,:), ZLATLON1 (:,:), ZCVR1 (:, :)
TYPE (SHUFFLE_t)  :: YLSHFL12
TYPE (HALO_t)     :: YLHALO1
TYPE (SHUFFLE4_t) :: YLSHFL21_4
TYPE (WEIGHTS4_t) :: YLWGHT21_4

INTEGER (KIND=JPIM) :: NGPTOT1, NGPTOT2
INTEGER (KIND=JPIM) :: JZS, JZS2, JZS_DX2, JZS_DY2, JZS_DXDY, JSAND, JCLAY
INTEGER (KIND=JPIM) :: JMIN_ZS, JMAX_ZS, JAVG_ZS, JSIL_ZS
INTEGER (KIND=JPIM) :: JSSO_DIR, JSSO_SLOPE, JSSO_ANIS, JSSO_STDEV
INTEGER (KIND=JPIM) :: JAOSIP, JAOSIM, JAOSJP, JAOSJM, JHO2IP, JHO2IM, JHO2JP, JHO2JM
INTEGER (KIND=JPIM) :: JLATGAUSS, JLONGAUSS, JLAT_G_XY, JLON_G_XY, JMESHGAUSS, JLONINF, JLONSUP, JLATINF, JLATSUP
INTEGER (KIND=JPIM) :: JFRAC_SEA, JFRAC_WATER, JFRAC_NATURE, JFRAC_TOWN
INTEGER (KIND=JPIM) :: JBATHY, JRUNOFFB, JWDRAIN
INTEGER (KIND=JPIM) :: IPASS
INTEGER (KIND=JPIM) :: IOFLD_INT, IOFLD_CAL, IOFLD_AUX, IOFLD_COV
INTEGER (KIND=JPIM) :: INFLD_INT, INFLD_CAL, INFLD_AUX, INFLD_COV
LOGICAL, ALLOCATABLE :: LCOVER (:)
LOGICAL :: LGARDEN, LWATER_TO_NATURE, LTOWN_TO_ROCK

NAMELIST / NAMGRID / YLGROPT
NAMELIST / NAMPGD / LGARDEN, LWATER_TO_NATURE, LTOWN_TO_ROCK


CHARACTER (LEN=32), PARAMETER :: CLDATA (4) = [ &
!                         01234567890123456789012345678901
&                        'relief                          ', &
&                        'CLAY                            ', &
&                        'SAND                            ', &
&                        'ecoclimap                       ' ]
INTEGER (KIND=JPIM), PARAMETER :: JOROG_DAT = 1, JCLAY_DAT = 2, JSAND_DAT = 3, JCOVR_DAT = 4

TYPE (HDR_t) :: YLHDRS (4), YLHDR
TYPE (FIELD_t), ALLOCATABLE :: YLFLDS (:)
INTEGER (KIND=JPIM) :: I, JFLD

REAL (KIND=JPRB) :: ZUNDEF

REAL (KIND=JPRB) :: ZHOOK_HANDLE



INTEGER(KIND=JPIM),DIMENSION(:), ALLOCATABLE :: dNglobalproc(:),eNglobalproc(:)

INTEGER (KIND=JPIM) :: ix,iy,ipt
CHARACTER (LEN=32) :: CLGRID


OPEN (4, FORM='FORMATTED')

NPRINTLEV = 0

CALL POSNAM (4, 'NAMPGD')
READ (4, NAMPGD)

CALL POSNAM (4, 'NAMGRID')
YLGROPT => YLGROPT2
READ (4, NAMGRID)

NULLIFY (YLGROPT)

CALL MPL_INIT (LDENV = .FALSE.)

MYPROC = MPL_MYRANK ()
NPROC  = MPL_NUMPROC

CALL SQUARE (NPROC, IA, IB)

NPRGPNS = IA
NPRGPEW = IB
NPRTRW  = IA
NPRTRV  = IB

CALL MPL_GROUPS_CREATE (NPRTRW, NPRTRV)



CALL SETUP_TRANS0 (KOUT=0, KERR=0, KPRINTLEV=NPRINTLEV, KMAX_RESOL=3,  &
                 & KPRGPNS=NPRGPNS, KPRGPEW=NPRGPEW, KPRTRW=NPRTRW,    &
                 & LDEQ_REGIONS=.TRUE.)

CALL MPL_BARRIER()

IF (LHOOK) CALL DR_HOOK ('AAINTERP',0,ZHOOK_HANDLE)


! Field setup

DO IPASS = 1, 2

  INFLD_INT = 0; INFLD_CAL = 0; INFLD_AUX = 0;

  CALL NEWFLD ('SFX.',  'ZS              ', 'INT', IPASS, JZS           )
  CALL NEWFLD ('SFX.',  'ZS2             ', 'AUX', IPASS, JZS2          )
  CALL NEWFLD ('SFX.',  'ZS_DX2          ', 'AUX', IPASS, JZS_DX2       )
  CALL NEWFLD ('SFX.',  'ZS_DY2          ', 'AUX', IPASS, JZS_DY2       )
  CALL NEWFLD ('SFX.',  'ZS_DXDY         ', 'AUX', IPASS, JZS_DXDY      )
  CALL NEWFLD ('SFX.',  'SAND            ', 'INT', IPASS, JSAND         )
  CALL NEWFLD ('SFX.',  'CLAY            ', 'INT', IPASS, JCLAY         )
  CALL NEWFLD ('SFX.',  'SSO_DIR         ', 'CAL', IPASS, JSSO_DIR      )
  CALL NEWFLD ('SFX.',  'SSO_SLOPE       ', 'CAL', IPASS, JSSO_SLOPE    )
  CALL NEWFLD ('SFX.',  'SSO_ANIS        ', 'CAL', IPASS, JSSO_ANIS     )
  CALL NEWFLD ('SFX.',  'SSO_STDEV       ', 'CAL', IPASS, JSSO_STDEV    )
  CALL NEWFLD ('SFX.',  'LATINF          ', 'CAL', IPASS, JLATINF       )
  CALL NEWFLD ('SFX.',  'LONINF          ', 'CAL', IPASS, JLONINF       )
  CALL NEWFLD ('SFX.',  'LATSUP          ', 'CAL', IPASS, JLATSUP       )
  CALL NEWFLD ('SFX.',  'LONSUP          ', 'CAL', IPASS, JLONSUP       )
  CALL NEWFLD ('SFX.',  'LATGAUSS        ', 'CAL', IPASS, JLATGAUSS     )
  CALL NEWFLD ('SFX.',  'LONGAUSS        ', 'CAL', IPASS, JLONGAUSS     )
  CALL NEWFLD ('SFX.',  'LAT_G_XY        ', 'CAL', IPASS, JLAT_G_XY     )
  CALL NEWFLD ('SFX.',  'LON_G_XY        ', 'CAL', IPASS, JLON_G_XY     )
  CALL NEWFLD ('SFX.',  'MESHGAUSS       ', 'CAL', IPASS, JMESHGAUSS    )
  
  CALL NEWFLD ('SFX.',  'AOSIP           ', 'CAL', IPASS, JAOSIP        )
  CALL NEWFLD ('SFX.',  'AOSIM           ', 'CAL', IPASS, JAOSIM        )
  CALL NEWFLD ('SFX.',  'AOSJP           ', 'CAL', IPASS, JAOSJP        )
  CALL NEWFLD ('SFX.',  'AOSJM           ', 'CAL', IPASS, JAOSJM        )
  CALL NEWFLD ('SFX.',  'HO2IP           ', 'CAL', IPASS, JHO2IP        )
  CALL NEWFLD ('SFX.',  'HO2IM           ', 'CAL', IPASS, JHO2IM        )
  CALL NEWFLD ('SFX.',  'HO2JP           ', 'CAL', IPASS, JHO2JP        )
  CALL NEWFLD ('SFX.',  'HO2JM           ', 'CAL', IPASS, JHO2JM        )

  CALL NEWFLD ('SFX.',  'FRAC_SEA        ', 'CAL', IPASS, JFRAC_SEA     )
  CALL NEWFLD ('SFX.',  'FRAC_TOWN       ', 'CAL', IPASS, JFRAC_TOWN    )
  CALL NEWFLD ('SFX.',  'FRAC_WATER      ', 'CAL', IPASS, JFRAC_WATER   )
  CALL NEWFLD ('SFX.',  'FRAC_NATURE     ', 'CAL', IPASS, JFRAC_NATURE  )

  CALL NEWFLD ('SFX.',  'BATHY           ', 'CAL', IPASS, JBATHY        )
  CALL NEWFLD ('SFX.',  'RUNOFFB         ', 'CAL', IPASS, JRUNOFFB      )
  CALL NEWFLD ('SFX.',  'WDRAIN          ', 'CAL', IPASS, JWDRAIN       )

  CALL NEWFLD ('SFX.',  'MIN_ZS          ', 'CAL', IPASS, JMIN_ZS       )
  CALL NEWFLD ('SFX.',  'MAX_ZS          ', 'CAL', IPASS, JMAX_ZS       )
  CALL NEWFLD ('SFX.',  'AVG_ZS          ', 'CAL', IPASS, JAVG_ZS       )
  CALL NEWFLD ('SFX.',  'SIL_ZS          ', 'CAL', IPASS, JSIL_ZS       )

  BLOCK
    INTEGER (KIND=JPIM) :: JCOVER, IRANK
    CHARACTER (LEN=16) :: CLSUFF
    IOFLD_COV = INFLD_INT + INFLD_AUX + INFLD_CAL
    INFLD_COV = 256
    DO JCOVER = 1, INFLD_COV
      WRITE (CLSUFF, '("COVER",I3.3)') JCOVER
      CALL NEWFLD ('SFX.', CLSUFF, 'CAL', IPASS, IRANK)
    ENDDO
  ENDBLOCK

  IF (IPASS == 1) THEN
    IOFLD_INT = 0
    IOFLD_AUX = IOFLD_INT + INFLD_INT
    IOFLD_CAL = IOFLD_AUX + INFLD_AUX
    ALLOCATE (YLFLDS (INFLD_INT + INFLD_AUX + INFLD_CAL))
  ENDIF

ENDDO

! Look at input files geometry

DO I = 1, 4
  CALL HDR (TRIM (CLDATA (I)), YLHDRS (I))
ENDDO

YLHDR = YLHDRS (1)

#ifdef UNDEF
DO I = 1, 4
  CALL COMPARE_GEOMHDR (YLHDR, YLHDRS (I))
ENDDO
#endif

! Create grids & distributions

!YLGROPT1%LSHIFTLON = .TRUE.
YLGROPT1%LZONAL    = .TRUE.   ! Force zonal distribution
YLGROPT1%LATLON    = .TRUE.
YLGROPT1%NDLON     = YLHDR%ICOLS
YLGROPT1%NDGLG     = YLHDR%IROWS
!YLGROPT1%RLONOFF   = YLHDR%ZWEST

CALL CREATE_GRID_DIST (YLGRID1, YLDIST1, YLGROPT1)
CALL CREATE_GRID_DIST (YLGRID2, YLDIST2, YLGROPT2)

allocate(dNglobalproc(YLGRID1%NGPTOTG))
allocate(eNglobalproc(YLGRID2%NGPTOTG))

CALL prepareDist (YLGRID1, YLDIST1, dNglobalproc)
CALL prepareDist (YLGRID2, YLDIST2, eNglobalproc)


IF (MYPROC == 1) THEN

BLOCK

INTEGER :: JLOC

DO JLOC = 1, SIZE (dNglobalproc)
WRITE (55, *) JLOC-1, dNglobalproc (JLOC)-1
ENDDO

DO JLOC = 1, SIZE (eNglobalproc)
WRITE (66, *) JLOC-1, eNglobalproc (JLOC)-1
ENDDO

ENDBLOCK

ENDIF



BLOCK

  TYPE (ATLAS_CONFIG) :: YLCFGR, YLCFDO
  REAL (KIND=JPRB) ZOFF
  CHARACTER (LEN=32) :: CLARG

  YLCFGR = ATLAS_CONFIG ()
  YLCFDO = ATLAS_CONFIG ()
  CALL YLCFGR%SET ("nx", YLGROPT1%NDLON)
  CALL YLCFGR%SET ("ny", YLGROPT1%NDGLG)
  CALL YLCFGR%SET ("type", "regular_lonlat")
  CALL YLCFDO%SET ("type", "rectangular")
  CALL YLCFDO%SET ("units", "degrees")

  CALL YLCFDO%SET ("xmin", +  0.0_JPRB)
  CALL YLCFDO%SET ("xmax", +360.0_JPRB)
  CALL YLCFDO%SET ("ymin", - 90.0_JPRB)
  CALL YLCFDO%SET ("ymax", + 90.0_JPRB)
  CALL YLCFGR%SET ("domain", YLCFDO)

  YLGRID1A = atlas_structuredgrid (YLCFGR)

  CALL YLCFGR%FINAL ()
  CALL YLCFDO%FINAL ()

! WRITE (CLGRID, '("S",I6.6,"x",I6.6)') YLGROPT1%NDLON, YLGROPT1%NDGLG
! YLGRID1A = atlas_structuredgrid (TRIM (CLGRID))

ENDBLOCK

YLDIST1A = atlas_griddistribution (dNglobalproc, part0=1)
YLGRID2AConf = atlas_Config ()
call YLGRID2AConf%set("type","rotated_schmidt")

YLGRID2A = atlas_reducedgaussiangrid (YLGROPT2%NLOEN (1:YLGROPT2%NDGLG),      &
                                 & [YLGROPT2%RLONCENT, YLGROPT2%RLATCENT], &
                                 & YLGROPT2%RSTRETCH)
YLDIST2A = atlas_griddistribution (eNglobalproc, part0=1)


IF (.TRUE.) THEN
BLOCK
  REAL (KIND=JPRB), ALLOCATABLE :: dRefLatLon (:,:), eRefLatLon (:,:)
  REAL (KIND=JPRB), ALLOCATABLE :: dTstLatLon (:,:), eTstLatLon (:,:)
  INTEGER (KIND=JPIM) :: IGPL
  
  ALLOCATE (dRefLatLon (YLGRID1%NGPTOTG, 2))
  ALLOCATE (eRefLatLon (YLGRID2%NGPTOTG, 2))
  ALLOCATE (dTstLatLon (YLGRID1%NGPTOTG, 2))
  ALLOCATE (eTstLatLon (YLGRID2%NGPTOTG, 2))

  PRINT *, " SIZE = ", YLGRID1A%SIZE ()
  PRINT *, " NGPTOTG = ", YLGRID1%NGPTOTG

  ipt=1
  DO iy = 1, INT (YLGRID1A%NY (), JPIM)
    DO ix = 1, INT (YLGRID1A%NX (iy), JPIM)
      dTstLatLon(ipt,2:1:-1) = YLGRID1A%LONLAT (ix, iy)
      ipt=ipt+1
    ENDDO
  ENDDO

  PRINT *, " SIZE = ", YLGRID2A%SIZE ()
  PRINT *, " NGPTOTG = ", YLGRID2%NGPTOTG
  ipt=1
  DO iy = 1, INT (YLGRID2A%NY (), JPIM)
    DO ix = 1, INT (YLGRID2A%NX (iy), JPIM)
      eTstLatLon(ipt,2:1:-1) = YLGRID2A%LONLAT (ix, iy)
      ipt=ipt+1
    ENDDO
  ENDDO
  
  CALL GET_MY_LATLON (YLGRID1, YLDIST1, dRefLatLon)
  CALL GET_MY_LATLON (YLGRID2, YLDIST2, eRefLatLon)
  
  dRefLatLon = dRefLatLon * RAD2DEG
  eRefLatLon = eRefLatLon * RAD2DEG

  WRITE (*, *) "----- d ------"
  DO IGPL = 1, SIZE (dRefLatLon, 1)
    WRITE (*, '(I8,2F12.4," | ",2F12.4)') IGPL, &
         & dRefLatLon (IGPL, 1), MODULO (dRefLatLon (IGPL, 2) + 180._JPRB, 360.0_JPRB) - 180.0_JPRB, &
         & dTstLatLon (IGPL, 1), dTstLatLon (IGPL, 2)
  ENDDO

  WRITE (*, *) "----- e ------"
  DO IGPL = 1, SIZE (eRefLatLon, 1)
    WRITE (*, '(I8,2F12.4," | ",2F12.4)') IGPL, &
         & eRefLatLon (IGPL, 1), eRefLatLon (IGPL, 2), &
         & eTstLatLon (IGPL, 1), eTstLatLon (IGPL, 2)
  ENDDO

ENDBLOCK
ENDIF


CALL DEMO_INTERPOLATIONA (YLDIST1,  YLGRID1,  YLDIST2,  YLGRID2,  &
                        & YLDIST1A, YLGRID1A, YLDIST2A, YLGRID2A, &
                        & .TRUE.)


#ifdef UNDEF


NGPTOT1 = YLDIST1%NGPTOTL (MYPROC)
NGPTOT2 = YLDIST2%NGPTOTL (MYPROC)

! Create a 1-point halo on grid #1 (for gradient computation)

CALL CREATE_HALO (YLGRID1, YLDIST1, 1, YLHALO1)

ALLOCATE (ZLATLON1 (YLHALO1%ISIZE+YLHALO1%IH_SIZE, 2))
ALLOCATE (ZFLD1 (YLHALO1%ISIZE+YLHALO1%IH_SIZE, INFLD_INT + INFLD_AUX))
ALLOCATE (ZFLD2 (NGPTOT2, INFLD_INT + INFLD_AUX + INFLD_CAL))
ALLOCATE (ZCVR1 (NGPTOT1, 1))

CALL GET_MY_LATLON (YLGRID1, YLDIST1, ZLATLON1)

! Read input data

CALL DIR (CLDATA (JOROG_DAT), YLHDRS (JOROG_DAT), ZFLD1 (:, JZS  ), YLGRID1, YLDIST1)
CALL DIR (CLDATA (JCLAY_DAT), YLHDRS (JCLAY_DAT), ZFLD1 (:, JCLAY), YLGRID1, YLDIST1)
CALL DIR (CLDATA (JSAND_DAT), YLHDRS (JSAND_DAT), ZFLD1 (:, JSAND), YLGRID1, YLDIST1)
CALL DIR (CLDATA (JCOVR_DAT), YLHDRS (JCOVR_DAT), ZCVR1 (:, 1),     YLGRID1, YLDIST1)

ZUNDEF = YLHDRS (JOROG_DAT)%RUNDEF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef block
type(atlas_functionspace_structuredcolumns) :: eFs,dFs
type(atlas_field) :: eField,dField
type(atlas_fieldset) :: eFieldSet,dFieldSet
real(KIND=JPRB), pointer   :: view(:,:)
type(atlas_Interpolation)             :: interpolation
type(atlas_Config)                    :: interpolation_config

real(KIND=JPRB), dimension(:),allocatable   :: eOut (:)
integer(KIND=JPIM):: jloop 


dFs=atlas_functionspace_structuredcolumns(YLGRID1A,YLDIST1A,halo=2)
eFs=atlas_functionspace_structuredcolumns(YLGRID2A,YLDIST2A,halo=2)


!write(*,*) dFs%size(),dFs%size_owned()
!write(*,*)  CALL dFs%sizeOwned()
!write(*,*)  CALL dFs%sizeHalo()
!write(*,*) 'j begin, j end', dFs%j_begin(),dFs%j_end()
!write(*,*) dFs%i_begin(dFs%j_begin()),dFs%i_end(dFs%j_end())

!write(*,*) 'j begin halo, j end halo', dFs%j_begin_halo(),dFs%j_end_halo()
!write(*,*) dFs%i_begin_halo(dFs%j_begin_halo()),dFs%i_end_halo(dFs%j_end_halo())

interpolation_config = atlas_Config()
!call interpolation_config%set("type","nearest-neighbour")
!call interpolation_config%set("type","structured-bilinear")
!call interpolation_config%set("halo","1")
!call interpolation_config%set("name","toto")
!call interpolation_config%set("type","finite-element")
call interpolation_config%set("type","structured-linear2D")
call interpolation_config%set("name","toto")
!

write(*,*) dFs%size()
write(*,*) 'j begin, j end', dFs%j_begin(),dFs%j_end()
write(*,*) dFs%i_begin(dFs%j_begin()),dFs%i_end(dFs%j_end())
interpolation = atlas_Interpolation(interpolation_config,dFs,eFs)
dFieldSet=atlas_fieldset("d")
eFieldSet=atlas_fieldset("e")

do jloop=1, 4
  dField = dFs%create_field(name=CLDATA(jloop),kind=ATLAS_REAL(JPRB),LEVELS=1)
  eField = eFs%create_field(name=CLDATA(jloop),kind=ATLAS_REAL(JPRB),LEVELS=1)

  CALL dField%data(view)
  view(1,:)=-1
  view(1,1:dFs%size_owned())=ZFLD1 (1:YLHALO1%ISIZE, jloop  )
  CALL dFs%halo_exchange(dField)
  CALL eField%data(view)
  view(:,:)=-1
  call dFieldSet%add(dField)
  call eFieldSet%add(eField)
  call interpolation%execute(dFieldSet,eFieldSet)
  CALL eFs%halo_exchange(eField)
end do 


CALL WFAATLAS ("ATLAS.T.fa",YLGRID2A,eFs,  &
        & ["SURF","SFX.","SFX.","SFX."], [0,0,0,0], ["GEOPOTENTIEL","COVER002","COVER003","COVER001"],eFieldSet)
!write(*,*) "shapes: ",shape(view),shape(ZFLD1 (:, JZS  ))
!write(*,*) "test: ",YLHALO1%ISIZE,YLHALO1%IH_SIZE



!write(*,*) view(1,dFs%size_owned()+1:dFs%size_owned()+10)


!write(*,*) view(1,dFs%size_owned()+1:dFs%size_owned()+10)




!call interpolation%execute(dField,eField)




!NOT NEEDED eGlobalFs=atlas_functionspace_structuredcolumns(YLGRID2A,eDist)
!eGlobalField =eFs%create_field(name="g",kind=ATLAS_REAL(JPRB),LEVELS=1,GLOBAL=.TRUE.)
!CALL eGlobalField%data(view)
!view(:,:)=-4
!CALL eFs%gather(eField,eGlobalField)
!write(*,*) 'field size', eField%size(),eGlobalField%size()
!write(*,*) 'oromean', sum(view(:,:)),eGlobalField%size()

!--------------------------------


!IF (MYPROC == 1) THEN

!  allocate(eOut(YLGRID2%NGPTOTG))
!  eOut(:)=view(1,:)
!
!  CALL WFA (eOut,  &
!        & YLGRID2%NGPTOTG, 1_JPIM, "ZFLD2.T.fa", &
!        & YLFLDS (IOFLD_CAL+1:IOFLD_CAL+INFLD_CAL)%CPREF, &
!        & YLFLDS (IOFLD_CAL+1:IOFLD_CAL+INFLD_CAL)%ILEVG, &
!        & YLFLDS (IOFLD_CAL+1:IOFLD_CAL+INFLD_CAL)%CSUFF, &
!        & YLGRID2, YLDIST2, PUNDEF=ZUNDEF)
!
!ENDIF


endblock
#endif
!CALL EXIT(0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! Create shuffle & weight functions

CALL CREATE_SHUFFLE (YLDIST1,YLDIST1A,YLGRID1,YLGRID1A,YLDIST2,YLDIST2A,YLGRID2,YLGRID2A, YLSHFL12,ylshflatlas)
CALL CREATE_SHUFFLE4 (YLDIST2, YLGRID2, YLDIST1, YLGRID1, YLSHFL21_4)
CALL CREATE_WEIGHTS4 (YLDIST2, YLGRID2, YLDIST1, YLGRID1, YLSHFL21_4, YLWGHT21_4)


! Create covers

BLOCK
  REAL (KIND=JPRB), ALLOCATABLE :: ZCVR2E (:,:)
  INTEGER (KIND=JPIM) :: ISUMRECV, II, ICOV, JLOC2, IOFF, ICNT
  ISUMRECV = SUM (YLSHFL12%YL_RECV%ISIZE)
  ALLOCATE (ZCVR2E (ISUMRECV, 1))

  CALL DO_SHUFFLE (YLDIST1, YLGRID1, YLDIST2, YLGRID2, YLSHFL12, ZCVR1, ZCVR2E)

!$OMP PARALLEL DO PRIVATE (JLOC2, IOFF, ICNT, II, ICOV)
  DO JLOC2 = 1, NGPTOT2

    IOFF = YLSHFL12%ILOCAL_OFF (JLOC2)
    ICNT = YLSHFL12%ILOCAL_CNT (JLOC2)

    ZFLD2 (JLOC2, IOFLD_COV+1:IOFLD_COV+INFLD_COV) = 0._JPRB

    DO II = IOFF+1, IOFF+ICNT
      ICOV = INT (ZCVR2E (II, 1))
      IF (ICOV < 1 .OR. ICOV > INFLD_COV) THEN
        CALL ABOR1 ('UNEXPECTED COVER VALUE')
      ENDIF
      ZFLD2 (JLOC2, IOFLD_COV+ICOV) = ZFLD2 (JLOC2, IOFLD_COV+ICOV) + 1
    ENDDO

    ZFLD2 (JLOC2, :) = ZFLD2 (JLOC2, :) / REAL (ICNT, JPRB)

  ENDDO
!$OMP END PARALLEL DO

ENDBLOCK

! Set missing values of orography

BLOCK
  WHERE (ZCVR1 (1:NGPTOT1, 1) == 1) ! Sea
    ZFLD1 (1:NGPTOT1, JZS) = ZUNDEF
  ENDWHERE
ENDBLOCK

! Get latlon and orography halos

CALL DO_HALO (YLHALO1, YLDIST1, ZFLD1 (:, JZS:JZS))
CALL DO_HALO (YLHALO1, YLDIST1, ZLATLON1)

! Compute orography gradient

BLOCK
  REAL (KIND=JPRB) :: ZSGRAD (NGPTOT1, 2)
  REAL (KIND=JPRB) :: ZFACTOR (NGPTOT1)
  REAL (KIND=JPRB) :: ZXYZ2 (NGPTOT2, 3)


  CALL GET_MY_XYZ (YLGRID2, YLDIST2, ZXYZ2)
  CALL GET_MY_FACTOR (YLGRID1, YLDIST1, ZFACTOR)
  CALL DO_GRADIENT (YLHALO1, YLDIST1, YLGRID1, ZLATLON1, ZFLD1 (:, JZS:JZS), &
                  & ZSGRAD (:, 1:1), ZSGRAD (:, 2:2), PUNDEF=[ZUNDEF], PFACTOR=ZFACTOR)


! Rotate orography gradient
  IF (YLGRID2%LROTATED .OR. YLGRID2%LSTRETCH) THEN
  BLOCK
    REAL (KIND=JPRB) :: ZGX1, ZGY1, U0U1, V0U1, U0V1, V0V1
    REAL (KIND=JPRB) :: ZUV1 (NGPTOT1, 4)
    REAL (KIND=JPRB) :: ZUV2 (NGPTOT2, 4)
    INTEGER (KIND=JPIM) :: ILOC1

    CALL GET_MY_UV (YLGRID2, YLDIST2, ZXYZ2, ZUV2)
    
    CALL DO_INTERPOLATION4 (YLDIST2, YLGRID2, YLDIST1, YLGRID1,      &
                          & ZUV2, ZUV1, YLWGHT21_4, YDSHFL12=YLSHFL21_4)

    DO ILOC1 = 1, NGPTOT1
      ZGX1 = ZSGRAD (ILOC1, 1); ZGY1 = ZSGRAD (ILOC1, 2)
      U0U1 = ZUV1 (ILOC1, 1); V0U1 = ZUV1 (ILOC1, 2)
      U0V1 = ZUV1 (ILOC1, 3); V0V1 = ZUV1 (ILOC1, 4)
      ZSGRAD (ILOC1, 1) = ZGX1 * U0U1 + ZGY1 * V0U1
      ZSGRAD (ILOC1, 2) = ZGX1 * U0V1 + ZGY1 * V0V1
    ENDDO
  ENDBLOCK
  ENDIF

  WHERE (ZFLD1 (1:NGPTOT1, JZS) /= ZUNDEF)
    ZFLD1 (:, JZS_DX2)  = ZSGRAD (:, 1) ** 2
    ZFLD1 (:, JZS_DY2)  = ZSGRAD (:, 2) ** 2
    ZFLD1 (:, JZS_DXDY) = ZSGRAD (:, 1) * ZSGRAD (:, 2)
    ZFLD1 (:, JZS2)     = ZFLD1 (:, JZS) ** 2
  ELSEWHERE
    ZFLD1 (:, JZS_DX2)  = ZUNDEF
    ZFLD1 (:, JZS_DY2)  = ZUNDEF
    ZFLD1 (:, JZS_DXDY) = ZUNDEF
    ZFLD1 (:, JZS2)     = ZUNDEF
  ENDWHERE

ENDBLOCK

! Interpolate fields

CALL DO_INTERPOLATION (YLDIST1, YLGRID1, YLDIST2, YLGRID2, ZFLD1 (:,1:INFLD_INT+INFLD_AUX), &
                     & ZFLD2 (:,1:INFLD_INT+INFLD_AUX), YDSHFL12=YLSHFL12, &
                     & PUNDEF=[(ZUNDEF, I = 1, INFLD_INT+INFLD_AUX)])

! Min, max, avg orography

BLOCK
  REAL (KIND=JPRB) :: ZZS1 (NGPTOT1, 2)

  ZZS1 (1:NGPTOT1, 1) = ZFLD1 (1:NGPTOT1, JZS)
  ZZS1 (1:NGPTOT1, 2) = ZFLD1 (1:NGPTOT1, JZS)

  IF (JMAX_ZS-JMIN_ZS /= 1) CALL ABOR1 ('MIN_ZS AND MAX_ZS NOT STORED TOGETHER')

  CALL DO_INTERPOLATION (YLDIST1, YLGRID1, YLDIST2, YLGRID2, ZZS1, &
                       & ZFLD2 (:,JMIN_ZS:JMAX_ZS), YDSHFL12=YLSHFL12, &
                       & PUNDEF=[ZUNDEF, ZUNDEF], CDTYPE=['MIN', 'MAX'])


  ZFLD2 (:, JAVG_ZS) = ZFLD2 (:, JZS)
  ZFLD2 (:, JSIL_ZS) = ZFLD2 (:, JZS)

ENDBLOCK


! Compute SSO parameters

BLOCK
  REAL (KIND=JPIM), PARAMETER :: XPI = RPI
  LOGICAL :: OSSO (NGPTOT2), OSSO_ANIS (NGPTOT2)
  REAL (KIND=JPRB), POINTER :: ZHXX (:), ZHYY (:), ZHXY (:)
  REAL (KIND=JPRB) :: ZK (NGPTOT2), ZL (NGPTOT2), ZM (NGPTOT2)
  REAL (KIND=JPRB), POINTER :: XSSO_DIR (:), XSSO_SLOPE (:), XSSO_ANIS (:), XSSO_STDEV (:)
  INTEGER (KIND=JPIM) :: II


  ZHXX       => ZFLD2 (:, JZS_DX2)
  ZHYY       => ZFLD2 (:, JZS_DY2)
  ZHXY       => ZFLD2 (:, JZS_DXDY)
  XSSO_DIR   => ZFLD2 (:, JSSO_DIR)
  XSSO_SLOPE => ZFLD2 (:, JSSO_SLOPE)
  XSSO_ANIS  => ZFLD2 (:, JSSO_ANIS)
  XSSO_STDEV => ZFLD2 (:, JSSO_STDEV)

  OSSO = ZFLD2 (:, JZS) /= ZUNDEF
  OSSO_ANIS = OSSO

  WHERE (OSSO (:))
    ZK(:)=0.5*(ZHXX(:)+ZHYY(:))
    ZL(:)=0.5*(ZHXX(:)-ZHYY(:))
    ZM(:)=     ZHXY(:)
  ELSE WHERE
    ZK(:) = ZUNDEF
    ZL(:) = ZUNDEF
    ZM(:) = ZUNDEF
    ZHXX = ZUNDEF
    ZHYY = ZUNDEF
    ZHXY = ZUNDEF
    XSSO_DIR   = ZUNDEF
    XSSO_SLOPE = ZUNDEF
    XSSO_ANIS  = ZUNDEF
  END WHERE
  !
  !*    8.     S.S.O. characteristics
  !            ----------------------
  !
  !*    8.1    S.S.O. direction of main axis
  !            -----------------------------
  !
  WHERE (OSSO(:))
    XSSO_DIR(:) = 0.5_JPRB * ATAN2 (ZM, ZL) * (180._JPRB / XPI)
  END WHERE
  !
  !*    8.2    S.S.O. slope
  !            ------------
  !
  WHERE (OSSO(:))
    XSSO_SLOPE(:) = SQRT( ZK+SQRT(ZL*ZL+ZM*ZM) )
  END WHERE
  !
  !*    8.3    S.S.O. anisotropy
  !            -----------------
  !
  WHERE (OSSO_ANIS(:) .AND. (ZK+SQRT(ZL*ZL+ZM*ZM)) >0. )
    XSSO_ANIS(:)=SQRT( MAX(ZK-SQRT(ZL*ZL+ZM*ZM),0.) / (ZK+SQRT(ZL*ZL+ZM*ZM)))
  END WHERE
  !
  WHERE (OSSO_ANIS(:) .AND. (ZK+SQRT(ZL*ZL+ZM*ZM))==0. )
    XSSO_ANIS(:)=1.
  END WHERE

! Orography standard deviation

  WHERE (OSSO (:))
    XSSO_STDEV = SQRT (ZFLD2 (:, JZS2) - ZFLD2 (:, JZS)**2)
  ELSE WHERE
    XSSO_STDEV = ZUNDEF
  END WHERE


ENDBLOCK

! Add enveloppe

WHERE (ZFLD2 (:, JZS) /= ZUNDEF)
  ZFLD2 (:, JZS) = ZFLD2 (:, JZS) + ZFLD2 (:, JSSO_STDEV) 
ENDWHERE

! Compute AOS parameters

BLOCK
! Orography of grid #2 interpolated on grid #1
  REAL (KIND=JPRB) :: ZS1INT4 (NGPTOT1, 1)     
! (Orography of grid #1) - (Orography of grid #2 interpolated on grid #1)
! we need a halo for gradient computation
  REAL (KIND=JPRB) :: ZS1MZS1INT4 (YLHALO1%ISIZE+YLHALO1%IH_SIZE, 1) 
  REAL (KIND=JPRB) :: ZDZXP (NGPTOT1, 1), ZDXP (NGPTOT1, 1)
  REAL (KIND=JPRB) :: ZDZXM (NGPTOT1, 1), ZDXM (NGPTOT1, 1)
  REAL (KIND=JPRB) :: ZDZYP (NGPTOT1, 1), ZDYP (NGPTOT1, 1)
  REAL (KIND=JPRB) :: ZDZYM (NGPTOT1, 1), ZDYM (NGPTOT1, 1)
  REAL (KIND=JPRB) :: ZAOS (NGPTOT1, 8)
  REAL (KIND=JPRB) :: ZTMP2 (NGPTOT2, 8)
  INTEGER (KIND=JPIM) :: JFLD
  
! Do 4-point interpolation from grid #2 to grid #1
  CALL DO_INTERPOLATION4 (YLDIST2, YLGRID2, YLDIST1, YLGRID1,      &
                        & ZFLD2 (:, JZS:JZS), ZS1INT4, YLWGHT21_4, &
                        & YDSHFL12=YLSHFL21_4, PUNDEF=[ZUNDEF])

! Difference of original orography with interpolated orography

  WHERE (ZFLD1 (1:NGPTOT1, JZS) /= ZUNDEF)
    ZS1MZS1INT4 (:, 1) = ZFLD1 (:, JZS) - ZS1INT4 (:, 1)
  ELSEWHERE
    ZS1MZS1INT4 (:, 1) = ZUNDEF
  ENDWHERE

! Get halo for orography difference

  CALL DO_HALO (YLHALO1, YLDIST1, ZS1MZS1INT4)

! Compute half-differences on grid #1

  CALL DO_HALFDIFF (YLHALO1, YLDIST1, YLGRID1, ZLATLON1, ZS1MZS1INT4,           &
                  & PDIFFXP=ZDZXP, PDIFFXM=ZDZXM, PDIFFYP=ZDZYP, PDIFFYM=ZDZYM, &
                  & PDXP=ZDXP,     PDXM=ZDXM,     PDYP=ZDYP,     PDYM=ZDYM,     &
                  & PUNDEF=[ZUNDEF])

  WHERE (ZFLD1 (1:NGPTOT1, JZS) /= ZUNDEF)
    ZAOS (:, 1) = + MAX (0._JPRB, ZDZXP (:, 1) / ZDXP (:, 1))    ! AOSIP
    ZAOS (:, 2) = - MIN (0._JPRB, ZDZXM (:, 1) / ZDXM (:, 1))    ! AOSIM
    ZAOS (:, 3) = + MAX (0._JPRB, ZDZYP (:, 1) / ZDYP (:, 1))    ! AOSJP
    ZAOS (:, 4) = - MIN (0._JPRB, ZDZYM (:, 1) / ZDYM (:, 1))    ! AOSJM
  ELSEWHERE
    ZAOS (:, 1) = ZUNDEF
    ZAOS (:, 2) = ZUNDEF
    ZAOS (:, 3) = ZUNDEF
    ZAOS (:, 4) = ZUNDEF
  ENDWHERE

  WHERE ((ZFLD1 (1:NGPTOT1, JZS) /= ZUNDEF) .AND. (ZDZXP (1:NGPTOT1, 1) > 0._JPRB)) 
    ZAOS (1:NGPTOT1, 5) = + ZDZXP (1:NGPTOT1, 1) * 0.5_JPRB       ! HO2IP
  ELSEWHERE
    ZAOS (1:NGPTOT1, 5) = ZUNDEF
  ENDWHERE

  WHERE ((ZFLD1 (1:NGPTOT1, JZS) /= ZUNDEF) .AND. (ZDZXM (1:NGPTOT1, 1) < 0._JPRB)) 
    ZAOS (1:NGPTOT1, 6) = - ZDZXM (1:NGPTOT1, 1) * 0.5_JPRB       ! HO2IM
  ELSEWHERE
    ZAOS (1:NGPTOT1, 6) = ZUNDEF
  ENDWHERE

  WHERE ((ZFLD1 (1:NGPTOT1, JZS) /= ZUNDEF) .AND. (ZDZYP (1:NGPTOT1, 1) > 0._JPRB)) 
    ZAOS (1:NGPTOT1, 7) = + ZDZYP (1:NGPTOT1, 1) * 0.5_JPRB       ! HO2JP
  ELSEWHERE
    ZAOS (1:NGPTOT1, 7) = ZUNDEF
  ENDWHERE

  WHERE ((ZFLD1 (1:NGPTOT1, JZS) /= ZUNDEF) .AND. (ZDZYM (1:NGPTOT1, 1) < 0._JPRB)) 
    ZAOS (1:NGPTOT1, 8) = - ZDZYM (1:NGPTOT1, 1) * 0.5_JPRB       ! HO2JM
  ELSEWHERE
    ZAOS (1:NGPTOT1, 8) = ZUNDEF
  ENDWHERE

  CALL DO_INTERPOLATION (YLDIST1, YLGRID1, YLDIST2, YLGRID2,              &
                       & ZAOS (1:NGPTOT1, 1:8),  ZFLD2 (:,JAOSIP:JHO2JM), &
                       & YDSHFL12=YLSHFL12, PUNDEF=[(ZUNDEF, I = 1, 8)])

  DO JFLD = JAOSIP, JHO2JM
    WHERE ((ZFLD2 (1:NGPTOT2, JZS) /= ZUNDEF) .AND. (ZFLD2 (1:NGPTOT2, JFLD) == ZUNDEF))
      ZFLD2 (1:NGPTOT2, JFLD) = 0._JPRB
    ENDWHERE

  ENDDO


ENDBLOCK

! Compute Gaussian grid parameters

BLOCK
  REAL (KIND=JPRB) :: ZLAT (0:YLGRID2%NDGLG+1)
  REAL (KIND=JPRB) :: ZCOORDYX2 (NGPTOT2, 2), ZFACTOR2 (NGPTOT2)
  REAL (KIND=JPRB) :: ZTEMP2 (NGPTOT2)
  REAL (KIND=JPRB) :: ZLATN, ZLATS, ZLOND, ZLAT2, ZLATD, ZLATLON2 (NGPTOT2, 2)
  INTEGER (KIND=JPIM) :: ILATLON2 (2, NGPTOT2)
  INTEGER (KIND=JPIM) :: JLAT2, JLOC2

  CALL GET_MY_FACTOR (YLGRID2, YLDIST2, ZFACTOR2)
  CALL GET_MY_LATLON (YLGRID2, YLDIST2, ZLATLON2)
  CALL GET_MY_JLATJLON (YLGRID2, YLDIST2, ILATLON2)
  CALL GET_MY_COORDYX (YLGRID2, YLDIST2, ZCOORDYX2)
  
  ZFLD2 (:, JLATGAUSS) = RAD2DEG * ZLATLON2 (:, 1)
  ZFLD2 (:, JLONGAUSS) = RAD2DEG * ZLATLON2 (:, 2)

  ZFLD2 (:, JLATINF) = RAD2DEG * ZLATLON2 (:, 1)
  ZFLD2 (:, JLONINF) = RAD2DEG * ZLATLON2 (:, 2)

  ZFLD2 (:, JLATSUP) = RAD2DEG * ZLATLON2 (:, 1)
  ZFLD2 (:, JLONSUP) = RAD2DEG * ZLATLON2 (:, 2)

  ZFLD2 (:, JLAT_G_XY) = RAD2DEG * ZCOORDYX2 (:, 1)
  ZFLD2 (:, JLON_G_XY) = RAD2DEG * ZCOORDYX2 (:, 2)

  ZLAT = [+RPI/2._JPRB, (ASIN (YLGRID2%RMU (JLAT2)), JLAT2 = 1, YLGRID2%NDGLG), -RPI/2._JPRB]

  DO JLOC2 = 1, NGPTOT2
    JLAT2 = ILATLON2 (1, JLOC2)

    ZLOND = R2PI / YLGRID2%NLOEN (JLAT2)

    ZLATN = ZLAT (JLAT2-1)
    ZLAT2 = ZLAT (JLAT2+0)
    ZLATS = ZLAT (JLAT2+1)

    IF (JLAT2 == 1) THEN
      ZLATD = (ZLATN - (ZLAT2 + ZLATS) / 2._JPRB)
    ELSEIF (JLAT2 == YLGRID2%NDGLG) THEN
      ZLATD = ((ZLATN + ZLAT2) / 2._JPRB - ZLATS)
    ELSE
      ZLATD = (ZLATN - ZLATS) / 2._JPRB
    ENDIF

    ZFLD2 (JLOC2, JMESHGAUSS) = RA * COS (ZLAT2) * ZLOND * RA * ZLATD / 2._JPRB 

    ZFLD2 (JLOC2, JMESHGAUSS) = ZFLD2 (JLOC2, JMESHGAUSS) / (ZFACTOR2 (JLOC2) ** 2)

  ENDDO

ENDBLOCK

! Compute fractions

BLOCK
  USE READCOVERS_MOD
  TYPE (COVERS_t) :: YLCOVERS
  REAL (KIND=JPRB) :: ZAVGCOVER (INFLD_COV)
  INTEGER (KIND=JPIM) :: JCOV

  CALL FIELDSTAT (YLDIST2, YLGRID2, ZFLD2 (:,IOFLD_COV+1:IOFLD_COV+INFLD_COV), PAVG=ZAVGCOVER)

  CALL READCOVERS (YLCOVERS, 1, LGARDEN, LWATER_TO_NATURE, LTOWN_TO_ROCK)

  ALLOCATE (LCOVER (YLCOVERS%JPCOVER))
  LCOVER = .FALSE.
  
  LCOVER (1:INFLD_COV) = ZAVGCOVER > 0._JPRB

  CALL MAKEFRAC (ZFLD2 (:, JFRAC_SEA),    ZFLD2 (:,IOFLD_COV+1:IOFLD_COV+INFLD_COV), YLCOVERS%XDATA_SEA,    LCOVER)
  CALL MAKEFRAC (ZFLD2 (:, JFRAC_NATURE), ZFLD2 (:,IOFLD_COV+1:IOFLD_COV+INFLD_COV), YLCOVERS%XDATA_NATURE, LCOVER)
  CALL MAKEFRAC (ZFLD2 (:, JFRAC_TOWN),   ZFLD2 (:,IOFLD_COV+1:IOFLD_COV+INFLD_COV), YLCOVERS%XDATA_TOWN,   LCOVER)
  CALL MAKEFRAC (ZFLD2 (:, JFRAC_WATER),  ZFLD2 (:,IOFLD_COV+1:IOFLD_COV+INFLD_COV), YLCOVERS%XDATA_WATER,  LCOVER)


ENDBLOCK

! Other

BLOCK
  WHERE (ZFLD2 (1:NGPTOT2, JFRAC_SEA) > 0._JPRB)
    ZFLD2 (:, JBATHY  ) = -300._JPRB
  ELSEWHERE
    ZFLD2 (:, JBATHY  ) = ZUNDEF
  ENDWHERE

  WHERE (ZFLD2 (1:NGPTOT2, JFRAC_NATURE) > 0._JPRB)
    ZFLD2 (:, JRUNOFFB) = 0.5_JPRB
    ZFLD2 (:, JWDRAIN ) = 0.0_JPRB
  ELSEWHERE
    ZFLD2 (:, JRUNOFFB) = ZUNDEF
    ZFLD2 (:, JWDRAIN ) = ZUNDEF
  ENDWHERE
ENDBLOCK

! Fix missing values: orography stuff, set to zero

BLOCK
  INTEGER (KIND=JPIM) :: JFLD, JIDX
  INTEGER (KIND=JPIM) :: IOROGIDX (17) 

  IOROGIDX = &
& [JZS, JSSO_DIR, JSSO_SLOPE, JSSO_ANIS, JSSO_STDEV, JAOSIP, JAOSIM, &
&  JAOSJP, JAOSJM, JHO2IP, JHO2IM, JHO2JP, JHO2JM, JMIN_ZS, JMAX_ZS, &
&  JAVG_ZS, JSIL_ZS]

  DO JFLD = 1, SIZE (IOROGIDX)
    JIDX = IOROGIDX (JFLD)
    WHERE (ZFLD2 (:, JIDX) == ZUNDEF)
      ZFLD2 (:, JIDX) = 0._JPRB
    ENDWHERE
  ENDDO



ENDBLOCK

! Fix missing values: sand & clay, use average values

BLOCK
  REAL (KIND=JPRB) :: ZAVGSANDCLAY (2)
  INTEGER (KIND=JPIM) :: ISACLIDX (2), JFLD, JIDX

  ISACLIDX = [JSAND, JCLAY]

  IF (JCLAY-JSAND /= 1) CALL ABOR1 ('SAND/CLAY NOT STORED TOGETHER')

  CALL FIELDSTAT (YLDIST2, YLGRID2, ZFLD2 (:,JSAND:JCLAY), &
                & PAVG=ZAVGSANDCLAY, PUNDEF=[ZUNDEF, ZUNDEF])

  DO JFLD = 1, 2
    JIDX = ISACLIDX (JFLD)
    WHERE ((ZFLD2 (:,JIDX) == ZUNDEF) .AND. (ZFLD2 (:, JFRAC_NATURE) > 0._JPRB))
      ZFLD2 (:,JIDX) = ZAVGSANDCLAY (JFLD)
    ENDWHERE 
    WHERE (ZFLD2 (:,JIDX) /= ZUNDEF) 
      ZFLD2 (:,JIDX) = ZFLD2 (:,JIDX) / 100._JPRB
    END WHERE
  ENDDO

ENDBLOCK

! Rank #1 writes metadata

BLOCK
  LOGICAL, PARAMETER :: LLCHECK = .TRUE.
  LOGICAL :: LLPRHDR
  INTEGER (KIND=JPIM) :: ILUN, JFLD, IREP, JRANK, ICNTUNDEF
  INTEGER (KIND=JPIM) :: IDIM_NATURE, IDIM_TOWN, IDIM_WATER, IDIM_SEA
  INTEGER (KIND=JPIM) :: NGPTOTG2
  INTEGER (KIND=JPIM) :: IERR
  TYPE (SFXFLDDESC) :: YLFLDDSC

  TYPE (SURFEX_FIELD_BUF_CACHE) :: YLSFXC

  YLSFXC%LRECORD = .FALSE.

  IDIM_NATURE = COUNT (ZFLD2 (1:NGPTOT2, JFRAC_NATURE) > 0._JPRB)
  IDIM_TOWN   = COUNT (ZFLD2 (1:NGPTOT2, JFRAC_TOWN  ) > 0._JPRB)
  IDIM_WATER  = COUNT (ZFLD2 (1:NGPTOT2, JFRAC_WATER ) > 0._JPRB)
  IDIM_SEA    = COUNT (ZFLD2 (1:NGPTOT2, JFRAC_SEA   ) > 0._JPRB)

  CALL MPL_ALLREDUCE (IDIM_NATURE, 'SUM', CDSTRING='AAINTERP:')
  CALL MPL_ALLREDUCE (IDIM_TOWN  , 'SUM', CDSTRING='AAINTERP:')
  CALL MPL_ALLREDUCE (IDIM_WATER , 'SUM', CDSTRING='AAINTERP:')
  CALL MPL_ALLREDUCE (IDIM_SEA   , 'SUM', CDSTRING='AAINTERP:')

  NGPTOTG2 = YLGRID2%NGPTOTG


  IF (MYPROC == 1) THEN

#define SET(a,b) CALL SURFEX_FIELD_BUF_ADD (YLSFXC, a, b)
SET (8_JPIM                    , 'VERSION'     ); SET (0_JPIM                    , 'BUG'         );
SET ('WATFLX'                  , 'WATER'       ); SET (.FALSE.                   , 'CH_EMIS'     );
SET (YLGRID2%RLATCENT * RAD2DEG, 'LAPO'        ); SET (YLGRID2%RLONCENT * RAD2DEG, 'LOPO'        );
SET (YLGRID2%RSTRETCH          , 'CODIL'       ); SET (YLGRID2%NLOEN             , 'NLOPA'       );
SET (YLGRID2%NDGLG             , 'NLATI'       ); SET (LCOVER                    , 'COVER_LIST'  );
SET (0_JPIM                    , 'COVER_PACKED'); SET (.FALSE.                   , 'CTI'         );
SET (.FALSE.                   , 'DATA_IRRIG'  ); SET (NGPTOTG2                  , 'DIM_FULL'    );
SET (IDIM_NATURE               , 'DIM_NATURE'  ); SET (IDIM_SEA                  , 'DIM_SEA'     );
SET (IDIM_TOWN                 , 'DIM_TOWN'    ); SET (IDIM_WATER                , 'DIM_WATER'   );
SET (0_JPIM                    , 'DUMMY_GR_NBR'); SET (.TRUE.                    , 'ECOCLIMAP'   );
SET (.FALSE.                   , 'GARDEN'      ); SET ('GAUSS'                   , 'GRID_TYPE'   );
SET (3_JPIM                    , 'GROUND_LAYER'); SET ('3-L'                     , 'ISBA'        );
SET (.FALSE.                   , 'L_GAMMAGV'   ); SET (.TRUE.                    , 'LCLIM_LAI'   );
                                                  SET (.FALSE.                   , 'L_ALBNIR_SOI');
SET (.FALSE.                   , 'L_ALBNIR_VEG'); SET (.FALSE.                   , 'L_ALBUV_SOI' );
SET (.FALSE.                   , 'L_ALBUV_VEG' ); SET (.FALSE.                   , 'L_ALBVIS_SOI');
SET (.FALSE.                   , 'L_ALBVIS_VEG'); SET (.FALSE.                   , 'L_BSLAI'     );
SET (.FALSE.                   , 'L_CE_NITRO'  ); SET (.FALSE.                   , 'L_CF_NITRO'  );
SET (.FALSE.                   , 'L_CNA_NITRO' ); SET (.FALSE.                   , 'L_CV'        );
SET (.FALSE.                   , 'L_DG'        ); SET (.FALSE.                   , 'L_DICE'      );
SET (.FALSE.                   , 'L_DMAX'      ); SET (.FALSE.                   , 'L_EMIS'      );
SET (.FALSE.                   , 'L_F2I'       ); SET (.FALSE.                   , 'L_GAMMA'     );
SET (.FALSE.                   , 'L_GC'        ); SET (.FALSE.                   , 'L_GMES'      );
SET (.FALSE.                   , 'L_GROUND_DPT'); SET (.FALSE.                   , 'L_H_TREE'    );
SET (.FALSE.                   , 'L_IRRIG'     ); SET (.FALSE.                   , 'L_LAI'       );
SET (.FALSE.                   , 'L_LAIMIN'    ); SET (.FALSE.                   , 'L_RE25'      );
SET (.FALSE.                   , 'L_RGL'       ); SET (.FALSE.                   , 'L_ROOTFRAC'  );
SET (.FALSE.                   , 'L_ROOT_DEPTH'); SET (.FALSE.                   , 'L_ROOT_EXT'  );
SET (.FALSE.                   , 'L_ROOT_LIN'  ); SET (.FALSE.                   , 'L_RSMIN'     );
SET (.FALSE.                   , 'L_SEFOLD'    ); SET (.FALSE.                   , 'L_STRESS'    );
SET (.FALSE.                   , 'L_VEG'       ); SET (.FALSE.                   , 'L_VEGTYPE'   );
SET (.FALSE.                   , 'L_WATSUP'    ); SET (.FALSE.                   , 'L_WRMAX_CF'  );
SET (.FALSE.                   , 'L_Z0'        ); SET (.FALSE.                   , 'L_Z0_O_Z0H'  );
SET ('ISBA'                    , 'NATURE'      ); SET (0_JPIM                    , 'NBIOMASS'    );
SET (.FALSE.                   , 'NO'          ); SET (1_JPIM                    , 'PATCH_NUMBER');
SET ('CH78'                    , 'PEDOTF'      ); SET (.FALSE.                   , 'PERMAFROST'  );
SET ('NON'                     , 'PHOTO'       ); SET (0.0_JPRB                  , 'RM_PATCH'    );
SET ('SEAFLX'                  , 'SEA'         ); SET (.FALSE.                   , 'SOCP'        );
SET (.FALSE.                   , 'SST_DATA'    ); SET ('PGD'                     , 'STORAGETYPE' );
SET ('NONE'                    , 'TOWN'        ); SET (.TRUE.                    , 'TOWN_TO_ROCK');
SET (.FALSE.                   , 'TR_ML'       ); SET (.FALSE.                   , 'WATER_TO_NAT');
SET (.FALSE.                   , 'L_GNDLITTER' ); SET (.FALSE.                   , 'L_H_VEG'     );
SET (.FALSE.                   , 'L_LAIGV'     ); SET (.FALSE.                   , 'L_RGLGV'     );
SET (.FALSE.                   , 'L_RSMINGV'   ); SET (.FALSE.                   , 'L_RTFRACGV'  );
SET (.FALSE.                   , 'L_RT_DEPTHGV'); SET (.FALSE.                   , 'L_RT_EXTGV'  );
SET (.FALSE.                   , 'L_WRMAX_CFGV'); SET (.FALSE.                   , 'L_Z0LITTER'  );
SET ([.FALSE.]                 , 'MEB_PATCH'   ); SET (.FALSE.                   , 'GWKEY'       );
#undef SET
  ENDIF

  LLPRHDR = .TRUE.
  DO JFLD = 1, SIZE (YLFLDS)
    JRANK = YLFLDS (JFLD)%IRANK

    IF (TRIM (YLFLDS (JFLD)%CTYPE) == 'AUX') CYCLE

    CALL SFXFLDDESC_LOOKUP (YLFLDS (JFLD)%CSUFF, YLFLDDSC, KERR=IERR)

    IF (IERR /= 0) THEN
      WRITE (0, *) " WARNING : UNKNOWN FIELD "//TRIM (YLFLDS (JFLD)%CSUFF)
    ENDIF

    IF (YLFLDDSC%CMASK == '') THEN
      WRITE (0, *) " WARNING : MASK IS NOT DEFINED FOR "//TRIM (YLFLDS (JFLD)%CSUFF)
    ENDIF

    IF (MYPROC == 1) THEN
      CALL SURFEX_FIELD_BUF_ADD (YLSFXC, ZFLD2 (:, JRANK),    &
                               & TRIM (YLFLDS (JFLD)%CSUFF),  &
                               & CDMASK=YLFLDDSC%CMASK)
    ENDIF


    IF (LLCHECK) THEN
! Check masks

      IF (YLFLDDSC%CMASK == 'NATURE') THEN
        ICNTUNDEF = COUNT ((ZFLD2 (:, JRANK) == ZUNDEF) .AND. (ZFLD2 (:, JFRAC_NATURE) > 0._JPRB))
      ELSEIF (YLFLDDSC%CMASK == 'SEA') THEN
        ICNTUNDEF = COUNT ((ZFLD2 (:, JRANK) == ZUNDEF) .AND. (ZFLD2 (:, JFRAC_SEA   ) > 0._JPRB))
      ELSEIF (YLFLDDSC%CMASK == 'WATER') THEN
        ICNTUNDEF = COUNT ((ZFLD2 (:, JRANK) == ZUNDEF) .AND. (ZFLD2 (:, JFRAC_WATER ) > 0._JPRB))
      ELSEIF (YLFLDDSC%CMASK == 'TOWN') THEN
        ICNTUNDEF = COUNT ((ZFLD2 (:, JRANK) == ZUNDEF) .AND. (ZFLD2 (:, JFRAC_TOWN  ) > 0._JPRB))
      ELSEIF (YLFLDDSC%CMASK == 'FULL') THEN
        ICNTUNDEF = COUNT (ZFLD2 (:, JRANK) == ZUNDEF)
      ELSE
        CALL ABOR1 ('UNKNOWN MASK '//TRIM (YLFLDDSC%CMASK))
      ENDIF

      CALL MPL_ALLREDUCE (ICNTUNDEF, 'SUM', CDSTRING='AAINTERP:')
  
      IF (ICNTUNDEF > 0 .AND. MYPROC == 1) THEN
        IF (LLPRHDR) THEN
          WRITE (0, *) "UNEXPECTED UNDEFINED VALUES..."
          WRITE (0, '(A16," ",A8," ",A10)') "FIELD NAME", "MASK", "UNDEFS"          
          LLPRHDR = .FALSE.
        ENDIF
        WRITE (0, '(A16," ",A8," ",I10)') YLFLDS (JFLD)%CSUFF, YLFLDDSC%CMASK, ICNTUNDEF
      ENDIF
    ENDIF

  ENDDO
    

  IF (MYPROC == 1) THEN
    ILUN = 77
    CALL OFA (ILUN, "ZFLD2.META.fa", YLGRID2)
    CALL SURFEX_FIELD_BUF_WRITE_MISC (YLSFXC, ILUN)
    CALL FAIRME (IREP, ILUN, 'KEEP')
    CALL SURFEX_FIELD_BUF_DEALLOC (YLSFXC)
  ENDIF

ENDBLOCK

! Write fields

CALL WFA (ZFLD2 (1:NGPTOT2,IOFLD_INT+1:IOFLD_INT+INFLD_INT),  &
        & NGPTOT2, INFLD_INT, "ZFLD2.INT.fa", &
        & YLFLDS (IOFLD_INT+1:IOFLD_INT+INFLD_INT)%CPREF, &
        & YLFLDS (IOFLD_INT+1:IOFLD_INT+INFLD_INT)%ILEVG, &
        & YLFLDS (IOFLD_INT+1:IOFLD_INT+INFLD_INT)%CSUFF, &
        & YLGRID2, YLDIST2, PUNDEF=ZUNDEF)

CALL WFA (ZFLD2 (1:NGPTOT2,IOFLD_CAL+1:IOFLD_CAL+INFLD_CAL),  &
        & NGPTOT2, INFLD_CAL, "ZFLD2.CAL.fa", &
        & YLFLDS (IOFLD_CAL+1:IOFLD_CAL+INFLD_CAL)%CPREF, &
        & YLFLDS (IOFLD_CAL+1:IOFLD_CAL+INFLD_CAL)%ILEVG, &
        & YLFLDS (IOFLD_CAL+1:IOFLD_CAL+INFLD_CAL)%CSUFF, &
        & YLGRID2, YLDIST2, PUNDEF=ZUNDEF)



BLOCK
  REAL (KIND=JPRB) :: ZXYZ2 (NGPTOT2, 3)
  REAL (KIND=JPRB) :: ZANGLE2 (NGPTOT2)
  REAL (KIND=JPRB) :: ZUV2 (NGPTOT2, 4)
  
  CALL GET_MY_XYZ (YLGRID2, YLDIST2, ZXYZ2)
  CALL GET_MY_UV (YLGRID2, YLDIST2, ZXYZ2, ZUV2)

  CALL GET_MY_ANGLE (YLGRID2, YLDIST2, ZXYZ2, ZANGLE2)

  CALL WFA (ZANGLE2, NGPTOT2, 1_JPIM, "ZFLD2.ANG.fa", &
          & ["SURF"], [0], [".ANGLE"], YLGRID2, YLDIST2)

  CALL WFA (ZUV2, NGPTOT2, 4_JPIM, "ZFLD2.UV.fa", &
          & ["SURF", "SURF", "SURF", "SURF"], [0, 0, 0, 0], &
          & [".UU", ".UV", ".VU", ".VV"], YLGRID2, YLDIST2)

ENDBLOCK

#endif

IF (LHOOK) CALL DR_HOOK ('AAINTERP',1,ZHOOK_HANDLE)


CALL MPL_END

CONTAINS

SUBROUTINE MAKEFRAC (PFRAC, PCOVER, PDATA, LDCOVER)

REAL (KIND=JPRB), INTENT (OUT) :: PFRAC (:)
REAL (KIND=JPRB), INTENT (IN)  :: PCOVER (:,:)
REAL (KIND=JPRB), INTENT (IN)  :: PDATA (:)
LOGICAL,          INTENT (IN)  :: LDCOVER (:)

INTEGER (KIND=JPIM) :: JLOC

!$OMP PARALLEL DO PRIVATE (JLOC)
DO JLOC = 1, SIZE (PFRAC)
  PFRAC (JLOC) = SUM (PCOVER (JLOC, :) * PDATA (:), MASK=LDCOVER) / SUM (PCOVER (JLOC, :), MASK=LDCOVER)
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE DIR (CDFILE, YDHDR, PDATA, YDGRID, YDDIST)

CHARACTER (LEN=*),   INTENT (IN)     :: CDFILE
TYPE (HDR_t),        INTENT (INOUT)  :: YDHDR
REAL (KIND=JPRB),    INTENT (OUT)    :: PDATA (:)
TYPE (GRID_t),       INTENT (IN)     :: YDGRID
TYPE (DIST_t),       INTENT (IN)     :: YDDIST

INTEGER*1, ALLOCATABLE :: IL1 (:)
INTEGER*2, ALLOCATABLE :: IL2 (:)
INTEGER (KIND=JPIM) :: IROW, ICOL, IOFF, JLAT, IROW1, IROW2
REAL (KIND=JPRB), ALLOCATABLE :: ZLAT (:)

IROW1 = YDDIST%NFRSTLAT (MYPROC)
IROW2 = YDDIST%NLSTLAT  (MYPROC)

OPEN (77, FILE=TRIM (CDFILE)//'.dir', FORM='UNFORMATTED', ACTION='READ', &
    & ACCESS='DIRECT', RECL=((YDHDR%ICOLS * YDHDR%IWIDTH) / 8))

SELECT CASE (YDHDR%IWIDTH)

  CASE ( 8)
  
    ALLOCATE (IL1 (YDHDR%ICOLS))

    IOFF = 0
    DO IROW = IROW1, IROW2
      READ (77, REC=IROW) IL1
      WHERE (IL1 < 0)
        PDATA (IOFF+1:IOFF+YDHDR%ICOLS) = IL1 + 256
      ELSEWHERE
        PDATA (IOFF+1:IOFF+YDHDR%ICOLS) = IL1 
      ENDWHERE
      IOFF = IOFF + YDHDR%ICOLS
    ENDDO

  CASE (16)

    ALLOCATE (IL2 (YDHDR%ICOLS))

    IOFF = 0
    DO IROW = IROW1, IROW2
      READ (77, REC=IROW) IL2
      PDATA (IOFF+1:IOFF+YDHDR%ICOLS) = IL2 
      IOFF = IOFF + YDHDR%ICOLS
    ENDDO

  CASE DEFAULT

    STOP

END SELECT

CLOSE (77)


WHERE (PDATA == YDHDR%RUNDEF)
  PDATA = HUGE (YDHDR%RUNDEF)
END WHERE
YDHDR%RUNDEF = HUGE (YDHDR%RUNDEF)


END SUBROUTINE 

SUBROUTINE HDR (CDFILE, YDHDR)

CHARACTER (LEN=*), INTENT (IN)  :: CDFILE
TYPE (HDR_t),      INTENT (OUT) :: YDHDR

CHARACTER (LEN=256) :: CLLINE

OPEN (77, FILE=TRIM (CDFILE)//'.hdr', FORM='FORMATTED')

READ (77, '(A)') CLLINE

YDHDR%CLCOMMENT = CLLINE

DO 
  READ (77, '(A256)', END=999) CLLINE
  IF (CLLINE (1:7) == 'nodata:'             ) READ (CLLINE ( 8:), *) YDHDR%RUNDEF
  IF (CLLINE (1:6) == 'north:'              ) READ (CLLINE ( 7:), *) YDHDR%ZNORTH
  IF (CLLINE (1:6) == 'south:'              ) READ (CLLINE ( 7:), *) YDHDR%ZSOUTH
  IF (CLLINE (1:5) == 'west:'               ) READ (CLLINE ( 6:), *) YDHDR%ZWEST
  IF (CLLINE (1:5) == 'east:'               ) READ (CLLINE ( 6:), *) YDHDR%ZEAST
  IF (CLLINE (1:5) == 'rows:'               ) READ (CLLINE ( 6:), *) YDHDR%IROWS
  IF (CLLINE (1:5) == 'cols:'               ) READ (CLLINE ( 6:), *) YDHDR%ICOLS
  IF (CLLINE (1:19) == 'recordtype: integer') READ (CLLINE (20:), *) YDHDR%IWIDTH
ENDDO

999 CONTINUE

CLOSE (77)

SELECT CASE (YDHDR%IWIDTH)

  CASE ( 8)
  
    IF (YDHDR%RUNDEF < 0) YDHDR%RUNDEF = YDHDR%RUNDEF + 256._JPRB

  CASE DEFAULT

END SELECT


YDHDR%IROWS = 40
YDHDR%ICOLS = 80


END SUBROUTINE 

SUBROUTINE COMPARE_GEOMHDR (YDHDR1, YDHDR2)

TYPE (HDR_t),  INTENT (IN) :: YDHDR1, YDHDR2

#include "abor1.intfb.h"

IF (YDHDR1%ZWEST  /= YDHDR2%ZWEST ) CALL ABOR1 ('COMPARE_GEOMHDR: DIFFERENT ZWEST ')   
IF (YDHDR1%ZEAST  /= YDHDR2%ZEAST ) CALL ABOR1 ('COMPARE_GEOMHDR: DIFFERENT ZEAST ')   
IF (YDHDR1%ZNORTH /= YDHDR2%ZNORTH) CALL ABOR1 ('COMPARE_GEOMHDR: DIFFERENT ZNORTH') 
IF (YDHDR1%ZSOUTH /= YDHDR2%ZSOUTH) CALL ABOR1 ('COMPARE_GEOMHDR: DIFFERENT ZSOUTH')
IF (YDHDR1%ICOLS  /= YDHDR2%ICOLS ) CALL ABOR1 ('COMPARE_GEOMHDR: DIFFERENT ICOLS ') 
IF (YDHDR1%IROWS  /= YDHDR2%IROWS ) CALL ABOR1 ('COMPARE_GEOMHDR: DIFFERENT IROWS ') 

END SUBROUTINE

SUBROUTINE SQUARE (KN, KA, KB)

INTEGER (KIND=JPIM) :: KN, KA, KB

KB = INT (SQRT (REAL (KN))) + 1

DO
  KA = KN / KB
  IF (KA * KB == KN) EXIT
  KB = KB - 1
ENDDO

END SUBROUTINE SQUARE

SUBROUTINE NEWFLD (CDPREF, CDSUFF, CDTYPE, KPASS, KRANK)

CHARACTER (LEN=*),   INTENT (IN)  :: CDPREF, CDSUFF, CDTYPE
INTEGER (KIND=JPIM), INTENT (IN)  :: KPASS
INTEGER (KIND=JPIM), INTENT (OUT) :: KRANK

#include "abor1.intfb.h"

SELECT CASE (CDTYPE)
  CASE ('INT')
    INFLD_INT = INFLD_INT + 1
    IF (KPASS > 1) THEN
      KRANK = IOFLD_INT + INFLD_INT
    ENDIF
  CASE ('CAL')
    INFLD_CAL = INFLD_CAL + 1
    IF (KPASS > 1) THEN
      KRANK = IOFLD_CAL + INFLD_CAL
    ENDIF
  CASE ('AUX')
    INFLD_AUX = INFLD_AUX + 1
    IF (KPASS > 1) THEN
      KRANK = IOFLD_AUX + INFLD_AUX
    ENDIF
  CASE DEFAULT
    CALL ABOR1 ('NEWFLD: UNEXPECTED CDTYPE')
END SELECT 

IF (KPASS > 1) THEN
  YLFLDS (KRANK)%CPREF = CDPREF
  YLFLDS (KRANK)%CSUFF = CDSUFF
  YLFLDS (KRANK)%CTYPE = CDTYPE
  YLFLDS (KRANK)%IRANK = KRANK
ENDIF

END SUBROUTINE

ELEMENTAL LOGICAL FUNCTION LLZZ (P)

REAL (KIND=JPRB), PARAMETER :: ZZ2 = 1.7E+308_JPRB, ZZ1 = 8.0E+307_JPRB

REAL (KIND=JPRB), INTENT (IN) :: P

LLZZ = (ZZ1 <= P) .AND. (P <= ZZ2)

END FUNCTION

END PROGRAM 

