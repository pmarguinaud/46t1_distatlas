! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FAIENO_MT64                                           &
&                     (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, &
&                      PCHAMP, LDCOSP, LDUNDF, PUNDF, LDREVERT)
USE FA_MOD, ONLY : FA_COM, FAGR1TAB
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!      Sous-programme d'ECRITURE d'un CHAMP HORIZONTAL sur un fichier
!     ARPEGE, avec reordonnement des coefficients spectraux si necessaire.
!       ( Integration par Ecriture d'un (Nouveau ?) Champ )
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KNUMER (Entree) ==> Numero de l'unite logique;
!                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
!                KNIVAU (Entree) ==> Niveau vertical eventuel;
!                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
!    ( Tableau ) PCHAMP (Entree) ==> Valeurs REELLES du champ a ecrire;
!                                    rangement modele.
!                LDCOSP (Entree) ==> Vrai si le champ est represente
!                                    par des coefficients spectraux.
!                LDUNDF (Entree) ==> Vrai si ce champ a des valeurs 
!                                    indefinies
!                PUNDF  (Entree) ==> Dans le cas ou LDUNDF est vrai,
!                                    valeur non definie
!
TYPE(FA_COM)           FA
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KNUMER                                 ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKB)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
REAL (KIND=JPDBLR)     PCHAMP     (*)                         ! IN   
LOGICAL                LDCOSP                                 ! IN   
LOGICAL,               OPTIONAL :: LDUNDF                     ! IN
REAL (KIND=JPDBLR),    OPTIONAL :: PUNDF                      ! IN
LOGICAL,               OPTIONAL :: LDREVERT                   ! IN
!
INTEGER (KIND=JPLIKB) IREP
INTEGER (KIND=JPLIKB) IRANG, INIMES, INGRIB
INTEGER (KIND=JPLIKB) ISMAX, IMSMAX
!
INTEGER (KIND=JPLIKB) IRANGC
!
LOGICAL LLVERF, LLRLFI, LLREORD
!
REAL (KIND=JPDBLR), ALLOCATABLE :: ZCHAMP (:)
LOGICAL               :: LLUNDF                   
REAL (KIND=JPDBLR)    :: ZUNDF                    
LOGICAL               :: LLREVERT
TYPE (FAGR1TAB)       :: YLGR1TAB
!
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!**
!     1.  -  CONTROLES ET INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FAIENO_MT',0,ZHOOK_HANDLE)

LLUNDF = .FALSE.
IF (PRESENT (LDUNDF  )) LLUNDF   = LDUNDF
ZUNDF  = 0._JPDBLR
IF (PRESENT (PUNDF   )) ZUNDF    = PUNDF 
LLREVERT = .TRUE.
IF (PRESENT (LDREVERT)) LLREVERT = LDREVERT


IREP=0
LLVERF=.FALSE.
LLRLFI=.FALSE.
!
CALL FANUMU_MT64                 &
&               (FA, KNUMER,IRANG)
!
IF (IRANG.EQ.0) THEN
  IREP=-51
  GOTO 1001
ENDIF
!
!         Verrouillage eventuel du fichier.
!
IF (FA%LFAMUL) CALL LFIVER_MT64                               &
&                              (FA%LFI, FA%FICHIER(IRANG)%VRFICH,'ON')
LLVERF=FA%LFAMUL
!
IRANGC=FA%FICHIER(IRANG)%NUCADR
INGRIB=FA%FICHIER(IRANG)%NFGRIB

LLREORD = LDCOSP .AND. (.NOT.(INGRIB==-1 .OR. INGRIB==3 .OR. FALGRA (INGRIB)))

IF (LLREORD) THEN
  ISMAX  = FA%CADRE(IRANGC)%NSMAX     
  IMSMAX = FA%CADRE(IRANGC)%NMSMAX     
  ALLOCATE (ZCHAMP (4 * (IMSMAX+1) * (ISMAX+1))) ! Assez grand
  CALL FAREOR_MT64 (FA, IREP, KNUMER, PCHAMP, ZCHAMP, .FALSE.)
  IF (IREP /= 0) GOTO 1001
  CALL FAIEN1_MT64 (FA, IREP, KNUMER, CDPREF, KNIVAU, CDSUFF, ZCHAMP, LDCOSP, &
                  & LLUNDF, ZUNDF, YLGR1TAB, LLREVERT)
  IF (IREP /= 0) GOTO 1001
  DEALLOCATE (ZCHAMP)
ELSE
  CALL FAIEN1_MT64 (FA, IREP, KNUMER, CDPREF, KNIVAU, CDSUFF, PCHAMP, LDCOSP, &
                  & LLUNDF, ZUNDF, YLGR1TAB, LLREVERT)
ENDIF

1001 CONTINUE
KREP=IREP
LLFATA=LLMOER (IREP,IRANG)
!
!        Deverrouillage eventuel du fichier.
!
IF (LLVERF) CALL LFIVER_MT64                                &
&                           (FA%LFI, FA%FICHIER(IRANG)%VRFICH,'OFF')
!
IF (LLFATA) THEN
  INIMES=2
ELSE
  INIMES=IXNVMS(IRANG)
ENDIF
!
IF (.NOT.LLFATA.AND.INIMES.NE.2)  THEN 
  IF (LHOOK) CALL DR_HOOK('FAIENO_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
CLNSPR='FAIENO'
!
WRITE (UNIT=CLMESS,FMT='(''KREP='',I5,'', KNUMER='',I3,         &
&       '', CDPREF='''''',A,'''''', KNIVAU='',I6,               &
&       '', CDSUFF='''''',A,'''''', LDCOSP= '',L1)')            &
&   KREP,KNUMER,TRIM (CDPREF),KNIVAU,TRIM (CDSUFF),LDCOSP
CALL FAIPAR_MT64                                       &
&               (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS, &
&                CLNSPR, '',LLRLFI)
!

IF (LHOOK) CALL DR_HOOK('FAIENO_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"
#include "facom2.ixnvms.h"
#include "falgra.h"

END SUBROUTINE FAIENO_MT64



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FAIENO64 (KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, PCHAMP, LDCOSP, &
                   & LDUNDF, PUNDF, LDREVERT)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                  FA_COM_DEFAULT_INIT,  &
&                  NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KNUMER                                 ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKB)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
REAL (KIND=JPDBLR)     PCHAMP     (*)                         ! IN   
LOGICAL                LDCOSP                                 ! IN   
LOGICAL,               OPTIONAL :: LDUNDF                     ! IN
REAL (KIND=JPDBLR),    OPTIONAL :: PUNDF                      ! IN
LOGICAL,               OPTIONAL :: LDREVERT                   ! IN

#include "faieno_mt64.h"

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAIENO_MT64 (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, PCHAMP, LDCOSP, &
                & LDUNDF, PUNDF, LDREVERT)

END SUBROUTINE FAIENO64

SUBROUTINE FAIENO                                          &
&           (KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, PCHAMP, &
&            LDCOSP, LDUNDF, PUNDF, LDREVERT)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                  FA_COM_DEFAULT_INIT,  &
&                  NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKM)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
REAL (KIND=JPDBLR)     PCHAMP     (*)                         ! IN   
LOGICAL                LDCOSP                                 ! IN   
LOGICAL,               OPTIONAL :: LDUNDF                     ! IN
REAL (KIND=JPDBLR),    OPTIONAL :: PUNDF                      ! IN
LOGICAL,               OPTIONAL :: LDREVERT                   ! IN

#include "faieno_mt.h"

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAIENO_MT (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, PCHAMP, LDCOSP, &
              & LDUNDF, PUNDF, LDREVERT)

END SUBROUTINE FAIENO

SUBROUTINE FAIENO_MT (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, PCHAMP, LDCOSP, &
                    & LDUNDF, PUNDF, LDREVERT)
USE FA_MOD, ONLY : FA_COM
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKM)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
REAL (KIND=JPDBLR)     PCHAMP     (*)                         ! IN   
LOGICAL                LDCOSP                                 ! IN   
LOGICAL,               OPTIONAL :: LDUNDF                     ! IN
REAL (KIND=JPDBLR),    OPTIONAL :: PUNDF                      ! IN
LOGICAL,               OPTIONAL :: LDREVERT                   ! IN

#include "faieno_mt64.h"

! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  INIVAU                                 ! IN   
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)
INIVAU     = INT (    KNIVAU, JPLIKB)

CALL FAIENO_MT64 (FA, IREP, INUMER, CDPREF, INIVAU, CDSUFF, PCHAMP, LDCOSP, &
                & LDUNDF, PUNDF, LDREVERT)

KREP       = INT (      IREP, JPLIKM)

END SUBROUTINE FAIENO_MT

!INTF KREP            OUT                               
!INTF KNUMER        IN                                  
!INTF CDPREF        IN                                  
!INTF KNIVAU        IN                                  
!INTF CDSUFF        IN                                  
!INTF PCHAMP        IN    DIMS=*                        
!INTF LDCOSP        IN                                  
!INTF LDUNDF        IN                                  
!INTF PUNDF         IN                                  

