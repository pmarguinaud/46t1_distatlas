! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FAIENC_MT64                                              &
&                     (FA,  KREP,   KNUMER, CDPREF, KNIVAU, CDSUFF, &
&                      PCHAMP, LDCOSP)
USE FA_MOD, ONLY : FA_COM, FAGR1TAB
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!      Sous-programme d'ECRITURE d'un CHAMP HORIZONTAL sur un fichier
!     ARPEGE.
!       ( Integration par Ecriture d'un (Nouveau ?) Champ )
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KNUMER (Entree) ==> Numero de l'unite logique;
!                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
!                KNIVAU (Entree) ==> Niveau vertical eventuel;
!                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
!    ( Tableau ) PCHAMP (Entree) ==> Valeurs REELLES du champ a ecrire;
!                LDCOSP (Entree) ==> Vrai si le champ est represente
!                                    par des coefficients spectraux.
!
!     Modifications
!     -------------
!
!    Avril 1998: Partie "codage" (paragraphe 3 du sous-programme)
!                demenagee dans un sous-programme a usage interne au
!                logiciel (FACINE). Le but est de pouvoir, sur machine
!                a memoire distribuee, separer codage (via FACOND) et
!                ecriture (via FAISAN) afin de paralleliser le codage.
!
!  Avril 2004, D. Paradis, DSI/DEV:
!
!    -Declaration IVALCO en ALLOCATABLE (gain memoire)
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KREP, KNUMER, KNIVAU
!
REAL (KIND=JPDBLR) PCHAMP (*), ZUNDF
!
CHARACTER CDPREF*(*), CDSUFF*(*)
!
LOGICAL LDCOSP, LLUNDF
TYPE (FAGR1TAB) YLGR1TAB
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FAIENC_MT',0,ZHOOK_HANDLE)

LLUNDF = .FALSE.
ZUNDF  = 0._JPDBLR

CALL FAIEN1_MT64 (FA,  KREP,   KNUMER, CDPREF, KNIVAU, CDSUFF, &
                & PCHAMP, LDCOSP, LLUNDF, ZUNDF, YLGR1TAB, .TRUE.)

IF (LHOOK) CALL DR_HOOK('FAIENC_MT',1,ZHOOK_HANDLE)

END SUBROUTINE FAIENC_MT64

! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FAIENC64                                        &
&           (KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, PCHAMP, &
&            LDCOSP)
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

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAIENC_MT64                                               &
&           (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, PCHAMP, &
&            LDCOSP)

END SUBROUTINE FAIENC64

SUBROUTINE FAIENC                                          &
&           (KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, PCHAMP, &
&            LDCOSP)
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

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FAIENC_MT                                                 &
&           (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, PCHAMP, &
&           LDCOSP)

END SUBROUTINE FAIENC

SUBROUTINE FAIENC_MT                                           &
&           (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, PCHAMP, &
&            LDCOSP)
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
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  INIVAU                                 ! IN   
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)
INIVAU     = INT (    KNIVAU, JPLIKB)

CALL FAIENC_MT64                                               &
&           (FA, IREP, INUMER, CDPREF, INIVAU, CDSUFF, PCHAMP, &
&            LDCOSP)

KREP       = INT (      IREP, JPLIKM)

END SUBROUTINE FAIENC_MT

!INTF KREP            OUT                               
!INTF KNUMER        IN                                  
!INTF CDPREF        IN                                  
!INTF KNIVAU        IN                                  
!INTF CDSUFF        IN                                  
!INTF PCHAMP        IN    DIMS=*                        
!INTF LDCOSP        IN                                  


