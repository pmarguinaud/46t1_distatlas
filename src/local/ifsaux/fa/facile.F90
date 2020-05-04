! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FACILE_MT64                                             &
&                     (FA,  KREP, KNUMER, CDPREF, KNIVAU, CDSUFF,  &
&                      KCHAMP, LDCOSP)
USE FA_MOD, ONLY : FA_COM, FAGR1TAB
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!      Sous-programme de LECTURE d'un CHAMP HORIZONTAL sur un fichier
!     ARPEGE.
!       ( Champ d'Interet en LEcture )
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KNUMER (Entree) ==> Numero de l'unite logique;
!                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
!                KNIVAU (Entree) ==> Niveau vertical eventuel;
!                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
!    ( Tableau ) KCHAMP (Sortie) ==> Valeurs REELLES du champ lu;
!                LDCOSP (Entree) ==> Vrai si le champ est represente
!                                    par des coefficients spectraux.
!     MODIF:
!     JM AUDOIN GMAP/EXT 10/05/95 intro de IVALC3 pour eviter ecrasement
!     D  PARADIS TTI/DEV 12/10/98 partie controle et decodage de l'article
!                                 demenagee dans un ss-prg a usage interne
!                                 du logiciel (FADECI).
!     D  PARADIS DSI/DEV 15/04/04 nettoyage code + declaration IVALCO en
!                                 ALLOCATABLE
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KREP, KNUMER, KNIVAU
!
CHARACTER CDPREF*(*), CDSUFF*(*)
!
INTEGER (KIND=JPLIKB) IUNDF
!
INTEGER (KIND=JPLIKB) KCHAMP (*)
!
TYPE (FAGR1TAB) YLGR1TAB
!
LOGICAL LDCOSP, LLUNDF

!**
!     1.  -  CONTROLES ET INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FACILE_MT',0,ZHOOK_HANDLE)

CALL FACIL1_MT64 (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, KCHAMP, LDCOSP, &
                & LLUNDF, IUNDF, YLGR1TAB, .TRUE.)


IF (LHOOK) CALL DR_HOOK('FACILE_MT',1,ZHOOK_HANDLE)

END SUBROUTINE FACILE_MT64

! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FACILE64                                        &
&           (KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, KCHAMP, &
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
INTEGER (KIND=JPLIKB)  KCHAMP     (*)                         !   OUT
LOGICAL                LDCOSP                                 ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FACILE_MT64                                               &
&           (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, KCHAMP, &
&            LDCOSP)

END SUBROUTINE FACILE64

SUBROUTINE FACILE                                          &
&           (KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, KCHAMP, &
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
INTEGER (KIND=JPLIKB)  KCHAMP     (*)                         !   OUT
LOGICAL                LDCOSP                                 ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FACILE_MT                                                 &
&           (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, KCHAMP, &
&            LDCOSP)

END SUBROUTINE FACILE

SUBROUTINE FACILE_MT                                           &
&           (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, KCHAMP, &
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
INTEGER (KIND=JPLIKB)  KCHAMP     (*)                         !   OUT
LOGICAL                LDCOSP                                 ! IN   
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  INIVAU                                 ! IN   
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)
INIVAU     = INT (    KNIVAU, JPLIKB)

CALL FACILE_MT64                                               &
&           (FA, IREP, INUMER, CDPREF, INIVAU, CDSUFF, KCHAMP, &
&            LDCOSP)

KREP       = INT (      IREP, JPLIKM)

END SUBROUTINE FACILE_MT

!INTF KREP            OUT                                                              
!INTF KNUMER        IN                                                                 
!INTF CDPREF        IN                                                                 
!INTF KNIVAU        IN                                                                 
!INTF CDSUFF        IN                                                                 
!INTF KCHAMP          OUT DIMS=*                         KIND=JPLIKB                   
!INTF LDCOSP        IN                                                                 

