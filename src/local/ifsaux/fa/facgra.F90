SUBROUTINE FACGRA_MT64                                              &
&                     (FA,  KREP,   KRANG,  CDPREF, KNIVAU, CDSUFF, &
&                      PCHAMP, LDCOSP, KVALCO, KLONGD,              &
&                      LDUNDF, PUNDF, LDREVERT)
USE FA_MOD, ONLY : FA_COM, JPNIIL, FACADR, FAFICH,                  &
                 & NGRIB2_GLO_SH, NGRIB2_GLO_GP, NGRIB2_LAM_GP,     &
                 & NGRIB2_LAM_BF, NGRIB2_LATLON, NGRIB1_LATLON,     &
                 & NUNDEF
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
USE GRIB_API_INTERFACE
USE GRIB_API
IMPLICIT NONE
!****
!      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
!      PREPARATION (codage GRIB_API) d'un CHAMP HORIZONTAL
!      destine a etre ecrit sur un fichier ARPEGE/ALADIN.
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KRANG  (Entree) ==> Rang de l'unite logique;
!                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
!                KNIVAU (Entree) ==> Niveau vertical eventuel;
!                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
!    ( Tableau ) PCHAMP (Entree) ==> Valeurs REELLES du champ a ecrire;
!                LDCOSP (Entree) ==> Vrai si le champ est represente
!                                    par des coefficients spectraux;
!    ( Tableau ) KVALCO (Sortie) ==> Donnees destinees a l'ecriture;
!                KLONGD (Entree/Sortie) 
!                                ==> Nombre de mots a ecrire;
!*
!
TYPE(FA_COM)   :: FA
INTEGER (KIND=JPLIKB) KREP, KRANG, KNIVAU, KLONGD, ILONGD
!
INTEGER (KIND=JPLIKB) KVALCO(*)
REAL (KIND=JPDBLR), TARGET :: PCHAMP(*)
REAL (KIND=JPDBLR) PUNDF, ZUNDF
!
LOGICAL LDCOSP, LDUNDF, LLFATA, LDREVERT
!
CHARACTER CDPREF*(*), CDSUFF*(*)
!
CHARACTER(LEN=FA%JPXNOM)   CLACTI 
CHARACTER(LEN=FA%JPLSPX)   CLNSPR
CHARACTER(LEN=FA%JPLMES)   CLMESS 
INTEGER (KIND=JPLIKB) :: INIMES, INUMER
CHARACTER, ALLOCATABLE :: CLGRIB (:)
INTEGER (KIND=JPKSIZE_T) :: ILGRIB
INTEGER (KIND=JPLIKM) :: IRET, IGRIBH
INTEGER (KIND=JPLIKB) :: INGRIB, INBITS
INTEGER :: J

!
REAL (KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FACGRA_MT',0,ZHOOK_HANDLE)

KREP = 0

CALL FACGRM_MT64 (FA, KREP, KRANG, CDPREF, KNIVAU, CDSUFF, PCHAMP, &
                & LDCOSP, IGRIBH, LDUNDF, PUNDF, 1_JPLIKB, LDREVERT)

IF (KREP /= 0) GOTO 1001

CALL IGRIB_GET_VALUE (IGRIBH, 'INGRIB', INGRIB)
CALL IGRIB_GET_VALUE (IGRIBH, 'INBITS', INBITS)

IF (FA%FICHIER(KRANG)%NCOGRIF(12) == 1) THEN
  CALL IGRIB_SET_VALUE (IGRIBH, 'setLocalDefinition', 0)
ENDIF

CALL IGRIB_GET_MESSAGE_SIZE (IGRIBH, ILGRIB)

! Make ILGRIB a multiple of JPLIKB

ILGRIB = (1+(ILGRIB-1) / JPLIKB) * JPLIKB

ALLOCATE (CLGRIB (ILGRIB))

CLGRIB (ILGRIB-7:ILGRIB) = '        '

CALL GRIB_COPY_MESSAGE (IGRIBH, CLGRIB, STATUS=IRET)

IF (IRET == GRIB_SUCCESS) THEN
  ILONGD = 3+ILGRIB/JPLIKB
  KVALCO (1) = INGRIB
  IF (LDCOSP) THEN
    KVALCO (2) = 1
  ELSE
    KVALCO (2) = 0
  ENDIF
  KVALCO (3) = INBITS
  IF ((KLONGD < ILONGD) .AND. (KLONGD > 0)) THEN
    KREP=-130
    GOTO 1001
  ELSE
    KLONGD = ILONGD
  ENDIF

  DO J = 4, ILONGD
    KVALCO (J) = TRANSFER (CLGRIB (1+(J-4)*8:(J-3)*8), KVALCO (J))
  ENDDO

ELSE
  KREP = IRET-1000
  RETURN
ENDIF

DEALLOCATE (CLGRIB)

CALL IGRIB_RELEASE (IGRIBH)

1001 CONTINUE
!
LLFATA=LLMOER (KREP,KRANG)
!
IF (FA%LFAMOP.OR.LLFATA) THEN
  INIMES=2
  CLNSPR='FACGRA'
  INUMER=JPNIIL
!
  WRITE (UNIT=CLMESS,FMT='(''KREP='',I5,'', KRANG='',I4,  &
&         '', CDPREF='''''',A,'''''', KNIVAU='',I6,       &
&         '', CDSUFF='''''',A,'''''', LDCOSP= '',L1,      &
&         '', KLONGD='',I10,'' < '',I10)')                &
&     KREP, KRANG, CDPREF(1:LEN_TRIM(CDPREF)), KNIVAU,    &
&     CDSUFF(1:LEN_TRIM(CDSUFF)), LDCOSP, KLONGD, ILONGD

  CALL FAIPAR_MT64                                        &
&                 (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                  CLNSPR,CLACTI,.FALSE.)
ENDIF

IF (LHOOK) CALL DR_HOOK('FACGRA_MT',1,ZHOOK_HANDLE)
!
CONTAINS

#include "facom2.llmoer.h"

END SUBROUTINE

