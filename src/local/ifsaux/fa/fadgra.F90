! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FADGRA_MT64                                          &
&                     (FA, KREP, KRANG, CDNOMA, KVALCO, KLONGA, &
&                      PCHAMP, LDCOSP, CDPREF, KNIVAU, CDSUFF,  &
&                      LDUNDF, PUNDF, LDREVERT)
USE FA_MOD, ONLY : FA_COM, JPNIIL, FAFICH, FACADR
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
USE GRIB_API_INTERFACE
USE GRIB_API
IMPLICIT NONE
!****
!      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KRANG  (Entree) ==> Rang de l'unite logique;
!                CDNOMA (Entree) ==> Nom d'article (prefabrique);
!    ( Tableau ) KVALCO (Entree) ==> Donnees issues de la lecture;
!                KLONGA (Entree) ==> Nombre de mots lus;
!    ( Tableau ) PCHAMP (Sortie) ==> Valeurs REELLES du champ lu;
!                LDCOSP (Entree) ==> Vrai si le champ est represente
!                                    par des coefficients spectraux;
!                CDPREF (Entree) ==> Prefixe au sens FA
!                KNIVAU (Entree) ==> Niveau au sens FA
!                CDSUFF (Entree) ==> Suffixe au sens FA
!                LDUNDF (Entree) ==> Si ce champ a des valeurs indefinies
!                                    alors inserer PUNDF sur les points
!                                    manquants
!                PUNDF  (Entree) ==> Dans le cas ou LDUNDF est vrai,
!                                    valeur non definie a inserer dans le champ
!                LDUNDF (Sortie) ==> Vrai si ce champ a des valeurs 
!                                    indefinies
!                PUNDF  (Sortie) ==> Dans le cas ou LDUNDF est vrai (en sortie),
!                                    valeur non definie a inserer dans le champ
!
!
TYPE(FA_COM)   :: FA
INTEGER (KIND=JPLIKB) KREP, KRANG, KLONGA, KNIVAU
!
INTEGER (KIND=JPLIKB), TARGET :: KVALCO(KLONGA)
REAL (KIND=JPDBLR) PCHAMP(*)
!
REAL (KIND=JPDBLR) PUNDF
!
LOGICAL LDCOSP, LDUNDF, LLUNDF, LLLTLN, LDREVERT
!
CHARACTER CDNOMA*(*), CDPREF*(*), CDSUFF*(*)
!
REAL (KIND=JPDBLR) ZUNDF
!
INTEGER (KIND=JPLIKB) ILCHAM
INTEGER (KIND=JPLIKB) INIMES
INTEGER (KIND=JPLIKB) INUMER
!
LOGICAL LLMLAM, LLCOSP, LLMGLO
!
TYPE (FAFICH), POINTER :: YLFICH
TYPE (FACADR), POINTER :: YLCADR
!
CHARACTER, ALLOCATABLE :: CLGRIB (:)
!
REAL (KIND=JPDBLR), ALLOCATABLE :: ZCHAMP (:)
!
INTEGER (KIND=JPLIKB) ILGRIB, IRANGC
INTEGER (KIND=JPLIKB) JLAT, JLON, JN, IDX, J
INTEGER (KIND=JPLIKM) IGRIBH, IRET, IBITMAP, INDATV, IBTMP
CHARACTER(LEN=FA%JPLSPX)   CLNSPR
CHARACTER(LEN=FA%JPLMES)   CLMESS 
CHARACTER(LEN=FA%JPXNOM)   CLNOMU
!
INTEGER (KIND=JPLIKB)    IMULTM, IMULTE
LOGICAL                  LLFATA
REAL (KIND=JPDBLR)       ZMULTI
INTEGER                  IEDITION, IPARAM
LOGICAL                  LLLOCSEC, LLGRIB1
INTEGER (KIND=JPLIKB)    IVERSI, INGRIB, ITRONC, ISUBTR, IDECSF
INTEGER (KIND=JPLIKB)    II
LOGICAL                  LLBUG_SH_DEC

!**
!     1.  -  CONTROLES ET INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FADGRA_MT',0,ZHOOK_HANDLE)

KREP = 0

YLFICH => FA%FICHIER(KRANG)
IRANGC = YLFICH%NUCADR
YLCADR => FA%CADRE(IRANGC)
!
LLMLAM = YLCADR%LIMLAM
LLLTLN = YLCADR%SINLAT(2) < 0 .AND. LLMLAM
LLMGLO = (.NOT. LLLTLN) .AND. (.NOT. LLMLAM)
!
LLCOSP = LDCOSP
!
INUMER=YLFICH%NULOGI
!
!**
!     2.  -  CONTROLE DES DONNEES DE L'ARTICLE
!-----------------------------------------------------------------------
!
IF ((.NOT. FALGRA (KVALCO(1))).OR.                     &
&   KVALCO(2).LT.0.OR.KVALCO(2).GT.1) THEN
  KREP=-91
  GOTO 1001
ELSE
  LLCOSP=KVALCO(2).EQ.1
ENDIF
!
IF ((LLCOSP.AND..NOT.LDCOSP).OR.(.NOT.LLCOSP.AND.LDCOSP)) THEN
  KREP=-92
  GOTO 1001
ENDIF
!
IF (LLCOSP) THEN
  IF (LLMLAM) THEN
    ILCHAM = YLCADR%NSFLAM
  ELSE    
    ILCHAM =(1+YLCADR%MTRONC)*(2+YLCADR%MTRONC)
  ENDIF   
ELSE
  ILCHAM=YLCADR%NVAPDG
ENDIF
!
!**
!     3.  -  DECODAGE GRIB_API DES DONNEES DE L'ARTICLE
!-----------------------------------------------------------------------
!

ILGRIB = (KLONGA-3)*8

ALLOCATE (CLGRIB (ILGRIB))
CLGRIB = TRANSFER (KVALCO (4:KLONGA), CLGRIB)
CALL GRIB_NEW_FROM_MESSAGE_CHAR (IGRIBH, CLGRIB, STATUS=IRET)
DEALLOCATE (CLGRIB)

CALL IGRIB_GET_VALUE (IGRIBH, 'editionNumber', IEDITION)
LLGRIB1 = IEDITION == 1

IF (LLGRIB1) THEN
  CALL IGRIB_SET_VALUE (IGRIBH, 'setLocalDefinition', 1)
  LLLOCSEC = .FALSE.
ELSE
! Scaling factor may be encoded in the message
  CALL IGRIB_IS_DEFINED (IGRIBH, 'grib2LocalSectionNumber', LLLOCSEC)
  IF (.NOT. LLLOCSEC) THEN
    CALL IGRIB_SET_VALUE (IGRIBH, 'grib2LocalSectionPresent', 1)
    CALL IGRIB_SET_VALUE (IGRIBH, 'grib2LocalSectionNumber', 1)
  ENDIF
ENDIF

IF (.NOT. LLLOCSEC) THEN
! Restore encoding parameters if needed
  IF (KNIVAU > 0) THEN
    CALL STRU (CDSUFF, CLNOMU)
  ELSE
    CALL STRU (CDNOMA, CLNOMU)
  ENDIF
  CALL IGRIB_SET_VALUE (IGRIBH, 'faFieldName', TRIM (CLNOMU), IRET)
ENDIF

IF (LLGRIB1) THEN
  CALL IGRIB_GET_VALUE (IGRIBH, 'indicatorOfParameter', IPARAM)
ELSE
  CALL IGRIB_GET_VALUE (IGRIBH, 'parameterNumber', IPARAM)
ENDIF

IF (IPARAM == 255) THEN
  WRITE (FA%NULOUT, '(" FADGRA: Field `",A,"'' is not &
       &declared in `faFieldName.def'' and has no encoded &
       &FMULTM and FMULTE")') TRIM (CDNOMA)
ENDIF

CALL IGRIB_GET_VALUE (IGRIBH, 'FMULTM', IMULTM, IRET)
IF (IRET /= 0) IMULTM = 1

CALL IGRIB_GET_VALUE (IGRIBH, 'FMULTE', IMULTE, IRET)
IF (IRET /= 0) IMULTE = 0

ZMULTI = REAL (IMULTM, JPDBLR) * 10._JPDBLR ** IMULTE

CALL IGRIB_GET_VALUE (IGRIBH, 'bitmapPresent', IBTMP)
IF (IBTMP == 0) THEN
! When there is not bitmap, numberOfDataPoints may be broken; in this case we
! use numberOfValues
  CALL IGRIB_GET_VALUE (IGRIBH, 'numberOfValues', INDATV) 
ELSE
! numberOfValues = number of non-missing values when a bitmap is present
  CALL IGRIB_GET_VALUE (IGRIBH, 'numberOfDataPoints', INDATV) 
ENDIF

! Basic check

IF (INDATV < ILCHAM) THEN
  KREP=-93
  GOTO 1001
ELSEIF (INDATV > ILCHAM) THEN
  KREP=-94
  GOTO 1001
ENDIF

IF (LLLTLN .AND. LDREVERT) THEN

  ALLOCATE (ZCHAMP (ILCHAM))
  CALL IGRIB_GET_VALUE (IGRIBH, 'values', ZCHAMP (1:ILCHAM))
  DO JLAT = 1, YLCADR%NLATIT
    DO JLON = 1, YLCADR%NXLOPA
      JN  = JLON+YLCADR%NXLOPA*(JLAT-1)
      IDX = JLON+YLCADR%NXLOPA*(YLCADR%NLATIT-JLAT)
      PCHAMP (JN) = ZCHAMP (IDX)
    ENDDO
  ENDDO
  DEALLOCATE (ZCHAMP)

ELSE
  CALL IGRIB_GET_VALUE (IGRIBH, 'values', PCHAMP (1:ILCHAM))
ENDIF

CALL IGRIB_GET_VALUE (IGRIBH, 'bitmapPresent', IBITMAP)
LLUNDF = IBITMAP /= 0

IF (LLUNDF) THEN
  CALL IGRIB_GET_VALUE (IGRIBH, 'missingValue',  ZUNDF)
ENDIF

INGRIB = KVALCO (1)
CALL GRIB_GET_API_VERSION (IVERSI)

LLBUG_SH_DEC = (INGRIB == 101) .AND. (IVERSI /= 11400) .AND. LDCOSP

IF (LLBUG_SH_DEC) THEN
  CALL IGRIB_GET_VALUE (IGRIBH, 'subSetJ',            ISUBTR)
  CALL IGRIB_GET_VALUE (IGRIBH, 'J',                  ITRONC)
  CALL IGRIB_GET_VALUE (IGRIBH, 'decimalScaleFactor', IDECSF)
  LLBUG_SH_DEC = (ISUBTR /= ITRONC) .AND. (IDECSF /= 0)
ENDIF

CALL IGRIB_RELEASE (IGRIBH)

IF (LLBUG_SH_DEC) CALL FIXBUG_SH_DEC 

!
! Facteur d'echelle eventuel
!
IF (ZMULTI /= REAL (1._4, JPDBLR)) THEN
  PCHAMP (1:ILCHAM) = PCHAMP (1:ILCHAM) / ZMULTI
  ZUNDF             = ZUNDF             / ZMULTI
ENDIF
!
IF (LDUNDF .AND. LLUNDF) THEN
  DO J = 1, ILCHAM
    IF (PCHAMP (J) == ZUNDF) THEN
      PCHAMP (J) = PUNDF
    ENDIF
  ENDDO
  ZUNDF = PUNDF
ENDIF
!
LDUNDF = LLUNDF
PUNDF  = ZUNDF

!**
!    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
!            VIA LE SOUS-PROGRAMME "FAIPAR" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
LLFATA=LLMOER (KREP,KRANG)
!
IF (FA%LFAMOP.OR.LLFATA) THEN
  INIMES=2
  CLNSPR='FADGRA'
  INUMER=YLFICH%NULOGI
!
  WRITE (UNIT=CLMESS,FMT='(''KREP='',I5,'', KRANG='',I4,  &
&         '', CDNOMA='''''',A,'''''', KLONGA= '',I8,      &
&         '', LDCOSP='',L1)')                             &
&     KREP, KRANG, CDNOMA, KLONGA, LDCOSP
  CALL FAIPAR_MT64                                        &
&                 (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                  CLNSPR,CDNOMA,.FALSE.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('FADGRA_MT',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE FIXBUG_SH_DEC 

INTEGER (KIND=JPLIKB) :: ISP, JN, JM, IM
REAL (KIND=JPDBLR) :: ZDEC

ZDEC = 10._JPDBLR ** (-IDECSF)

DO JN = 0, ITRONC
  DO JM = -JN, JN
    IM = ABS (JM)
    IF (JM < 0) THEN
      ISP = FA%CADRE(IRANGC)%NDIM0GG (IM) + (JN - IM) * 2 + 1
    ELSE
      ISP = FA%CADRE(IRANGC)%NDIM0GG (IM) + (JN - IM) * 2
    ENDIF
    IF ((JN <= ISUBTR) .AND. (JM <= ISUBTR)) THEN
      PCHAMP (ISP) = PCHAMP (ISP) * ZDEC
    ENDIF
  ENDDO
ENDDO

END SUBROUTINE

#include "facom2.llmoer.h"
#include "falgra.h"

SUBROUTINE STRU (CDS, CDU)
CHARACTER (LEN=*) :: CDS, CDU
INTEGER (KIND=JPLIKB) :: J

DO J = 1, LEN (CDU)
  CDU (J:J) = ' '
ENDDO

DO J = 1, LEN_TRIM (CDS)
  IF (CDS (J:J) == ' ') THEN
    CDU (J:J) = '_'
  ELSE
    CDU (J:J) = CDS (J:J)
  ENDIF
ENDDO

END SUBROUTINE STRU

END SUBROUTINE

!INTF KREP            OUT                                                              
!INTF KRANG         IN                                                                 
!INTF CDNOMA        IN                                                                 
!INTF KVALCO        IN    DIMS=*                         KIND=JPLIKB                   
!INTF KLONGA        IN                                                                 
!INTF PCHAMP          OUT DIMS=*                                                       
!INTF LDCOSP        IN                                                                 
!INTF LDUNDF          OUT                                                              
!INTF PUNDF           OUT                                                              

