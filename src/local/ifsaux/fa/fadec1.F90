! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FADEC1_MT64                                                &
&                     (FA,  KREP,   KNUMER, CDPREF, KNIVAU, CDSUFF,   &
&                      LDCOSP, CDNOMA, KLNOMA, KVALCO, KLONGD,        &
&                      PCHAMP, LDUNDF, PUNDF, YDGR1TAB)
USE FA_MOD, ONLY : FA_COM, FAGR1TAB
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!      Sous-programme de controle et de DECODAGE d'un CHAMP HORIZONTAL
!      venant d'etre lu sur un fichier ARPEGE/ALADIN.
!       ( DECOdage de donnees )
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KNUMER (Entree) ==> Numero de l'unite logique;
!                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
!                KNIVAU (Entree) ==> Niveau vertical eventuel;
!                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
!                LDCOSP (Entree) ==> Vrai si le champ est represente
!                                    par des coefficients spectraux;
!                CDNOMA (Sortie) ==> Nom de l'article-champ lu;
!                KLNOMA (Sortie) ==> Nombre de caracteres utiles dans
!                                    CDNOMA;
!    ( Tableau ) KVALCO (Entree) ==> Donnees issues de la lecture;
!                KLONGD (Entree) ==> Nombre de valeurs (mots de 64 bits
!                                    en principe) lues;
!    ( Tableau ) PCHAMP (Sortie) ==> Valeurs REELLES du champ lu.
!                LDUNDF (Sortie) ==> Vrai si ce champ a des valeurs 
!                                    indefinies
!                PUNDF  (Sortie) ==> Dans le cas ou LDUNDF est vrai,
!                                    valeur non definie
!
!    Remarques:
!
!    - KVALCO est type entier, et doit avoir une longueur
!      suffisante pour stocker les donnees codees. Le dimensionnement
!      "tous terrains" est (2+ILCHAM), qui permet le cas echeant de
!      stocker un champ a pleine resolution sans codage effectif.
!      (ILCHAM est le nombre de valeurs du champ a decoder)
!
!    - CDNOMA doit avoir au moins FA%JPXNOM caracteres.
!
!
TYPE(FA_COM)   :: FA
TYPE(FAGR1TAB) :: YDGR1TAB
INTEGER (KIND=JPLIKB) KREP, KNUMER, KNIVAU, KLNOMA, KLONGD
!
!
INTEGER (KIND=JPLIKB) IREP, ILPRFU, ILSUFU, ILNOMU
INTEGER (KIND=JPLIKB) IRANG, INIMES
INTEGER (KIND=JPLIKB) ILPREF, ILSUFF, ILCDNO, IRANGC, IVALC1
INTEGER (KIND=JPLIKB) IB1PAR (FA%JPLB1P)
!
REAL (KIND=JPDBLR) PCHAMP (*), PUNDF
INTEGER (KIND=JPLIKB) KVALCO(*)
!
LOGICAL LLVERF, LLRLFI, LDCOSP, LLNOMU, LDUNDF
!
CHARACTER CDPREF*(*), CDSUFF*(*), CDNOMA*(*)
CHARACTER CLPREF*(FA%JPXNOM), CLSUFF*(FA%JPXSUF)
!
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!**
!     1.  -  CONTROLES ET INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FADEC1_MT',0,ZHOOK_HANDLE)

LLVERF=.FALSE.
LLRLFI=.FALSE.
LLNOMU=.FALSE.
ILPRFU=INT (LEN (CDPREF), JPLIKB)
ILSUFU=INT (LEN (CDSUFF), JPLIKB)
ILCDNO=INT (LEN (CDNOMA), JPLIKB)
KLNOMA=0
CALL FANUMU_MT64                 &
&               (FA, KNUMER,IRANG)
!
IF (IRANG.EQ.0) THEN
  IREP=-51
  GOTO 1001
ELSEIF (ILCDNO.LT.FA%JPXNOM) THEN
  IREP=-65
  GOTO 1001
ELSE
  CDNOMA=' '
ENDIF
!
!         Verrouillage eventuel du fichier.
!
IF (FA%LFAMUL) CALL LFIVER_MT64                               &
&                              (FA%LFI, FA%FICHIER(IRANG)%VRFICH,'ON')
LLVERF=FA%LFAMUL
!
IF (FA%FICHIER(IRANG)%LCREAF) THEN
  IREP=-85
  GOTO 1001
ENDIF
!**
!     2.  -  FABRICATION DU NOM D'ARTICLE VIA LE SOUS-PROGRAMME "FANFAR"
!            ( controles de CDPREF, KNIVAU, CDSUFF inclus )
!-----------------------------------------------------------------------
!
CALL FANFAR_MT64                                             &
&               (FA, IREP,IRANG,CDPREF,KNIVAU,CDSUFF,CDNOMA, &
&                IB1PAR(6),ILPRFU,ILSUFU,ILNOMU)
IF (IREP.NE.0) GOTO 1001
LLNOMU=.TRUE.
KLNOMA=ILNOMU
!**
!     3.  -  CONTROLE ET DECODAGE DE L'ARTICLE DEJA LU SUR LE FICHIER.
!-----------------------------------------------------------------------
!
!
!  Controle de l'homogeneite du type de rangement des coeff. spectraux
!  parmi les champs lus/ecrits: ces champs compactes avec
!  FA%NIGRIB=-1 ou 3 doivent etre ranges comme dans le modele ("verticalement"
!  soit selon des colonnes JM=cst consecutives) et contrairement si compactes
!  avec FA%NIGRIB= 0,1 ou 2.
! 
IRANGC=FA%FICHIER(IRANG)%NUCADR
IVALC1=KVALCO(1)
IF (LDCOSP) THEN
  IF (IVALC1.EQ.-1.OR.IVALC1.EQ.3) THEN
    FA%FICHIER(IRANG)%NRASVE=FA%FICHIER(IRANG)%NRASVE+1
    IF (FA%FICHIER(IRANG)%NRASVE.EQ.1.AND.FA%FICHIER(IRANG)%NRASHO.GT.0) THEN
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
      WRITE(FA%NULOUT,*)' FADEC1 :  WARNING !!!!!           '
      WRITE(FA%NULOUT,*)' Un champ de coeff. spectraux avec'
      WRITE(FA%NULOUT,*)                            &
&      ' rangement type modele va etre lu alors que'
      WRITE(FA%NULOUT,*)                                 &
&      ' d''autres champs spect. ont un rangt different.'
      WRITE(FA%NULOUT,*)                                 &
&      ' ***  Prenez en compte cette heterogeneite!  ***'
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
    ENDIF
  ELSEIF (IVALC1.GE.0.AND.IVALC1.LE.2) THEN
    FA%FICHIER(IRANG)%NRASHO=FA%FICHIER(IRANG)%NRASHO+1
    IF (FA%FICHIER(IRANG)%NRASHO.EQ.1.AND.FA%FICHIER(IRANG)%NRASVE.GT.0) THEN
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
      WRITE(FA%NULOUT,*)' FADEC1 :  WARNING !!!!!           '
      WRITE(FA%NULOUT,*)' Un champ de coeff. spectraux avec'
      WRITE(FA%NULOUT,*)                                &
&      ' rangement autre que celui du modele va etre lu'
      WRITE(FA%NULOUT,*)                                &
&      ' alors que d''autres champs ont le rangt modele'
      WRITE(FA%NULOUT,*)                                 &
&      ' ***  Prenez en compte cette heterogeneite!  ***'
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
    ENDIF
  ENDIF
ENDIF
!
IF (FALGRA (IVALC1)) THEN
! Cas d'un champ gribe avec GRIB_API
  CALL FADGRA_MT64                                                 &
&                 (FA,  IREP,   IRANG,  CDNOMA(1:ILNOMU), KVALCO,  &
&                  KLONGD, PCHAMP, LDCOSP, CDPREF, KNIVAU, CDSUFF, &
&                  LDUNDF, PUNDF, .FALSE.)
ELSEIF (IVALC1.EQ.3) THEN
! Cas d'un champ gribe avec GRIBEX
  CALL FADECX_MT64                                                 &
&                 (FA,  IREP,   IRANG,  CDNOMA(1:ILNOMU), KVALCO,  &
&                  KLONGD, PCHAMP, LDCOSP, CDPREF, KNIVAU, CDSUFF, &
&                  LDUNDF, PUNDF, YDGR1TAB)
ELSEIF (IVALC1.EQ.4) THEN
  CALL FADCPL_MT64                                                 &
&                 (FA,  IREP,   IRANG,  CDNOMA(1:ILNOMU), KVALCO,  &
&                  KLONGD, PCHAMP, LDCOSP, LDUNDF, PUNDF)
ELSE
  CALL FADECI_MT64                                                 &
&                 (FA,  IREP,   IRANG,  CDNOMA(1:ILNOMU), KVALCO,  &
&                  KLONGD, PCHAMP, LDCOSP )
ENDIF
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
!            VIA LE SOUS-PROGRAMME "FAIPAR" .
!-----------------------------------------------------------------------
!
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
  IF (LHOOK) CALL DR_HOOK('FADEC1_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
CLNSPR='FADEC1'
!
IF (ILPRFU.GE.1) THEN
  ILPREF=MIN (ILPRFU,INT (LEN (CLPREF), JPLIKB))
  CLPREF(1:ILPREF)=CDPREF(1:ILPREF)
ELSE
  ILPREF=8
  CLPREF(1:ILPREF)=FA%CHAINC(:ILPREF)
ENDIF
!
IF (ILSUFU.GE.1) THEN
  ILSUFF=MIN (ILSUFU,INT (LEN (CLSUFF), JPLIKB))
  CLSUFF(1:ILSUFF)=CDSUFF(1:ILSUFF)
ELSE
  ILSUFF=8
  CLSUFF(1:ILSUFF)=FA%CHAINC(:ILSUFF)
ENDIF
!
IF (.NOT.LLNOMU) THEN
  ILNOMU=MIN (ILPREF,FA%NCPCAD)
  CDNOMA(1:ILNOMU)=CLPREF(1:ILPREF)
ENDIF
!
WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,         &
&       '', CDPREF='''''',A,'''''', KNIVAU='',I6,               &
&       '', CDSUFF='''''',A,'''''', LDCOSP= '',L1)')            &
&   KREP,KNUMER,CLPREF(1:ILPREF),KNIVAU,CLSUFF(1:ILSUFF),LDCOSP
CALL FAIPAR_MT64                                       &
&               (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS, &
&                CLNSPR,CDNOMA(1:ILNOMU),LLRLFI)
!
IF (LHOOK) CALL DR_HOOK('FADEC1_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"
#include "facom2.ixnvms.h"
#include "falgra.h"

END SUBROUTINE FADEC1_MT64
