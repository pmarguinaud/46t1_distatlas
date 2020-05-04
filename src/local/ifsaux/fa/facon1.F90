! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FACON1_MT64                                              &
&                     (FA,  KREP,   KNUMER, CDPREF, KNIVAU, CDSUFF, &
&                      PCHAMP, LDCOSP, CDNOMA, KLNOMA, KVALCO,      &
&                      KLONGD, LDUNDF, PUNDF, YDGR1TAB)
USE FA_MOD, ONLY : FA_COM, FAGR1TAB
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!      Sous-programme de CODAGE d'un CHAMP HORIZONTAL destine a etre
!      ecrit sur un fichier ARPEGE/ALADIN.
!       ( COdage de (Nouvelles ?) Donnees )
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KNUMER (Entree) ==> Numero de l'unite logique;
!                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
!                KNIVAU (Entree) ==> Niveau vertical eventuel;
!                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
!    ( Tableau ) PCHAMP (Entree) ==> Valeurs REELLES du champ a ecrire;
!                LDCOSP (Entree) ==> Vrai si le champ est represente
!                                    par des coefficients spectraux;
!                CDNOMA (Sortie) ==> Nom de l'article-champ a ecrire;
!                KLNOMA (Sortie) ==> Nombre de caracteres utiles dans
!                                    CDNOMA;
!    ( Tableau ) KVALCO (Sortie) ==> Donnees destinees a l'ecriture;
!                KLONGD (Sortie) ==> Nombre de valeurs (mots de 64 bits
!                                    en principe) a ecrire.
!                LDUNDF (Entree) ==> Vrai si ce champ a des valeurs 
!                                    indefinies
!                PUNDF  (Entree) ==> Dans le cas ou LDUNDF est vrai,
!                                    valeur non definie
!
!    Remarques:
!
!    - KVALCO doit avoir une longueur
!      suffisante pour stocker les donnees codees. Le dimensionnement
!      "tous terrains" est (2+ILCHAM), qui permet le cas echeant de
!      stocker un champ a pleine resolution sans codage effectif.
!      (ILCHAM est le nombre de valeurs du champ a ecrire)
!
!    - CDNOMA doit avoir au moins FA%JPXNOM caracteres.
!
!
TYPE(FA_COM)   :: FA
TYPE(FAGR1TAB) :: YDGR1TAB
INTEGER (KIND=JPLIKB) KREP, KNUMER, KNIVAU, KLNOMA, KLONGD
!
REAL (KIND=JPDBLR) PCHAMP (*), PUNDF
INTEGER (KIND=JPLIKB) KVALCO (*)
!
CHARACTER CDPREF*(*), CDSUFF*(*), CDNOMA*(*)
!
INTEGER (KIND=JPLIKB) IREP, ILPRFU, ILSUFU, ILNOMU
INTEGER (KIND=JPLIKB) IRANG, INIMES
INTEGER (KIND=JPLIKB) ILPREF, ILSUFF, ILCDNO, IRANGC
INTEGER (KIND=JPLIKB) IB1PAR (FA%JPLB1P), INGRIB
!
LOGICAL LLVERF, LLRLFI, LDCOSP, LLNOMU, LLNOPA, LDUNDF
!
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
IF (LHOOK) CALL DR_HOOK('FACON1_MT',0,ZHOOK_HANDLE)
LLVERF=.FALSE.
LLRLFI=.FALSE.
LLNOMU=.FALSE.
LLNOPA=.FALSE.
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
&               IB1PAR(6),ILPRFU,ILSUFU,ILNOMU)
IF (IREP.NE.0) GOTO 1001
LLNOMU=.TRUE.
KLNOMA=ILNOMU
!**
!     3.  -  FABRICATION DE L'ARTICLE A ECRIRE SUR LE FICHIER.
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
IF (LDCOSP) THEN
  IF (FA%FICHIER(IRANG)%NFGRIB.EQ.-1.OR.FA%FICHIER(IRANG)%NFGRIB.EQ.3) THEN
    FA%FICHIER(IRANG)%NRASVE=FA%FICHIER(IRANG)%NRASVE+1
    IF (FA%FICHIER(IRANG)%NRASVE.EQ.1.AND.FA%FICHIER(IRANG)%NRASHO.GT.0) THEN
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
      WRITE(FA%NULOUT,*)' FACON1 :  WARNING !!!!!           '
      WRITE(FA%NULOUT,*)' Un champ de coeff. spectraux avec'
      WRITE(FA%NULOUT,*)                               &
&      ' rangement type modele va etre ecrit alors que'
      WRITE(FA%NULOUT,*)                                &
&      ' d''autres champs spec. ont un rangt different.'
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
    ENDIF
  ELSEIF (FA%FICHIER(IRANG)%NFGRIB.GE.0.AND.FA%FICHIER(IRANG)%NFGRIB.LE.2) THEN
    FA%FICHIER(IRANG)%NRASHO=FA%FICHIER(IRANG)%NRASHO+1
    IF (FA%FICHIER(IRANG)%NRASHO.EQ.1.AND.FA%FICHIER(IRANG)%NRASVE.GT.0) THEN
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
      WRITE(FA%NULOUT,*)                    &
&      ' FACON1 :  WARNING !!!!!           '
      WRITE(FA%NULOUT,*)                   &
&      ' Un champ de coeff. spectraux avec'
      WRITE(FA%NULOUT,*)                               &
&      ' rangt autre que celui du modele va etre ecrit'
      WRITE(FA%NULOUT,*)                                &
&      ' alors que d''autres champs ont le rangt modele'
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
    ENDIF
  ENDIF
ENDIF
!
500 CONTINUE
!
IF (FA%FICHIER(IRANG)%NFGRIB.EQ.3) THEN
! Cas d'un champ qu'il faut "griber" avec GRIBEX
  CALL FACODX_MT64                                            &
&                 (FA,  IREP, IRANG, CDPREF, KNIVAU, CDSUFF,  &
&                  PCHAMP, LDCOSP, KVALCO, KLONGD,            &
&                  LDUNDF, PUNDF, YDGR1TAB)
!
! Cas particulier de l'erreur GRIBEX num 710: OUTPUT ARRAY TOO SMALL
! On s'en sert pour detecter un probleme de compactage lie a ce que
! le champ compacte + les descripteurs prennent plus de place que le
! champ non compacte...
! On sort donc du compactage (FACODX) pour demander un codage sans
! compactage (FACINE) avec rangement des valeurs selon le modele:
! FA%NFGRIB=-1.
!
  IF (IREP==-1710) THEN
    IREP = 0
    FA%FICHIER(IRANG)%NFGRIB = -1
    LLNOPA = .TRUE.
    GOTO 500
  ENDIF
ELSEIF (FALGRA (FA%FICHIER(IRANG)%NFGRIB)) THEN
! Cas d'un champ qu'il faut "griber" avec GRIB_API
  IF (LDCOSP .AND. (FALGRA_SP (FA%FICHIER(IRANG)%NFGRIB) == 102)) THEN
    INGRIB = FA%FICHIER(IRANG)%NFGRIB
    FA%FICHIER(IRANG)%NFGRIB = 2_JPLIKB
    CALL FACINE_MT64                                             &
&                   (FA,  IREP, IRANG, CDNOMA(1:ILNOMU), PCHAMP, &
&                    LDCOSP, KVALCO, KLONGD, IB1PAR,             &
&                    LDUNDF, PUNDF)
    FA%FICHIER(IRANG)%NFGRIB = INGRIB
  ELSE
    CALL FACGRA_MT64 (FA,  IREP, IRANG, CDPREF, KNIVAU, CDSUFF,  &
                    & PCHAMP(1), LDCOSP, KVALCO, KLONGD,         &
                    & LDUNDF, PUNDF, .TRUE.)
  ENDIF
ELSEIF (FA%FICHIER(IRANG)%NFGRIB.EQ.4) THEN
  CALL FACCPL_MT64                                            &
&                 (FA,  IREP, IRANG, CDPREF, KNIVAU, CDSUFF,  &
&                  PCHAMP, LDCOSP, KVALCO, KLONGD, IB1PAR)
ELSE
  CALL FACINE_MT64                                             &
&                 (FA,  IREP, IRANG, CDNOMA(1:ILNOMU), PCHAMP, &
&                  LDCOSP, KVALCO, KLONGD, IB1PAR,             &
&                  LDUNDF, PUNDF)
  IF (LLNOPA) FA%FICHIER(IRANG)%NFGRIB = 3
!  Le codage num 3 avait ete demande mais se revelait etre
!  plus gourmand en place que le num -1: on avait donc force
!  l'absence de compactage (-1). On revient maintenant au codage
!  num 3 pour ce cadre IRANG et les eventuels codages suivants.
!
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
  IF (LHOOK) CALL DR_HOOK('FACON1_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
CLNSPR='FACON1'
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
IF (LHOOK) CALL DR_HOOK('FACON1_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"
#include "facom2.ixnvms.h"
#include "falgra.h"

END SUBROUTINE FACON1_MT64

