! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FAIEN1_MT64                                              &
&                     (FA,  KREP,   KNUMER, CDPREF, KNIVAU, CDSUFF, &
&                      PCHAMP, LDCOSP, LDUNDF, PUNDF, YDGR1TAB, LDREVERT)
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
!                LDUNDF (Entree) ==> Vrai si ce champ a des valeurs 
!                                    indefinies
!                PUNDF  (Entree) ==> Dans le cas ou LDUNDF est vrai,
!                                    valeur non definie
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
TYPE(FA_COM)   :: FA
TYPE(FAGR1TAB) :: YDGR1TAB
INTEGER (KIND=JPLIKB) KREP, KNUMER, KNIVAU
!
REAL (KIND=JPDBLR) PCHAMP (*)
!
CHARACTER CDPREF*(*), CDSUFF*(*)
!
INTEGER (KIND=JPLIKB) IREP, ILPRFU, ILSUFU, ILNOMU
INTEGER (KIND=JPLIKB) ILONGA, IRANG, INIMES
INTEGER (KIND=JPLIKB) ILPREF, ILSUFF
!
INTEGER (KIND=JPLIKB), ALLOCATABLE :: IVALCO(:)
INTEGER (KIND=JPLIKB) IB1PAR (FA%JPLB1P)
!
INTEGER (KIND=JPLIKB) IVALC1, IRANGC, ILCHAM, INGRIB, IPFAOS
!
LOGICAL LLVERF, LLRLFI, LDCOSP, LLNOMU, LLMLAM, LLNOPA, LDUNDF, LDREVERT
!
REAL (KIND=JPDBLR) :: PUNDF
!
CHARACTER CLPREF*(FA%JPXNOM), CLSUFF*(FA%JPXSUF)
!
CHARACTER(LEN=FA%JPXNOM) CLNOMA
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!**
!     1.  -  CONTROLES ET INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FAIEN1_MT',0,ZHOOK_HANDLE)
LLVERF=.FALSE.
LLRLFI=.FALSE.
LLNOMU=.FALSE.
LLNOPA=.FALSE.
ILPRFU=INT (LEN (CDPREF), JPLIKB)
ILSUFU=INT (LEN (CDSUFF), JPLIKB)
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
IF (FA%FICHIER(IRANG)%LCREAF) THEN
  IREP=-85
  GOTO 1001
ENDIF
!**
!     2.  -  FABRICATION DU NOM D'ARTICLE VIA LE SOUS-PROGRAMME "FANFAR"
!            ( controles de CDPREF, KNIVAU, CDSUFF inclus )
!-----------------------------------------------------------------------
!
CALL FANFAR_MT64                                           &
&               (FA, IREP,IRANG,CDPREF,KNIVAU,CDSUFF,CLNOMA, &
&             IB1PAR(6), ILPRFU,ILSUFU,ILNOMU)
IF (IREP.NE.0) GOTO 1001
LLNOMU=.TRUE.
!**
!     3.  -  CALCUL D'UN MAJORANT POUR LA LONGUEUR DE L'ARTICLE (mots)
!            ( on va prendre le nombre de valeurs du champ +2 :
!              l'absence de compactage est un majorant et les 2 mots
!              correspondent a l'enrobage FA dans ce cas )
!-----------------------------------------------------------------------
!
IVALC1=FA%FICHIER(IRANG)%NFGRIB
IRANGC=FA%FICHIER(IRANG)%NUCADR
LLMLAM=FA%CADRE(IRANGC)%LIMLAM
IF (LDCOSP) THEN
  IF (LLMLAM) THEN
    ILCHAM=FA%CADRE(IRANGC)%NSFLAM
  ELSE
    IF (IVALC1.EQ.-1 .OR. IVALC1.EQ.3) THEN
      ILCHAM=(1+FA%CADRE(IRANGC)%MTRONC)*(2+FA%CADRE(IRANGC)%MTRONC)
    ELSE
      ILCHAM=(1+FA%CADRE(IRANGC)%MTRONC)**2
    ENDIF
  ENDIF
ELSE
  ILCHAM=FA%CADRE(IRANGC)%NVAPDG
ENDIF
!

CALL FASGRA_MT64 (FA, IREP, FA%CADRE(IRANGC)%CNOMCA, IPFAOS)

IF (IREP.NE.0) GOTO 1001

ILONGA = ILCHAM+IPFAOS

ALLOCATE (IVALCO (ILONGA))
IVALCO = 0
!**
!     4.  -  FABRICATION DE L'ARTICLE A ECRIRE SUR LE FICHIER.
!-----------------------------------------------------------------------
!
!  Controle de l'homogeneite du type de rangement de coeff. spectraux
!  parmi les champs lus/ecrits: ces champs compactes avec
!  FA%NIGRIB=-1 ou 3 doivent etre ranges comme dans le modele ("verticalement"
!  soit selon des colonnes JM=cst consecutives) et contrairement si compactes
!  avec FA%NIGRIB= 0,1 ou 2.
!
IRANGC=FA%FICHIER(IRANG)%NUCADR
IF (LDCOSP) THEN
  IF (FA%FICHIER(IRANG)%NFGRIB.EQ.-1 .OR. FA%FICHIER(IRANG)%NFGRIB.EQ.3) THEN
    FA%FICHIER(IRANG)%NRASVE=FA%FICHIER(IRANG)%NRASVE+1
    IF (FA%FICHIER(IRANG)%NRASVE.EQ.1 .AND. FA%FICHIER(IRANG)%NRASHO.GT.0) THEN
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
      WRITE(FA%NULOUT,*)' FAIEN1 :  WARNING !!!!!           '
      WRITE(FA%NULOUT,*)' Un champ de coeff. spectraux avec'
      WRITE(FA%NULOUT,*)                               &
&      ' rangement type modele va etre ecrit alors que'
      WRITE(FA%NULOUT,*)                               &
&      ' les autres champs ont un rangement different.'
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
    ENDIF
  ELSEIF (FA%FICHIER(IRANG)%NFGRIB.GE.0 .AND. FA%FICHIER(IRANG)%NFGRIB.LE.2) THEN
    FA%FICHIER(IRANG)%NRASHO=FA%FICHIER(IRANG)%NRASHO+1
    IF (FA%FICHIER(IRANG)%NRASHO.EQ.1 .AND. FA%FICHIER(IRANG)%NRASVE.GT.0) THEN
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
      WRITE(FA%NULOUT,*)' FAIEN1 :  WARNING !!!!!           '
      WRITE(FA%NULOUT,*)' Un champ de coeff. spectraux avec'
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
&                  PCHAMP(1), LDCOSP, IVALCO, ILONGA,         &
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
&                   (FA,  IREP, IRANG, CLNOMA(1:ILNOMU), PCHAMP, &
&                    LDCOSP, IVALCO, ILONGA, IB1PAR,             &
&                    LDUNDF, PUNDF)
    FA%FICHIER(IRANG)%NFGRIB = INGRIB
  ELSE
    CALL FACGRA_MT64 (FA,  IREP, IRANG, CDPREF, KNIVAU, CDSUFF,  &
                    & PCHAMP(1), LDCOSP, IVALCO, ILONGA,         &
                    & LDUNDF, PUNDF, LDREVERT)
  ENDIF
ELSEIF (FA%FICHIER(IRANG)%NFGRIB.EQ.4) THEN
  CALL FACCPL_MT64                                            &
&                 (FA,  IREP, IRANG, CDPREF, KNIVAU, CDSUFF,  &
&                  PCHAMP(1), LDCOSP, IVALCO, ILONGA, IB1PAR)
ELSE
  CALL FACINE_MT64                                             &
&                 (FA,  IREP, IRANG, CLNOMA(1:ILNOMU), PCHAMP, &
&                  LDCOSP, IVALCO, ILONGA, IB1PAR,             &
&                  LDUNDF, PUNDF)
  IF (LLNOPA) FA%FICHIER(IRANG)%NFGRIB = 3
!  Le codage num 3 avait ete demande mais se revelait etre
!  plus gourmand en place que le num -1: on avait donc force
!  l'absence de compactage (-1). On revient maintenant au codage
!  num 3 pour ce cadre IRANG et les eventuels codages suivants.
!
ENDIF
IF (IREP.NE.0) GOTO 1001
!**
!     5.  -  ECRITURE DE L'ARTICLE "CHAMP" SUR LE FICHIER.
!-----------------------------------------------------------------------
!
CALL FAISAN_MT64 (FA, IREP, KNUMER, CLNOMA(1:ILNOMU), IVALCO, ILONGA)
LLRLFI=IREP.NE.0
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
!            VIA LE SOUS-PROGRAMME "FAIPAR" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
IF (ALLOCATED( IVALCO )) DEALLOCATE ( IVALCO )
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
  IF (LHOOK) CALL DR_HOOK('FAIEN1_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
CLNSPR='FAIEN1'
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
  CLNOMA(1:ILNOMU)=CLPREF(1:ILPREF)
ENDIF
!
WRITE (UNIT=CLMESS,FMT='(''KREP='',I5,'', KNUMER='',I3,        &
&       '', CDPREF='''''',A,'''''', KNIVAU='',I6,               &
&       '', CDSUFF='''''',A,'''''', LDCOSP= '',L1)')            &
&   KREP,KNUMER,CLPREF(1:ILPREF),KNIVAU,CLSUFF(1:ILSUFF),LDCOSP
CALL FAIPAR_MT64                                     &
&               (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS, &
&                CLNSPR, CLNOMA(1:ILNOMU),LLRLFI)
!
IF (LHOOK) CALL DR_HOOK('FAIEN1_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"
#include "facom2.ixnvms.h"
#include "falgra.h"

END SUBROUTINE FAIEN1_MT64
