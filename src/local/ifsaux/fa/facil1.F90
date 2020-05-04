! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FACIL1_MT64                                             &
&                     (FA,  KREP, KNUMER, CDPREF, KNIVAU, CDSUFF,  &
&                      PCHAMP, LDCOSP, LDUNDF, PUNDF, YDGR1TAB, LDREVERT)
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
!    ( Tableau ) PCHAMP (Sortie) ==> Valeurs REELLES du champ lu;
!                LDCOSP (Entree) ==> Vrai si le champ est represente
!                                    par des coefficients spectraux.
!                LDUNDF (Sortie) ==> Vrai si ce champ a des valeurs 
!                                    indefinies
!                PUNDF  (Sortie) ==> Dans le cas ou LDUNDF est vrai,
!                                    valeur non definie
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
TYPE(FAGR1TAB) :: YDGR1TAB
INTEGER (KIND=JPLIKB) KREP, KNUMER, KNIVAU
!
INTEGER (KIND=JPLIKB) IREP, ILPRFU, ILSUFU, ILNOMU
INTEGER (KIND=JPLIKB) ILONGA, IRANG, INIMES
INTEGER (KIND=JPLIKB) ILPREF, ILSUFF, IPOSEX, IRANGC
!
REAL (KIND=JPDBLR) PCHAMP (*)
REAL (KIND=JPRB) PUNDF
INTEGER (KIND=JPLIKB), ALLOCATABLE :: IVALCO(:)
INTEGER (KIND=JPLIKB) IB1PAR (FA%JPLB1P)
!
LOGICAL LLVERF, LLRLFI, LDCOSP, LLNOMU, LDUNDF, LDREVERT
!
CHARACTER CDPREF*(*), CDSUFF*(*)
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
IF (LHOOK) CALL DR_HOOK('FACIL1_MT',0,ZHOOK_HANDLE)

LLVERF=.FALSE.
LLRLFI=.FALSE.
LLNOMU=.FALSE.
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
CALL FANFAR_MT64                                             &
&               (FA, IREP,IRANG,CDPREF,KNIVAU,CDSUFF,CLNOMA, &
&                IB1PAR(6),ILPRFU,ILSUFU,ILNOMU)
IF (IREP.NE.0) GOTO 1001
LLNOMU=.TRUE.
!**
!     3.  -  LECTURE DE L'ARTICLE SUR LE FICHIER
!-----------------------------------------------------------------------
!
CALL LFINFO_MT64                                       &
&               (FA%LFI, IREP,KNUMER,CLNOMA(1:ILNOMU), &
&                ILONGA,IPOSEX)
!
IF (IREP.NE.0) THEN
  LLRLFI=.TRUE.
  GOTO 1001
ELSEIF (ILONGA.EQ.0) THEN
  IREP=-89
  GOTO 1001
ELSEIF (ILONGA.GT.FA%JPXCHA+2) THEN
  IREP=-90
  GOTO 1001
ENDIF
!
ALLOCATE (IVALCO (ILONGA))
CALL LFILEC_MT64                                       &
&               (FA%LFI, IREP,KNUMER,CLNOMA(1:ILNOMU), &
&               IVALCO,ILONGA)
LLRLFI=IREP.NE.0
IF (LLRLFI) GOTO 1001
!
!**
!     4.  -  CONTROLES ET DECODAGE DE L'ARTICLE
!----------------------------------------------
!
!  Controle de l'homogeneite du type de rangement de coeff. spectraux
!  parmi les champs lus/ecrits: ces champs compactes avec
!  FA%NIGRIB=-1 ou 3 doivent etre ranges comme dans le modele ("verticalement"
!  soit selon des colonnes JM=cst consecutives) et contrairement si compactes
!  avec FA%NIGRIB= 0,1 ou 2.
! 
IRANGC=FA%FICHIER(IRANG)%NUCADR
IF (LDCOSP) THEN
  IF (IVALCO(1).EQ.-1.OR.IVALCO(1).EQ.3) THEN
    FA%FICHIER(IRANG)%NRASVE=FA%FICHIER(IRANG)%NRASVE+1
    IF (FA%FICHIER(IRANG)%NRASVE.EQ.1.AND.FA%FICHIER(IRANG)%NRASHO.GT.0) THEN
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
      WRITE(FA%NULOUT,*)' FACIL1 :  WARNING !!!!!           '
      WRITE(FA%NULOUT,*)' Un champ de coeff. spectraux avec'
      WRITE(FA%NULOUT,*)                            &
&      ' rangement type modele va etre lu alors que'
      WRITE(FA%NULOUT,*)                                &
&      ' d''autres champs spec. ont un rangt different.'
      WRITE(FA%NULOUT,*)                                 &
&      ' ***  Prenez en compte cette heterogeneite!  ***'
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
    ENDIF
  ELSEIF (IVALCO(1).GE.0.AND.IVALCO(1).LE.2) THEN
    FA%FICHIER(IRANG)%NRASHO=FA%FICHIER(IRANG)%NRASHO+1
    IF (FA%FICHIER(IRANG)%NRASHO.EQ.1.AND.FA%FICHIER(IRANG)%NRASVE.GT.0) THEN
      WRITE(FA%NULOUT,*)                                 &
&      '------------------------------------------------'
      WRITE(FA%NULOUT,*)' FACIL1 :  WARNING !!!!!           '
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
IF (FALGRA (IVALCO(1))) THEN
! Cas d'un champ gribe avec GRIB_API
  CALL FADGRA_MT64                             &
&                 (FA, IREP,IRANG,CLNOMA,      &
&                  IVALCO,ILONGA,PCHAMP,LDCOSP,&
&                  CDPREF, KNIVAU, CDSUFF,     &
&                  LDUNDF, PUNDF, LDREVERT)
ELSEIF (IVALCO(1).EQ.3) THEN
! Cas d'un champ gribe avec GRIBEX
  CALL FADECX_MT64                             &
&                 (FA, IREP,IRANG,CLNOMA,      &
&                  IVALCO,ILONGA,PCHAMP,LDCOSP,&
&                  CDPREF, KNIVAU, CDSUFF,     &
&                  LDUNDF, PUNDF, YDGR1TAB)
ELSEIF (IVALCO(1).EQ.4) THEN
  CALL FADCPL_MT64                             &
&                 (FA, IREP,IRANG,CLNOMA,      &
&                  IVALCO,ILONGA,PCHAMP,LDCOSP,&
&                  LDUNDF, PUNDF)
ELSE
  CALL FADECI_MT64                             &
&                 (FA, IREP,IRANG,CLNOMA,      &
&                  IVALCO,ILONGA,PCHAMP,LDCOSP)
ENDIF
!
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

IF (LLFATA) THEN
  INIMES=2
ELSE
  INIMES=IXNVMS(IRANG)
ENDIF
!
IF (.NOT.LLFATA.AND.INIMES.NE.2)  THEN 
  IF (LHOOK) CALL DR_HOOK('FACIL1_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
CLNSPR='FACIL1'
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
IF (LHOOK) CALL DR_HOOK('FACIL1_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"
#include "facom2.ixnvms.h"
#include "falgra.h"

END SUBROUTINE FACIL1_MT64

