! Feb-2013 P. Marguinaud Use JNGEOM & JNEXPL parameters
!                        Reallocate cadre when redefinition happens
! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FACADI_MT64                                           &
&                   (FA, KREP, CDNOMC, KTYPTR, PSLAPO, PCLOPO,   &
&                    PSLOPO,                                     &
&                    PCODIL, KTRONC, KNLATI, KNXLON, KNLOPA,     &
&                    KNOZPA, PSINLA, KNIVER, PREFER, PAHYBR,     &
&                    PBHYBR, LDMODC, LDREDF, KPHASE, KRANGC,     &
&                    KLNOMC, KGARDE)
USE FA_MOD, ONLY : FA_COM, NEW_CADRE, FREE_CADRE, JPNIIL, JNGEOM, JNEXPL
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!        Sous-programme A USAGE INTERNE AU LOGICIEL. Fait la plupart
!     des controles en vue de Definir un CADre, voire le redefinir.
!        En mode multi-taches, il doit y avoir verrouillage global
!     de la zone d'appel au sous-programme.
!**
!        Arguments : KREP   ==> Code-reponse du sous-programme;
!                    CDNOMC ==> Nom symbolique du cadre;
!  (tous d'Entree)   KTYPTR ==> Type de transformation horizontale;
!   sauf KRANGC      PSLAPO ==> Sinus de la latitude du pole d'interet;
!     et KLNOMC)     PCLOPO ==> Cosinus " " longitude "   "       "   ;
!                    PSLOPO ==> Sinus   " " longitude "   "       "   ;
!                    PCODIL ==> Coefficient de dilatation;
!                    KTRONC ==> Troncature;
!                    KNLATI ==> Nombre de latitudes (de pole a pole);
!                    KNXLON ==> Nombre maxi de longitudes par parallele;
!         (Tableau)  KNLOPA ==> Nombre de longitudes par parallele;
!                               (du pole nord vers l'equateur seulement)
!         (Tableau)  KNOZPA ==> Nombre d'onde zonal maxi par parallele;
!                               (du pole nord vers l'equateur seulement)
!         (Tableau)  PSINLA ==> Sinus des latitudes de l'hemisphere nord
!                               (du pole nord vers l'equateur seulement)
!                    KNIVER ==> Nombre de niveaux verticaux;
!                    PREFER ==> Pression de reference (facteur multipli-
!                               catif de la premiere fonction de la
!                               coordonnee hybride)
!         (Tableau)  PAHYBR ==> Valeurs de la fonction "A" de la coordo-
!                               nnee hybride AUX LIMITES DE COUCHES;
!         (Tableau)  PBHYBR ==> Valeurs de la fonction "B" de la coordo-
!                               nnee hybride AUX LIMITES DE COUCHES;
!                    LDMODC ==> Vrai s'il y a modification d'un cadre
!                               deja defini au prealable;
!                    LDREDF ==> Vrai s'il y a redefinition d'un cadre
!                               au sens large du terme (avec ou sans
!                               modification).
!                    KPHASE ==> Indique quelle(s) phase(s) du sous-prog.
!                               on doit executer:
!                               0 ==> Toutes,
!                               1 ==> Controle des variables simples,
!                               2 ==> Controle des tableaux,
!                               3 ==> Definition du cadre seule.
!        (Sortie)    KRANGC ==> Rang du cadre dans les tables.
! (Sortie si phase 1,KLNOMC ==> Longueur en caracteres du nom de cadre.
!  Entree sinon ! )
!                    KGARDE ==> Option de conservation du cadre
!                               apres la fermeture du dernier fichier
!                               qui s'y rattache. A noter que lors dans
!                               le cas d'une definition dynamique de
!                               cadre (appel par FAITOU, avec KGARDE=1),
!                               une redefinition de cadre n'est toleree
!                               qu'a l'identique.
!
!     N.B. :    En mode multi-taches, si l'on appelle le sous-programme
!            avec KPHASE=0 ou KPHASE=3, on doit verrouiller dans le
!            programme appelant l'appel au sous-programme.
!               Par ailleurs, LDMODC et LDREDF ne sont definis que si
!            KPHASE=0 ou KPHASE=3.
!*
!        La "redefinition" d'un cadre est possible a l'une de ces
!     conditions:
!
!     - le cadre a ete defini, mais n'a aucun fichier qui s'y rattache;
!     - le cadre defini a au moins un fichier qui s'y rattache, et les
!       nouveaux parametres de definition sont identiques a ceux deja
!       definis (a l'exception de l'option de conservation).
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KTYPTR, KTRONC, KNLATI
INTEGER (KIND=JPLIKB) KNXLON, KNIVER, KREP, KPHASE
INTEGER (KIND=JPLIKB) KRANGC, KLNOMC, KGARDE
!
INTEGER (KIND=JPLIKB) KNLOPA (FA%JPXPAH), KNOZPA (FA%JPXIND)
!
REAL (KIND=JPDBLR) PSLAPO, PCLOPO, PSLOPO, PCODIL, PREFER
!
REAL (KIND=JPDBLR) PSINLA (FA%JPXGEO), PAHYBR (0:KNIVER)
REAL (KIND=JPDBLR) PBHYBR (0:KNIVER)
REAL (KIND=JPDBLR),PARAMETER ::  ZEPS=1.E-15_JPDBLR
!
CHARACTER CDNOMC*(*)
!
LOGICAL LDREDF, LDMODC
!
INTEGER (KIND=JPLIKB) INPAHE
INTEGER (KIND=JPLIKB) ILCDNO, J, IPREC, ICOMPT, IMSMAX
INTEGER (KIND=JPLIKB) ISFLAM, JL, IK, INIMES, INUMER, ILNOMC
!
INTEGER (KIND=JPLIKB) IESN0 (0:FA%JPXTRO)
INTEGER (KIND=JPLIKB) IKNTMP(0:FA%JPXTRO)
INTEGER (KIND=JPLIKB) IKMTMP(0:FA%JPXTRO)
INTEGER (KIND=JPLIKB) ICPL4N(0:FA%JPXTRO)
!
REAL (KIND=JPDBLR) ZMIN, ZPMIN, ZPMAX, ZPMINP, ZPMAXP
!
LOGICAL LLMLAM
CHARACTER(LEN=FA%JPXNOM) CLACTI 
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!**
!     0.  -  AIGUILLAGE EN FONCTION DE *KPHASE*.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FACADI_MT',0,ZHOOK_HANDLE)

CLACTI=''
KREP=0
LDREDF=.FALSE.
LDMODC=.FALSE.
!
IF (KTYPTR .LE. 0 ) THEN
   LLMLAM = .TRUE.
ELSE
   LLMLAM = .FALSE.
ENDIF
!
INPAHE=(1+KNLATI)/2
!
IF (KPHASE.EQ.2) THEN
  GOTO 200
ELSEIF (KPHASE.EQ.3) THEN
  GOTO 300
ELSEIF (KPHASE.LT.0.OR.KPHASE.GT.3) THEN
  KREP=-66
  GOTO 1001
ENDIF
!**
!     1.  -  CONTROLE DES VARIABLES SIMPLES (SYNTAXE ET COHERENCE).
!            (sauf pression de reference)
!-----------------------------------------------------------------------
!
ILCDNO=INT (LEN (CDNOMC), JPLIKB)
KLNOMC=1
!
IF (ILCDNO.LE.0) THEN
  KREP=-65
  GOTO 1001
ELSEIF (CDNOMC.EQ.' ') THEN
  KREP=-68
  GOTO 1001
ELSEIF (KGARDE.LT.0.OR.KGARDE.GT.2) THEN
  KREP=-66
  GOTO 1001
ENDIF
!
DO J=ILCDNO,1,-1
!
IF (CDNOMC(J:J).NE.' ') THEN
  KLNOMC=J
  GOTO 102
ENDIF
!
ENDDO
!
102 CONTINUE
!
IF (KLNOMC.GT.FA%NCPCAD) THEN
  KREP=-65
  GOTO 1001
ENDIF
!
IF (KTRONC.LE.0.OR.KTRONC.GT.FA%NXTRON) THEN
  KREP=-70
  GOTO 1001
ELSEIF (KNLATI.LE.0.OR.KNLATI.GT.FA%NXLATI) THEN
  KREP=-71
  GOTO 1001
ELSEIF (KNIVER.LE.0.OR.KNIVER.GT.FA%NXNIVV) THEN
  KREP=-72
  GOTO 1001
ELSEIF (KNXLON.LE.0.OR.KNXLON.GT.FA%NXLONG) THEN
  KREP=-83
  GOTO 1001
ENDIF

IF (LLMLAM) THEN
!        IF (-2*KTYPTR+1.GT.KNXLON) THEN
!          KREP=-115
!          GOTO 1001
!        ELSEIF (2*KTRONC+1.GT.KNLATI) THEN
!          KREP=-116
!          GOTO 1001
!        ENDIF
ELSE
  IF (PCODIL.LT.1._JPDBLR) THEN
    KREP=-73
    GOTO 1001
  ELSEIF (KTYPTR.LE.0.OR.KTYPTR.GT.FA%NTYPTX) THEN
    KREP=-109
    GOTO 1001
  ELSEIF (MAX(ABS(PSLAPO),ABS(PCLOPO),ABS(PSLOPO)) &
&          .GT.1._JPDBLR) THEN
    KREP=-100
    GOTO 1001
  ELSEIF (ABS (1._JPDBLR-(PCLOPO**2+PSLOPO**2)) &
&          .GT.1.E-5_JPDBLR) THEN
    KREP=-101
    GOTO 1001
  ELSEIF (2*KTRONC+1.GT.KNXLON) THEN
    KREP=-84
    GOTO 1001
  ELSEIF (2*KTRONC+1.GT.4*(KNLATI/2)) THEN
!
!       Le test ci-dessus est "dur" car il fait l'hypothese que,
!     dans le cas ou KNLATI est impair, la grille comporte les poles.
!
    KREP=-79
    GOTO 1001
  ENDIF
ENDIF
!
IF (KPHASE.EQ.1) THEN
  GOTO 1001
ENDIF
!**
!     2.  -  CONTROLE DES TABLEAUX (SYNTAXE ET COHERENCE).
!            (et de la pression de reference)
!-----------------------------------------------------------------------
!
200 CONTINUE
!
!
IF (PREFER.LT.0._JPDBLR.OR.                     &
&    PREFER.GT.REAL (10*FA%MPRESX, JPDBLR)) THEN
  KREP=-108
  GOTO 1001
ENDIF
!
!     No Mount Everest test
!
IF (.FALSE.) THEN
DO J=0,KNIVER
IPREC=MAX (0_JPLIKB ,J-1)
ZMIN=MIN (PAHYBR(J),PBHYBR(J))
ZPMIN=PREFER*PAHYBR(J)+FA%SPSMIN*PBHYBR(J)
ZPMAX=PREFER*PAHYBR(J)+FA%SPSMAX*PBHYBR(J)
ZPMINP=PREFER*PAHYBR(IPREC)+FA%SPSMIN*PBHYBR(IPREC)
ZPMAXP=PREFER*PAHYBR(IPREC)+FA%SPSMAX*PBHYBR(IPREC)
!
IF (ZMIN.LT.0._JPDBLR.OR.PBHYBR(J).GT.1._JPDBLR) THEN
  KREP=-80
  GOTO 1001
ELSEIF (J.NE.0.AND.(PBHYBR(J).LT.PBHYBR(IPREC).OR.         &
&                 ZPMIN.LE.ZPMINP.OR.ZPMAX.LE.ZPMAXP)) THEN
  KREP=-81
  GOTO 1001
ENDIF
!
ENDDO
ENDIF ! No Mount Everest test
!
IF (.NOT.LLMLAM) THEN
!
   DO J=1,INPAHE
   IPREC=MAX (1_JPLIKB ,J-1)
!
   IF (KNLOPA(J).LE.0.OR.KNLOPA(J).GT.KNXLON) THEN
     KREP=-74
     GOTO 1001
   ELSEIF (KNLOPA(J).LT.KNLOPA(IPREC)) THEN
     KREP=-75
     GOTO 1001
   ELSEIF (KNOZPA(J).LT.0.OR.KNOZPA(J).GT.KTRONC) THEN
     KREP=-76
     GOTO 1001
   ELSEIF (KNOZPA(J).LT.KNOZPA(IPREC)) THEN
     KREP=-77
     GOTO 1001
   ELSEIF ((2*KNOZPA(J)+1).GT.KNLOPA(J)) THEN
     KREP=-78
     GOTO 1001
   ELSEIF (ABS (PSINLA(J)).GT.1._JPDBLR) THEN
     KREP=-102
     GOTO 1001
   ELSEIF (PSINLA(J).GE.PSINLA(IPREC).AND.J.NE.1) THEN
     KREP=-103
     GOTO 1001
   ENDIF
!
   ENDDO
!
ELSE
!
!        *****  ERROR HANDLING FOR LAM CASE
!
   IF (ABS(KNLOPA(2)).GT.1) THEN
     KREP=-117
     GOTO 1001
   ELSEIF (KNLOPA(3).LE.0.OR.KNLOPA(3).GT.KNXLON) THEN
     KREP=-118
     GOTO 1001
   ELSEIF (KNLOPA(4).LT.KNLOPA(3).OR.KNLOPA(4).GT.KNXLON) THEN
     KREP=-119
     GOTO 1001
   ELSEIF (KNLOPA(5).LE.0.OR.KNLOPA(5).GT.KNLATI) THEN
     KREP=-120
     GOTO 1001
   ELSEIF (KNLOPA(6).LE.KNLOPA(5).OR.KNLOPA(6).GT.KNLATI) THEN
     KREP=-121
     GOTO 1001
   ELSEIF (2*KNLOPA(7).GT.(KNLOPA(4)-KNLOPA(3))) THEN
     KREP=-122
     GOTO 1001
   ELSEIF (2*KNLOPA(8).GT.(KNLOPA(6)-KNLOPA(5))) THEN
     KREP=-123
     GOTO 1001
   ENDIF
!
ENDIF
!
IF (KPHASE.EQ.2) GOTO 1001
!**
!     3.  -  CONTROLES LIES A LA DEFINITION DU CADRE PROPREMENT DITE.
!-----------------------------------------------------------------------
!
300 CONTINUE
!
!        Le nom de cadre specifie est-il deja defini ?
!
CALL FANUCA_MT64                          &
&               (FA, CDNOMC,KRANGC,.FALSE.)
LDREDF=KRANGC.NE.0
IF (LDREDF) GOTO 500
!
!        En arrivant ici, il s'agit donc d'un nouveau cadre.
!
IF (FA%NCADEF.GE.FA%JPNXCA) THEN
!
!        Trop de cadres deja definis pour en stocker un de plus.
!
  KREP=-56
  GOTO 1001
ENDIF
!
!       Recherche d'un emplacement disponible dans les tables de cadres,
!     lequel devrait en bonne logique exister...
!
DO J=1,FA%JPNXCA
!
IF (FA%CADRE(J)%CNOMCA.EQ.' ') THEN
  KRANGC=J
  GOTO 303
ENDIF
!
ENDDO
!
KREP=-66
GOTO 1001
!
303 CONTINUE
!
!           Nouveau cadre, mise a jour des tables partagees de cadres.
!
FA%NCADEF=FA%NCADEF+1
FA%NCAIND(FA%NCADEF)=KRANGC

400 CONTINUE

CALL NEW_CADRE (FA%CADRE(KRANGC), KTYPTR, KNLATI, KTRONC, KNIVER)

FA%CADRE(KRANGC)%CNOMCA=CDNOMC
FA%CADRE(KRANGC)%NLCCAD=KLNOMC
!**
!     4.  -  STOCKAGE DES PARAMETRES DU CADRE (NOUVEAU, OU REDEFINI).
!-----------------------------------------------------------------------
!
FA%CADRE(KRANGC)%NULCAD=0
FA%CADRE(KRANGC)%NTYPTR=KTYPTR
FA%CADRE(KRANGC)%MTRONC=KTRONC
FA%CADRE(KRANGC)%NNIVER=KNIVER
FA%CADRE(KRANGC)%NLATIT=KNLATI
FA%CADRE(KRANGC)%NXLOPA=KNXLON
FA%CADRE(KRANGC)%SSLAPO=PSLAPO
FA%CADRE(KRANGC)%SCLOPO=PCLOPO
FA%CADRE(KRANGC)%SSLOPO=PSLOPO
FA%CADRE(KRANGC)%SCODIL=PCODIL
FA%CADRE(KRANGC)%SPREFE=PREFER
!
FA%CADRE(KRANGC)%LIMLAM=LLMLAM
FA%CADRE(KRANGC)%NSFLAM=0
!
IF (.NOT.LDREDF.OR.KGARDE.NE.1) FA%CADRE(KRANGC)%NGARDE=KGARDE
!
IF (.NOT.LLMLAM) THEN
   ICOMPT=0
!
   DO J=1,INPAHE
   ICOMPT=ICOMPT+KNLOPA(J)
   FA%CADRE(KRANGC)%NLOPAR(J)=KNLOPA(J)
   FA%CADRE(KRANGC)%NOZPAR(J)=KNOZPA(J)
   FA%CADRE(KRANGC)%SINLAT(J)=PSINLA(J)
   ENDDO
!
   IF (KNLATI.EQ.2*INPAHE) THEN
     FA%CADRE(KRANGC)%NVAPDG=ICOMPT*2
   ELSE
     FA%CADRE(KRANGC)%NVAPDG=ICOMPT*2-KNLOPA(INPAHE)
   ENDIF
!
ELSE
! *****  CALCULATION OF KNOZPA(), THEN ALSO SETTING OF FACOM1-TABLES  *****
!
   IMSMAX = -KTYPTR
   ISFLAM = 0
   CALL ELLIPS64  (KTRONC,IMSMAX,IKNTMP,IKMTMP)
!DP      CALL ELLIPS(IMSMAX,KTRONC,IKNTMP,IKMTMP)
!
! Initialisation de FA%NOMPAR (du module FAMODU)
!
   FA%CADRE(KRANGC)%NOMPAR(2) = 0
   DO JL=0,IMSMAX
     FA%CADRE(KRANGC)%NOMPAR(2*JL+3) = FA%CADRE(KRANGC)%NOMPAR(2*JL+2) + 1
     FA%CADRE(KRANGC)%NOMPAR(2*JL+4) = FA%CADRE(KRANGC)%NOMPAR(2*JL+3) &
&                             + 4*(IKNTMP(JL)+1) -1
   ENDDO
   FA%CADRE(KRANGC)%NOMPAR(1) = KTRONC
   FA%CADRE(KRANGC)%NOMPAR(2) = IMSMAX
!
   DO JL=0,KTRONC
      IK=IKMTMP(JL)
!DP         IK=IKNTMP(JL)
      ICPL4N(JL)=4*(IK+1)
      ISFLAM = ISFLAM + 4*(IK+1)
   ENDDO
!
   IESN0(0)=1
!
   DO J=1,KTRONC
      IESN0(J)=IESN0(J-1)+ICPL4N(J-1)
   ENDDO
!
! -----  NOW SETTING OF TABLES  -----
   DO J=1,JNEXPL
      FA%CADRE(KRANGC)%NLOPAR(J)=KNLOPA(J)
   ENDDO
   FA%CADRE(KRANGC)%SINLAT = 0._JPDBLR
   DO J=1,JNGEOM
      FA%CADRE(KRANGC)%SINLAT(J)=PSINLA(J)
   ENDDO
   FA%CADRE(KRANGC)%NOZPAR(1)=KTRONC
   FA%CADRE(KRANGC)%NOZPAR(2)=IMSMAX
!
   DO J=0,KTRONC
      FA%CADRE(KRANGC)%NOZPAR(2*J+3)=IESN0(J)
      FA%CADRE(KRANGC)%NOZPAR(2*J+4)=IESN0(J)+ICPL4N(J)-1
   ENDDO
  
   IF (FA%CADRE(KRANGC)%NOZPAR(2*KTRONC+4).NE. &
&       FA%CADRE(KRANGC)%NOMPAR(2*IMSMAX+4))    &
&   THEN
     KREP=-127
     GOTO 1001
   ENDIF
!
   FA%CADRE(KRANGC)%NSFLAM=ISFLAM
!
! *****  DETERMINATION OF FA%NVAPDG()  *****
!
   FA%CADRE(KRANGC)%NVAPDG=KNLATI*KNXLON
!
ENDIF
!
DO J=0,KNIVER
FA%CADRE(KRANGC)%SFOHYB(1,J)=PAHYBR(J)
FA%CADRE(KRANGC)%SFOHYB(2,J)=PBHYBR(J)
ENDDO
!
GOTO 1001
!**
!     5.  -  TENTATIVE DE REDEFINITION D'UN CADRE. CONTROLES AD HOC.
!-----------------------------------------------------------------------
!
500 CONTINUE
!
IF (FA%CADRE(KRANGC)%MTRONC.NE.KTRONC)        GOTO 505
IF (FA%CADRE(KRANGC)%NNIVER.NE.KNIVER)        GOTO 505
IF (FA%CADRE(KRANGC)%NLATIT.NE.KNLATI)        GOTO 505
IF (FA%CADRE(KRANGC)%NXLOPA.NE.KNXLON)        GOTO 505
IF (FA%CADRE(KRANGC)%NTYPTR.NE.KTYPTR)        GOTO 505
IF (ABS(REAL (FA%CADRE(KRANGC)%SSLAPO, JPDBLR)-REAL (PSLAPO, JPDBLR))>ZEPS) GOTO 505
IF (ABS(REAL (FA%CADRE(KRANGC)%SCLOPO, JPDBLR)-REAL (PCLOPO, JPDBLR))>ZEPS) GOTO 505
IF (ABS(REAL (FA%CADRE(KRANGC)%SSLOPO, JPDBLR)-REAL (PSLOPO, JPDBLR))>ZEPS) GOTO 505
IF (ABS(REAL (FA%CADRE(KRANGC)%SCODIL, JPDBLR)-REAL (PCODIL, JPDBLR))>ZEPS) GOTO 505
IF (ABS(REAL (FA%CADRE(KRANGC)%SPREFE, JPDBLR)-REAL (PREFER, JPDBLR))>ZEPS) GOTO 505
!
IF (.NOT.LLMLAM) THEN
   DO J=1,INPAHE
     IF (FA%CADRE(KRANGC)%NLOPAR(J).NE.KNLOPA(J))        GOTO 505
     IF (FA%CADRE(KRANGC)%NOZPAR(J).NE.KNOZPA(J))        GOTO 505
     IF (ABS(REAL (FA%CADRE(KRANGC)%SINLAT(J), JPDBLR) - REAL (PSINLA(J), JPDBLR))>ZEPS) GOTO 505
   ENDDO
ELSE
   DO J=1,JNEXPL
   IF (FA%CADRE(KRANGC)%NLOPAR(J).NE.KNLOPA(J)) GOTO 505
   ENDDO
   DO J=1,JNGEOM
   IF (ABS(REAL (FA%CADRE(KRANGC)%SINLAT(J), JPDBLR)-REAL (PSINLA(J), JPDBLR))>ZEPS) GOTO 505
   ENDDO
ENDIF
!
DO J=0,KNIVER
IF (ABS(REAL (FA%CADRE(KRANGC)%SFOHYB(1,J), JPDBLR)-REAL (PAHYBR(J), JPDBLR))>ZEPS) GOTO 505
IF (ABS(REAL (FA%CADRE(KRANGC)%SFOHYB(2,J), JPDBLR)-REAL (PBHYBR(J), JPDBLR))>ZEPS) GOTO 505
ENDDO
!
!        Si on arrive ici, il y a redefinition a l'identique,
!     du moins pour les parametres numeriques.
!        L'option de conservation du cadre peut, elle, etre modifiee
!     dans le cas d'une definition non dynamique.
!
IF (KGARDE.NE.1) FA%CADRE(KRANGC)%NGARDE=KGARDE
GOTO 1001
!
505 CONTINUE
LDMODC=.TRUE.
!
!        Il y a donc redefinition avec changement de parametre(s),
!     ce qui n'est possible que s'il n'y a pas de fichier rattache,
!     et s'il ne s'agit pas d'une definition dynamique de cadre
!     (appel par FAITOU avec KGARDE=1).
!
IF (KGARDE.EQ.1) THEN
  KREP=-58
ELSEIF (FA%CADRE(KRANGC)%NULCAD.NE.0) THEN
  KREP=-59
ELSE
  CALL FREE_CADRE (FA%CADRE(KRANGC))
  GOTO 400
ENDIF
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE EVENTUELLE,
!            VIA LE sous-programme "FAIPAR" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
!
LLFATA=KREP.NE.0.AND.FA%NRFAGA.NE.2
!
IF (FA%LFAMOP.OR.LLFATA) THEN
  INIMES=2
  CLNSPR='FACADI'
  INUMER=JPNIIL
!
  IF (KREP.EQ.-65.AND.ILCDNO.LE.0) THEN
    ILNOMC=8
    CLACTI(1:ILNOMC)=FA%CHAINC(:ILNOMC)
  ELSE
    ILNOMC=MIN (KLNOMC,FA%NCPCAD,INT (LEN (CLACTI), JPLIKB))
    CLACTI(1:ILNOMC)=CDNOMC(1:ILNOMC)
  ENDIF
!
  WRITE (UNIT=CLMESS,FMT='(''ARGUM.SIMPLES='',I4,'','''''',A, &
&         '''''''',4('','',F7.4),4('','',I4),'','',F10.3,      &
&         2('','',L1),2('','',I2),'','',I3,'','',I1)')         &
&  KREP,CLACTI(1:ILNOMC),PSLAPO,PCLOPO,PSLOPO,PCODIL,          &
&  KTRONC,KNLATI,KNXLON,KNIVER,PREFER,LDMODC,LDREDF,KPHASE,    &
&  KRANGC,KLNOMC,KGARDE
  CALL FAIPAR_MT64                                      &
&                 (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                  CLNSPR,CLACTI(1:ILNOMC),.FALSE.)
ELSEIF (KTRONC.LE.FA%NSTROI.AND.(KPHASE.EQ.0.OR.KPHASE.EQ.1)) THEN
  INIMES=1
  CLNSPR='FACADI'
  INUMER=JPNIIL
  ILNOMC=MIN (KLNOMC,FA%NCPCAD)
  WRITE (UNIT=CLMESS,                                              &
&         FMT='(''TRONCATURE ('',I2,'') INFERIEURE '',              &
& ''OU EGALE A LA SOUS-TRONCATURE "NON COMPACTEE" IMPLICITE ('',I2, &
& ''), CADRE '''''',A,'''''''')') KTRONC,FA%NSTROI,CDNOMC(1:ILNOMC)
  CALL FAIPAR_MT64                                      &
&                 (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                  CLNSPR,CLACTI,.FALSE.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('FACADI_MT',1,ZHOOK_HANDLE)
END SUBROUTINE FACADI_MT64

