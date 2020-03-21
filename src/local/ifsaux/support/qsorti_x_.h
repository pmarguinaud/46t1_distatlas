      
      USE PARKIND1

      IMPLICIT NONE

      INTEGER(KIND=JPKIND) :: N
      INTEGER(KIND=JPKIND) :: ORD(N)
      INTEGER(KIND=JPKIND) :: I,IP,IQ,IX,IZ,L,L1,NDEEP,P
      INTEGER(KIND=JPKIND) :: POPLST(2,20)
      INTEGER(KIND=JPKIND) :: Q,U,U1,YP

      INTEGER(KIND=JPKIND) :: A(W,N) 
      INTEGER(KIND=JPKIND) :: X (W), XX (W), Z (W), ZZ (W), Y (W)
      NDEEP=0
      U1=N
      L1=1
      DO 1  I=1,N
    1 ORD(I)=I
    2 IF (U1.LE.L1) RETURN

    3 L=L1
      U=U1
    4 P=L
      Q=U

      X=A(1:W,ORD(P))
      Z=A(1:W,ORD(Q))
      IF (LE (X,Z)) GO TO 5
      Y=X
      X=Z
      Z=Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
    5 IF (U-L.LE.1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q

    6 P=P+1
      IF (P.GE.Q) GO TO 7
      X=A(1:W,ORD(P))
      IF (GE (X,XX)) GO TO 8
      GO TO 6
    7 P=Q-1
      GO TO 13

    8 Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=A(1:W,ORD(Q))
      IF (LE (Z,ZZ)) GO TO 10
      GO TO 8
    9 Q=P
      P=P-1
      Z=X
      X=A(1:W,ORD(P))

   10 IF (LE (X,Z)) GO TO 11
      Y=X
      X=Z
      Z=Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
   11 IF (LE (X,XX)) GO TO 12
      XX=X
      IX=P
   12 IF (GE (Z,ZZ)) GO TO 6
      ZZ=Z
      IZ=Q
      GO TO 6

   13 CONTINUE
      IF (.NOT.(P.NE.IX.AND.NE (X,XX))) GO TO 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
   14 CONTINUE
      IF (.NOT.(Q.NE.IZ.AND.NE (Z,ZZ))) GO TO 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
   15 CONTINUE
      IF (U-Q.LE.P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
   16 U1=U
      L1=Q+1
      U=P-1
   17 CONTINUE
      IF (U1.LE.L1) GO TO 18

      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
   18 IF (U.GT.L) GO TO 4

      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18

      CONTAINS


