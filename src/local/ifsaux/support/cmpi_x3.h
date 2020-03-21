      LOGICAL FUNCTION LE (X, Y)
      INTEGER (KIND=JPKIND) :: X (3), Y (3)

      LE = (X (1) <  Y (1)) .OR. &
     & ((X (1) == Y (1)) .AND. ((X (2) <  Y (2)) .OR. &
     & ((X (2) == Y (2)) .AND.  (X (3) <= Y (3)))))

      END FUNCTION LE

      LOGICAL FUNCTION GE (X, Y)
      INTEGER (KIND=JPKIND) :: X (3), Y (3)

      GE = (X (1) >  Y (1)) .OR. &
     & ((X (1) == Y (1)) .AND. ((X (2) >  Y (2)) .OR. &
     & ((X (2) == Y (2)) .AND.  (X (3) >= Y (3)))))

      END FUNCTION GE

      LOGICAL FUNCTION NE (X, Y)
      INTEGER (KIND=JPKIND) :: X (3), Y (3)

      NE = (X (1) /= Y (1)) .OR. &
         & (X (2) /= Y (2)) .OR. &
         & (X (3) /= Y (3)) 

      END FUNCTION NE


