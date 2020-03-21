      LOGICAL FUNCTION LE (X, Y)
      INTEGER (KIND=JPKIND) :: X (2), Y (2)

      LE = (X (1) <  Y (1)) .OR. &
     & ((X (1) == Y (1)) .AND. (X (2) <= Y (2)))

      END FUNCTION LE

      LOGICAL FUNCTION GE (X, Y)
      INTEGER (KIND=JPKIND) :: X (2), Y (2)

      GE = (X (1) >  Y (1)) .OR. &
     & ((X (1) == Y (1)) .AND. (X (2) >= Y (2)))

      END FUNCTION GE

      LOGICAL FUNCTION NE (X, Y)
      INTEGER (KIND=JPKIND) :: X (2), Y (2)

      NE = (X (1) /= Y (1)) .OR. &
         & (X (2) /= Y (2))

      END FUNCTION NE


