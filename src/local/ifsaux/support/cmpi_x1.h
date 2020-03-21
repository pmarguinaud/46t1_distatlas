      LOGICAL FUNCTION LE (X, Y)
      INTEGER (KIND=JPKIND) :: X (1), Y (1)

      LE = (X (1) <= Y (1))

      END FUNCTION LE

      LOGICAL FUNCTION GE (X, Y)
      INTEGER (KIND=JPKIND) :: X (1), Y (1)

      GE = (X (1) >= Y (1))

      END FUNCTION GE

      LOGICAL FUNCTION NE (X, Y)
      INTEGER (KIND=JPKIND) :: X (1), Y (1)

      NE = (X (1) /= Y (1))

      END FUNCTION NE


