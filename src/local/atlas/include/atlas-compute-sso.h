INTERFACE
SUBROUTINE ATLAS_COMPUTE_SSO (YDFSSC, PUNDEF, PZS, PZS2, PZS_DXDX, PZS_DYDY, PZS_DXDY, &
                            & PSSO_DIR, PSSO_SLOPE, PSSO_ANIS, PSSO_STDEV)
USE ATLAS_MODULE  
USE PARKIND1, ONLY : JPIM, JPRB
IMPLICIT NONE
TYPE (ATLAS_FUNCTIONSPACE_STRUCTUREDCOLUMNS), INTENT (IN) :: YDFSSC
REAL (KIND=JPRB), INTENT (IN)  :: PUNDEF
REAL (KIND=JPRB), INTENT (IN)  :: PZS (:), PZS2 (:), PZS_DXDX (:), PZS_DYDY (:), PZS_DXDY (:)
REAL (KIND=JPRB), INTENT (OUT) :: PSSO_DIR (:), PSSO_SLOPE (:), PSSO_ANIS (:), PSSO_STDEV (:)
END SUBROUTINE 
END INTERFACE