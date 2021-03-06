MODULE ATLAS_FMT_SOIL_AND_VEG

USE ATLAS_FMT, ONLY : ATLAS_FMT_t

USE PARKIND1, ONLY : JPRB, JPIM
USE ATLAS_MODULE

#include "atlas-abort.h"

IMPLICIT NONE

TYPE, EXTENDS (ATLAS_FMT_t) :: ATLAS_FMT_SOIL_AND_VEG_t

CONTAINS

  PROCEDURE :: OPEN
  PROCEDURE :: CLOSE
  PROCEDURE :: GRID
  PROCEDURE :: READ    
  PROCEDURE :: WRITE   

END TYPE

PRIVATE

PUBLIC :: ATLAS_FMT_SOIL_AND_VEG_t

CONTAINS

TYPE (ATLAS_STRUCTUREDGRID) FUNCTION GRID (THIS, CDFILE) RESULT (YLGRID)
CLASS (ATLAS_FMT_SOIL_AND_VEG_t), INTENT (IN) :: THIS
CHARACTER (LEN=*),                INTENT (IN) :: CDFILE

INTEGER (KIND=JPIM) :: ILEN, ILUN, ILAT, ILON

ILUN = 77_JPIM

OPEN (UNIT=ILUN, FILE=TRIM (CDFILE), FORM='UNFORMATTED', STATUS='OLD', ACCESS='DIRECT', RECL=4)

READ (ILUN, REC=1) ILEN

CLOSE (ILUN)

IF (MODULO (ILEN, 8)) THEN
  CALL ATLAS_ABORT ('UNEXPECTED GEOMETRY')
ENDIF

ILEN = ILEN / 8

ILON = NINT (SQRT (REAL (2 * ILEN, JPRB)), JPIM)

IF (ILON * ILON /= 2 * ILEN) THEN
  CALL ATLAS_ABORT ('UNEXPECTED GEOMETRY')
ENDIF

ILAT = ILON / 2

BLOCK
  TYPE (ATLAS_CONFIG) :: YLCFGR, YLCFDO
  REAL (KIND=JPRB) :: ZLONW, ZLATS, ZLONE, ZLATN, ZLOND, ZLATD

  ZLOND = 360._JPRB / REAL (ILON, JPRB)
  ZLATD = 180._JPRB / REAL (ILAT, JPRB)
  ZLONW =   0._JPRB + ZLOND / 2._JPRB
  ZLONE = 360._JPRB - ZLOND / 2._JPRB
  ZLATS = -90._JPRB + ZLATD / 2._JPRB
  ZLATN = +90._JPRB - ZLATD / 2._JPRB

  YLCFGR = ATLAS_CONFIG ()
  YLCFDO = ATLAS_CONFIG ()

  CALL YLCFGR%SET ("nx", ILON)
  CALL YLCFGR%SET ("ny", ILAT)
  CALL YLCFGR%SET ("type", "shifted_lonlat")
  CALL YLCFDO%SET ("type", "global") 
  CALL YLCFDO%SET ("xmin", ZLONW)
  CALL YLCFDO%SET ("xmax", ZLONE+ZLOND)
  CALL YLCFDO%SET ("west", ZLOND / 2._JPRB)

  CALL YLCFDO%SET ("ymin", ZLATS)
  CALL YLCFDO%SET ("ymax", ZLATN)
  CALL YLCFDO%SET ("units", "degrees")
  CALL YLCFGR%SET ("domain", YLCFDO)

  YLGRID = ATLAS_STRUCTUREDGRID (YLCFGR)

  CALL YLCFGR%FINAL ()
  CALL YLCFDO%FINAL ()
ENDBLOCK

CALL YLGRID%RETURN ()

END FUNCTION

SUBROUTINE OPEN (THIS, CDFILE, YDGRID)
CLASS (ATLAS_FMT_SOIL_AND_VEG_t),   INTENT (INOUT) :: THIS
CHARACTER (LEN=*),                  INTENT (IN)    :: CDFILE
CLASS (ATLAS_STRUCTUREDGRID),       INTENT (IN)    :: YDGRID
END SUBROUTINE OPEN

SUBROUTINE CLOSE (THIS)
CLASS (ATLAS_FMT_SOIL_AND_VEG_t),  INTENT (INOUT) :: THIS
END SUBROUTINE


SUBROUTINE READ (THIS, CDNAME, PFLDG, LDUNDEF, PUNDEF)
CLASS (ATLAS_FMT_SOIL_AND_VEG_t),  INTENT (IN)    :: THIS
CHARACTER (LEN=*),                 INTENT (IN)    :: CDNAME
REAL (KIND=JPRB),                  INTENT (OUT)   :: PFLDG (:)
LOGICAL,                           INTENT (INOUT) :: LDUNDEF
REAL (KIND=JPRB),                  INTENT (INOUT) :: PUNDEF

INTEGER (KIND=JPIM) :: ILUN

ILUN = 77_JPIM

OPEN (ILUN, FILE=TRIM (CDNAME), FORM='UNFORMATTED', STATUS='OLD')
READ (ILUN) PFLDG
CLOSE (ILUN)

WHERE (PFLDG == -9999._JPRB)
  PFLDG = PUNDEF
ENDWHERE

END SUBROUTINE

SUBROUTINE WRITE (THIS, CDNAME, PFLDG, LDUNDEF, PUNDEF)
CLASS (ATLAS_FMT_SOIL_AND_VEG_t), INTENT (IN)    :: THIS
CHARACTER (LEN=*),                INTENT (IN)    :: CDNAME
REAL (KIND=JPRB),                 INTENT (IN)    :: PFLDG (:)
LOGICAL,                          INTENT (INOUT) :: LDUNDEF
REAL (KIND=JPRB),                 INTENT (INOUT) :: PUNDEF

INTEGER (KIND=JPIM) :: ILUN

ILUN = 77_JPIM

OPEN (ILUN, FILE=TRIM (CDNAME), FORM='UNFORMATTED', STATUS='NEW')
WRITE (ILUN) PFLDG
CLOSE (ILUN)

END SUBROUTINE

END MODULE

