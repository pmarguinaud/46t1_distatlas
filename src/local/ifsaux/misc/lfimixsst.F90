SUBROUTINE LFIMIXSST

USE PARKIND1, ONLY : JPRB, JPIM
USE XRD_GETOPTIONS
USE LFI_PRECISION
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT

IMPLICIT NONE

#include "faieno.h"
#include "facilo.h"
#include "abor1.intfb.h"

CHARACTER (LEN=*), PARAMETER :: CLNOMC1 = 'c1', CLNOMC2 = 'c2'
INTEGER, PARAMETER :: ILUN1 = 77, ILUN2 = 88
CHARACTER (LEN=128) :: CLF1, CLF2
REAL (KIND=JPRB), ALLOCATABLE :: ZCHAMP1 (:), ZCHAMP2 (:)
INTEGER :: INBARP, INBARI, IREP, ILCHAM1, ILCHAM2
INTEGER :: INGRIB, INBPDG, INBCSP, ISTRON, IPUILA, IDMOPL, INBITS
LOGICAL :: LLEXIS, LLCOSP
REAL (KIND=JPRB) :: ZUNDEF
LOGICAL :: LLUNDEF


CALL INITOPTIONS (KOPTMIN=1)
CALL GETOPTION ('--fa-file-1', CLF1, MND = .TRUE., USE = "Take SST from this file")
CALL GETOPTION ('--fa-file-2', CLF2, MND = .TRUE., USE = "Update SST in this file")
CALL CHECKOPTIONS ()

INBARI=0
CALL FAITOU (IREP, ILUN1, .TRUE., CLF1, 'OLD', .TRUE., .FALSE., 0_JPIM, INBARP, INBARI, CLNOMC1)
INBARI=0
CALL FAITOU (IREP, ILUN2, .TRUE., CLF2, 'OLD', .TRUE., .FALSE., 0_JPIM, INBARP, INBARI, CLNOMC2)

CALL FATCHA (IREP, CLNOMC1, .FALSE., ILCHAM1)
CALL FATCHA (IREP, CLNOMC2, .FALSE., ILCHAM2)

IF (ILCHAM1 /= ILCHAM2) CALL ABOR1 ('SIZE MISMATCH')

ALLOCATE (ZCHAMP1 (ILCHAM1), ZCHAMP2 (ILCHAM2))

LLUNDEF = .TRUE.; ZUNDEF = HUGE (ZUNDEF)
CALL FACILO (IREP, ILUN1, 'SFX.', 0, 'SST', ZCHAMP1, .FALSE., LDUNDF=LLUNDEF, PUNDF=ZUNDEF)
LLUNDEF = .TRUE.; ZUNDEF = HUGE (ZUNDEF)
CALL FACILO (IREP, ILUN2, 'SFX.', 0, 'SST', ZCHAMP2, .FALSE., LDUNDF=LLUNDEF, PUNDF=ZUNDEF)

CALL FAVEUR (IREP, ILUN2, INGRIB, INBPDG, INBCSP, ISTRON, IPUILA, IDMOPL)
CALL FANION (IREP, ILUN2, 'SFX.', 0, 'SST', LLEXIS, LLCOSP, INGRIB, INBITS, ISTRON, IPUILA)
CALL FAGOTE (IREP, ILUN2, INGRIB, INBITS, INBCSP, ISTRON, IPUILA, IDMOPL)

WRITE (*, *) COUNT (ZCHAMP2 == ZUNDEF)
WRITE (*, *) COUNT (.NOT. ((ZCHAMP1 == ZUNDEF) .EQV. (ZCHAMP2 == ZUNDEF)))
WRITE (*, *) COUNT ((ZCHAMP1 /= ZUNDEF) .AND. (ZCHAMP2 == ZUNDEF))
WRITE (*, *) COUNT ((ZCHAMP2 /= ZUNDEF) .AND. (ZCHAMP1 == ZUNDEF))

WHERE (ZCHAMP1 (:) == ZUNDEF)
  ZCHAMP1 (:) = ZCHAMP2 (:)
ENDWHERE

WRITE (*, *) COUNT (ZCHAMP1 == ZUNDEF)

CALL FAIENO (IREP, ILUN2, 'SFX.', 0, 'SST', ZCHAMP1, .FALSE., LDUNDF=.TRUE., PUNDF=ZUNDEF)

CALL FAIRME (IREP, ILUN1, 'KEEP')
CALL FAIRME (IREP, ILUN2, 'KEEP')


END SUBROUTINE LFIMIXSST

