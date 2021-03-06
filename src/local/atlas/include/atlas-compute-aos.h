INTERFACE
FUNCTION ATLAS_COMPUTE_AOS (YDFSSC1, YDFSSC2, YDOROG1, YDOROG2, YDINTEA, YDINTE4, LDOPENMP) RESULT (YLFLDSAOS2)

USE PARKIND1, ONLY : JPIM, JPRB
USE INTERPOLATIONA_MOD
USE INTERPOLATION4_MOD
USE ATLAS_MODULE  
IMPLICIT NONE

TYPE (ATLAS_FUNCTIONSPACE_STRUCTUREDCOLUMNS), INTENT (IN) :: YDFSSC1
TYPE (ATLAS_FUNCTIONSPACE_STRUCTUREDCOLUMNS), INTENT (IN) :: YDFSSC2
TYPE (ATLAS_FIELD),                           INTENT (IN) :: YDOROG1
TYPE (ATLAS_FIELD),                           INTENT (IN) :: YDOROG2
TYPE (INTERPOLATIONA),                        INTENT (IN) :: YDINTEA
TYPE (INTERPOLATION4),                        INTENT (IN) :: YDINTE4
LOGICAL,                                      INTENT (IN) :: LDOPENMP

TYPE (ATLAS_FIELDSET) :: YLFLDSAOS2
END FUNCTION
END INTERFACE

