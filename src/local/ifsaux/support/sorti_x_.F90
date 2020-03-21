#define JPKIND JPIM

#define W 1
#define QSORT QSORTI4X1
SUBROUTINE SORTI4X1 (N, ORD, A)
#include "sorti_x_.h"
#include "cmpi_x1.h"
END SUBROUTINE 
#undef QSORT
#undef W

#define W 2
#define QSORT QSORTI4X2
SUBROUTINE SORTI4X2 (N, ORD, A)
#include "sorti_x_.h"
#include "cmpi_x2.h"
END SUBROUTINE 
#undef QSORT
#undef W

#define W 3
#define QSORT QSORTI4X3
SUBROUTINE SORTI4X3 (N, ORD, A)
#include "sorti_x_.h"
#include "cmpi_x3.h"
END SUBROUTINE
#undef QSORT
#undef W

#undef JPKIND

#define JPKIND JPIB

#define W 1
#define QSORT QSORTI8X1
SUBROUTINE SORTI8X1 (N, ORD, A)
#include "sorti_x_.h"
#include "cmpi_x1.h"
END SUBROUTINE 
#undef QSORT
#undef W

#define W 2
#define QSORT QSORTI8X2
SUBROUTINE SORTI8X2 (N, ORD, A)
#include "sorti_x_.h"
#include "cmpi_x2.h"
END SUBROUTINE 
#undef QSORT
#undef W

#define W 3
#define QSORT QSORTI8X3
SUBROUTINE SORTI8X3 (N, ORD, A)
#include "sorti_x_.h"
#include "cmpi_x3.h"
END SUBROUTINE
#undef QSORT
#undef W

#undef JPKIND
