#define JPKIND JPIM

#define W 1
SUBROUTINE QSORTI4X1 (N,ORD,A)
#include "qsorti_x_.h"
#include "cmpi_x1.h"
END SUBROUTINE 
#undef W

#define W 2
SUBROUTINE QSORTI4X2 (N,ORD,A)
#include "qsorti_x_.h"
#include "cmpi_x2.h"
END SUBROUTINE 
#undef W

#define W 3
SUBROUTINE QSORTI4X3 (N,ORD,A)
#include "qsorti_x_.h"
#include "cmpi_x3.h"
END SUBROUTINE 
#undef W

#undef JPKIND

#define JPKIND JPIB

#define W 1
SUBROUTINE QSORTI8X1 (N,ORD,A)
#include "qsorti_x_.h"
#include "cmpi_x1.h"
END SUBROUTINE 
#undef W

#define W 2
SUBROUTINE QSORTI8X2 (N,ORD,A)
#include "qsorti_x_.h"
#include "cmpi_x2.h"
END SUBROUTINE 
#undef W

#define W 3
SUBROUTINE QSORTI8X3 (N,ORD,A)
#include "qsorti_x_.h"
#include "cmpi_x3.h"
END SUBROUTINE 
#undef W

#undef JPKIND
