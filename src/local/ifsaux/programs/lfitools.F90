PROGRAM LFITOOLS

! PURPOSE.
! --------
!  Driver to various programs using FA or LFI software

! INTERFACE.
! ----------
!  Usage : $(LFITOOLS) 'name' [arg-1] [arg-1] ... [arg-N]
!          where arg-x are the arguments of the program 'name'
!          When invoked without argument, it returns the
!          usage and the available programs.
! METHOD.
! -------
!  I don't know yet 

! AUTHOR.
! -------
!  RYAD EL KHATIB *METEO-FRANCE*
!  ORIGINAL : 29-Mar-2012

!  MODIFICATIONS.
!  P.Marguinaud : 11-09-2012 : More programs
!  ------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IARGUMENTS=0 ! Number of arguments in command line
INTEGER(KIND=JPIM) :: IARGC
CHARACTER(LEN=64)  :: CLNAME=' '  ! argument program name
CHARACTER(LEN=256) :: CLMAIN=' '  ! this program name


IARGUMENTS = IARGC()

IF (IARGUMENTS == 0) THEN

  CALL GETARG(0,CLMAIN)
  CALL USAGE

ELSE

  CALL GETARG(1,CLNAME)
  CALL LC (CLNAME)

  SELECT CASE (TRIM(CLNAME))
    CASE ('--help','help')
      CALL USAGE
    CASE ('lficfm')
      CALL LFICFM
    CASE ('faop')
      CALL FAOP
    CASE ('fainterp')
!     CALL FAINTERP
    CASE ('facplspec')
!     CALL FACPLSPEC
    CASE ('fahis2cpl')
!     CALL FAHIS2CPL
    CASE ('faconvcpl')
      CALL FACONVCPL
    CASE ('lfistress')
!     CALL LFISTRESS
    CASE ('lfi_alt_merge')
      CALL LFI_ALT_MERGE
    CASE ('lfi_alt_copy')
      CALL LFI_ALT_COPY
    CASE ('lfi_alt_size')
      CALL LFI_ALT_SIZE
    CASE ('lfi_alt_pack')
      CALL LFI_ALT_PACK
    CASE ('lfi_alt_remove')
      CALL LFI_ALT_REMOVE
    CASE ('lfi_alt_index')
      CALL LFI_ALT_INDEX
    CASE ('faidx')
    CASE ('testfa','TESTFA')
      CALL TESTFA
    CASE ('tstlfi','test')
      CALL TSTLFI
    CASE ('faempty')
      CALL FAEMPTY
    CASE ('facat')
      CALL FACAT
    CASE ('fadate')
      CALL FADATE
    CASE ('datefa')
      CALL DATEFA
    CASE ('lficat','cat')
      CALL LFICAT
    CASE ('lfisplit','split')
      CALL LFISPLIT
    CASE ('lfidiff','diff')
      CALL LFIDIFF
    CASE ('falist','FALIST')
      CALL FALIST
    CASE ('lfilist','list')
      CALL LFILIST
!   CASE ('testfagrib')
!     CALL TESTFAGRIB
    CASE ('extractgrib')
      CALL EXTRACTGRIB
    CASE ('lfifactm')
      CALL LFIFACTM
    CASE ('faprogrid')
      CALL FAPROGRID
    CASE ('faconvgrib')
      CALL FACONVGRIB
    CASE ('fadiff')
      CALL FADIFF
!   CASE ('fastat')
!     CALL FASTAT
!   CASE ('lfitestread')
!     CALL LFITESTREAD
!   CASE ('lfitestwrite')
!     CALL LFITESTWRITE
    CASE ('lfixxx')
      CALL LFIXXX
    CASE DEFAULT
      PRINT *, 'Unknown program '//TRIM(CLNAME)
      CALL USAGE
  END SELECT

ENDIF

STOP

CONTAINS

SUBROUTINE USAGE
PRINT *, 'Usage: '//TRIM(CLMAIN)//' name [arg-1] [arg-1] ... [arg-N]'
PRINT *, 'where <name> is the name of a tool'
PRINT *, 'and [arg-1] [arg-1] ... [arg-N] are the arguments of <name>'
PRINT *, ' '
PRINT *, 'Existing programs today :'
PRINT *, '  faidx           Create an index for FA files'
PRINT *, '  testfa          Interactive tool to test FA software'
PRINT *, '  tstlfi          Interactive tool to test LFI software'
PRINT *, '  faempty         Create empty FA file'
PRINT *, '  facat           Concatenate several FA files'
PRINT *, '  datefa          Display the date of a FA file' 
PRINT *, '  lficat          Display the catalog of a LFI file'
PRINT *, '  lfisplit        Split a LFI file into several ones'
PRINT *, '  lfilist         List all fields in a LFI file'
PRINT *, '  lfidiff         Make the difference between two LFI files'
PRINT *, ' '
PRINT *, '  --help to print this help again'

END SUBROUTINE USAGE

SUBROUTINE LC (CLSTR)

CHARACTER(*), INTENT(IN OUT) :: CLSTR
INTEGER :: I

DO I = 1, LEN(CLSTR)
  SELECT CASE(CLSTR(I:I))
    CASE("A":"Z")
      CLSTR(I:I) = ACHAR(IACHAR(CLSTR(I:I))+32)
  END SELECT
END DO  

END SUBROUTINE LC

END PROGRAM LFITOOLS

