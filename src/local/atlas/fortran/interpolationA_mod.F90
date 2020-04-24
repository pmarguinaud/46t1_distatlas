#include "atlas/atlas_f.h"

module interpolationA_mod

use fckit_owned_object_module, only : fckit_owned_object
use atlas_functionspace_StructuredColumns_module
use atlas_fieldset_module, only : atlas_FieldSet
use atlas_GridDistribution_module

implicit none

private :: fckit_owned_object

public :: interpolationA

private

type, extends (fckit_owned_object) :: interpolationA
contains
  procedure, public :: interpolate
  procedure, public :: shuffle
  procedure, public :: getcnt
  procedure, public :: getoff
  procedure, private :: getlen
end type

interface interpolationA
  module procedure interpolationA__ctor
end interface

interface

function interpolationA__new (dist1, fs1, dist2, fs2) bind (C, name="interpolationA__new")
  use iso_c_binding, only : c_ptr
  type (c_ptr) :: interpolationA__new
  type (c_ptr), value :: dist1
  type (c_ptr), value :: fs1
  type (c_ptr), value :: dist2
  type (c_ptr), value :: fs2
end function

function interpolationA__interpolate (this, pgp1) bind (C, name="interpolationA__interpolate")
  use iso_c_binding, only : c_ptr
  type (c_ptr) :: interpolationA__interpolate
  type (c_ptr), value :: this
  type (c_ptr), value :: pgp1
end function

function interpolationA__shuffle (this, pgp1) bind (C, name="interpolationA__shuffle")
  use iso_c_binding, only : c_ptr
  type (c_ptr) :: interpolationA__shuffle
  type (c_ptr), value :: this
  type (c_ptr), value :: pgp1
end function

function interpolationA__getlen (this) bind (C, name="interpolationA__getlen")
  use iso_c_binding, only : c_ptr
  integer :: interpolationA__getlen
  type (c_ptr), value :: this
end function

subroutine interpolationA__getcnt (this, cnt) bind (C, name="interpolationA__getcnt")
  use iso_c_binding, only : c_ptr, c_int
  type (c_ptr), value :: this
  integer (c_int), dimension(*) :: cnt
end subroutine

subroutine interpolationA__getoff (this, off) bind (C, name="interpolationA__getoff")
  use iso_c_binding, only : c_ptr, c_int
  type (c_ptr), value :: this
  integer (c_int), dimension(*) :: off
end subroutine

end interface

contains

function interpolationA__ctor (dist1, fs1, dist2, fs2) result (this)
  type (interpolationA) :: this
  type (atlas_GridDistribution),                intent(in) :: dist1
  type (atlas_functionspace_StructuredColumns), intent(in) :: fs1
  type (atlas_GridDistribution),                intent(in) :: dist2
  type (atlas_functionspace_StructuredColumns), intent(in) :: fs2
  call this%reset_c_ptr (interpolationA__new (dist1%CPTR_PGIBUG_A, fs1%CPTR_PGIBUG_A, &
                       &                      dist2%CPTR_PGIBUG_A, fs2%CPTR_PGIBUG_A))
  call this%return ()
end function

function interpolate (this, pgp1) result (pgp2)
  class (interpolationA), intent (in) :: this
  type (atlas_FieldSet), intent(in) :: pgp1  
  type (atlas_FieldSet) :: pgp2
  pgp2 = atlas_FieldSet (interpolationA__interpolate (this%CPTR_PGIBUG_A, pgp1%CPTR_PGIBUG_A)) 
  call pgp2%detach () ! Required here, because this comes from a temporary atlas::FieldSet object
                      ! and the implementation object had its count increased to
                      ! avoid deletion
  call pgp2%return ()
end function

function shuffle (this, pgp1) result (pgp2)
  class (interpolationA), intent (in) :: this
  type (atlas_FieldSet), intent(in) :: pgp1  
  type (atlas_FieldSet) :: pgp2
  pgp2 = atlas_FieldSet (interpolationA__shuffle (this%CPTR_PGIBUG_A, pgp1%CPTR_PGIBUG_A)) 
  call pgp2%detach () 
  call pgp2%return ()
end function

function getlen (this)
  class (interpolationA), intent (in) :: this
  integer :: getlen
  getlen = interpolationA__getlen (this%CPTR_PGIBUG_A)
end function

function getcnt (this) result (cnt)
  use iso_c_binding, only : c_int
  class (interpolationA), intent (in) :: this
  integer (kind=c_int), allocatable :: cnt (:)
  allocate (cnt (this%getlen ()))
  call interpolationA__getcnt (this%CPTR_PGIBUG_A, cnt)
end function

function getoff (this) result (off)
  use iso_c_binding, only : c_int
  class (interpolationA), intent (in) :: this
  integer (kind=c_int), allocatable :: off (:)
  allocate (off (this%getlen ()))
  call interpolationA__getoff (this%CPTR_PGIBUG_A, off)
end function

end module

