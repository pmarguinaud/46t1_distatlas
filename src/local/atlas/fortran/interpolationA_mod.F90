#include "atlas/atlas_f.h"

module interpolationA_mod

use iso_c_binding, only : c_int
use fckit_owned_object_module, only : fckit_owned_object
use atlas_functionspace_StructuredColumns_module
use atlas_fieldset_module, only : atlas_FieldSet
use atlas_field_module, only : atlas_Field
use atlas_GridDistribution_module

implicit none

private :: fckit_owned_object


public :: interpolationA

private

type, extends (fckit_owned_object) :: interpolationA
contains
  generic :: interpolate => interpolate_field, interpolate_fieldset
  procedure, private :: interpolate_field
  procedure, private :: interpolate_fieldset
  generic :: shuffle => shuffle_field, shuffle_fieldset
  procedure, private :: shuffle_field
  procedure, private :: shuffle_fieldset
  procedure, public :: getcnt
  procedure, public :: getoff
  procedure, private :: getlen
  procedure, public :: opt_avg 
  procedure, public :: opt_min
  procedure, public :: opt_max
  procedure, public :: opt_sum
end type

interface interpolationA
  module procedure interpolationA__ctor
end interface

interface

function interpolationA__new (dist1, fs1, dist2, fs2, ldopenmp) &
        & bind (C, name="interpolationA__new")
  use iso_c_binding, only : c_ptr, c_int
  type (c_ptr) :: interpolationA__new
  type (c_ptr), value :: dist1
  type (c_ptr), value :: fs1
  type (c_ptr), value :: dist2
  type (c_ptr), value :: fs2
  integer (c_int), value :: ldopenmp
end function

function interpolationA__interpolate (this, pgp1, opt) &
        & bind (C, name="interpolationA__interpolate")
  use iso_c_binding, only : c_ptr, c_int
  type (c_ptr) :: interpolationA__interpolate
  type (c_ptr), value :: this
  type (c_ptr), value :: pgp1
  integer (c_int), value :: opt
end function

function interpolationA__shuffle (this, pgp1) &
        & bind (C, name="interpolationA__shuffle")
  use iso_c_binding, only : c_ptr
  type (c_ptr) :: interpolationA__shuffle
  type (c_ptr), value :: this
  type (c_ptr), value :: pgp1
end function

function interpolationA__getlen (this) &
       & bind (C, name="interpolationA__getlen")
  use iso_c_binding, only : c_ptr
  integer :: interpolationA__getlen
  type (c_ptr), value :: this
end function

subroutine interpolationA__getcnt (this, cnt) &
          & bind (C, name="interpolationA__getcnt")
  use iso_c_binding, only : c_ptr, c_int
  type (c_ptr), value :: this
  integer (c_int), dimension(*) :: cnt
end subroutine

subroutine interpolationA__getoff (this, off) &
         & bind (C, name="interpolationA__getoff")
  use iso_c_binding, only : c_ptr, c_int
  type (c_ptr), value :: this
  integer (c_int), dimension(*) :: off
end subroutine

end interface

contains

function interpolationA__ctor (dist1, fs1, dist2, fs2, ldopenmp) result (this)
  type (interpolationA) :: this
  type (atlas_GridDistribution),                intent(in) :: dist1
  type (atlas_functionspace_StructuredColumns), intent(in) :: fs1
  type (atlas_GridDistribution),                intent(in) :: dist2
  type (atlas_functionspace_StructuredColumns), intent(in) :: fs2
  logical, optional,                            intent(in) :: ldopenmp
  integer (c_int) :: llopenmp
  llopenmp = 1
  if (present (ldopenmp)) then
    if (.not. ldopenmp) llopenmp = 0
  endif
  call this%reset_c_ptr (interpolationA__new (dist1%CPTR_PGIBUG_A, fs1%CPTR_PGIBUG_A, &
                       &                      dist2%CPTR_PGIBUG_A, fs2%CPTR_PGIBUG_A, llopenmp))
  call this%return ()
end function

function interpolate_fieldset (this, pgp1, opt) result (pgp2)
  use iso_c_binding, only : c_int
  class (interpolationA), intent (in) :: this
  type (atlas_FieldSet), intent(in) :: pgp1  
  integer (c_int), optional, intent (in) :: opt
  type (atlas_FieldSet) :: pgp2
  integer (c_int) :: opt_
  opt_ = 0
  if (present (opt)) opt_ = opt
  pgp2 = atlas_FieldSet (interpolationA__interpolate (this%CPTR_PGIBUG_A, pgp1%CPTR_PGIBUG_A, opt_)) 
  call pgp2%detach () ! Required here, because this comes from a temporary atlas::FieldSet object
                      ! and the implementation object had its count increased to
                      ! avoid deletion
  call pgp2%return ()
end function

function interpolate_field (this, f1, opt) result (f2)
  use iso_c_binding, only : c_int
  class (interpolationA), intent (in) :: this
  type (atlas_Field), intent(in) :: f1  
  integer (c_int), optional, intent (in) :: opt
  type (atlas_Field) :: f2
  type (atlas_FieldSet) :: pgp1, pgp2
  integer (c_int) :: opt_
  opt_ = 0
  if (present (opt)) opt_ = opt
  pgp1 = atlas_FieldSet ()
  call pgp1%add (f1)
  pgp2 = atlas_FieldSet (interpolationA__interpolate (this%CPTR_PGIBUG_A, pgp1%CPTR_PGIBUG_A, opt_)) 
  f2 = pgp2%field (1)
  call pgp2%detach () 
  call pgp2%final ()
  call pgp1%final ()
  call f2%return ()
end function

function shuffle_fieldset (this, pgp1) result (pgp2)
  class (interpolationA), intent (in) :: this
  type (atlas_FieldSet), intent(in) :: pgp1  
  type (atlas_FieldSet) :: pgp2
  pgp2 = atlas_FieldSet (interpolationA__shuffle (this%CPTR_PGIBUG_A, pgp1%CPTR_PGIBUG_A)) 
  call pgp2%detach () 
  call pgp2%return ()
end function

function shuffle_field (this, f1) result (f2)
  class (interpolationA), intent (in) :: this
  type (atlas_Field), intent(in) :: f1  
  type (atlas_Field) :: f2
  type (atlas_FieldSet) :: pgp1, pgp2
  pgp1 = atlas_FieldSet ()
  call pgp1%add (f1)
  pgp2 = atlas_FieldSet (interpolationA__shuffle (this%CPTR_PGIBUG_A, pgp1%CPTR_PGIBUG_A)) 
  f2 = pgp2%field (1)
  call pgp2%detach () 
  call pgp2%final ()
  call pgp1%final ()
  call f2%return ()
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

integer (c_int) function opt_avg (this)
  class (interpolationA), intent (in) :: this
  opt_avg = 0
end function

integer (c_int) function opt_min (this)
  class (interpolationA), intent (in) :: this
  opt_min = 1
end function

integer (c_int) function opt_max (this)
  class (interpolationA), intent (in) :: this
  opt_max = 2
end function

integer (c_int) function opt_sum (this)
  class (interpolationA), intent (in) :: this
  opt_sum = 3
end function


end module

