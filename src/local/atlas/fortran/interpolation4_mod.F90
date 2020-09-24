#include "atlas/atlas_f.h"

module interpolation4_mod

use iso_c_binding, only : c_int
use fckit_owned_object_module, only : fckit_owned_object
use atlas_functionspace_StructuredColumns_module
use atlas_fieldset_module, only : atlas_FieldSet
use atlas_field_module, only : atlas_Field
use atlas_GridDistribution_module

implicit none

private :: fckit_owned_object

public :: interpolation4

private

type, extends (fckit_owned_object) :: interpolation4
contains
  generic :: interpolate => interpolate_field, interpolate_fieldset
  procedure, private :: interpolate_field, interpolate_fieldset
end type

interface interpolation4
  module procedure interpolation4__ctor
end interface


interface
function interpolation4__new (dist1, fs1, dist2, fs2, ldopenmp) bind (C, name="interpolation4__new")
    use iso_c_binding, only: c_ptr, c_int
    type (c_ptr) :: interpolation4__new
    type (c_ptr), value :: dist1
    type (c_ptr), value :: fs1
    type (c_ptr), value :: dist2
    type (c_ptr), value :: fs2
    integer (c_int), value :: ldopenmp
end function

function interpolation4__interpolate (this, pgp1) bind (C, name="interpolation4__interpolate")
    use iso_c_binding, only: c_ptr
    type (c_ptr) :: interpolation4__interpolate
    type (c_ptr), value :: this
    type (c_ptr), value :: pgp1
end function
end interface

contains

function interpolation4__ctor (dist1, fs1, dist2, fs2, ldopenmp) result (this)
  type (interpolation4) :: this
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
  call this%reset_c_ptr (interpolation4__new (dist1%CPTR_PGIBUG_A, fs1%CPTR_PGIBUG_A, &
                       &                      dist2%CPTR_PGIBUG_A, fs2%CPTR_PGIBUG_A, llopenmp))
  call this%return ()
end function

function interpolate_fieldset (this, pgp1) result (pgp2)
  class (interpolation4), intent (in) :: this
  type (atlas_FieldSet), intent(in) :: pgp1  
  type (atlas_FieldSet) :: pgp2
  pgp2 = atlas_FieldSet (interpolation4__interpolate (this%CPTR_PGIBUG_A, pgp1%CPTR_PGIBUG_A)) 
  call pgp2%detach () ! Required here, because this comes from a temporary atlas::FieldSet object
                      ! and the implementation object had its count increased to
                      ! avoid deletion
  call pgp2%return ()
end function

function interpolate_field (this, f1) result (f2)
  use iso_c_binding, only : c_int
  class (interpolation4), intent (in) :: this
  type (atlas_Field), intent(in) :: f1  
  type (atlas_Field) :: f2
  type (atlas_FieldSet) :: pgp1, pgp2
  pgp1 = atlas_FieldSet ()
  call pgp1%add (f1)
  pgp2 = atlas_FieldSet (interpolation4__interpolate (this%CPTR_PGIBUG_A, pgp1%CPTR_PGIBUG_A)) 
  f2 = pgp2%field (1)
  call pgp2%detach () 
  call pgp2%final ()
  call pgp1%final ()
  call f2%return ()
end function

end module

