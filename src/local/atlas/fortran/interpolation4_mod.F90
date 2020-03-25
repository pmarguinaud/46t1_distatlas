#include "atlas/atlas_f.h"

module interpolation4_mod

use fckit_owned_object_module, only : fckit_owned_object
use atlas_functionspace_StructuredColumns_module
use atlas_fieldset_module, only : atlas_FieldSet
use atlas_GridDistribution_module

implicit none

private :: fckit_owned_object

public :: interpolation4

private

type, extends (fckit_owned_object) :: interpolation4
contains
  procedure, public :: interpolate
end type

interface interpolation4
  module procedure interpolation4__ctor
end interface


interface
function interpolation4__new (dist1, fs1, dist2, fs2) bind (C, name="interpolation4__new")
    use iso_c_binding, only: c_ptr
    type (c_ptr) :: interpolation4__new
    type (c_ptr), value :: dist1
    type (c_ptr), value :: fs1
    type (c_ptr), value :: dist2
    type (c_ptr), value :: fs2
end function

function interpolation4__interpolate (this, pgp1) bind (C, name="interpolation4__interpolate")
    use iso_c_binding, only: c_ptr
    type (c_ptr) :: interpolation4__interpolate
    type (c_ptr), value :: this
    type (c_ptr), value :: pgp1
end function
end interface

contains

function interpolation4__ctor (dist1, fs1, dist2, fs2) result (this)
  type (interpolation4) :: this
  type (atlas_GridDistribution),                intent(in) :: dist1
  type (atlas_functionspace_StructuredColumns), intent(in) :: fs1
  type (atlas_GridDistribution),                intent(in) :: dist2
  type (atlas_functionspace_StructuredColumns), intent(in) :: fs2
  call this%reset_c_ptr (interpolation4__new (dist1%CPTR_PGIBUG_A, fs1%CPTR_PGIBUG_A, &
                       &                      dist2%CPTR_PGIBUG_A, fs2%CPTR_PGIBUG_A))
  call this%return ()
end function

function interpolate (this, pgp1) result (pgp2)
  class (interpolation4), intent (in) :: this
  type (atlas_FieldSet), intent(in) :: pgp1  
  type (atlas_FieldSet) :: pgp2
  pgp2 = atlas_FieldSet (interpolation4__interpolate (this%CPTR_PGIBUG_A, pgp1%CPTR_PGIBUG_A)) 
  call pgp2%detach () ! Required here, because this comes from a temporary atlas::FieldSet object
                      ! and the implementation object had its count increased to
                      ! avoid deletion
  call pgp2%return ()
end function

end module

