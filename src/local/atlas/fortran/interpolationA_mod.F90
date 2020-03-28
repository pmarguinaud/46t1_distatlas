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
end type

interface interpolationA
  module procedure interpolationA__ctor
end interface


interface
function interpolationA__new (dist1, fs1, dist2, fs2) bind (C, name="interpolationA__new")
    use iso_c_binding, only: c_ptr
    type (c_ptr) :: interpolationA__new
    type (c_ptr), value :: dist1
    type (c_ptr), value :: fs1
    type (c_ptr), value :: dist2
    type (c_ptr), value :: fs2
end function

function interpolationA__interpolate (this, pgp1) bind (C, name="interpolationA__interpolate")
    use iso_c_binding, only: c_ptr
    type (c_ptr) :: interpolationA__interpolate
    type (c_ptr), value :: this
    type (c_ptr), value :: pgp1
end function
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

end module

