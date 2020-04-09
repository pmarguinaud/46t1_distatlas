#include "atlas/atlas_f.h"

module gradient_mod

use fckit_owned_object_module, only : fckit_owned_object
use atlas_functionspace_StructuredColumns_module
use atlas_fieldset_module, only : atlas_FieldSet
use atlas_GridDistribution_module

implicit none

private :: fckit_owned_object

public :: gradient

private

interface
function gradient__ (pgp) bind (C, name="gradient__")
    use iso_c_binding, only: c_ptr
    type (c_ptr) :: gradient__
    type (c_ptr), value :: pgp
end function
end interface

contains

function gradient (pgp) result (pgpg)
  type (atlas_FieldSet), intent(in) :: pgp  
  type (atlas_FieldSet) :: pgpg
  pgpg = atlas_FieldSet (gradient__ (pgp%CPTR_PGIBUG_A)) 
  call pgpg%detach () ! Required here, because this comes from a temporary atlas::FieldSet object
                      ! and the implementation object had its count increased to
                      ! avoid deletion
  call pgpg%return ()
end function

end module

