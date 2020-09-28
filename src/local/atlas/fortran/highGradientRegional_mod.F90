#include "atlas/atlas_f.h"

module highGradientRegional_mod

use fckit_owned_object_module, only : fckit_owned_object
use atlas_functionspace_StructuredColumns_module
use atlas_fieldset_module, only : atlas_FieldSet
use atlas_Field_module, only : atlas_Field
use atlas_Grid_module, only : atlas_StructuredGrid
use atlas_GridDistribution_module

implicit none

private :: fckit_owned_object

public :: highGradientRegional

interface highGradientRegional
  module procedure :: highGradientRegional_field, highGradientRegional_fieldset
end interface 

private

interface
function highGradientRegional__ (fs, pgp, order) bind (C, name="highGradientRegional__")
    use iso_c_binding, only: c_ptr, c_int
    type (c_ptr) :: highGradientRegional__
    type (c_ptr), value :: fs
    type (c_ptr), value :: pgp
    integer (c_int), value :: order
end function
end interface

contains

function highGradientRegional_fieldset (fs, pgp, order) result (pgpg)
  use iso_c_binding, only: c_int
  type (atlas_functionspace_StructuredColumns), intent(in) :: fs
  type (atlas_FieldSet), intent(in) :: pgp  
  integer (c_int),       intent(in) :: order
  type (atlas_FieldSet) :: pgpg
  pgpg = atlas_FieldSet (highGradientRegional__ (fs%CPTR_PGIBUG_A, pgp%CPTR_PGIBUG_A, order)) 
  call pgpg%detach () ! Required here, because this comes from a temporary atlas::FieldSet object
                      ! and the implementation object had its count increased to
                      ! avoid deletion
  call pgpg%return ()
end function

function highGradientRegional_field (fs, f, order) result (pgpg)
  use iso_c_binding, only: c_int
  type (atlas_functionspace_StructuredColumns), intent(in) :: fs
  type (atlas_Field), intent(in) :: f
  integer (c_int),    intent(in) :: order
  type (atlas_FieldSet) :: pgp, pgpg
  pgp = atlas_FieldSet ()
  call pgp%add (f)
  pgpg = atlas_FieldSet (highGradientRegional__ (fs%CPTR_PGIBUG_A, pgp%CPTR_PGIBUG_A, order)) 
  call pgpg%detach () 
  call pgpg%return ()
  call pgp%final ()
end function

end module

