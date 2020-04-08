#include "atlas/atlas_f.h"

module do_halfdiff_mod

use atlas_module

implicit none


public :: do_halfdiff

private

interface
function do_halfdiff__do_halfdiff (fs, field) bind (C,name="do_halfdiff__do_halfdiff")
    use iso_c_binding, only: c_ptr
    type (c_ptr) :: do_halfdiff__do_halfdiff
    type (c_ptr), value :: fs
    type (c_ptr), value :: field
end function
end interface


contains

function do_halfdiff(fs,field) result(fieldset)
  use atlas_functionspace_c_binding
  class(atlas_functionspace_StructuredColumns), intent(in) :: fs
  type(atlas_Field),intent(in) :: field
  type(atlas_FieldSet) :: fieldset
  fieldset = atlas_FieldSet(do_halfdiff__do_halfdiff(fs%CPTR_PGIBUG_A,field%CPTR_PGIBUG_A))
  call fieldset%detach()
  call fieldset%return()
end function

end module

