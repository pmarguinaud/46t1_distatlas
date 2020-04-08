#pragma once

#include "atlas/functionspace.h"
#include "atlas/field.h"


atlas::FieldSet do_halfdiff (const atlas::functionspace::StructuredColumns & fs, const atlas::Field & field);

extern "C"
{
atlas::field::FieldSetImpl * do_halfdiff__do_halfdiff (atlas::functionspace::detail::StructuredColumns *, atlas::field::FieldImpl *);
};

