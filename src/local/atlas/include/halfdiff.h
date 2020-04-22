#pragma once

#include "atlas/functionspace.h"
#include "atlas/field.h"


atlas::FieldSet halfdiff (const atlas::functionspace::StructuredColumns & fs, const atlas::FieldSet & pgp);

extern "C"
{
atlas::field::FieldSetImpl * halfdiff__ (atlas::functionspace::detail::StructuredColumns *, atlas::field::FieldSetImpl *);
};

