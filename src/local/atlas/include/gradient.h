#pragma once

#include "atlas/grid.h"
#include "atlas/functionspace.h"

atlas::FieldSet
gradient (const atlas::functionspace::StructuredColumns & fs, const atlas::FieldSet & pgp);


extern "C"
{
atlas::field::FieldSetImpl * gradient__
(const atlas::functionspace::detail::StructuredColumns * fs, atlas::field::FieldSetImpl * pgp);
};
