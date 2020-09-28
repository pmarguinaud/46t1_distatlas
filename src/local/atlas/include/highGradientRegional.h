#pragma once

#include "atlas/grid.h"
#include "atlas/functionspace.h"

template <typename T, int O>
atlas::FieldSet
highGradientRegional (const atlas::functionspace::StructuredColumns & fs, const atlas::FieldSet & pgp);

extern "C"
{
atlas::field::FieldSetImpl * highGradientRegional__
(const atlas::functionspace::detail::StructuredColumns * fs, atlas::field::FieldSetImpl * pgp, const int order);
};

