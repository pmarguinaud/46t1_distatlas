#pragma once

#include "atlas/grid.h"
#include "atlas/functionspace.h"

void 
rotate (const atlas::functionspace::StructuredColumns & fs, 
        const atlas::StructuredGrid & grid, atlas::FieldSet & pgp);

atlas::FieldSet
gradient (const atlas::functionspace::StructuredColumns & fs, const atlas::FieldSet & pgp);


extern "C"
{
void rotate__
(const atlas::functionspace::detail::StructuredColumns * fs, 
 const atlas::grid::detail::grid::Structured * grid, atlas::field::FieldSetImpl * pgp);
atlas::field::FieldSetImpl * gradient__
(const atlas::functionspace::detail::StructuredColumns * fs, atlas::field::FieldSetImpl * pgp);
};
