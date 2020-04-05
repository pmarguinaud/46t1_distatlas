#pragma once

#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <omp.h>
#include <stdlib.h>  

#include "quicksort.h"


template <typename I, typename C>
void ompsort (I b, I e, C cmp)
{
  using tt = typename I::value_type;

  size_t n = e - b;
  int nt = omp_get_max_threads ();
  int nr = nt;
  std::vector<size_t> off (nr + 1, 0);
  off[nr] = n;

  // Sort sections
#pragma omp parallel 
  {
    int it = omp_get_thread_num ();
    size_t i1 = ((it + 0) * n) / nt;
    size_t i2 = ((it + 1) * n) / nt;
    i2 = std::min (i2, n);
//  std::sort (b + i1, b + i2, cmp);
    quicksort (b + i1, b + i2, cmp);
    off[it] = i1;
  }

  // Merge adjacents sections
  while (off.size () > 2)
    {

#pragma omp parallel for
      for (int ir = 0; ir < off.size () - 2; ir += 2)
        {
          auto cpr = [&] (size_t i1, size_t i2) 
          {
            std::vector<tt> o;
            o.reserve (i2-i1); 
            std::copy (b + i1, b + i2, 
                       std::back_inserter (o));
            return o;
          };

          std::vector<tt> ord1 = cpr (off[ir+0], off[ir+1]);
          std::vector<tt> ord2 = cpr (off[ir+1], off[ir+2]);
  
          std::merge (ord1.begin (), ord1.end (), ord2.begin (), ord2.end (), b + off[ir], cmp);
        }
                    

      for (int ir = 0; ir < off.size () / 2; ir++)
        off[ir] = off[2*ir];

      int kr = off.size () - 1;
      if (kr % 2 == 1)
        {
          off[kr/2+1] = off[kr];
          off.resize (kr/2+2);
        }
      else
        {
          off[kr/2+0] = off[kr];
          off.resize (kr/2+1);
        }

    }

}

