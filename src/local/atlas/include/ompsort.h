#pragma once

#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <omp.h>
#include <stdlib.h>  


template <typename T, typename I, typename C>
void ompsort (std::vector<I> & ord, const std::vector<T> & vec, C cmp)
{
  I n = vec.size ();
  I nt = omp_get_max_threads ();
  I nr = nt;
  std::vector<I> off (nr + 1, 0);
  off[nr] = n;


  ord.resize (n);
  for (I i = 0; i < ord.size (); i++)
    ord[i] = i;

  // Sort sections
#pragma omp parallel 
  {
    int it = omp_get_thread_num ();
    I i1 = ((it + 0) * n) / nt;
    I i2 = ((it + 1) * n) / nt;
    i2 = std::min (i2, n);
    std::sort (ord.begin () + i1, ord.begin () + i2, cmp);
    off[it] = i1;
  }

  // Merge adjacents sections
  while (off.size () > 2)
    {

#pragma omp parallel for
      for (int ir = 0; ir < off.size () - 1; ir += 2)
        {
          auto cpr = [&] (I i1, I i2) 
          {
            std::vector<I> o;
            o.reserve (i2-i1); 
            std::copy (ord.begin () + i1, ord.begin () + i2, 
                       std::back_inserter (o));
            return o;
          };

          std::vector<I> ord1 = cpr (off[ir+0], off[ir+1]);
          std::vector<I> ord2 = cpr (off[ir+1], off[ir+2]);
  
          std::merge (ord1.begin (), ord1.end (), ord2.begin (), ord2.end (), ord.begin () + off[ir], cmp);
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

template <typename T, typename I>
void ompsort (std::vector<I> & ord, const std::vector<T> & vec)
{
  auto cmp = [&] (const I o1, const I o2) 
    { return vec[o1] < vec[o2]; };
  ompsort (ord, vec, cmp);
}
