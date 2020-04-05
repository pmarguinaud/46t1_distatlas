
#pragma once

#include <utility>
#include <algorithm>
#include <stdlib.h>

template <typename I, typename C>
void quicksort (I b, I e, C cmp)
{
  class stack_t
  {
  public:
    stack_t () { top = &arr[0]; }
    void push (I lo, I hi) { top->lo = lo; top->hi = hi; top++; }
    void pop (I & lo, I & hi) { top--; lo = top->lo; hi = top->hi; }
    size_t size () const { return top - &arr[0]; }
  private:
    struct node_t { I lo; I hi; };
    node_t arr[64];
    node_t * top;
  };

  const int max = 4;

  if (e - b == 0)
    return;

  if (e - b > max)
    {
      stack_t stack;

      I lo = b, hi = e - 1;
      stack.push (e, e);

      while (stack.size ())
        {
          I le, ri, mi = lo + (hi - lo) / 2;

          if (cmp (*mi, *lo))
            std::swap (*mi, *lo);
          if (cmp (*hi, *mi))
            std::swap (*mi, *hi);
          else
            goto jump_over;

          if (cmp (*mi, *lo))
            std::swap (*mi, *lo);

jump_over:
          le = lo + 1;
          ri = hi - 1;

          do
            {
              while (cmp (*le, *mi)) le++;
              while (cmp (*mi, *ri)) ri--;
              if (le < ri)
                {
                  std::swap (*le, *ri);
                  if (mi == le)
                    mi = ri;
                  else if (mi == ri)
                    mi = le;
                  le++;
                  ri--;
                }
              else if (le == ri)
                {
                  le++;
                  ri--;
                  break;
                }
            }
          while (le <= ri);

          if ((ri - lo) <= max)
            {
              if ((hi - le) <= max)
                stack.pop (lo, hi);
              else
                lo = le;
            }
          else if ((hi - le) <= max)
            {
              hi = ri;
            }
          else if ((ri - lo) > (hi - le))
            {
              stack.push (lo, ri);
              lo = le;
            }
          else
            {
              stack.push (le, hi);
              hi = ri;
            }
        }
    }

  I end = e - 1, tmp = b, thr = std::min (end, b + max), run;
  for (run = tmp + 1; run <= thr; run++)
    if (cmp (*run, *tmp))
      tmp = run;
  if (tmp != b)
    std::swap (*tmp, *b);
  run = b + 1;
  while (++run <= end)
    {
      tmp = run - 1;
      while (cmp (*run, *tmp))
        tmp--;
      tmp++;
      if (tmp != run)
        {
          I x1 = tmp, x2 = run;
          typename I::value_type v2 = *x2;
          for (I x = x2; x > x1; x--)
            *x = *(x - 1);
          *x1 = v2;
        }
    }

}


