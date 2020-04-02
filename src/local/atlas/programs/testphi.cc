
#include "ompsort.h"


int main (int argc, char * argv[]) 
{
  const int n = 100;
  std::vector<int> vec (n);  
  std::vector<int> ord (n);

  srand (0);

  for (int i = 0; i < n; i++)
    vec[i] = rand () % n; 

  auto cmp = [&] (const int o1, const int o2) 
    { return vec[o1] > vec[o2]; };

  ompsort (ord, vec, cmp);
  
  for (int i = 0; i < vec.size (); i++)
    printf (" %8d > %8d\n", i, vec[ord[i]]);


  return 0;
}

