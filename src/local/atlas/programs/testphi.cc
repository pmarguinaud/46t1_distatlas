#include <stdlib.h>
#include "ompsort.h"


int main (int argc, char * argv[]) 
{
  const int n = atoi (argv[1]);
  std::vector<int> vec (n);  
  std::vector<int> ord (n);

  srand (0);

  for (int i = 0; i < n; i++)
    {
      ord[i] = i;
      vec[i] = rand () % n; 
    }

#ifdef UNDEF
  for (int i = 0; i < n; i++)
    printf (" %8d > %8d\n", i, vec[i]);
  printf ("----\n");
#endif

  auto cmp = [&] (const int o1, const int o2) 
    { return vec[o1] < vec[o2]; };

  ompsort (ord.begin (), ord.end (), cmp);
  

  for (int i = 0; i < vec.size ()-1; i++)
    {
      int j = i + 1;
      if (cmp (ord[j], ord[i]))
        printf (" %8d > %8d, %8d\n", i, vec[ord[i]], vec[ord[j]]);
    }

  return 0;
}

