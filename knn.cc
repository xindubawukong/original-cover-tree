#include "cover_tree.h"

// Compute the k nearest neighbors

int main(int argc, char *argv[])
{
  int k = atoi(argv[1]);
  v_array<point> set_of_points = parse_points(fopen(argv[2],"r"));
  v_array<point> set_of_queries = parse_points(fopen(argv[3],"r"));

  node top = batch_create(set_of_points);
  node top_query = batch_create(set_of_queries);
  
  v_array<v_array<point> > res;
  k_nearest_neighbor(top,top_query,res,k);
  
  printf("Printing results\n");
  for (int i = 0; i < res.index; i++)
    {
      for (int j = 0; j<res[i].index; j++)
	print(res[i][j]);
      printf("\n");
    }
  printf("results printed\n");
}
