#include "cover_tree.h"

using namespace original_cover_tree;

// Compute the k nearest neighbors

int main(int argc, char *argv[])
{
  int k = atoi(argv[1]);
  FILE* input = fopen(argv[2],"r");
  FILE* queries = fopen(argv[3],"r");
  v_array<point> set_of_points = parse_points(input);
  v_array<point> set_of_queries = parse_points(queries);
  fclose(input);
  fclose(queries);

  node top = batch_create(set_of_points, distance);
  node top_query = batch_create(set_of_queries, distance);
  
  v_array<v_array<point> > res;
  k_nearest_neighbor(top,top_query,res,k, distance);
  
  printf("Printing results\n");
  for (int i = 0; i < res.index; i++)
    {
      for (int j = 0; j<res[i].index; j++)
	print(res[i][j]);
      printf("\n");
    }
  printf("results printed\n");
  free_results(res);
  free_tree(top);
  free_tree(top_query);
  for (int i = 0; i < set_of_points.index; i++) free(const_cast<void*>(set_of_points[i]));
  for (int i = 0; i < set_of_queries.index; i++) free(const_cast<void*>(set_of_queries[i]));
  free(set_of_points.elements);
  free(set_of_queries.elements);
}
