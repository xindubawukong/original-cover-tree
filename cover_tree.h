#pragma once
#include<string.h>
#include<math.h>
#include<list>
#include<stdlib.h>
#include<assert.h>
#include "point.h"

namespace original_cover_tree {

struct node {
  point p;
  double max_dist;  // The maximum distance to any grandchild.
  double parent_dist; // The distance to the parent.
  node* children;
  int num_children; // The number of children.
  int scale; // Essentially, an upper bound on the distance to any child.
};

void print(int depth, node &top_node);

//construction
node new_leaf(const point &p);
node batch_create(v_array<point> points, distance_function metric, const void* context = NULL);
// These functions free allocations, but never the caller-owned points.
void free_tree(node &top_node);
void free_results(v_array<v_array<point> > &results);
//node insert(point, node *top_node); // not yet implemented
//void remove(point, node *top_node); // not yet implemented
//query
void k_nearest_neighbor(const node &tope_node, const node &query, 
			v_array<v_array<point> > &results, int k,
                        distance_function metric, const void* context = NULL);
void epsilon_nearest_neighbor(const node &tope_node, const node &query, 
			      v_array<v_array<point> > &results, double epsilon,
                              distance_function metric, const void* context = NULL);
void unequal_nearest_neighbor(const node &tope_node, const node &query, 
			      v_array<v_array<point> > &results,
                              distance_function metric, const void* context = NULL);
//information gathering
int height_dist(const node top_node,v_array<int> &heights);
void breadth_dist(const node top_node,v_array<int> &breadths);
void depth_dist(int top_scale, const node top_node,v_array<int> &depths);

} // namespace original_cover_tree
