#pragma once
#include "stack.h"
#include<stdio.h>
#include<stdlib.h>

namespace original_cover_tree {

typedef const void* point;
// Exact when <= upper_bound; otherwise any value > upper_bound is valid.
typedef double (*distance_function)(point, point, double upper_bound, const void* context);

double complete_distance(point v1, point v2);
double distance(point v1, point v2, double upper_bound, const void* context);
v_array<point > parse_points(FILE *input);
void print(point &p);

} // namespace original_cover_tree
