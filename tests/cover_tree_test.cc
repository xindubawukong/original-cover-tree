#include "cover_tree.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace original_cover_tree;
using coordinates = std::array<double, 3>;
struct metric_state { double factor; bool manhattan; };

double numeric_distance(point a, point b, double, const void* context)
{
  const auto &x = *static_cast<const coordinates*>(a), &y = *static_cast<const coordinates*>(b);
  const auto &state = *static_cast<const metric_state*>(context);
  double result = 0.;
  for (int i = 0; i < 3; i++)
    result = state.manhattan ? result + std::abs(x[i] - y[i]) : std::hypot(result, x[i] - y[i]);
  return state.factor * result;
}
void require(bool condition) { if (!condition) throw std::runtime_error("correctness check failed"); }

v_array<point> references(const std::vector<coordinates> &points)
{
  v_array<point> result;
  for (const auto &p : points) push(result, static_cast<point>(&p));
  return result;
}

void check(const node &tree, const std::vector<coordinates> &points,
           const std::vector<coordinates> &queries, const metric_state &metric, int k, bool batch)
{
  auto query_points = references(queries);
  node query_tree = batch ? batch_create(query_points, numeric_distance, &metric) : new_leaf(query_points[0]);
  v_array<v_array<point> > rows;
  k_nearest_neighbor(tree, query_tree, rows, k, numeric_distance, &metric);
  require(rows.index == (batch ? static_cast<int>(queries.size()) : 1));
  std::vector<bool> seen(queries.size());
  for (int i = 0; i < rows.index; i++) {
    auto &row = rows[i];
    int query_id = 0;
    while (query_id < query_points.index && query_points[query_id] != row[0]) ++query_id;
    require(query_id < query_points.index && !seen[query_id]);
    seen[query_id] = true;
    std::vector<double> expected, actual;
    for (const auto &p : points) expected.push_back(numeric_distance(row[0], &p, 0., &metric));
    for (int j = 1; j < row.index; j++) {
      require(std::any_of(points.begin(), points.end(), [&](const coordinates &p) { return &p == row[j]; }));
      for (int previous = 1; previous < j; previous++) require(row[previous] != row[j]);
      actual.push_back(numeric_distance(row[0], row[j], 0., &metric));
    }
    std::sort(expected.begin(), expected.end());
    std::sort(actual.begin(), actual.end());
    require(actual.size() >= static_cast<std::size_t>(k));
    for (int j = 0; j < k; j++)
      require(std::abs(actual[j] - expected[j]) <= 1e-12 * std::max(1e-200, expected[j]));
    for (double d : actual) require(d <= expected[k - 1] + 1e-12 * std::max(1e-200, expected[k - 1]));
  }
  free_results(rows);
  free_tree(query_tree);
  free(query_points.elements);
}

void run_case(const std::vector<coordinates> &points, const std::vector<coordinates> &queries,
              const metric_state &metric)
{
  auto refs = references(points);
  node tree = batch_create(refs, numeric_distance, &metric);
  for (int k : {1, std::min(7, refs.index), refs.index}) {
    check(tree, points, queries, metric, k, false);
    check(tree, points, queries, metric, k, true);
  }
  free_tree(tree);
  free(refs.elements);
}

double hamming(point a, point b, double, const void* context)
{
  const auto &x = *static_cast<const std::string*>(a), &y = *static_cast<const std::string*>(b);
  double d = 0.;
  for (std::size_t i = 0; i < x.size(); i++) d += x[i] != y[i];
  return d * *static_cast<const double*>(context);
}

int main()
{
  const metric_state euclidean = {1., false}, manhattan = {2.5, true};
  std::vector<coordinates> points(97), queries(11);
  unsigned random = 42;
  for (auto *dataset : {&points, &queries})
    for (auto &p : *dataset) for (auto &x : p) {
      random = 1664525u * random + 1013904223u;
      x = static_cast<double>(random) / 4294967296.;
    }
  run_case(points, queries, euclidean);
  run_case(points, queries, manhattan);
  run_case({{{1., 2., 3.}}}, queries, euclidean);
  run_case({{{0., 0., 0.}}, {{1e-20, 0., 0.}}, {{3e-20, 0., 0.}}},
           {{{2e-20, 0., 0.}}, {{0., 0., 0.}}}, euclidean);
  run_case({{{1., 0., 0.}}, {{1. + 1e-12, 0., 0.}}, {{1. + 3e-12, 0., 0.}}},
           {{{1. + 1.1e-12, 0., 0.}}}, euclidean);
  std::vector<coordinates> scales;
  for (int exponent = -200; exponent <= 200; exponent += 10)
    scales.push_back({{std::pow(10., exponent), 0., 0.}});
  run_case(scales, {{{1e-199, 0., 0.}}, {{0., 0., 0.}}, {{1e100, 0., 0.}}}, euclidean);

  // Calls with different k and different metric contexts must not share state.
  auto refs = references(points);
  node trees[] = {batch_create(refs, numeric_distance, &euclidean), batch_create(refs, numeric_distance, &manhattan)};
  std::array<std::exception_ptr, 4> failures;
  std::array<std::thread, 4> workers;
  for (int i = 0; i < 4; i++) workers[i] = std::thread([&, i] {
    try {
      for (int repeat = 0; repeat < 8; repeat++)
        check(trees[i % 2], points, queries, i % 2 ? manhattan : euclidean, i + 1, true);
    } catch (...) { failures[i] = std::current_exception(); }
  });
  for (auto &worker : workers) worker.join();
  for (auto failure : failures) if (failure) std::rethrow_exception(failure);
  for (auto &tree : trees) free_tree(tree);
  free(refs.elements);

  std::vector<coordinates> duplicates(65537, {{1., 2., 3.}});
  refs = references(duplicates);
  node tree = batch_create(refs, numeric_distance, &euclidean);
  v_array<v_array<point> > rows;
  k_nearest_neighbor(tree, new_leaf(&duplicates[0]), rows, 10, numeric_distance, &euclidean);
  require(rows.index == 1 && rows[0].index == 65538);
  free_results(rows);
  free_tree(tree);
  free(refs.elements);

  const std::string strings[] = {"aaaa", "aabb", "bbbb", "aabb"};
  const std::string query = "aaab";
  const double factor = 3.5;
  v_array<point> words;
  for (const auto &s : strings) push(words, static_cast<point>(&s));
  tree = batch_create(words, hamming, &factor);
  k_nearest_neighbor(tree, new_leaf(&query), rows, 1, hamming, &factor);
  require(rows.index == 1 && rows[0].index == 4);
  for (int i = 1; i < rows[0].index; i++) require(hamming(&query, rows[0][i], 0., &factor) == factor);
  free_results(rows);
  epsilon_nearest_neighbor(tree, new_leaf(&query), rows, factor, hamming, &factor);
  require(rows.index == 1 && rows[0].index == 4);
  free_results(rows);
  unequal_nearest_neighbor(tree, new_leaf(&strings[0]), rows, hamming, &factor);
  require(rows.index == 1 && rows[0].index == 4);
  free_results(rows);
  free_tree(tree);
  free(words.elements);

  v_array<point> empty;
  tree = batch_create(empty, numeric_distance, &euclidean);
  k_nearest_neighbor(tree, new_leaf(&queries[0]), rows, 1, numeric_distance, &euclidean);
  require(rows.index == 0);
  free_tree(tree);
  free_tree(tree);
  free_results(rows);
  free_results(rows);
  std::cout << "native cover tree checks passed\n";
}
