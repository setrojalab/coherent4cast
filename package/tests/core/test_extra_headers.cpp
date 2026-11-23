#include <librecurcast/fibre.hpp>
#include <librecurcast/continuous_target.hpp>
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
  // Integrating two ordered event times over a common window of length 0.6
  // gives triangular volume 0.6^2/2. Repeated records of one mark do not
  // create a second labelled history. Distinct recurrent marks permit both.
  const std::vector<double> nodes{0.5}, weights{1.0};
  const auto repeated = recurcast::coarsening_fibre(
      {{0.2, 0.8, 1, false}, {0.2, 0.8, 1, false}}, 1.0, {}, nodes, weights);
  const auto distinct = recurcast::coarsening_fibre(
      {{0.2, 0.8, 1, false}, {0.2, 0.8, 2, false}}, 1.0, {}, nodes, weights);
  const auto stopped = recurcast::coarsening_fibre(
      {{0.2, 0.8, 1, false}, {0.2, 0.8, 0, true}}, 1.0, {}, nodes, weights);
  assert(std::abs(repeated.volume - 0.18) < 1e-14);
  assert(std::abs(distinct.volume - 0.36) < 1e-14);
  assert(std::abs(stopped.volume - 0.18) < 1e-14);
  // A singular information matrix identifies only the common direction.
  // Its Moore-Penrose solution is (1,1), not an arbitrary null-space shift.
  const auto common = recurcast::ct_psolve({1.0, 1.0, 1.0, 1.0}, {2.0, 2.0});
  assert(std::abs(common[0] - 1.0) < 1e-12);
  assert(std::abs(common[1] - 1.0) < 1e-12);
  const auto boundary = recurcast::ct_psolve({2.0, 0.0, 0.0, 0.0}, {6.0, 0.0});
  assert(std::abs(boundary[0] - 3.0) < 1e-12);
  assert(boundary[1] == 0.0);
  std::cout << "Fibre-volume and singular-information checks passed\n";
}
