#include <libcoherent4cast/rng.hpp>
#include <array>
#include <iostream>
#include <stdexcept>
int main() {
  libcoherent4cast::Pcg32 rng(42U, 54U);
  const std::array<std::uint32_t, 6> expected = {
    0xa15c02b7U, 0x7b47f409U, 0xba1d3330U,
    0x83d2f293U, 0xbfa4784bU, 0xcbed606eU};
  for (const auto value : expected) {
    if (rng.next() != value) throw std::runtime_error("PCG32 reference stream mismatch");
  }
  libcoherent4cast::Pcg32 stream_a(123U, 1U), stream_b(123U, 2U);
  if (stream_a.next() == stream_b.next()) throw std::runtime_error("task streams coincide at first draw");
  for (int j = 0; j < 10000; ++j) {
    const double u = rng.uniform();
    if (!(u > 0.0 && u < 1.0)) throw std::runtime_error("uniform outside open unit interval");
    if (!(rng.exponential(2.0) > 0.0)) throw std::runtime_error("invalid exponential draw");
  }
  bool caught = false;
  try { rng.exponential(0.0); } catch (const std::invalid_argument&) { caught = true; }
  if (!caught) throw std::runtime_error("zero rate not rejected");
  std::cout << "PCG32 reference stream and numeric guards passed\n";
  return 0;
}
