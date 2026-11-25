#ifndef LIBcoherent4cast_RNG_HPP
#define LIBcoherent4cast_RNG_HPP
#include <cstdint>
#include <cmath>
#include <stdexcept>

namespace libcoherent4cast {
// PCG XSH-RR, 64-bit state and 32-bit output. The seed and stream are explicit.
class Pcg32 {
 private:
  std::uint64_t state_ = 0U;
  std::uint64_t increment_ = 1U;
 public:
  Pcg32(std::uint64_t seed, std::uint64_t stream) {
    increment_ = (stream << 1U) | 1U;
    next(); state_ += seed; next();
  }
  std::uint32_t next() noexcept {
    const std::uint64_t previous = state_;
    state_ = previous * UINT64_C(6364136223846793005) + increment_;
    const auto shifted = static_cast<std::uint32_t>(((previous >> 18U) ^ previous) >> 27U);
    const auto rotation = static_cast<std::uint32_t>(previous >> 59U);
    return (shifted >> rotation) | (shifted << ((0U - rotation) & 31U));
  }
  double uniform() noexcept {
    return (static_cast<double>(next()) + 0.5) / 4294967296.0;
  }
  double exponential(double rate) {
    if (!(rate > 0.0) || !std::isfinite(rate))
      throw std::invalid_argument("rate must be positive and finite");
    return -std::log(uniform()) / rate;
  }
};
}  // namespace libcoherent4cast
#endif
