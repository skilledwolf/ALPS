// SPDX-License-Identifier: MIT
#include <alps/alea.h>
#include <array>
#include <cmath>
#include <stdexcept>
#include <valarray>
#include <vector>

namespace {
void close(double actual, double expected) {
  if (!std::isfinite(actual) || std::abs(actual - expected) > 1e-12)
    throw std::runtime_error("binning differs from direct sample statistics");
}
}

int main() {
  alps::SimpleBinning<std::valarray<double>> vector;
  std::array<alps::SimpleBinning<double>, 3> scalars;
  std::vector<std::array<double, 3>> samples;
  // Check partial bins as well as every power-of-two growth boundary.
  for (int n = 0; n < 513; ++n) {
    std::array<double, 3> sample;
    for (int j = 0; j < 3; ++j) {
      sample[j] = ((n*n*(j+1) + n*(j+3)) % 101 - 50) / 8.0;
      scalars[j] << sample[j];
    }
    samples.push_back(sample);
    vector << std::valarray<double>(sample.data(), sample.size());
    if (vector.count() != samples.size()) throw std::runtime_error("incorrect sample count");
    for (int j = 0; j < 3; ++j) {
      for (std::size_t width = 1, level = 0; width <= samples.size(); width *= 2, ++level) {
        const auto count = samples.size() / width;
        double sum = 0, squares = 0;
        for (std::size_t bin = 0; bin < count; ++bin) {
          double mean = 0;
          for (std::size_t k = 0; k < width; ++k) mean += samples[bin*width+k][j];
          mean /= width;
          sum += mean;
          squares += mean*mean;
        }
        const double mean = sum / count;
        const double variance = squares / count - mean*mean;
        close(vector.binmean_element(j, level), mean);
        close(vector.binvariance_element(j, level), variance);
        if (level == 0) {
          close(vector.mean()[j], mean);
          close(scalars[j].mean(), mean);
          if (count > 1) {
            close(vector.variance()[j], variance*count/(count-1));
            close(scalars[j].variance(), variance*count/(count-1));
          }
        }
      }
    }
  }
  bool rejected = false;
  try { vector << std::valarray<double>(4); }
  catch (const std::runtime_error&) { rejected = true; }
  if (!rejected || vector.count() != samples.size())
    throw std::runtime_error("size mismatch must preserve accumulated samples");
  vector.reset();
  if (vector.count() != 0 || vector.size() != 0) throw std::runtime_error("reset failed");
  vector << std::valarray<double>{1, 2};
  if (vector.size() != 2 || vector.count() != 1) throw std::runtime_error("resize after reset failed");
  close(vector.mean()[0], 1);
  close(vector.mean()[1], 2);
}
