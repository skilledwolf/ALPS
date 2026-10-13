// SPDX-License-Identifier: MIT
// Counts the random numbers a scheduler run draws while it is constructed.
// Every scheduler run seeds the same kind of engine from SEED, so a run's next
// number after construction is found in the sequence of a run with the same
// parameters that has drawn nothing; its position is the number of draws.
#ifndef ALPS_QMC_TESTS_CONSTRUCTION_DRAWS_H
#define ALPS_QMC_TESTS_CONSTRUCTION_DRAWS_H

#include <alps/scheduler/montecarlo.h>

#include <functional>
#include <vector>

namespace construction_draws_detail {

class reference_run : public alps::scheduler::MCRun {
public:
    explicit reference_run(const alps::Parameters& p)
        : alps::scheduler::MCRun(alps::ProcessList(1, alps::Process(0)), p, 0) {}
    void dostep() {}
    double work_done() const { return 0.; }
    double draw() { return random_01(); }
};

} // namespace construction_draws_detail

// The number of random numbers drawn before next() is called, or -1 if it is
// more than horizon. next() draws the run's next number from its engine.
inline int construction_draws(const alps::Parameters& p, const std::function<double()>& next,
                              int horizon = 256) {
    construction_draws_detail::reference_run reference(p);
    std::vector<double> sequence(horizon + 1);
    for (double& x : sequence)
        x = reference.draw();
    const double first = next();
    const double second = next();
    for (int d = 0; d < horizon; ++d)
        if (sequence[d] == first && sequence[d + 1] == second)
            return d;
    return -1;
}

#endif // ALPS_QMC_TESTS_CONSTRUCTION_DRAWS_H
