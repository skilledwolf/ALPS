// SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <boost/filesystem/operations.hpp>
#include <array>
#include <complex>
#include <csignal>
#include <stdexcept>
#include <vector>

#if !defined(BOOST_MSVC) && !defined(ALPS_NGS_NO_SIGNALS)
namespace {
void caller_signal_handler(int) {}
}
#endif

int main() {
    const char* filename = "hdf5-component-contract.h5";
    const std::vector<std::complex<double>> expected{{1., 2.}, {3., -4.}};
#if !defined(BOOST_MSVC) && !defined(ALPS_NGS_NO_SIGNALS)
    const std::array<int, 2> signals{{SIGSEGV, SIGBUS}};
    std::array<struct sigaction, 2> previous{};
    struct sigaction action{};
    action.sa_handler = caller_signal_handler;
    sigemptyset(&action.sa_mask);
    action.sa_flags = SA_RESTART;
    for (std::size_t i = 0; i < signals.size(); ++i)
        if (sigaction(signals[i], &action, &previous[i]) != 0)
            throw std::runtime_error("Could not install caller signal handler");
    auto check_signal_handlers = [&] {
        for (int number : signals) {
            struct sigaction observed{};
            if (sigaction(number, nullptr, &observed) != 0
                || observed.sa_handler != caller_signal_handler
                || !(observed.sa_flags & SA_RESTART))
                throw std::runtime_error("HDF5 changed the caller's signal handler");
        }
    };
    // Failed opens must not install process handlers either.
    const char* missing_file = "hdf5-missing-signal-contract.h5";
    boost::filesystem::remove(missing_file);
    try {
        alps::hdf5::archive missing(missing_file, "r");
        throw std::runtime_error("Missing archive unexpectedly opened");
    } catch (const alps::hdf5::archive_not_found&) {
    }
    check_signal_handlers();
#endif
    {
        alps::hdf5::archive writer(filename, "w");
#if !defined(BOOST_MSVC) && !defined(ALPS_NGS_NO_SIGNALS)
        check_signal_handlers();
#endif
        writer["/values"] << expected;
        writer["/values/@description"] << std::string("complex vector");
        // Multiple handles must share one archive context, including typed I/O.
        alps::hdf5::archive reader(filename, "r");
        std::vector<std::complex<double>> actual;
        std::string description;
        reader["/values"] >> actual;
        reader["/values/@description"] >> description;
        if (actual != expected || description != "complex vector")
            throw std::runtime_error("HDF5 component round trip failed");
        bool caught = false;
        try {
            int missing = 0;
            reader["/missing"] >> missing;
        } catch (const alps::hdf5::path_not_found&) {
            caught = true;
        }
        if (!caught)
            throw std::runtime_error("HDF5 component exception identity failed");
    }
#if !defined(BOOST_MSVC) && !defined(ALPS_NGS_NO_SIGNALS)
    check_signal_handlers();
    for (std::size_t i = 0; i < signals.size(); ++i)
        if (sigaction(signals[i], &previous[i], nullptr) != 0)
            throw std::runtime_error("Could not restore caller signal handler");
#endif
    boost::filesystem::remove(filename);
}
