// SPDX-License-Identifier: MIT
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <stdexcept>
int main() {
    // The lattice-model applications and tutorials hand scalar
    // dictionaries to the legacy lattice and model libraries.
    // This internal bridge has a live caller; legacy file ingress is retired.
    alps::params p;
    p["LATTICE"]="chain lattice"; p["L"]=10; p["T"]=2.25;
    const auto legacy=alps::make_deprecated_parameters(p);
    if(std::string(legacy["LATTICE"])!="chain lattice" ||
       std::string(legacy["L"])!="10" || std::string(legacy["T"])!="2.25")
        throw std::runtime_error("typed parameters to model-library bridge failed");
}
