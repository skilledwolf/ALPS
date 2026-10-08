// SPDX-License-Identifier: MIT
#include <alps/parameter.h>
#include <stdexcept>
int main() {
    // The lattice and model libraries evaluate text symbols; typed
    // parameters reach them through this one conversion.
    alps::params p;
    p["LATTICE"]="chain lattice"; p["L"]=10; p["T"]=2.25;
    const alps::Parameters legacy(p);
    if(std::string(legacy["LATTICE"])!="chain lattice" ||
       std::string(legacy["L"])!="10" || std::string(legacy["T"])!="2.25")
        throw std::runtime_error("typed parameters to model-library bridge failed");
}
