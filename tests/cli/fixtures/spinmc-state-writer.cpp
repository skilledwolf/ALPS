// SPDX-License-Identifier: MIT
// Protocol fixture producer, not a complete released spinmc application.
// Field order comes from ALPS v3.0.0 Worker::save_worker,
// AbstractSpinSim::save, SpinSim::save and the moment ODump operators.
#include <alps/osiris/xdrdump.h>
#include <cstdint>
#include <string>

int main(int argc, char** argv) {
    if (argc != 2) return 1;
    for (std::string model : {"Ising", "Potts", "XY", "Heisenberg", "O4"}) {
        alps::OXDRFileDump dump(std::string(argv[1]) + "/spinmc-" + model + ".xdr");
        dump << int32_t(3) << int32_t(0) << int32_t(400);
        dump << uint64_t(0x100000005) << double(1.125) << uint64_t(0x100000003);
        dump << uint32_t(2); // std::vector<Moment> length
        if (model == "Ising") dump << false << true;
        else if (model == "Potts") dump << uint32_t(0) << uint32_t(2);
        else {
            int dim = model == "XY" ? 2 : model == "Heisenberg" ? 3 : 4;
            // TinyVector stores components directly, without an inner length.
            for (int site = 0; site < 2; ++site)
                for (int i = 0; i < dim; ++i)
                    dump << double(i == site ? (site ? -1 : 1) : 0);
        }
    }
}
