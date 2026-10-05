// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/fortran/fortran_wrapper.h>
#include <alps/hdf5/archive.hpp>
#include <boost/filesystem.hpp>
#include <iostream>

int main() {
    const auto path = boost::filesystem::temp_directory_path()/boost::filesystem::unique_path("fortran-%%%%-%%%%.h5");
    try {
        alps::params p;
        p["L"]=2; p["TEMPERATURE"]=3.; p["SWEEPS"]=100; p["THERMALIZATION"]=3;
        p["SEED"]=42; p["RNG"]="mt19937";
        alps::fortran_wrapper live(p,8,0), reference(p,8,0);
        for (int i=0;i<19;++i) { live.update(); reference.update(); }
        static_cast<alps::mcbase const&>(live).save(path);
        {
            alps::hdf5::archive ar(path,"a");
            ar["/simulation/realizations/0/clones/0/checkpoint/fortran/fields/1/value"] << std::vector<int>(4,0);
        }
        bool rejected=false;
        try { static_cast<alps::mcbase&>(live).load(path); }
        catch (std::exception const&) { rejected=true; }
        if (!rejected) throw std::runtime_error("Invalid physical state accepted");
        // This also exercises destruction of the rejected staged Fortran state.
        for (int i=0;i<37;++i) { live.update(); reference.update(); }
        if (live.completed_sweeps()!=reference.completed_sweeps() || live.get_random()()!=reference.get_random()())
            throw std::runtime_error("Rejected checkpoint changed progress or RNG");
        for (auto const& name : {"Energy", "Magnetization"})
            if (!(live.measurement(name)->result()==reference.measurement(name)->result()))
                throw std::runtime_error("Rejected checkpoint changed physical state or measurements");
        boost::filesystem::remove(path);
        return 0;
    } catch (std::exception const& e) {
        boost::filesystem::remove(path);
        std::cerr << e.what() << '\n'; return 1;
    }
}
