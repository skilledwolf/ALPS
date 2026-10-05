#ifdef _WIN32
#include <windows.h>
#endif
#include <boost/spirit/include/classic_core.hpp>
#include <alps/parapack/integer_range.h>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <alps/utility/encode.hpp>
#include <alps/mcbase.hpp>
#include <alps/parser/xslt_path.h>
#include <alps/osiris/xdrdump.h>
#include <alps/numeric/functional.hpp>
#include <boost/filesystem/operations.hpp>
#include <vector>

class simulation : public alps::mcbase {
public:
    simulation() : alps::mcbase(alps::params{}) {
        measurements.emplace("samples", std::make_shared<alps::alea::batch_acc<double>>(1, 8));
    }
    void update() override { ++steps; }
    void measure() override { *measurements.at("samples") << double(steps); }
    double fraction_completed() const override { return double(steps)/3; }
private:
    int steps = 0;
};

int main(int argc, char** argv) {
    if (argc == 2) {
        for (const char* resource : {"lattices.xml", "models.xml", "ALPS.xsl"}) {
            if (!boost::filesystem::equivalent(alps::search_xml_library_path(resource),
                                               boost::filesystem::path(argv[1]) / resource))
                return 1;
        }
        return 0;
    }
    const std::string archive_name = "a/path with spaces";
    if (alps::hdf5_name_decode(alps::hdf5_name_encode(archive_name)) != archive_name)
        return 1;
    alps::params parameters;
    parameters["count"] = 3;
    const alps::integer_range<int> range("[2:5]");
    const std::vector<double> expected{1.0, 2.0, 3.0};
    {
        alps::OXDRFileDump output("sdk-contract.xdr");
        output << 42;
    }
    int checkpoint_value = 0;
    {
        alps::IXDRFileDump input("sdk-contract.xdr");
        input >> checkpoint_value;
    }
    {
        alps::hdf5::archive output("sdk-contract.h5", "w");
        output["/values"] << expected;
    }
    std::vector<double> actual;
    alps::hdf5::archive input("sdk-contract.h5", "r");
    input["/values"] >> actual;
    simulation sim;
    sim.run([] { return false; });
    const auto result = sim.collect_results().at("samples");
    return actual == expected && checkpoint_value == 42 && int(parameters["count"]) == 3
        && (range.min)() == 2 && (range.max)() == 5
        && result.count() == 3 && result.mean()(0) == 2.0 ? 0 : 1;
}
