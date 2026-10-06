#ifdef _WIN32
#include <windows.h>
#endif
#include <boost/spirit/include/classic_core.hpp>
#include <alps/parapack/integer_range.h>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <alps/utility/encode.hpp>
#include <alps/mc/driver.hpp>
#include <alps/mc/replica_exchange.hpp>
#include <alps/parser/xslt_path.h>
#include <alps/osiris/xdrdump.h>
#include <alps/numeric/functional.hpp>
#include <boost/filesystem/operations.hpp>
#include <vector>

class simulation : public alps::mcbase {
public:
    simulation(alps::params const& p={},size_t bins=8,size_t chain=0) : alps::mcbase(p,chain) {
        alps::mc::add_measurement(*this,"samples",1,bins);
    }
    static alps::params checkpoint_parameters(alps::params p) { return p; }
    void update() override { ++steps; }
    void measure() override { alps::mc::record(*this,"samples",double(steps)); }
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
    // The installed driver must compile without exposing toml++ or private
    // application headers, and its dynamic schema helper must link transitively.
    const auto schema=alps::mc::schema({},"[parameters.SWEEPS]\ntype='int64'\n",{});
    if (schema.find("SWEEPS")==std::string::npos) return 1;
    // Generic MC clients (including Fortran) need no XML graph/model library.
    alps::run_configuration run;
    run.execution["chains"]=2; run.execution["seed"]=42;
    run.execution["bins"]=8; run.execution["rng"]="mt19937";
    auto chains=alps::mc::prepare_chains<simulation>(run,
        [](alps::params&,alps::run_configuration const&){},alps::mc::parallel{});
    if (chains.size()!=2) return 1;
    alps::params replicas;
    replicas["INVERSE_TEMPERATURE_SET"]=std::vector<double>{.5,1.};
    replicas["SWEEPS"]=3; replicas["THERMALIZATION"]=0;
    alps::mc::replica_exchange<double> exchange(replicas,0,0.);
    exchange.step([](auto const&,auto const&,bool){},[]{return std::vector<double>{-1.,1.};},
                  [](double energy,double beta){return -beta*energy;},[](size_t,char const*,double){});
    if (exchange.production_sweeps()!=1) return 1;
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
    const auto result = std::get<alps::alea::batch_result<double>>(sim.collect_results().at("samples"));
    return actual == expected && checkpoint_value == 42 && int(parameters["count"]) == 3
        && (range.min)() == 2 && (range.max)() == 5
        && result.count() == 3 && result.mean()(0) == 2.0 ? 0 : 1;
}
