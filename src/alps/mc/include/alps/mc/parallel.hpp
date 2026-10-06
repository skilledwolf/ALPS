// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/config.h>
#include <alps/hdf5/archive.hpp>
#include <alps/params.hpp>
#include <alps/run_config.hpp>
#include <array>
#include <set>
#include <boost/filesystem.hpp>
#include <fstream>
#include <memory>
#include <sstream>
#include <climits>
#ifdef ALPS_HAVE_MPI
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#endif

namespace alps::mc {
// Independent chains retain global IDs and seeds regardless of process count.
struct parallel {
#ifdef ALPS_HAVE_MPI
    // MPI thread support the program requires. A group whose chains run
    // threads beside MPI calls on the main thread raises this to funneled.
    static constexpr boost::mpi::threading::level threading=boost::mpi::threading::single;
    boost::mpi::communicator world;
    int rank() const { return world.rank(); }
    int size() const { return world.size(); }
#else
    int rank() const { return 0; }
    int size() const { return 1; }
#endif
    void configure(alps::run_configuration const&) {}
    bool owns(size_t id) const { return id%size()==size_t(rank()); }
    template<class Simulation>
    auto make(alps::params const& p,size_t bins,size_t id) const {
        return std::make_unique<Simulation>(p,bins,id);
    }
    // Independent chains stop locally between scheduled collective checks.
    // A spatial worker's group instead agrees before every collective sweep.
    bool stopped(bool local) const { return local; }
    template<class Runs> void verify(Runs const&,bool) const {}
    bool any(bool local) const {
#ifdef ALPS_HAVE_MPI
        return boost::mpi::all_reduce(world,local,std::logical_or<bool>());
#else
        return local;
#endif
    }
    template<class F> void checked(F&& operation) const {
        std::string message;
        try { operation(); }
        catch (std::exception const& error) { message=*error.what() ? error.what() : "Native MC failure"; }
        catch (...) { message="Unknown native MC failure"; }
#ifdef ALPS_HAVE_MPI
        int failed=boost::mpi::all_reduce(world,message.empty()?size():rank(),boost::mpi::minimum<int>());
        if (failed<size()) boost::mpi::broadcast(world,message,failed);
#endif
        if (!message.empty()) throw std::runtime_error(message);
    }
    template<class Save,class Load> void transport(Save const& save,Load const& load) const {
#ifdef ALPS_HAVE_MPI
        if (size()==1) return;
        // Transport the native checkpoint bytes, not a second state schema.
        // A private spool also works when ranks do not share a filesystem.
        struct spool {
            boost::filesystem::path directory;
            ~spool() { boost::system::error_code error; if (!directory.empty()) boost::filesystem::remove_all(directory,error); }
        } temporary;
        std::string bytes;
        checked([&] {
            auto directory=boost::filesystem::temp_directory_path()/boost::filesystem::unique_path("alps-mpi-%%%%-%%%%-%%%%");
            if (!boost::filesystem::create_directory(directory)) throw std::runtime_error("Cannot create MPI checkpoint spool");
            temporary.directory=std::move(directory);
            auto path=(temporary.directory/"chains.h5").string();
            {
                alps::hdf5::archive ar(path,"w");
                save(ar);
            }
            std::ifstream file(path,std::ios::binary);
            if (!file) throw std::runtime_error("Cannot read MPI checkpoint spool");
            bytes.assign(std::istreambuf_iterator<char>(file),{});
            if (bytes.size()>size_t(INT_MAX)/size()-1024) throw std::overflow_error("MPI checkpoint payload exceeds transport capacity");
        });
        std::vector<std::string> states;
        boost::mpi::gather(world,bytes,states,0);
        checked([&] {
            if (rank()!=0) return;
            auto path=(temporary.directory/"chains.h5").string();
            for (int owner=1;owner<size();++owner) {
                {
                    std::ofstream file(path,std::ios::binary|std::ios::trunc);
                    file.exceptions(std::ios::badbit|std::ios::failbit);
                    file.write(states[owner].data(),states[owner].size());
                    file.close();
                }
                alps::hdf5::archive ar(path);
                load(ar,owner);
            }
        });
#else
        (void)save;(void)load;
#endif
    }
    template<class Chains> void synchronize(Chains& chains) const {
        transport([&](auto& ar) {
            for (size_t id=0;id<chains.size();++id) if (owns(id)) {
                ar.set_context("/clones/"+std::to_string(id));chains[id]->save(ar);
            }
        },[&](auto& ar,int owner) {
            for (size_t id=owner;id<chains.size();id+=size()) {
                ar.set_context("/clones/"+std::to_string(id));chains[id]->load(ar);
            }
        });
    }
};

// Collective workers must agree before entering any physical communication.
// Construction and checkpoint validation remain local until this preflight.
struct collective : parallel {
    bool owns(size_t) const {return true;}
    bool stopped(bool local) const {return any(local);}
    template<class Runs> void verify(Runs const& runs,bool validate=false) const {
#ifdef ALPS_HAVE_MPI
        if (size()==1) return;
        std::vector<std::string> local,expected;
        checked([&] {
            local.push_back(validate ? "validate" : "execute");
            for (auto const& run:runs) local.push_back(alps::format_run_configuration(run));
            if (rank()==0) expected=local;
        });
        boost::mpi::broadcast(world,expected,0);
        checked([&] {
            if (local!=expected) throw std::invalid_argument("Collective ranks require identical run configurations");
        });
        // Matching paths alone do not guarantee matching rank-local input
        // files. Compare bytes once before any physical collective, keeping
        // the existing schemas and bounded memory even for large checkpoints.
        std::set<std::string> inputs;
        checked([&] {for (auto const& run:runs) for (auto const& [key,value]:run.input)
            for (auto const& path:alps::run_paths(value)) inputs.insert(path);});
        for (auto const& path:inputs) {
            std::ifstream file(path,std::ios::binary);
            checked([&] {if (!file) throw std::runtime_error("Cannot read collective input: "+path);});
            std::array<char,65536> actual{},reference{};
            bool different=false;
            for (;;) {
                int count=0;
                if (rank()==0) {file.read(reference.data(),reference.size());count=int(file.gcount());}
                boost::mpi::broadcast(world,count,0);
                if (!count) break;
                boost::mpi::broadcast(world,reference.data(),count,0);
                if (rank()!=0) {
                    file.read(actual.data(),count);
                    different=different || file.gcount()!=count || !std::equal(actual.begin(),actual.begin()+count,reference.begin());
                }
            }
            different=different || file.bad() || file.peek()!=std::char_traits<char>::eof();
            checked([&] {if (different) throw std::invalid_argument("Collective ranks require identical input files: "+path);});
        }
#else
        (void)runs;(void)validate;
#endif
    }
};
}
