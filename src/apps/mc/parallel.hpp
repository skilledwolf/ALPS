// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/config.h>
#include <alps/hdf5/archive.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <sstream>
#include <climits>
#ifdef ALPS_HAVE_MPI
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#endif

namespace native_mc {
// Independent chains retain global IDs and seeds regardless of process count.
struct parallel {
#ifdef ALPS_HAVE_MPI
    boost::mpi::communicator world;
    int rank() const { return world.rank(); }
    int size() const { return world.size(); }
#else
    int rank() const { return 0; }
    int size() const { return 1; }
#endif
    bool owns(size_t id) const { return id%size()==size_t(rank()); }
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
    template<class Chains> void synchronize(Chains& chains) const {
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
                for (size_t id=0;id<chains.size();++id) if (owns(id)) {
                    ar.set_context("/clones/"+std::to_string(id));
                    chains[id]->save(ar);
                }
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
                for (size_t id=owner;id<chains.size();id+=size()) {
                    ar.set_context("/clones/"+std::to_string(id));
                    chains[id]->load(ar);
                }
            }
        });
#else
        (void)chains;
#endif
    }
};
}
