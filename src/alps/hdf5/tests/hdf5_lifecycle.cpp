// SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>

#include <hdf5.h>

#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, char const* message) {
    if (!condition) throw std::runtime_error(message);
}

template<class Operation> void rejects(Operation operation) {
    bool rejected = false;
    try { operation(); }
    catch (alps::hdf5::archive_error const&) { rejected = true; }
    require(rejected, "invalid archive operation was accepted");
}

int value(alps::hdf5::archive& ar, std::string const& path) {
    int result = 0;
    ar.read(path, result);
    return result;
}

void modes(std::string const& filename) {
    {
        alps::hdf5::archive ar(filename, "a");
        ar.write("/first", 1);
    }
    {
        alps::hdf5::archive ar(filename, "a");
        ar.write("/second", 2);
        require(value(ar, "/first") == 1, "append lost existing contents");
    }
    for (auto mode : {"", "rw", "aw", "rr", "aa", "ww", "rc", "al", "wm", "c", "l", "m", "x"})
        rejects([&] { alps::hdf5::archive ar(filename, mode); });
    {
        alps::hdf5::archive ar(filename);
        require(value(ar, "/first") == 1 && value(ar, "/second") == 2,
                "invalid mode modified the file");
    }
    {
        alps::hdf5::archive ar(filename, "w");
        require(ar.list_children("/").empty(), "write mode did not truncate at open");
        ar.write("/replacement", 3);
        ar.close();
        require(!ar.is_open(), "explicit close left the archive open");
        rejects([&] { ar.write("/closed", 4); });
    }
    {
        alps::hdf5::archive ar(filename);
        require(!ar.is_data("/first") && !ar.is_data("/second") && value(ar, "/replacement") == 3,
                "write mode preserved old datasets");
    }
}

void independent_opens(std::string const& filename) {
    {
        alps::hdf5::archive writer(filename, "a");
        alps::hdf5::archive reader(filename, "r");
        rejects([&] { reader.write("/forbidden", 4); });
        rejects([&] { reader.delete_data("/replacement"); });
        rejects([&] { reader.create_group("/forbidden"); });
        writer.close();
        require(value(reader, "/replacement") == 3, "closing another open invalidated a reader");
    }
    {
        alps::hdf5::archive reader(filename, "r");
        // HDF5 refuses an incompatible writable open while a read-only file
        // is open. ALPS must propagate that error instead of promoting it.
        rejects([&] { alps::hdf5::archive writer(filename, "a"); });
        rejects([&] { alps::hdf5::archive writer(filename, "w"); });
        rejects([&] { reader.write("/forbidden", 4); });
        require(value(reader, "/replacement") == 3, "failed open changed existing contents");
    }
}

void copies(std::string const& filename) {
    {
        alps::hdf5::archive first(filename, "a");
        {
            alps::hdf5::archive temporary(first);
            temporary.write("/copy", 5);
        }
        require(first.is_open() && value(first, "/copy") == 5, "destroying a view closed its owner");
        alps::hdf5::archive second(first);
        first.close();
        require(!first.is_open() && !second.is_open(), "explicit close left a view open");
        rejects([&] { second.write("/closed", 6); });
        second.close();
    }
    {
        alps::hdf5::archive first(filename, "r");
        alps::hdf5::archive second(first);
        require(value(second, "/copy") == 5, "read-only copy lost its file");
        rejects([&] { second.write("/forbidden", 4); });
        second.close();
        require(!first.is_open(), "closing a read-only view left another open");
    }
}

void checkpoint_copies(std::string const& filename) {
    std::unique_ptr<alps::hdf5::archive> escaped;
    alps::hdf5::save_checkpoint(filename, [&](auto& ar) {
        ar.write("/checkpoint", 6);
        escaped = std::make_unique<alps::hdf5::archive>(ar);
    });
    require(!escaped->is_open(), "escaped checkpoint copy remained open");
    rejects([&] { escaped->write("/late", 7); });
    escaped.reset();
    bool rejected = false;
    try {
        alps::hdf5::save_checkpoint(filename, [&](auto& ar) {
            ar.write("/checkpoint", 7);
            escaped = std::make_unique<alps::hdf5::archive>(ar);
            throw std::runtime_error("checkpoint failure");
        });
    } catch (std::runtime_error const& error) {
        require(std::string(error.what()) == "checkpoint failure", "checkpoint lost the original exception");
        rejected = true;
    }
    require(rejected && !escaped->is_open(), "failed checkpoint retained an open copy");
    rejects([&] { int result; escaped->read("/checkpoint", result); });
    escaped.reset();
    alps::hdf5::archive ar(filename);
    require(value(ar, "/checkpoint") == 6, "failed checkpoint modified its destination");
}
}

int main() {
    std::string const filename = "test_hdf5_lifecycle.h5";
    std::filesystem::remove(filename);
    try {
        auto const objects_before = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
        modes(filename);
        independent_opens(filename);
        copies(filename);
        checkpoint_copies(filename);
        alps::hdf5::update_archive(filename,[](auto& ar){ ar.write("/derived",8); });
        bool failed=false;
        try {
            alps::hdf5::update_archive(filename,[](auto& ar){
                ar.write("/checkpoint",99); throw std::runtime_error("evaluation failed");
            });
        } catch (std::runtime_error const& error) {
            failed=std::string(error.what())=="evaluation failed";
        }
        require(failed,"archive update lost the original exception");
        {
            alps::hdf5::archive ar(filename);
            require(value(ar,"/checkpoint")==6 && value(ar,"/derived")==8,
                    "archive update lost existing data or published a failed update");
        }
        require(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL) == objects_before, "file handles leaked");
        std::filesystem::remove(filename);
        return 0;
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        std::filesystem::remove(filename);
        return 1;
    }
}
