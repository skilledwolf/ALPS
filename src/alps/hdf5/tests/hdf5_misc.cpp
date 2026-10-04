/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/hdf5.hpp>

#include <boost/filesystem.hpp>

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

static_assert(!std::is_copy_assignable<alps::hdf5::archive>::value,
              "archive assignment would bypass registry ownership");

class my_class {
    public:
        my_class(double v = 0): d(v) {}
        void save(alps::hdf5::archive & ar) const {
            using alps::make_pvp;
            ar << make_pvp("value", d);
        }
        void load(alps::hdf5::archive & ar) { 
            using alps::make_pvp;
            ar >> make_pvp("value", d); 
        }
    private:
        double d;
};

namespace {

    struct context_failure : std::runtime_error {
        context_failure() : std::runtime_error("save/load hook failed") {}
    };

    void require_context(alps::hdf5::archive const & ar, std::string const & expected) {
        if (ar.get_context() != expected)
            throw std::runtime_error("archive context was not restored: " + ar.get_context());
    }

    template<typename Error, typename Operation>
    void expect_failure(alps::hdf5::archive & ar, std::string const & context, Operation operation) {
        try {
            operation();
        } catch (Error const &) {
            require_context(ar, context);
            return;
        }
        throw std::runtime_error("save/load operation did not throw");
    }

    struct context_value {
        int value = 0;
        bool fail = false;

        void save(alps::hdf5::archive & ar) const {
            ar.write("value", value);
            if (fail) {
                ar.set_context("changed/by/hook");
                throw context_failure();
            }
        }

        void load(alps::hdf5::archive & ar) {
            ar.read("value", value);
            if (fail) {
                ar.set_context("changed/by/hook");
                throw context_failure();
            }
        }
    };

    struct recovering_parent {
        context_value child{17, true};
        int after = 19;

        void save(alps::hdf5::archive & ar) const {
            std::string const context = ar.get_context();
            try {
                ar << alps::make_pvp("child", child);
            } catch (context_failure const &) {}
            require_context(ar, context);
            ar.write("after", after);
        }

        void load(alps::hdf5::archive & ar) {
            std::string const context = ar.get_context();
            try {
                ar >> alps::make_pvp("child", child);
            } catch (context_failure const &) {}
            require_context(ar, context);
            ar.read("after", after);
        }
    };

    void test_context_restoration() {
        using alps::make_pvp;
        std::string const context = "/context/base";
        context_value const ordinary{23, false};
        std::vector<context_value> const sequence{{29, false}, {31, false}};
        std::map<std::string, context_value> const mapping{{"first", {37, false}}, {"second", {41, false}}};

        {
            alps::hdf5::archive ar("data.h5", "w");
            ar.set_context(context);
            ar << make_pvp("ordinary", ordinary);
            require_context(ar, context);
            ar << make_pvp("/context/absolute", ordinary);
            require_context(ar, context);
            ar << make_pvp("sequence", sequence) << make_pvp("mapping", mapping);
            require_context(ar, context);

            context_value const failing{43, true};
            expect_failure<context_failure>(ar, context, [&] { ar << make_pvp("failed", failing); });
            expect_failure<context_failure>(ar, context, [&] { ar << make_pvp("/context/failed_absolute", failing); });
            std::vector<context_value> const failed_sequence{{47, false}, failing};
            expect_failure<context_failure>(ar, context, [&] { ar << make_pvp("failed_sequence", failed_sequence); });
            std::map<std::string, context_value> const failed_mapping{{"first", {53, false}}, {"second", failing}};
            expect_failure<context_failure>(ar, context, [&] { ar << make_pvp("failed_mapping", failed_mapping); });

            recovering_parent parent;
            ar << make_pvp("parent", parent);
            require_context(ar, context);
            ar.write("after_failures", 59);
            if (!ar.is_data("/context/base/after_failures") || !ar.is_data("/context/base/parent/after"))
                throw std::runtime_error("subsequent writes used the failed hook's context");
        }

        {
            alps::hdf5::archive ar("data.h5", "r");
            ar.set_context(context);
            my_class missing;
            expect_failure<alps::hdf5::path_not_found>(ar, context, [&] { ar >> make_pvp("missing", missing); });
            expect_failure<alps::hdf5::archive_error>(ar, context, [&] { ar << make_pvp("read_only", ordinary); });

            context_value value{0, true};
            expect_failure<context_failure>(ar, context, [&] { ar >> make_pvp("ordinary", value); });
            expect_failure<context_failure>(ar, context, [&] { ar >> make_pvp("/context/absolute", value); });
            std::vector<context_value> loaded_sequence{{0, false}, {0, true}};
            expect_failure<context_failure>(ar, context, [&] { ar >> make_pvp("sequence", loaded_sequence); });
            std::map<std::string, context_value> loaded_mapping{{"second", {0, true}}};
            expect_failure<context_failure>(ar, context, [&] { ar >> make_pvp("mapping", loaded_mapping); });

            recovering_parent parent;
            parent.after = 0;
            ar >> make_pvp("parent", parent);
            require_context(ar, context);
            if (parent.child.value != 17 || parent.after != 19)
                throw std::runtime_error("nested hook recovery changed stored values");
            value.fail = false;
            ar >> make_pvp("ordinary", value);
            require_context(ar, context);
            if (value.value != ordinary.value)
                throw std::runtime_error("subsequent read used the failed hook's context");
            loaded_sequence.clear();
            loaded_mapping.clear();
            ar >> make_pvp("sequence", loaded_sequence) >> make_pvp("mapping", loaded_mapping);
            require_context(ar, context);
            if (loaded_sequence.size() != 2 || loaded_sequence[0].value != 29 || loaded_sequence[1].value != 31
                || loaded_mapping.size() != 2 || loaded_mapping.at("first").value != 37 || loaded_mapping.at("second").value != 41)
                throw std::runtime_error("container object roundtrip changed stored values");
        }
    }

}

int main () {

    if (boost::filesystem::exists(boost::filesystem::path("data.h5")))
        boost::filesystem::remove(boost::filesystem::path("data.h5"));

    {
        alps::hdf5::archive ar("data.h5", "w");
        ar << alps::make_pvp("/value", 42);
    }
    
    {
        alps::hdf5::archive ar("data.h5");
        int i;
        ar >> alps::make_pvp("/value", i);
    }

    {
        alps::hdf5::archive ar("data.h5");
        std::string s;
        expect_failure<alps::hdf5::wrong_type>(ar, "/", [&] { ar >> alps::make_pvp("/value", s); });
    }

    {
        alps::hdf5::archive ar("data.h5", "w");
        std::vector<double> vec(5, 42);
        ar << alps::make_pvp("/path/2/vec", vec);
    }
    
    {
        std::vector<double> vec;
        // fill the vector
        alps::hdf5::archive ar("data.h5");
        ar >> alps::make_pvp("/path/2/vec", vec);
    }

    {
        std::string str("foobar");
        alps::hdf5::archive ar("data.h5", "w");
        ar << alps::make_pvp("/foo/bar", str);
    }
    
    {
        alps::hdf5::archive ar("data.h5");
        std::string str;
        ar >> alps::make_pvp("/foo/bar", str);
    }

    {
        long *d = new long[17]{};
        alps::hdf5::archive ar("data.h5", "w");
        ar << alps::make_pvp("/c/array", d, 17);
        delete[] d;
    }

    {
        alps::hdf5::archive ar("data.h5");
        std::size_t size = ar.extent("/c/array")[0];
        long *d = new long[size];
        ar >> alps::make_pvp("/c/array", d, size);
        delete[] d;
    }

    {
        {
                my_class c(42);
                alps::hdf5::archive ar("data.h5", "w");
                ar << alps::make_pvp("/my/class", c);
        }
        {
                my_class c;
                alps::hdf5::archive ar("data.h5");
                ar >> alps::make_pvp("/my/class", c);
        }
    }

    {
        alps::hdf5::archive ar("data.h5", "w"); 
        // the parent of an attribute must exist
        ar.create_group("/foo");
        ar << alps::make_pvp("/foo/@bar", std::string("hello"));
    }

    {
        alps::hdf5::archive ar("data.h5");
        std::string str;
        ar >> alps::make_pvp("/foo/@bar", str);
    }

    try {
        test_context_restoration();
    } catch (std::exception const & error) {
        std::cerr << error.what() << '\n';
        return 1;
    }

    boost::filesystem::remove(boost::filesystem::path("data.h5"));
    return 0;
}
