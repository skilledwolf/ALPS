// Derived from ALPSCore params at 7146b9e1f017938a94e5dae35d88467cc5ba7969.
// Copyright (C) 1998-2018 ALPS Collaboration; modifications (C) 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
// See ALPSCore-LICENSE.txt for the original permission notice.
#include <algorithm>
#include <alps/dictionary.hpp>
#ifdef ALPS_HAVE_MPI
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/serialization/complex.hpp>
#endif
namespace alps::params_ns {
dictionary::map_type::const_iterator dictionary::find_nonempty_(const std::string &key) const {
    auto it = map_.find(key);
    return it != map_.end() && !it->second.empty() ? it : map_.end();
}
dictionary::value_type &dictionary::operator[](const std::string &key) {
    auto it = map_.lower_bound(key);
    if (it == map_.end() || map_.key_comp()(key, it->first))
        it = map_.insert(it, map_type::value_type(key, value_type(key)));
    return it->second;
}
const dictionary::value_type &dictionary::operator[](const std::string &key) const {
    auto it = map_.find(key);
    if (it == map_.end())
        throw exception::uninitialized_value(key, "value is not set");
    return it->second;
}
bool dictionary::equals(const dictionary &rhs) const {
    return size() == rhs.size() &&
           std::equal(begin(), end(), rhs.begin(), [](const auto &a, const auto &b) {
               return a.first == b.first && a.second.equals(b.second);
           });
}
std::ostream &operator<<(std::ostream &out, const dictionary &values) {
    for (const auto &entry : values)
        out << entry.first << " = " << entry.second << '\n';
    return out;
}
std::ostream &operator<<(std::ostream &out, const dict_value &value) {
    value.apply_visitor([&](const auto &x) {
        using T = std::decay_t<decltype(x)>;
        if constexpr (std::is_same_v<T, dict_value::None>)
            out << "[NONE]";
        else if constexpr (detail::vector_type<T>::value) {
            out << '[';
            for (std::size_t i = 0; i < x.size(); ++i) {
                if (i)
                    out << ", ";
                out << x[i];
            }
            out << ']';
        } else
            out << x;
    });
    return out;
}
#ifdef ALPS_HAVE_MPI
void dictionary::broadcast(const boost::mpi::communicator &comm, int root) {
    boost::mpi::broadcast(comm, map_, root);
}
#endif
} // namespace alps::params_ns
