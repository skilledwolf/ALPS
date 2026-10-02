// Derived from ALPSCore params at 7146b9e1f017938a94e5dae35d88467cc5ba7969.
// Copyright (C) 1998-2018 ALPS Collaboration; modifications (C) 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
// See ALPSCore-LICENSE.txt in this module for the original permission notice.

#ifndef ALPS_PARAMS_DICTIONARY_HPP_e15039548f43464996cad06f9c8a3220
#define ALPS_PARAMS_DICTIONARY_HPP_e15039548f43464996cad06f9c8a3220

#include <alps/config.h>
#include <alps/params_export.h>
#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>
#include <map>
#ifdef ALPS_HAVE_MPI
#include <boost/mpi/communicator.hpp>
#endif
#include "./params/dict_value.hpp"

namespace alps {
namespace params_ns {

/// Python-like dictionary
class ALPS_PARAMS_DECL dictionary {
  public:
    typedef dict_value value_type;

  private:
    typedef std::map<std::string, value_type> map_type;
    map_type map_;
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive &ar, unsigned) { ar & map_; }

  public:
    typedef map_type::const_iterator const_iterator;

    /// Const-iterator to the beginning of the contained map
    const_iterator begin() const { return map_.begin(); }

    /// Const-iterator to the end of the contained map
    const_iterator end() const { return map_.end(); }

    /// Virtual destructor to make dictionary inheritable
    virtual ~dictionary() {}

    /// True if the cdictionary does not contain elements (even empty ones)
    bool empty() const { return map_.empty(); }

    /// Size of the dictionary (including empty elements)
    std::size_t size() const { return map_.size(); }

    /// Erase an element if it exists
    void erase(const std::string &key) { map_.erase(key); }

    /// Access with intent to assign
    value_type &operator[](const std::string &key);

    /// Read-only access
    const value_type &operator[](const std::string &key) const;

    /// Obtain read-only iterator to a name
    const_iterator find(const std::string &key) const { return map_.find(key); }

  private:
    /// Check if the key exists and has a value; return the iterator
    map_type::const_iterator find_nonempty_(const std::string &key) const;

  public:
    /// Check if a key exists and has a value (without creating the key)
    bool exists(const std::string &key) const { return find_nonempty_(key) != map_.end(); }

    /// Check if a key exists and has a value of a particular type (without creating the key)
    template <typename T> bool exists(const std::string &key) const {
        map_type::const_iterator it = find_nonempty_(key);
        return it != map_.end() && (it->second).isType<T>();
    }

    /// Swap the dictionaries
    friend void swap(dictionary &d1, dictionary &d2) {
        using std::swap;
        swap(d1.map_, d2.map_);
    }

    /// Compare two dictionaries (true if all entries are of the same type and value)
    bool equals(const dictionary &rhs) const;

    /// Save the dictionary to an archive
    void save(alps::hdf5::archive &ar) const;

    /// Load the dictionary from an archive
    void load(alps::hdf5::archive &ar);

    template <class T> T value_or(const std::string &key, const T &fallback) const {
        auto it = find_nonempty_(key);
        return it == end() ? fallback : it->second.template as<T>();
    }
    std::string value_or(const std::string &key, const char *fallback) const {
        return value_or(key, std::string(fallback));
    }

    friend ALPS_PARAMS_DECL std::ostream &operator<<(std::ostream &, const dictionary &);

#ifdef ALPS_HAVE_MPI
    /// Broadcast the dictionary
    void broadcast(const boost::mpi::communicator &comm, int root);
#endif
};

inline bool operator==(const dictionary &lhs, const dictionary &rhs) { return lhs.equals(rhs); }

inline bool operator!=(const dictionary &lhs, const dictionary &rhs) { return !(lhs == rhs); }

/// Const-access visitor to a value by an iterator
/** @param visitor A functor that should be callable as `R result=visitor(bound_value_const_ref)`
    @param it Iterator to the dictionary
    @tparam F The functor type; must define typename `F::result_type`.
*/
template <typename F>
inline typename F::result_type apply_visitor(F &visitor, dictionary::const_iterator it) {
    return it->second.apply_visitor(visitor);
}

/// Const-access visitor to a value by an iterator
/** @param visitor A functor that should be callable as `R result=visitor(bound_value_const_ref)`
    @param it Iterator to the dictionary
    @tparam F The functor type; must define typename `F::result_type`.
*/
template <typename F>
inline typename F::result_type apply_visitor(const F &visitor, dictionary::const_iterator it) {
    return it->second.apply_visitor(visitor);
}

} // namespace params_ns
} // namespace alps

#endif /* ALPS_PARAMS_DICTIONARY_HPP_e15039548f43464996cad06f9c8a3220 */
