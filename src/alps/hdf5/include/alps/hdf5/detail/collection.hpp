// SPDX-License-Identifier: MIT
// Bespoke preflight for the owned-children replacement contract adapted from
// ALPS reconciliation 81e7ad00d. Native HDF5 names and links need validation
// before the ALPS logical-path resolver or typed deletion APIs see them.
#ifndef ALPS_HDF5_DETAIL_COLLECTION_HPP
#define ALPS_HDF5_DETAIL_COLLECTION_HPP

#include <alps/hdf5/archive.hpp>
#include <algorithm>
#include <string>
#include <vector>

namespace alps::hdf5::detail {

inline std::string collection_child_path(const std::string &parent,
                                         const std::string &name) {
    if (name.empty() || name == "." || name == ".."
        || name.find_first_of("/@") != std::string::npos
        || name.find('\0') != std::string::npos)
        throw invalid_path("collection child is not a flat archive name: " + name + ALPS_STACKTRACE);
    return parent + (parent == "/" ? "" : "/") + name;
}

inline std::vector<std::string> collection_children(archive &ar, const std::string &root) {
    auto children = ar.list_children(root);
    // Validate every child before any deletion, including names that would
    // resolve to attributes or ancestors and links that typed deletion skips.
    for (const auto &name : children) {
        const auto path = collection_child_path(root, name);
        if (!ar.is_group(path) && !ar.is_data(path))
            throw invalid_path("unsupported collection child: " + path + ALPS_STACKTRACE);
    }
    return children;
}

inline std::vector<std::string> preflight_collection_replace(archive &ar,
                                                            const std::string &path) {
    const auto root = ar.complete_path(path);
    if (root.find('@') != std::string::npos || root.find('\0') != std::string::npos)
        throw invalid_path("collections require a group path" + ALPS_STACKTRACE);
    std::string parent = "/";
    for (std::size_t begin = 1; begin < root.size();) {
        const auto end = root.find('/', begin);
        const auto name = root.substr(begin, end == std::string::npos ? end : end - begin);
        const auto current = collection_child_path(parent, name);
        if (!ar.is_group(current)) {
            // Replacing the selected dataset is supported; traversing one is
            // not. Listing the parent distinguishes absence from dangling links
            // and named datatypes that is_group/is_data both report as false.
            const auto siblings = ar.list_children(parent);
            const bool exists = std::find(siblings.begin(), siblings.end(), name) != siblings.end();
            if (exists && (end != std::string::npos || !ar.is_data(current)))
                throw invalid_path("unsupported collection target or ancestor: " + current + ALPS_STACKTRACE);
            return {};
        }
        parent = current;
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return collection_children(ar, root);
}

} // namespace alps::hdf5::detail
#endif
