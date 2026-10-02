// Derived from ALPSCore params at 7146b9e1f017938a94e5dae35d88467cc5ba7969.
// Copyright (C) 1998-2018 ALPS Collaboration; modifications (C) 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
// See ALPSCore-LICENSE.txt in this module for the original permission notice.
#pragma once
#include <stdexcept>
#include <string>
namespace alps::params_ns::exception {
struct exception_base : std::runtime_error {
    exception_base(const std::string &key, const std::string &reason)
        : std::runtime_error("Parameter '" + key + "': " + reason), name_(key) {}
    const std::string &name() const noexcept { return name_; }

  private:
    std::string name_;
};
struct uninitialized_value : exception_base {
    using exception_base::exception_base;
};
struct type_mismatch : exception_base {
    using exception_base::exception_base;
};
struct value_mismatch : exception_base {
    using exception_base::exception_base;
};
} // namespace alps::params_ns::exception
