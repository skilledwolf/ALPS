/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2015 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "parameter_conversion.hpp"

namespace {
struct xml_writer {
  static const char* extension() { return ".xml"; }
  static const char* name() { return "XML"; }
  static void write(const std::string& filename, const alps::Parameters& params) {
    alps::oxstream task(boost::filesystem::path((filename).c_str()));
    task << alps::header("UTF-8")
         << alps::stylesheet(alps::xslt_path("ALPS.xsl"));
    task << alps::start_tag("SIMULATION")
         << alps::xml_namespace("xsi",
                                "http://www.w3.org/2001/XMLSchema-instance")
         << alps::attribute("xsi:noNamespaceSchemaLocation",
                            "http://xml.comp-phys.org/2002/10/QMCXML.xsd");
    task << params;
    task << alps::end_tag("SIMULATION");
  }
};
}

int main(int argc, char** argv) {
  return parameter_conversion::run<xml_writer>(argc, argv);
}
