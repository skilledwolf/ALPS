/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: convert2xml.C 3523 2009-12-12 05:52:24Z troyer $ */

#include <alps/config.h>

#include <alps/hdf5.hpp>

#include <alps/scheduler/convert.h>
#include <alps/osiris/xdrdump.h>
#include <alps/parser/xslt_path.h>
#include <alps/scheduler/types.h>
#include <alps/scheduler/diag.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>

#include <fstream>
#include <stdexcept>

namespace alps {

void convert_spectrum(const std::string& inname) 
{
  boost::filesystem::path p(inname);
  ProcessList nowhere;
  scheduler::DiagTask<double> sim(nowhere,p);
  sim.checkpoint(p,true);
}

// Monte Carlo runs moved to native HDF5 results; their XDR and XML
// checkpoints are no longer read.
[[noreturn]] void reject_monte_carlo(const std::string& inname)
{
  boost::throw_exception(std::runtime_error(inname + " is neither a parameter file nor a "
    "spectrum task; convert2xml no longer converts Monte Carlo runs or checkpoints. "
    "Convert released results with alps-hdf5-convert --alea-results, or finish the "
    "run with ALPS 3.0."));
}

void convert_xml(const std::string& inname)
{
  bool is_spectrum=false;
  std::string h5name = inname.substr(0, inname.find_last_of('.')) + ".h5";
  if (boost::filesystem::exists(boost::filesystem::path(h5name))) 
  {
    hdf5::archive ar(h5name);
    if (ar.is_group("/spectrum"))
      is_spectrum=true;
  }
  if (!is_spectrum)
    reject_monte_carlo(inname);
  convert_spectrum(inname);
}



void convert_params(const std::string& inname)
{
  ParameterList list;
  {
    std::ifstream in(inname.c_str());
    in >> list;
  }

  std::string basename = boost::filesystem::path(inname).filename().string();
  std::cout << "Converting parameter file " << inname << " to "
            <<  basename+".in.xml" << std::endl;

  oxstream out(boost::filesystem::path((basename+".in.xml").c_str()));
  out << header("UTF-8")
      << stylesheet(xslt_path("ALPS.xsl"))
      << start_tag("JOB")
      << xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
      << attribute("xsi:noNamespaceSchemaLocation",
                         "http://xml.comp-phys.org/2003/8/job.xsd")
      << start_tag("OUTPUT")
      << attribute("file", basename+".out.xml")
      << end_tag("OUTPUT");

  for (unsigned int i = 0; i < list.size(); ++i) {
    std::string taskname =
      basename+".task"+boost::lexical_cast<std::string,int>(i+1);
    out << start_tag("TASK") << attribute("status","new")
        << start_tag("INPUT")
        << attribute("file", taskname + ".in.xml")
        << end_tag("INPUT")
        << start_tag("OUTPUT")
        << attribute("file", taskname + ".out.xml")
        << end_tag("OUTPUT")
        << end_tag("TASK");
    //      out << "    <CPUS min=\"1\">\n";
    oxstream task(boost::filesystem::path((taskname+".in.xml").c_str()));
    task << header("UTF-8")
         << stylesheet(xslt_path("ALPS.xsl"));
    task << start_tag("SIMULATION")
         << xml_namespace("xsi",
                                "http://www.w3.org/2001/XMLSchema-instance")
         << attribute("xsi:noNamespaceSchemaLocation",
                            "http://xml.comp-phys.org/2002/10/QMCXML.xsd");
    task << list[i];
    task << end_tag("SIMULATION");
  }

  out << end_tag("JOB");
}

std::string convert2xml(std::string const& inname)
{
    IXDRFileDump dump=IXDRFileDump(boost::filesystem::path(inname));
    int type;
    dump >> type;
    switch (type) {
    case scheduler::MCDump_scheduler:
    case scheduler::MCDump_task:
    case scheduler::MCDump_run:
      reject_monte_carlo(inname);
    default:
      {
        bool isxml=false;
        {
          std::ifstream is(inname.c_str());
          char c1=is.get();
          char c2=is.get();
          isxml = (c1=='<' && c2 =='?');
        }
        if (isxml)
          convert_xml(inname);
        else
          convert_params(inname);
      }
    }
  return inname+".in.xml";
}

} // end namespace
