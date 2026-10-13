/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// How the legacy scheduler names the checkpoint files of a task's runs. The
// search was moved unchanged from WorkerTask::write_xml_body (workertask.C).
//
// A run is named when its task is first checkpointed and keeps that name in
// later checkpoints, also under another task file name. Each run has an XDR
// checkpoint and an HDF5 checkpoint whose name adds ".h5". Names are
// relative to the directory of the task file.

#ifndef ALPS_SCHEDULER_RUN_FILES_H
#define ALPS_SCHEDULER_RUN_FILES_H

#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>

namespace alps {
namespace scheduler {

/// The XDR checkpoint name for a run without one, of a task checkpointed as
/// task_file: the task file name up to its last '.', followed by ".run<j>"
/// with the smallest j >= 1 that is not among the names already taken.
inline std::string run_checkpoint_name(const boost::filesystem::path& task_file,
                                       const std::vector<std::string>& taken)
{
  int j=0;
  bool found=false;
  std::string name;
  do {
    found = false;
    name =task_file.filename().string();
    name = name.substr(0, name.find_last_of('.'));
    name+= ".run" + boost::lexical_cast<std::string,int>(j+1);
    for (unsigned int k=0;k<taken.size();++k)
    if(taken[k]==name)
      found=true;
    j++;
  } while (found);
  return name;
}

} // namespace scheduler
} // namespace alps

#endif // ALPS_SCHEDULER_RUN_FILES_H
