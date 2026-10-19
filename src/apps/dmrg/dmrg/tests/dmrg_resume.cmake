# Resuming an interrupted dmrg run. A run with CHECKPOINT_DIRECTORY is
# killed outright during its regular sweeps (execute_process's TIMEOUT, as a
# batch system's time limit would); its checkpoint names where it continues,
# and the run resumed there with START_SWEEP, START_DIR and START_ITER writes
# the energies, iteration rows and measurements of an uninterrupted run
# without a checkpoint, bit for bit (dmrg_compare_results). A checkpointed
# run that completes leaves its directory empty. Resume input that the
# checkpoint does not support is rejected and leaves the checkpoint usable.
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}/xml")
file(COPY "${xml_dir}/lattices.xml" "${xml_dir}/models.xml"
  "${xml_dir}/ALPS.xsl" DESTINATION "${test_dir}/xml")
set(ENV{ALPS_XML_PATH} "${test_dir}/xml")

set(model [=[
LATTICE="open chain lattice"
MODEL="spin"
CONSERVED_QUANTUMNUMBERS="Sz"
Sz_total=0
L=40
J=1
SWEEPS=6
MAXSTATES=100
MEASURE_LOCAL[Sz]="Sz"
MEASURE_CORRELATIONS[SzSz]="Sz:Sz"
]=])
set(checkpoint "${test_dir}/checkpoint")
set(keep "CHECKPOINT_DIRECTORY=\"${checkpoint}\"\n")

# Runs the model with the parameter lines `extra` in ${test_dir}/${name};
# sets ${name}_status and ${name}_output. TIMEOUT <seconds> kills the run.
function(dmrg_run name extra)
  cmake_parse_arguments(RUN "" "TIMEOUT" "" ${ARGN})
  set(dir "${test_dir}/${name}")
  file(REMOVE_RECURSE "${dir}")
  file(MAKE_DIRECTORY "${dir}")
  file(WRITE "${dir}/parameters" "${model}${extra}")
  execute_process(COMMAND "${parameter2xml}" parameters
    WORKING_DIRECTORY "${dir}" OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
  set(limit)
  if(RUN_TIMEOUT)
    set(limit TIMEOUT ${RUN_TIMEOUT})
  endif()
  execute_process(COMMAND "${dmrg}" --write-xml parameters.in.xml
    WORKING_DIRECTORY "${dir}" ${limit}
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE output)
  set(${name}_status "${status}" PARENT_SCOPE)
  set(${name}_output "${output}" PARENT_SCOPE)
endfunction()

function(expect_success name)
  if(NOT ${name}_status EQUAL 0)
    message(FATAL_ERROR "${name}: the run failed (${${name}_status}):\n${${name}_output}")
  endif()
  execute_process(COMMAND "${compare}" "${test_dir}/reference/parameters.task1.out.h5"
    "${test_dir}/${name}/parameters.task1.out.h5"
    RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE out)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "${name}: the results differ from the uninterrupted run:\n${out}")
  endif()
  file(GLOB left "${checkpoint}/*")
  if(left)
    message(FATAL_ERROR "${name}: the checkpoint directory was not emptied: ${left}")
  endif()
endfunction()

# The scheduler reports a task's error and goes on, so a rejected run is
# one that wrote no result.
function(expect_rejected name pattern)
  if(EXISTS "${test_dir}/${name}/parameters.task1.out.h5")
    message(FATAL_ERROR "${name}: accepted; expected an error matching '${pattern}'")
  endif()
  if(NOT ${name}_output MATCHES "${pattern}")
    message(FATAL_ERROR "${name}: expected an error matching '${pattern}':\n${${name}_output}")
  endif()
endfunction()

dmrg_run(reference "")
if(NOT reference_status EQUAL 0)
  message(FATAL_ERROR "reference: the run failed:\n${reference_output}")
endif()
file(MAKE_DIRECTORY "${checkpoint}")
dmrg_run(complete "${keep}")
expect_success(complete)

dmrg_run(no_directory "START_SWEEP=2\n")
expect_rejected(no_directory "set CHECKPOINT_DIRECTORY")
dmrg_run(nothing "${keep}START_SWEEP=2\n")
expect_rejected(nothing "holds no checkpoint to resume")
dmrg_run(bad_direction "${keep}START_SWEEP=2\nSTART_DIR=2\n")
expect_rejected(bad_direction "START_DIR can assume the values 0")

# Kill the run ever later until it is stopped in its regular sweeps.
set(stage "")
foreach(seconds 0.25 0.5 1 2 4 8)
  file(REMOVE_RECURSE "${checkpoint}")
  file(MAKE_DIRECTORY "${checkpoint}")
  dmrg_run(interrupted "${keep}" TIMEOUT ${seconds})
  if(interrupted_status EQUAL 0)
    message(FATAL_ERROR "the run completed within ${seconds} s; no interruption landed in its sweeps")
  endif()
  set(stage "")
  if(EXISTS "${checkpoint}/checkpoint.txt")
    file(STRINGS "${checkpoint}/checkpoint.txt" stage REGEX "^stage ")
  endif()
  if(stage STREQUAL "stage sweeps")
    set(killed_after ${seconds})
    break()
  endif()
endforeach()
if(NOT stage STREQUAL "stage sweeps")
  message(FATAL_ERROR "no run was interrupted in its regular sweeps (last checkpoint: '${stage}')")
endif()
file(STRINGS "${checkpoint}/checkpoint.txt" next REGEX "^next ")
if(NOT next MATCHES "^next ([0-9]+) ([01]) ([0-9]+)$")
  message(FATAL_ERROR "unexpected checkpoint point: '${next}'")
endif()
set(sweep ${CMAKE_MATCH_1})
set(direction ${CMAKE_MATCH_2})
set(iteration ${CMAKE_MATCH_3})
message(STATUS "killed after ${killed_after} s; the checkpoint continues at sweep ${sweep}, "
  "direction ${direction}, iteration ${iteration}")

if(sweep GREATER 1)
  math(EXPR other "${sweep} - 1")
else()
  math(EXPR other "${sweep} + 1")
endif()
dmrg_run(wrong_point "${keep}START_SWEEP=${other}\nSTART_DIR=${direction}\nSTART_ITER=${iteration}\n")
expect_rejected(wrong_point "checkpoint continues at sweep ${sweep}, ")
dmrg_run(fresh_start "${keep}")
expect_rejected(fresh_start "holds the checkpoint of an interrupted solve")

dmrg_run(resumed "${keep}START_SWEEP=${sweep}\nSTART_DIR=${direction}\nSTART_ITER=${iteration}\n")
expect_success(resumed)
