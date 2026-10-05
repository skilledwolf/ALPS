! Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
! USE ISO_C_BINDING before IMPLICIT NONE; include this file after it.
integer(c_int), parameter :: ALPS_CHAR=0, ALPS_INT=1, ALPS_INT64=2, ALPS_REAL=3, ALPS_DOUBLE_PRECISION=4
interface
  function alps_completed_sweeps(caller) bind(C) result(value)
    import
    type(c_ptr), value :: caller
    integer(c_int64_t) :: value
  end function
  function alps_get_context(caller) bind(C) result(value)
    import
    type(c_ptr) :: value
    type(c_ptr), value :: caller
  end function
  function alps_failed(caller) bind(C) result(value)
    import
    logical(c_bool) :: value
    type(c_ptr), value :: caller
  end function
  function alps_random(caller) bind(C) result(value)
    import
    real(c_double) :: value
    type(c_ptr), value :: caller
  end function
  function alps_parameter_defined(caller, name) bind(C) result(value)
    import
    logical(c_bool) :: value
    type(c_ptr), value :: caller
    character(c_char) :: name(*)
  end function
  subroutine alps_set_context(caller, context) bind(C)
    import
    type(c_ptr), value :: caller, context
  end subroutine
  subroutine alps_fail(caller, message) bind(C)
    import
    type(c_ptr), value :: caller
    character(c_char) :: message(*)
  end subroutine
  subroutine alps_get_parameter(caller, data, name, kind, width) bind(C)
    import
    type(c_ptr), value :: caller, data
    character(c_char) :: name(*)
    integer(c_int), value :: kind
    integer(c_size_t), value :: width
  end subroutine
  subroutine alps_init_observable(caller, count, name) bind(C)
    import
    type(c_ptr), value :: caller
    integer(c_size_t), value :: count
    character(c_char) :: name(*)
  end subroutine
  subroutine alps_accumulate_observable(caller, data, count, kind, name) bind(C)
    import
    type(c_ptr), value :: caller, data
    integer(c_size_t), value :: count
    integer(c_int), value :: kind
    character(c_char) :: name(*)
  end subroutine
  subroutine alps_dump(caller, data, count, kind, width) bind(C)
    import
    type(c_ptr), value :: caller, data
    integer(c_size_t), value :: count
    integer(c_int), value :: kind
    integer(c_size_t), value :: width
  end subroutine
  subroutine alps_restore(caller, data, count, kind, width) bind(C)
    import
    type(c_ptr), value :: caller, data
    integer(c_size_t), value :: count
    integer(c_int), value :: kind
    integer(c_size_t), value :: width
  end subroutine
end interface
