! Copyright (C) 2011 Synge Todo; 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
module hello_state
  use iso_c_binding
  implicit none
  include 'alps/fortran/alps_fortran.h'
  type, bind(C) :: state
    integer(c_int) :: finished = 0
    real(c_double) :: x
    integer(c_int) :: y
    character(c_char) :: world(32)
    real(c_float) :: fractions(4)
    integer(c_int64_t) :: wide(2)
    character(c_char) :: words(8,2)
  end type
contains
  function instance(caller) result(s)
    type(c_ptr), value :: caller
    type(state), pointer :: s
    call c_f_pointer(alps_get_context(caller), s)
  end function
end module

subroutine alps_init(caller) bind(C)
  use hello_state
  implicit none
  type(c_ptr), value :: caller
  type(state), pointer :: s
  allocate(s)
  call alps_set_context(caller, c_loc(s))
  call alps_get_parameter(caller, c_loc(s%x), 'X'//c_null_char, ALPS_DOUBLE_PRECISION, 0_c_size_t)
  call alps_get_parameter(caller, c_loc(s%y), 'Y'//c_null_char, ALPS_INT, 0_c_size_t)
  call alps_get_parameter(caller, c_loc(s%world), 'WORLD'//c_null_char, ALPS_CHAR, 32_c_size_t)
  call alps_get_parameter(caller, c_loc(s%fractions(1)), 'X'//c_null_char, ALPS_REAL, 0_c_size_t)
  call alps_get_parameter(caller, c_loc(s%wide(1)), 'Y'//c_null_char, ALPS_INT64, 0_c_size_t)
  if (alps_failed(caller)) return
  s%fractions = s%fractions(1)*[1,2,3,4]
  s%wide(2) = 1234567890123_c_int64_t
  s%words(:,1) = ['o','n','e',' ',' ',' ',' ',' ']
  s%words(:,2) = ['T','W','O',' ',' ',' ',' ',' ']
  print *, 'Hello ',s%world,' X=',s%x,' Y=',s%y,' Z defined: ',alps_parameter_defined(caller,'Z'//c_null_char)
end subroutine

subroutine alps_run(caller) bind(C)
  use hello_state
  implicit none
  type(c_ptr), value :: caller
  type(state), pointer :: s
  s => instance(caller)
  call alps_accumulate_observable(caller, c_loc(s%fractions), 4_c_size_t, ALPS_REAL, 'Values'//c_null_char)
  s%finished = 1
end subroutine

subroutine alps_init_observables(caller) bind(C)
  use hello_state
  implicit none
  type(c_ptr), value :: caller
  call alps_init_observable(caller, 4_c_size_t, 'Values'//c_null_char)
end subroutine

subroutine alps_progress(progress,caller) bind(C)
  use hello_state
  implicit none
  type(c_ptr), value :: caller
  real(c_double), intent(out) :: progress
  type(state), pointer :: s
  s => instance(caller)
  progress=s%finished
end subroutine

subroutine alps_is_thermalized(thermalized,caller) bind(C)
  use iso_c_binding
  implicit none
  type(c_ptr), value :: caller
  integer(c_int), intent(out) :: thermalized
  thermalized=1
end subroutine

subroutine alps_save(caller) bind(C)
  use hello_state
  implicit none
  type(c_ptr), value :: caller
  type(state), pointer :: s
  s => instance(caller)
  call alps_dump(caller, c_loc(s%finished), 1_c_size_t, ALPS_INT, 0_c_size_t)
  call alps_dump(caller, c_loc(s%x), 1_c_size_t, ALPS_DOUBLE_PRECISION, 0_c_size_t)
  call alps_dump(caller, c_loc(s%y), 1_c_size_t, ALPS_INT, 0_c_size_t)
  call alps_dump(caller, c_loc(s%world), 1_c_size_t, ALPS_CHAR, 32_c_size_t)
  call alps_dump(caller, c_loc(s%fractions), 4_c_size_t, ALPS_REAL, 0_c_size_t)
  call alps_dump(caller, c_loc(s%wide), 2_c_size_t, ALPS_INT64, 0_c_size_t)
  call alps_dump(caller, c_loc(s%words), 2_c_size_t, ALPS_CHAR, 8_c_size_t)
end subroutine

subroutine alps_load(caller) bind(C)
  use hello_state
  implicit none
  type(c_ptr), value :: caller
  type(state), pointer :: s
  s => instance(caller)
  call alps_restore(caller, c_loc(s%finished), 1_c_size_t, ALPS_INT, 0_c_size_t)
  call alps_restore(caller, c_loc(s%x), 1_c_size_t, ALPS_DOUBLE_PRECISION, 0_c_size_t)
  call alps_restore(caller, c_loc(s%y), 1_c_size_t, ALPS_INT, 0_c_size_t)
  call alps_restore(caller, c_loc(s%world), 1_c_size_t, ALPS_CHAR, 32_c_size_t)
  call alps_restore(caller, c_loc(s%fractions), 4_c_size_t, ALPS_REAL, 0_c_size_t)
  call alps_restore(caller, c_loc(s%wide), 2_c_size_t, ALPS_INT64, 0_c_size_t)
  call alps_restore(caller, c_loc(s%words), 2_c_size_t, ALPS_CHAR, 8_c_size_t)
  if (alps_failed(caller)) return
  if (s%finished/=alps_completed_sweeps(caller) .or. s%finished<0 .or. s%finished>1) &
    call alps_fail(caller, 'Invalid hello checkpoint progress'//c_null_char)
end subroutine

subroutine alps_finalize(caller) bind(C)
  use hello_state
  implicit none
  type(c_ptr), value :: caller
  type(state), pointer :: s
  s => instance(caller)
  if (associated(s)) deallocate(s)
  call alps_set_context(caller, c_null_ptr)
end subroutine
