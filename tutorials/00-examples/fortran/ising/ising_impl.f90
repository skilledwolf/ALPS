! Copyright (C) 2011 Synge Todo; 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
module ising_state
  use iso_c_binding
  implicit none
  include 'alps/fortran/alps_fortran.h'
  type state
    integer(c_int) :: length = 0
    integer(c_int64_t) :: steps = 0, warmup = 0, sweeps = 0
    real(c_double) :: temperature = 0, probability(-4:4)
    integer(c_int), pointer :: spin(:,:) => null()
  end type
contains
  function instance(caller) result(s)
    type(c_ptr), value :: caller
    type(state), pointer :: s
    call c_f_pointer(alps_get_context(caller), s)
  end function
end module

subroutine alps_init(caller) bind(C)
  use ising_state
  implicit none
  type(c_ptr), value :: caller
  type(state), pointer :: s
  integer :: m
  allocate(s)
  call alps_set_context(caller, c_loc(s))
  call alps_get_parameter(caller, c_loc(s%length), 'L'//c_null_char, ALPS_INT, 0_c_size_t)
  call alps_get_parameter(caller, c_loc(s%temperature), 'TEMPERATURE'//c_null_char, ALPS_DOUBLE_PRECISION, 0_c_size_t)
  call alps_get_parameter(caller, c_loc(s%warmup), 'THERMALIZATION'//c_null_char, ALPS_INT64, 0_c_size_t)
  call alps_get_parameter(caller, c_loc(s%sweeps), 'SWEEPS'//c_null_char, ALPS_INT64, 0_c_size_t)
  if (alps_failed(caller)) return
  if (s%length<2 .or. s%temperature<=0 .or. s%warmup<0 .or. s%sweeps<=0) then
    call alps_fail(caller, 'Invalid Ising parameters'//c_null_char)
    return
  end if
  allocate(s%spin(s%length,s%length))
  s%spin = 1
  do m=-4,4
    s%probability(m) = (1.0_c_double+tanh(real(m,c_double)/s%temperature))/2
  end do
end subroutine

subroutine alps_init_observables(caller) bind(C)
  use ising_state
  implicit none
  type(c_ptr), value :: caller
  call alps_init_observable(caller, 1_c_size_t, 'Energy'//c_null_char)
  call alps_init_observable(caller, 1_c_size_t, 'Magnetization'//c_null_char)
end subroutine

subroutine alps_run(caller) bind(C)
  use ising_state
  implicit none
  type(c_ptr), value :: caller
  type(state), pointer :: s
  integer :: i,j,l,m
  real(c_double), target :: energy, magnetization
  s => instance(caller)
  l = s%length
  do j=1,l
    do i=1,l
      m = s%spin(modulo(i-2,l)+1,j)+s%spin(modulo(i,l)+1,j) &
        + s%spin(i,modulo(j-2,l)+1)+s%spin(i,modulo(j,l)+1)
      s%spin(i,j) = -1
      if (alps_random(caller)<s%probability(m)) s%spin(i,j) = 1
    end do
  end do
  energy = 0
  magnetization = 0
  do j=1,l
    do i=1,l
      energy = energy-s%spin(i,j)*(s%spin(modulo(i,l)+1,j)+s%spin(i,modulo(j,l)+1))
      magnetization = magnetization+s%spin(i,j)
    end do
  end do
  energy = energy/(real(l,c_double)*l)
  magnetization = magnetization/(real(l,c_double)*l)
  call alps_accumulate_observable(caller, c_loc(energy), 1_c_size_t, ALPS_DOUBLE_PRECISION, 'Energy'//c_null_char)
  call alps_accumulate_observable(caller, c_loc(magnetization), 1_c_size_t, ALPS_DOUBLE_PRECISION, &
    'Magnetization'//c_null_char)
  s%steps = s%steps+1
end subroutine

subroutine alps_progress(progress, caller) bind(C)
  use ising_state
  implicit none
  type(c_ptr), value :: caller
  real(c_double), intent(out) :: progress
  type(state), pointer :: s
  s => instance(caller)
  progress = real(s%steps,c_double)/(real(s%warmup,c_double)+s%sweeps)
end subroutine

subroutine alps_is_thermalized(thermalized, caller) bind(C)
  use ising_state
  implicit none
  type(c_ptr), value :: caller
  integer(c_int), intent(out) :: thermalized
  type(state), pointer :: s
  s => instance(caller)
  thermalized = 0
  if (s%steps>=s%warmup) thermalized = 1
end subroutine

subroutine alps_save(caller) bind(C)
  use ising_state
  implicit none
  type(c_ptr), value :: caller
  type(state), pointer :: s
  s => instance(caller)
  call alps_dump(caller, c_loc(s%steps), 1_c_size_t, ALPS_INT64, 0_c_size_t)
  call alps_dump(caller, c_loc(s%spin), size(s%spin,kind=c_size_t), ALPS_INT, 0_c_size_t)
end subroutine

subroutine alps_load(caller) bind(C)
  use ising_state
  implicit none
  type(c_ptr), value :: caller
  type(state), pointer :: s
  s => instance(caller)
  call alps_restore(caller, c_loc(s%steps), 1_c_size_t, ALPS_INT64, 0_c_size_t)
  call alps_restore(caller, c_loc(s%spin), size(s%spin,kind=c_size_t), ALPS_INT, 0_c_size_t)
  if (alps_failed(caller)) return
  if (s%steps/=alps_completed_sweeps(caller) .or. any(s%spin/=1 .and. s%spin/=-1)) &
    call alps_fail(caller, 'Invalid Ising checkpoint state'//c_null_char)
end subroutine

subroutine alps_finalize(caller) bind(C)
  use ising_state
  implicit none
  type(c_ptr), value :: caller
  type(state), pointer :: s
  s => instance(caller)
  if (.not.associated(s)) return
  if (associated(s%spin)) deallocate(s%spin)
  deallocate(s)
  call alps_set_context(caller, c_null_ptr)
end subroutine
