// Copyright (C) 2011 Synge Todo; 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
// Application callbacks use BIND(C), with TYPE(C_PTR), VALUE :: caller.
extern "C" {
void alps_init(void* caller);
void alps_init_observables(void* caller);
void alps_run(void* caller);
void alps_progress(double* progress, void* caller);
void alps_is_thermalized(int* thermalized, void* caller);
void alps_finalize(void* caller);
void alps_save(void* caller);
void alps_load(void* caller);
}
