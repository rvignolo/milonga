/*------------ -------------- -------- --- ----- ---   --       -            -
 * milonga's method of characteristics protoypes
 *
 *  Copyright (C) 2016--2018 ramiro vignolo
 *
 *  This file is part of milonga.
 *
 *  milonga is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  milonga is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with wasora.  If not, see <http://www.gnu.org/licenses/>.
 *------------------- ------------  ----    --------  --     -       -         -
 */

#include "../milonga.h"

#define angular_index(t,d,p,g)     ((t)* 2 * milonga.moc_solver.polar_quadrature->n_polar_2 * milonga.groups + \
                                    (d)    * milonga.moc_solver.polar_quadrature->n_polar_2 * milonga.groups + \
                                    (p)                                                     * milonga.groups + \
                                    (g))
#define reduced_angular_index(p,g) ((p)                                                     * milonga.groups + \
                                    (g))

extern int moc_volumes_problem_init(void);
extern int moc_volumes_allocate_solver(void);
extern int moc_volumes_compute_total_weights(tracks_t *);
extern int moc_volumes_results_fill_args(function_t *);
extern int moc_volumes_solve_eigen_problem(int);
extern int moc_volumes_define_exp_evaluator(void);
extern int moc_volumes_compute_exp_table(void);
//extern int moc_volumes_verify_exponential_approx(void);
extern int moc_volumes_set_uniform_phi(double);
extern int moc_volumes_update_prev_phi(void);
extern int moc_volumes_set_uniform_start_boundary_psi(double);
extern int moc_volumes_update_boundary_psi(void);
extern int moc_volumes_normalize_fluxes(void);
extern int moc_volumes_compute_q(void);
extern int moc_volumes_compute_phi(void);
extern int moc_volumes_tally_phi(segment_t *, double *, int);
extern int moc_volumes_set_start_boundary_psi(track_t *, double *, int);
extern int moc_volumes_add_source_to_phi(void);
extern int moc_volumes_compute_keff(void);
extern int moc_volumes_compute_residual(double *);
extern int moc_volumes_results_fill_flux(void);
extern int moc_volumes_normalize_flux(void);
extern int moc_volumes_results_fill_power(void);
extern int moc_volumes_problem_free(void);

extern double moc_volumes_compute_linear_exponential(double, int);
extern double moc_volumes_compute_intrinsic_exponential(double, int);
extern double moc_volumes_cell_integral(cell_t *, expr_t *);

extern int wasora_dummy_true_cmp(void *, void *);
extern int wasora_dummy_equal_cmp(void *, void *);
extern int wasora_dummy_false_cmp(void *, void *);
