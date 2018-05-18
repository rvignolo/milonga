/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's method of characteristics routines
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

#include "moc_volumes.h"

int moc_volumes_problem_init(void) {
  
  int g;
  
  if (milonga.moc_solver.tracks == NULL) {
    wasora_push_error_message("no tracking performed");
    return WASORA_RUNTIME_ERROR;
  }
  
  // nos aseguramos de usar la malla asociada al ray tracing dado en MILONGA_PROBLEM
  if ((wasora_mesh.main_mesh = milonga.moc_solver.tracks->mesh) == NULL) {
    wasora_push_error_message("no mesh found");
    return WASORA_RUNTIME_ERROR;
  }
  
  if (wasora_mesh.main_mesh->structured == 0) {
    if (wasora_mesh.main_mesh->cell == NULL) {
      wasora_call(mesh_element2cell(wasora_mesh.main_mesh));
    }
    if (wasora_mesh.main_mesh->cell[0].ifaces == NULL) {
      // moc no hace uso de esto por ahora, pero por completitud lo ponemos (quizas en esquema de aceleracion sirva)
      wasora_call(mesh_find_neighbors(wasora_mesh.main_mesh));
    }
    
    wasora_call(mesh_compute_coords(wasora_mesh.main_mesh));
    // moc no hace uso de esto por ahora, pero por completitud lo ponemos (quizas en esquema de aceleracion sirva)
    wasora_call(mesh_fill_neighbors(wasora_mesh.main_mesh));
  } else {
    wasora_mesh_struct_init_rectangular_for_cells(wasora_mesh.main_mesh);
  }
  
  milonga.spatial_unknowns = wasora_mesh.main_mesh->n_cells;
  milonga.problem_size = milonga.spatial_unknowns * milonga.groups;
  
  wasora_call(moc_volumes_allocate_solver());
  wasora_call(moc_volumes_compute_total_weights(milonga.moc_solver.tracks));
  
  wasora_var(wasora_mesh.vars.cells) = (double)wasora_mesh.main_mesh->n_cells;
  wasora_var(wasora_mesh.vars.nodes) = (double)wasora_mesh.main_mesh->n_nodes;
  wasora_var(wasora_mesh.vars.elements) = (double)wasora_mesh.main_mesh->n_elements;
  wasora_mesh.main_mesh->data_type = data_type_element;
  
  for (g = 0; g < milonga.groups; g++) {
    wasora_call(moc_volumes_results_fill_args(milonga.functions.phi[g]));
  }
  wasora_call(moc_volumes_results_fill_args(milonga.functions.pow));
  
  wasora_call(mesh_cell_indexes(wasora_mesh.main_mesh, milonga.groups));
  
  return WASORA_RUNTIME_OK;  
}

int moc_volumes_allocate_solver(void) {
  
  int scalar_size, angular_size;
  
  // primero nos fijamos a cual cuadratura polar tenemos que apuntar:
  //   si nos dieron una nueva en MILONGA_PROBLEM apuntamos a ella
  //   si no nos dieron ninguna, conservamos la cargada por defecto en un ray tracing (no hacemos nada)
  if (milonga.moc_solver.polar_quadrature != NULL) {
    milonga.moc_solver.tracks->quadrature->polar = milonga.moc_solver.polar_quadrature;
  }
  
  scalar_size = milonga.problem_size;
  angular_size = milonga.moc_solver.tracks->n_total_tracks * 2 * milonga.moc_solver.tracks->quadrature->polar->n_polar_2 * milonga.groups;
  
  milonga.moc_solver.scalar_size = scalar_size;
  milonga.moc_solver.angular_size = angular_size;
  
  milonga.moc_solver.phi = calloc(scalar_size, sizeof(double));
  milonga.moc_solver.old_phi = calloc(scalar_size, sizeof(double));
  
  milonga.moc_solver.reduced_source = calloc(scalar_size, sizeof(double));
  
  milonga.moc_solver.boundary_psi = calloc(angular_size, sizeof(double));
  milonga.moc_solver.start_boundary_psi = calloc(angular_size, sizeof(double));
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_compute_total_weights(tracks_t *tracks) {
  
  int i, j;
  
  tracks->quadrature->w_total = calloc(tracks->quadrature->azimuthal->n_azim_2, sizeof(double *));
  for (i = 0; i < tracks->quadrature->azimuthal->n_azim_2; i++) {
    // allocamos solo n_polar_2 porque hay simetria polar y simplemente multiplicaremos por 2 al hacer el tally del scalar flux
    tracks->quadrature->w_total[i] = calloc(tracks->quadrature->polar->n_polar_2, sizeof(double));
  }
  
  for (i = 0; i < tracks->quadrature->azimuthal->n_azim_2; i++) {
    for(j = 0; j < tracks->quadrature->polar->n_polar_2; j++) {
      tracks->quadrature->w_total[i][j]  = 4.0 * M_PI;
      tracks->quadrature->w_total[i][j] *= tracks->quadrature->azimuthal->w[i];
      tracks->quadrature->w_total[i][j] *= tracks->quadrature->polar->w[j];
      tracks->quadrature->w_total[i][j] *= tracks->quadrature->effective_spacing[i];
      tracks->quadrature->w_total[i][j] *= tracks->quadrature->polar->sin_theta[j];
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_results_fill_args(function_t *function) {

  // estructurado o no, tenemos data
  function->data_size = milonga.spatial_unknowns;
  
  // pero tambien variables por si queremos hacer cuentitas
  function->var_argument = calloc(3, sizeof(var_t *));
  function->var_argument = wasora_mesh.vars.arr_x;

  function->data_argument = wasora_mesh.main_mesh->cells_argument;
  function->data_value = calloc(function->data_size, sizeof(double));
  
  // en volumes finitos decimos que la funcion es tipo mesh cell
  function->type = type_pointwise_mesh_cell;
  function->multidim_threshold = DEFAULT_MULTIDIM_INTERPOLATION_THRESHOLD;
  function->mesh = wasora_mesh.main_mesh;

  return WASORA_RUNTIME_OK;
}

int moc_volumes_solve_eigen_problem(int verbose) {
  
  int iter;
  
  double residual;
  
  // informacion del solver
  milonga.moc_solver.max_moc_iterations = (int)wasora_value(milonga.vars.moc_n_iter);
  milonga.moc_solver.max_moc_residual = wasora_value(milonga.vars.moc_rel_error);
  
  if(verbose) {
    printf("\n");
    printf("# milonga's moc volumes scheme\n");
    printf("\n");
    printf("# moc solver parameters:\n");
    printf("#   mesh name is '%s'\n", milonga.moc_solver.tracks->mesh->name);
    printf("#   ray tracing name is '%s'\n", milonga.moc_solver.tracks->name);
    if (milonga.moc_solver.polar_quadrature != NULL) {
      printf("#   polar quadrature name is '%s'\n", milonga.moc_solver.polar_quadrature->name);
    } else {
      printf("#   polar quadrature is set to milonga's default: tabuchi-yamamoto with 6 polar angles\n");
    }
    printf("#   maximum number of iterations is set to %d\n", milonga.moc_solver.max_moc_iterations);
    printf("#   maximum residual error set to %e\n", milonga.moc_solver.max_moc_residual);
    printf("\n");
    printf("# iter \t residual\n");
  }
  
  // esto va aca porque pueden cambiar los tau (sigmaT_g) en diferentes static steps
  // pero va a ir modificando sobre ray tracings posiblemente ya segmentados en pasos previos
  wasora_call(moc_volumes_define_exp_evaluator());
  
  // comenzamos a resolver el problema con algunas inicializaciones
  wasora_var_value(milonga.vars.keff) = 1.0;
  wasora_call(moc_volumes_set_uniform_phi(1.0));
  wasora_call(moc_volumes_update_old_phi());
  wasora_call(moc_volumes_set_uniform_start_boundary_psi(0.0));
  wasora_call(moc_volumes_update_boundary_psi());
  
  for (iter = 0; iter < milonga.moc_solver.max_moc_iterations; iter++) {
    
    wasora_call(moc_volumes_normalize_fluxes());
    wasora_call(moc_volumes_compute_q());
    wasora_call(moc_volumes_compute_phi());
    wasora_call(moc_volumes_compute_keff());
    wasora_call(moc_volumes_compute_residual(&residual));
    wasora_call(moc_volumes_update_old_phi());
    
    if (verbose) {
      printf("%d \t %e\n", iter, residual);
    }
    
    // por lo menos tenemos que iterar tres veces:
    //    i = 0: las boundary condition no estan (a menos que sea todo vaccum)
    //    i = 1: aparecen las boundary pero computar el residual con el paso anterior es fruta
    //    i = 2: podemos calcular un correcto residual (las boundary estan en los dos pasos que se usan para comparar)
    if (iter > 1 && residual < milonga.moc_solver.max_moc_residual) {
      break;
    }
  }
  
  // TODO: encapsular de alguna forma los warnings de wasora
  if (iter == milonga.moc_solver.max_moc_iterations) {
    printf("warning: moc calculation did not converge");
  }
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_define_exp_evaluator(void) {
  
  if (milonga.moc_solver.exp_evaluator != NULL) {
    free(milonga.moc_solver.exp_evaluator->table);
    free(milonga.moc_solver.exp_evaluator);
  }
  milonga.moc_solver.exp_evaluator = calloc(1, sizeof(exp_evaluator_t));
  
  // computamos los caminos opticos y ademas obtenemos el maximo q por ahi lo necesitemos
  wasora_call(track_compute_segments_tau(milonga.moc_solver.tracks));
  
  if (milonga.moc_solver.do_not_linearize_exp == 0) {
    
    milonga.moc_solver.exp_evaluator->max_err_allowed = EXP_PRECISION;
    milonga.moc_solver.exp_evaluator->max_tau_allowed = MAX_TAU;
    
    // segmentamos los segmentos si superamos el maximo tau permitido y actualizamos el valor max de tau sobre el ray tracing
    if (milonga.moc_solver.tracks->max_tau > milonga.moc_solver.exp_evaluator->max_tau_allowed) {
      wasora_call(track_segmentize_segments(milonga.moc_solver.tracks, milonga.moc_solver.exp_evaluator->max_tau_allowed));
      wasora_call(track_compute_segments_tau(milonga.moc_solver.tracks));
    }
    
    // generamos la tabla de interpolacion y apuntamos a la funcion a utilizar
    wasora_call(moc_volumes_compute_exp_table());
    milonga.moc_solver.exp_evaluator->compute_exponential = moc_volumes_compute_linear_exponential;
    
    // verifico si esta bien la tabla
    //wasora_call(moc_volumes_verify_exponential_approx());
    
  } else {
    // sino, usamos la evalucion intrínseca
    milonga.moc_solver.exp_evaluator->compute_exponential = moc_volumes_compute_intrinsic_exponential;
  }
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_compute_exp_table(void) {
  
  int i, p;
  int n_intervals;
  
  double sin_theta;
  double tau;
  double exponential;
  double q_pn, b_pn;
  
  // apuntamos por comodidad
  tracks_t *tracks = milonga.moc_solver.tracks;
  exp_evaluator_t *exp_evaluator = milonga.moc_solver.exp_evaluator;
  polar_quadrature_t *polar_quad = milonga.moc_solver.tracks->quadrature->polar;
  
  // el numero de intervalos depende del error maximo que querramos
  n_intervals = round(tracks->max_tau / sqrt(8.0 * exp_evaluator->max_err_allowed));
  // la longitud de cada intervalo sera
  exp_evaluator->delta = tracks->max_tau / n_intervals;
  // pero voy a agregar un intervalo para evaluar correctamente a tau_max
  n_intervals++;;
  // en la tabla final tenemos una pendiente y ordenada para cada angulo polar que recorramos en un transport sweep
  exp_evaluator->size = 2 * polar_quad->n_polar_2 * n_intervals;
  
  exp_evaluator->table = calloc(exp_evaluator->size, sizeof(double));
  
  for ( i = 0; i < n_intervals; i++) {
    for (p = 0; p < polar_quad->n_polar_2; p++) {
      
      sin_theta = polar_quad->sin_theta[p];
      
      // en el primer intervalo hacemos que la recta pase por el primer punto (para no superar exp(0) = 1.0)
      if (i == 0) {
        tau = i * exp_evaluator->delta;
        
        // pero en los otros tomo como punto el centro del intervalo
      } else {
        tau = (i + 0.5) * exp_evaluator->delta;
      }
      
      // hacemos las cuentitas necesarias
      exponential = exp(- tau / sin_theta);
      q_pn = - exponential / sin_theta;
      b_pn = exponential * (1 + tau / sin_theta);
      
      // y guardamos en la tabla
      exp_evaluator->table[2 * polar_quad->n_polar_2 * i + 2 * p] = q_pn;
      exp_evaluator->table[2 * polar_quad->n_polar_2 * i + 2 * p + 1] = b_pn;
    }
  }
  
  return WASORA_RUNTIME_OK;
}

double moc_volumes_compute_linear_exponential(double tau, int polar) {
  
  int n, n_polar_2;
  
  double exponential;
  
  exp_evaluator_t *exp_evaluator = milonga.moc_solver.exp_evaluator;
  
  n = floor(tau / exp_evaluator->delta);
  n_polar_2 = milonga.moc_solver.tracks->quadrature->polar->n_polar_2;
  
  exponential = exp_evaluator->table[2 * n_polar_2 * n + 2 * polar] * tau + exp_evaluator->table[2 * n_polar_2 * n + 2 * polar + 1];
  
  return 1.0 - exponential;
}

/*
int moc_volumes_verify_exponential_approx(void) {
  
  int i, p;
  int n_intervals;
  
  double tau, exact, approx;
  
  tracks_t *tracks = milonga.moc_solver.tracks;
  exp_evaluator_t *exp_evaluator = milonga.moc_solver.exp_evaluator;
  polar_quadrature_t *polar_quad = milonga.moc_solver.tracks->quadrature->polar_quad;
  
  FILE *gnuplot_approx_exp_points = NULL;
  static int file = 0;
  char filename[64];
  
  if (gnuplot_approx_exp_points == NULL) {
    file++;
    sprintf(filename, "gnuplot_approx_exp_points.dat-%d", file);
    gnuplot_approx_exp_points = fopen(filename, "w+");
  }
  
  n_intervals = round(tracks->max_tau / sqrt(8.0 * exp_evaluator->max_err_allowed));
  
  fprintf(gnuplot_approx_exp_points, "# max tau es %e\n", tracks->max_tau);
  
  // chekeo de las aproximaciones
  for (i = 0; i < (2*n_intervals+1); i++) {
    tau = i * exp_evaluator->delta / 2;
    fprintf(gnuplot_approx_exp_points, "%e\t", tau);
    for (p = 0; p < polar_quad->n_polar_2; p++) {
      
      exact = moc_volumes_compute_intrinsic_exponential(tau, p);
      approx = moc_volumes_compute_linear_exponential(tau, p);
      
      fprintf(gnuplot_approx_exp_points, "%e\t%e\t", exact, approx);
    }
    fprintf(gnuplot_approx_exp_points, "\n");
  }
  
  fclose(gnuplot_approx_exp_points);
  
  return WASORA_RUNTIME_OK;
}
*/

double moc_volumes_compute_intrinsic_exponential(double tau, int polar) {
  
  double exponential;
  double sin_theta;
  
  sin_theta = milonga.moc_solver.tracks->quadrature->polar->sin_theta[polar];
  
  exponential = exp(- tau / sin_theta);
  
  return 1.0 - exponential;
}

int moc_volumes_set_uniform_phi(double uniform_phi) {
  
  int i, g;
  
  cell_t *cell;
  
  for (i = 0; i < milonga.spatial_unknowns; i++) {
    for (g = 0; g < milonga.groups; g++) {
      cell = &wasora_mesh.main_mesh->cell[i];
      milonga.moc_solver.phi[cell->index[g]] = uniform_phi;
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_update_old_phi(void) {
  
  int i, g;
  
  cell_t *cell;
  
  for (i = 0; i < milonga.spatial_unknowns; i++) {
    cell = &wasora_mesh.main_mesh->cell[i];
    for (g = 0; g < milonga.groups; g++) {
      milonga.moc_solver.old_phi[cell->index[g]] = milonga.moc_solver.phi[cell->index[g]];
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_set_uniform_start_boundary_psi(double uniform_psi) {
  
  int t, d, p, g;
  int index;
  
  for (t = 0; t < milonga.moc_solver.tracks->n_total_tracks; t++) {
    for (d = 0; d < 2; d++) {
      for (p = 0; p < milonga.moc_solver.tracks->quadrature->polar->n_polar_2; p++) {
        for (g = 0; g < milonga.groups; g++) {
          index = angular_index(t,d,p,g);
          milonga.moc_solver.start_boundary_psi[index] = uniform_psi;
        }
      }
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_update_boundary_psi(void) {
  
  int t, d, p, g;
  int index;
  
  for (t = 0; t < milonga.moc_solver.tracks->n_total_tracks; t++) {
    for (d = 0; d < 2; d++) {
      for (p = 0; p < milonga.moc_solver.tracks->quadrature->polar->n_polar_2; p++) {
        for (g = 0; g < milonga.groups; g++) {
          index = angular_index(t,d,p,g);
          milonga.moc_solver.boundary_psi[index] = milonga.moc_solver.start_boundary_psi[index];
        }
      }
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_normalize_fluxes(void) {
  
  int i, g, g_prime;
  int t, d, p;
  int index;
  
  double total_fission_source;
  double normalization_factor;
  
  cell_t *cell;
  xs_t *material_xs;
  
  // conviene hacer sumas pairwise! hablarlo con jeremy
  // la fuente de fision total: integre en g (por ello el chi suma uno) y me quedan estas sumatorias
  total_fission_source = 0.0;
  for (i = 0; i < milonga.spatial_unknowns; i++) {
    
    cell = &wasora_mesh.main_mesh->cell[i];
    material_xs = (xs_t *) (cell->element->physical_entity->material->ext);
    
    for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
      total_fission_source += material_xs->xs_values.nuSigmaF[g_prime] * milonga.moc_solver.phi[cell->index[g_prime]] * cell->volume;
      //total_fission_source += moc_volumes_cell_integral(cell, material_xs->nuSigmaF[g_prime]) * milonga.moc_solver.phi[cell->index[g_prime]];
    }
  }
  
  normalization_factor = 1.0 / total_fission_source;
  
  for (i = 0; i < milonga.spatial_unknowns; i++) {
    for (g = 0; g < milonga.groups; g++) {
      cell = &wasora_mesh.main_mesh->cell[i];
      milonga.moc_solver.phi[cell->index[g]] *= normalization_factor;
      milonga.moc_solver.old_phi[cell->index[g]] *= normalization_factor;
    }
  }
  
  for (t = 0; t < milonga.moc_solver.tracks->n_total_tracks; t++) {
    for (d = 0; d < 2; d++) {
      for (p = 0; p < milonga.moc_solver.tracks->quadrature->polar->n_polar_2; p++) {
        for (g = 0; g < milonga.groups; g++) {
          index = angular_index(t,d,p,g);
          milonga.moc_solver.boundary_psi[index] *= normalization_factor;
          milonga.moc_solver.start_boundary_psi[index] *= normalization_factor;
        }
      }
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_compute_q(void) {
  
  int i;
  int g, g_prime;
  
  double sigma_t;
  
  cell_t *cell;
  xs_t *material_xs;
  
  for (i = 0; i < milonga.spatial_unknowns; i++) {
    
    cell = &wasora_mesh.main_mesh->cell[i];
    material_xs = (xs_t *) (cell->element->physical_entity->material->ext);
    
    for (g = 0; g < milonga.groups; g++) {
      
      milonga.moc_solver.reduced_source[cell->index[g]] = 0.0;
      for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
        // conviene hacer sumas pairwise! hablarlo con jeremy
        milonga.moc_solver.reduced_source[cell->index[g]] += material_xs->xs_values.SigmaS0[g_prime][g] * milonga.moc_solver.phi[cell->index[g_prime]];
        //milonga.moc_solver.reduced_source[cell->index[g]] += wasora_evaluate_expression(material_xs->SigmaS0[g_prime][g]) * milonga.moc_solver.phi[cell->index[g_prime]];
        
        milonga.moc_solver.reduced_source[cell->index[g]] += (1 / wasora_var_value(milonga.vars.keff)) * gsl_vector_get(wasora_value_ptr(milonga.vectors.chi), g) * material_xs->xs_values.nuSigmaF[g_prime] * milonga.moc_solver.phi[cell->index[g_prime]];
        //milonga.moc_solver.reduced_source[cell->index[g]] += (1 / wasora_var_value(milonga.vars.keff)) * gsl_vector_get(wasora_value_ptr(milonga.vectors.chi), g) * wasora_evaluate_expression(material_xs->nuSigmaF[g_prime]) * milonga.moc_solver.phi[cell->index[g_prime]];
      }
      
      //computo sigma total
      if (material_xs->SigmaT[g]->n_tokens != 0) {
        sigma_t = wasora_evaluate_expression(material_xs->SigmaT[g]);
      } else {
        sigma_t = wasora_evaluate_expression(material_xs->SigmaA[g]);
        for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
          sigma_t += wasora_evaluate_expression(material_xs->SigmaS0[g][g_prime]);
        }
      }
      
      milonga.moc_solver.reduced_source[cell->index[g]] /= (4.0 * M_PI * material_xs->xs_values.SigmaT[g]);
      //milonga.moc_solver.reduced_source[cell->index[g]] /= (4.0 * M_PI * sigma_t);
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_compute_phi(void) {
  
  int i, j;
  
  double *boundary_psi;
  
  track_t *track;
  segment_list_item_t *associated_segment;
  
  moc_volumes_set_uniform_phi(0.0);
  moc_volumes_update_boundary_psi();
  
  // loopeamos en todos los tracks
  for (i = 0; i < milonga.moc_solver.tracks->n_azim_2; i++) {
    for (j = 0; j < milonga.moc_solver.tracks->n_tracks[i]; j++) {
      
      track = &milonga.moc_solver.tracks->track[i][j];
      
      // nos movemos en direccion fwd del track
      // tomamos el flujo angular en la frontera y en la direccion fwd (d = 0)
      boundary_psi = &milonga.moc_solver.boundary_psi[angular_index(track->id,0,0,0)];
      // sumamos las contribuciones al flujo escalar recorriendo cada track en sentido fwd
      LL_FOREACH(track->associated_segments, associated_segment) {
        wasora_call(moc_volumes_tally_phi(associated_segment->segment, boundary_psi, track->azim_index));
      }
      
      // segun la condicion de contorno, transfiero los flujos angulares fwd obtenidos
      // (para cada indice polar y grupo) al correspondiente track
      wasora_call(moc_volumes_set_start_boundary_psi(track, boundary_psi, 1));
      
      // nos movemos en direccion bwd del track (ahora usar n_azim_2 en el loop "i" se traduce a n_azim)
      // tomamos el flujo angular en la frontera y en la direccion bwd (d = 1)
      boundary_psi = &milonga.moc_solver.boundary_psi[angular_index(track->id,1,0,0)];
      // primero revierto el sentido de linked list con los segmentos para movernos en sentido bwd
      LL_SORT(track->associated_segments, wasora_dummy_true_cmp);
      // sumamos las contribuciones al flujo escalar recorriendo cada track en sentido bwd
      LL_FOREACH(track->associated_segments, associated_segment) {
        wasora_call(moc_volumes_tally_phi(associated_segment->segment, boundary_psi, track->azim_index));
      }
      
      // acomodo nuevamente el orden de los segmentos
      LL_SORT(track->associated_segments, wasora_dummy_true_cmp);
      
      // segun la condicion de contorno, transfiero los flujos angulares bwd obtenidos
      // (para cada indice polar y grupo) al correspondiente track
      wasora_call(moc_volumes_set_start_boundary_psi(track, boundary_psi, 0));
    }
  }
  
  // y ahora falta sumar la reduced_source
  wasora_call(moc_volumes_add_source_to_phi());
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_tally_phi(segment_t *segment, double *boundary_psi, int azim_index) {
  
  int p, g;
  
  double exponential, delta_psi;
  
  for (g = 0; g < milonga.groups; g++) {
    for (p = 0; p < milonga.moc_solver.tracks->quadrature->polar->n_polar_2; p++) {
      
      exponential = milonga.moc_solver.exp_evaluator->compute_exponential(segment->tau[g], p);
      delta_psi = (boundary_psi[reduced_angular_index(p,g)] - milonga.moc_solver.reduced_source[segment->element->cell->index[g]]) * exponential;
      // como solo barremos n_polar_2, tenemos que multiplicar por 2.0 para tener en cuenta esos terminos que faltan
      milonga.moc_solver.phi[segment->element->cell->index[g]] += 2.0 * milonga.moc_solver.tracks->quadrature->w_total[azim_index][p] * delta_psi;
      boundary_psi[reduced_angular_index(p,g)] -= delta_psi;
      
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_set_start_boundary_psi(track_t *current_track, double *boundary_psi, int fwd) {
  
  int p, g;
  int flag;
  int next_track_id;
  int next_track_dir;
  
  track_t *next_track;
  
  // si nos movemos fwd
  if (fwd) {
    next_track = current_track->next_track_fwd;
    next_track_dir = current_track->dir_next_track_fwd;
    // miramos si es vacio
    flag = !(current_track->boundary_end->bc_type_phys == BC_VACUUM || current_track->boundary_end->bc_type_phys == BC_NULL);
    
    // si nos movemos bwd
  } else {
    next_track = current_track->next_track_bwd;
    next_track_dir = current_track->dir_next_track_bwd;
    // miramos si es vacio
    flag = !(current_track->boundary_start->bc_type_phys == BC_VACUUM || current_track->boundary_start->bc_type_phys == BC_NULL);
  }
  
  next_track_id = next_track->id;
  
  for (p = 0; p < milonga.moc_solver.tracks->quadrature->polar->n_polar_2; p++) {
    for (g = 0; g < milonga.groups; g++) {
      milonga.moc_solver.start_boundary_psi[angular_index(next_track_id,next_track_dir,p,g)] = boundary_psi[reduced_angular_index(p,g)] * flag;
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_add_source_to_phi(void) {
  
  int i;
  int g; 
  //int g_prime;
  
  //double xi;
  
  cell_t *cell;
  xs_t *material_xs;
  
  for (i = 0; i < milonga.spatial_unknowns; i++) {
    
    cell = &wasora_mesh.main_mesh->cell[i];
    material_xs = (xs_t *) (cell->element->physical_entity->material->ext);
    
    if (cell->volume == 0) {
      cell->volume = cell->element->type->element_volume(cell->element);
    }
    
    for (g = 0; g < milonga.groups; g++) {
      
/*
      if (material_xs->SigmaT[g]->n_tokens != 0) {
        xi = moc_volumes_cell_integral(cell, material_xs->SigmaT[g]);
      } else {
        xi = moc_volumes_cell_integral(cell, material_xs->SigmaA[g]);
        for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
          xi += moc_volumes_cell_integral(cell, material_xs->SigmaS0[g][g_prime]);
        }
      }
*/
      
      milonga.moc_solver.phi[cell->index[g]] /= (material_xs->xs_values.SigmaT[g] * cell->volume);
      //milonga.moc_solver.phi[cell->index[g]] /= xi;
      milonga.moc_solver.phi[cell->index[g]] += 4.0 * M_PI * milonga.moc_solver.reduced_source[cell->index[g]];
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_compute_keff(void) {
  
  int i, g_prime;
  
  double total_fission_source;
  
  cell_t *cell;
  xs_t *material_xs;
  
  // conviene hacer sumas pairwise! hablarlo con jeremy
  total_fission_source = 0.0;
  for (i = 0; i < milonga.spatial_unknowns; i++) {
    
    cell = &wasora_mesh.main_mesh->cell[i];
    material_xs = (xs_t *) (cell->element->physical_entity->material->ext);
    
    for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
      total_fission_source += material_xs->xs_values.nuSigmaF[g_prime] * milonga.moc_solver.phi[cell->index[g_prime]] * cell->volume;
      //total_fission_source += moc_volumes_cell_integral(cell, material_xs->nuSigmaF[g_prime]) * milonga.moc_solver.phi[cell->index[g_prime]];
    }
  }
  
  // dado que el phi ya lo normalice antes (con lo que seria el phi_old de esta iteracion), con hacer esto alcanza
  wasora_var_value(milonga.vars.keff) *= total_fission_source;
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_compute_residual(double *residual) {
  
  // TODO: el residual lo hago respecto a fuentes totales, pero estaria bueno agregar
  // todas las fuentes y que se computen todos y todos se cumplan (fission y flux)
  
  int i;
  int g, g_prime;
  
  double nu_sigma_f, sigma_s;
  double new_qi;
  double old_qi;
  
  cell_t *cell;
  xs_t *material_xs;
  
  // calculamos la fuente total (integrada en energia g) para cada celda i (fuente new y old)
  *residual = 0;
  for (i = 0; i < milonga.spatial_unknowns; i++) {
    
    cell = &wasora_mesh.main_mesh->cell[i];
    material_xs = (xs_t *) (cell->element->physical_entity->material->ext);
    
    new_qi = 0;
    old_qi = 0;
    
    // fuente total de fision en cada region (los chi los integro en g y suman 1.0 asi q chau)
    for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
      nu_sigma_f = material_xs->xs_values.nuSigmaF[g];
      //nu_sigma_f = wasora_evaluate_expression(material_xs->nuSigmaF[g_prime]);
      new_qi += nu_sigma_f * milonga.moc_solver.phi[cell->index[g_prime]];
      old_qi += nu_sigma_f * milonga.moc_solver.old_phi[cell->index[g_prime]];
    }
    new_qi /= wasora_var_value(milonga.vars.keff);
    old_qi /= wasora_var_value(milonga.vars.keff);
    
    // fuente total de scattering
    for (g = 0; g < milonga.groups; g++) {
      for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
        sigma_s = material_xs->xs_values.SigmaS0[g_prime][g];
        //sigma_s = wasora_evaluate_expression(material_xs->SigmaS0[g_prime][g]);
        new_qi += sigma_s * milonga.moc_solver.phi[cell->index[g_prime]];
        old_qi += sigma_s * milonga.moc_solver.old_phi[cell->index[g_prime]];
      }
    }
    
    // puede ocurrir que no haya fuente en la celda y que no sea vacio (?)
    if (old_qi > 0)
      *residual += pow(((new_qi - old_qi) / old_qi), 2);
  }
  
  *residual = sqrt((*residual) / milonga.spatial_unknowns);
  
  return WASORA_RUNTIME_OK;
}

int moc_volumes_results_fill_flux(void) {

  int i, g;

  for (i = 0; i < wasora_mesh.main_mesh->n_cells; i++) {
    for (g = 0; g < milonga.groups; g++) {
      milonga.functions.phi[g]->data_value[i] = milonga.moc_solver.phi[wasora_mesh.main_mesh->cell[i].index[g]];
    }
  }

  return WASORA_RUNTIME_OK;
}

int moc_volumes_normalize_flux(void) {

  int i, g;
  double factor;
  double num = 0;
  double den = 0;

  if (wasora_var(milonga.vars.power) == 0) {
    // calculamos el factor de normalizacion factor = num/den
    for (i = 0; i < wasora_mesh.main_mesh->n_cells; i++) {
      if (wasora_mesh.main_mesh->cell[i].volume == 0) {
        wasora_mesh.main_mesh->cell[i].volume = wasora_mesh.main_mesh->cell[i].element->type->element_volume(wasora_mesh.main_mesh->cell[i].element);
      }
      num += wasora_mesh.main_mesh->cell[i].volume;
      for (g = 0; g < milonga.groups; g++) {
        milonga.functions.phi[g]->data_value[i] = milonga.moc_solver.phi[wasora_mesh.main_mesh->cell[i].index[g]];
        den += wasora_mesh.main_mesh->cell[i].volume * milonga.functions.phi[g]->data_value[i];
      }
    }

  } else {

    xs_t *xs;

    num = wasora_var(milonga.vars.power);
    for (i = 0; i < wasora_mesh.main_mesh->n_cells; i++) {
      if (wasora_mesh.main_mesh->cell[i].element->physical_entity != NULL &&
          wasora_mesh.main_mesh->cell[i].element->physical_entity->material != NULL &&
          (xs = (xs_t *)wasora_mesh.main_mesh->cell[i].element->physical_entity->material->ext) != NULL) {

        for (g = 0; g < milonga.groups; g++) {
          milonga.functions.phi[g]->data_value[i] = milonga.moc_solver.phi[wasora_mesh.main_mesh->cell[i].index[g]];
          den += moc_volumes_cell_integral(&wasora_mesh.main_mesh->cell[i], xs->eSigmaF[g]) * milonga.functions.phi[g]->data_value[i];
        }
      }
    }

    if (den == 0) {
      wasora_push_error_message("power setpoint was given but eSigmaF is identically zero");
      return WASORA_RUNTIME_ERROR;
    }
  }


  factor = num/den;

  for (g = 0; g < milonga.groups; g++) {
    for (i = 0; i < milonga.spatial_unknowns; i++) {
      milonga.functions.phi[g]->data_value[i] *= factor;
    }
  }

  return WASORA_RUNTIME_OK;
}

int moc_volumes_results_fill_power(void) {

  int i, g;
  xs_t *xs;
  
  for (i = 0; i < wasora_mesh.main_mesh->n_cells; i++) {

    wasora_var(wasora_mesh.vars.x) = wasora_mesh.main_mesh->cell[i].x[0];
    wasora_var(wasora_mesh.vars.y) = wasora_mesh.main_mesh->cell[i].x[1];
    wasora_var(wasora_mesh.vars.z) = wasora_mesh.main_mesh->cell[i].x[2];
    
    milonga.functions.pow->data_value[i] = 0;

    if (wasora_mesh.main_mesh->cell[i].element->physical_entity != NULL &&
        wasora_mesh.main_mesh->cell[i].element->physical_entity->material != NULL &&
        (xs = (xs_t *)wasora_mesh.main_mesh->cell[i].element->physical_entity->material->ext) != NULL) {    

      for (g = 0; g < milonga.groups; g++) {
        milonga.functions.pow->data_value[i] += wasora_evaluate_expression(xs->eSigmaF[g]) * milonga.functions.phi[g]->data_value[i];
      }
    }
  }

  return WASORA_RUNTIME_OK;

}

double moc_volumes_cell_integral(cell_t *cell, expr_t *f) {

  if (f == NULL) {
    return 0;
  }

  // si todavia no calculamos el volumen de la celda, lo hacemos ahora
  if (cell->volume == 0) {
    cell->volume = cell->element->type->element_volume(cell->element);
  }

  wasora_value(wasora_mesh.vars.x) = cell->x[0];
  wasora_value(wasora_mesh.vars.y) = cell->x[1];
  wasora_value(wasora_mesh.vars.z) = cell->x[2];
  
  return wasora_evaluate_expression(f) * cell->volume;
}

int moc_volumes_problem_free(void) {
  
  int g;
  
  if (wasora_mesh.main_mesh != NULL && wasora_mesh.main_mesh->n_cells != 0) {
    if (milonga.functions.phi != NULL) {
      for (g = 0; g < milonga.groups; g++) {       
        free(milonga.functions.phi[g]->data_value);
        milonga.functions.phi[g]->data_argument = NULL;
        milonga.functions.phi[g]->data_value = NULL;
        milonga.functions.phi[g]->var_argument = NULL;
      }
      free(milonga.functions.pow->data_value);
      milonga.functions.pow->data_argument = NULL;
      milonga.functions.pow->data_value = NULL;
      milonga.functions.pow->var_argument = NULL;
    }

    mesh_free(wasora_mesh.main_mesh);
  }
   
  return WASORA_RUNTIME_OK;
}

int wasora_dummy_true_cmp(void *x, void *y) {
  return 1;
}

int wasora_dummy_equal_cmp(void *x, void *y) {
  return 0;
}

int wasora_dummy_false_cmp(void *x, void *y) {
  return -1;
}