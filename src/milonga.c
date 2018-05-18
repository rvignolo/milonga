/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga plugin for wasora
 *
 *  Copyright (C) 2010--2017 jeremy theler
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
#include <sys/time.h>
#include <sys/utsname.h>
#include <sys/resource.h>
#include <unistd.h>

#include "milonga.h"

#define time_checkpoint(which) \
  petsc_call(PetscTime(&wall.which)); \
  cpu.which = milonga_get_cpu_time();

#undef  __FUNCT__
#define __FUNCT__ "milonga_instruction_step"
int milonga_instruction_step(void *arg) {

  milonga_step_t *milonga_step = (milonga_step_t *)arg;
  milonga_times_t wall;
  milonga_times_t cpu;
  PetscStageLog     stageLog;  
  struct rusage resource_usage;


  double spectrum;
  int (*user_provided_eigensolver)(Mat, Mat, Vec, PetscScalar *);
  int (*user_provided_linearsolver)(Mat, Vec, Vec);
  int i;

  PetscFunctionBegin;

  if (!milonga.initialized) {
    petsc_call(PetscLogStagePush(milonga.petsc_stage_init));
    petsc_call(PetscLogEventBegin(milonga.petsc_event_init, 0, 0, 0, 0));
    time_checkpoint(init_begin);

    wasora_call(milonga.problem_init());
    wasora_call(milonga_read_boundaries());

    // chequeamos que el espectro de fision este normalizado
    spectrum = 0;
    for (i = 0; i < milonga.vectors.chi->size; i++) {
      spectrum += gsl_vector_get(wasora_value_ptr(milonga.vectors.chi), i);
    }
    if (gsl_fcmp(spectrum, 1.0, wasora_var(wasora_mesh.vars.eps)) != 0) {
      wasora_push_error_message("vector 'chi' with the fission spectrum is not normalized to one (sum of elements is %.6f)", spectrum);
      return WASORA_RUNTIME_ERROR;
    }

    wasora_var(milonga.vars.unknowns) = (double)milonga.problem_size;
    milonga.initialized = 1;

    petsc_call(PetscLogEventEnd(milonga.petsc_event_init, 0, 0, 0, 0));
    petsc_call(PetscLogStagePop());
    time_checkpoint(init_end);
  }
  

  // ------------------------------------
  // build
  // ------------------------------------
  if (milonga_step->do_not_build == 0) {
    petsc_call(PetscLogStagePush(milonga.petsc_stage_build));
    petsc_call(PetscLogEventBegin(milonga.petsc_event_build, 0, 0, 0, 0));
    time_checkpoint(build_begin);

    // esto es para problemas con static_steps != 0, si hacemos matsetzero usamos mucha memoria
    if (wasora_var_value(wasora.special_vars.step_static) != 1) {
      petsc_call(milonga_free_global_matrices());
      petsc_call(milonga_allocate_global_matrices(0, 0, 0));
    }
    wasora_call(milonga.matrices_build());
    
    PetscLogEventEnd(milonga.petsc_event_build, 0, 0, 0, 0);
    PetscLogStagePop();
    time_checkpoint(build_end);
  }

  // ------------------------------------
  // solve
  // ------------------------------------
  if (milonga_step->do_not_solve == 0) {
    petsc_call(PetscLogStagePush(milonga.petsc_stage_solve));
    petsc_call(PetscLogEventBegin(milonga.petsc_event_solve, 0, 0, 0, 0));
    time_checkpoint(solve_begin);
    
    if (milonga.has_fission != 0 && milonga.has_sources == 0) {
  
      // problema de autovalores
      if (milonga.user_provided_eigensolver != NULL) {
        user_provided_eigensolver = (void *)milonga.user_provided_eigensolver->routine;   // para cambiar los argumentos
        wasora_call(user_provided_eigensolver(milonga.R, milonga.F, milonga.phi, wasora_value_ptr(milonga.vars.keff)));
      } else {
        wasora_call(milonga_solve_eigen_slepc(milonga.R, milonga.F, milonga.phi, wasora_value_ptr(milonga.vars.keff)));
      }
      
    
      // normalizamos el autovector
      wasora_call(milonga.results_fill_flux());
      wasora_call(milonga.normalize_flux());
    
    } else if (milonga.has_sources != 0) {

      // problema lineal
      // calculamos una nueva R = (R-F)
      if (milonga.has_fission) {
        MatAXPY(milonga.R, -1.0, milonga.F, SUBSET_NONZERO_PATTERN);
      }

      if (milonga.user_provided_linearsolver != NULL) {
        user_provided_linearsolver = (void *)milonga.user_provided_linearsolver->routine;   // para cambiar los argumentos
        wasora_call(user_provided_linearsolver(milonga.R, milonga.S, milonga.phi));
      } else {
        wasora_call(milonga_solve_linear_petsc(milonga.R, milonga.S, milonga.phi));
      }
      wasora_call(milonga.results_fill_flux());
    }
  
    wasora_call(milonga.results_fill_power());

  
    // tomamos el resultado final como el guess inicial del paso siguiente
  /*
    VecCopy(problem.phi, problem.guess);
    VecCopy(problem.guess, problem.phi);
    VecAssemblyBegin(problem.guess);
    VecAssemblyEnd(problem.guess);
  */
  
    petsc_call(PetscLogEventEnd(milonga.petsc_event_solve, 0, 0, 0, 0));
    petsc_call(PetscLogStagePop());
    time_checkpoint(solve_end);
  }


/*
 * matt dijo que habia que hacer esto pero no camino
  petsc_call(PetscLogEventGetPerfInfo(milonga.petsc_stage_init, milonga.petsc_event_init, &info_init));
  petsc_call(PetscLogEventGetPerfInfo(milonga.petsc_stage_build, milonga.petsc_event_build, &info_build));
  petsc_call(PetscLogEventGetPerfInfo(milonga.petsc_stage_solve, milonga.petsc_event_solve, &info_solve));
 */

  PetscLogGetStageLog(&stageLog);
  
  
  
  wasora_var(milonga.vars.time_petsc_ini) = stageLog->stageInfo[milonga.petsc_stage_init].perfInfo.time;
  wasora_var(milonga.vars.time_wall_ini)  += wall.init_end  - wall.init_begin;
  wasora_var(milonga.vars.time_cpu_ini)   += cpu.init_end   - cpu.init_begin;

  if (milonga_step->do_not_build == 0) {
    wasora_var(milonga.vars.time_petsc_build) = stageLog->stageInfo[milonga.petsc_stage_build].perfInfo.time;
    wasora_var(milonga.vars.time_wall_build)  += wall.build_end  - wall.build_begin;
    wasora_var(milonga.vars.time_cpu_build)   += cpu.build_end   - cpu.build_begin;
  }
  
  if (milonga_step->do_not_solve == 0) {
    wasora_var(milonga.vars.time_petsc_solve) = stageLog->stageInfo[milonga.petsc_stage_solve].perfInfo.time;
    wasora_var(milonga.vars.time_wall_solve)  += wall.solve_end  - wall.solve_begin;
    wasora_var(milonga.vars.time_cpu_solve)   += cpu.solve_end   - cpu.solve_begin;
  }
  
  wasora_var(milonga.vars.time_petsc_total) = wasora_var(milonga.vars.time_petsc_ini) + wasora_var(milonga.vars.time_petsc_build) + wasora_var(milonga.vars.time_petsc_solve);
  wasora_var(milonga.vars.time_wall_total)  = wasora_var(milonga.vars.time_wall_ini)  + wasora_var(milonga.vars.time_wall_build)  + wasora_var(milonga.vars.time_wall_solve);
  wasora_var(milonga.vars.time_cpu_total)   = wasora_var(milonga.vars.time_cpu_ini)   + wasora_var(milonga.vars.time_cpu_build)   + wasora_var(milonga.vars.time_cpu_solve);

  wasora_value(milonga.vars.available_memory) = sysconf(_SC_PHYS_PAGES)*sysconf(_SC_PAGESIZE);
  getrusage(RUSAGE_SELF, &resource_usage);
  wasora_value(milonga.vars.memory_usage_global) = (double)(1024.0*resource_usage.ru_maxrss);
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

int milonga_instruction_moc_step(void *arg) {
    
  milonga_step_t *milonga_step = (milonga_step_t *)arg;
  int i;
  
  double spectrum;
  
  if (!milonga.initialized) {
    
    wasora_call(milonga.problem_init());
    wasora_call(milonga_read_boundaries());
    
    // dependiendo de las CC del problema, unimos los tracks
    wasora_call(tracks_set_tracks_boundary_conditions(milonga.moc_solver.tracks));
    wasora_call(tracks_set_tracks_next_tracks(milonga.moc_solver.tracks));
    
    // chequeamos que el espectro de fision este normalizado
    spectrum = 0;
    for (i = 0; i < milonga.vectors.chi->size; i++) {
      spectrum += gsl_vector_get(wasora_value_ptr(milonga.vectors.chi), i);
    }
    if (gsl_fcmp(spectrum, 1.0, wasora_var(wasora_mesh.vars.eps)) != 0) {
      wasora_push_error_message("vector 'chi' with the fission spectrum is not normalized to one (sum of elements is %.6f)", spectrum);
      return WASORA_RUNTIME_ERROR;
    }
    
    wasora_var(milonga.vars.unknowns) = (double)milonga.problem_size;
    milonga.initialized = 1;
  }
  
  // rellenamos las xs
  wasora_call(milonga_resolve_xs_values());
  
  // problema de autovalores
  if (milonga.has_fission != 0 && milonga.has_sources == 0) {
    
    // mandamos a calcular
    wasora_call(moc_volumes_solve_eigen_problem(milonga_step->verbose));
    
    // normalizamos el autovector
    wasora_call(milonga.results_fill_flux());
    wasora_call(milonga.normalize_flux());
    
  // problema lineal
  }  else if (milonga.has_sources != 0) {
    wasora_push_error_message("linear problem not yet supported for moc formulation");
    return WASORA_RUNTIME_ERROR;
  } else if (milonga.has_fission == 0 && milonga.has_sources == 0) {
    wasora_push_error_message("sources (independent & fission) are identically zero through the domain");
    return WASORA_RUNTIME_ERROR;
  }
  
  wasora_call(milonga.results_fill_power());
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "milonga_resolve_xs_expr"
void milonga_resolve_xs_expr(material_t *material, char *xsname, expr_t **expr, int g, int g_prime) {

  property_data_t *property_data;
  char property_name[32];

  // empezamos con cero
  *expr = &milonga.xs_zero;

  if (g_prime == -1) {

    // es una XS simple (i.e. no scattering)
    // si tenemos un solo grupo intentamos sin subindice
    if (milonga.groups == 1) {
      sprintf(property_name, "%s", xsname);
      HASH_FIND_STR(material->property_datums, property_name, property_data);
      if (property_data != NULL) {
        *expr = &property_data->expr;
        return;
      }
    }

    // probamos sin underscore
    sprintf(property_name, "%s%d", xsname, g+1);
    HASH_FIND_STR(material->property_datums, property_name, property_data);
    if (property_data != NULL) {
      *expr = &property_data->expr;
      return;
    }

    // probamos con underscore
    sprintf(property_name, "%s_%d", xsname, g+1);
    HASH_FIND_STR(material->property_datums, property_name, property_data);
    if (property_data != NULL) {
      *expr = &property_data->expr;
      return;
    }

  } else {

    // es una XS compuesta (scattering)

    // probamos sin nada  
    if (milonga.groups == 1) {
      sprintf(property_name, "%s", xsname);
      HASH_FIND_STR(material->property_datums, property_name, property_data);
      if (property_data != NULL) {
        *expr = &property_data->expr;
        return;
      }
    }
      
    // probamos sin underscore con punto
    sprintf(property_name, "%s%d.%d", xsname, g+1, g_prime+1);
    HASH_FIND_STR(material->property_datums, property_name, property_data);
    if (property_data != NULL) {
      *expr = &property_data->expr;
      return;
    }

    // probamos sin underscore con flecha
    // no vale tener un menos ni un mayor en el nombre!
/*    
    sprintf(property_name, "%s%d->%d", xsname, g+1, g_prime+1);
    HASH_FIND_STR(material->property_datums, property_name, property_data);
    if (property_data != NULL) {
      *expr = &property_data->expr;
      return;
    }
 */

    // probamos con underscore con punto
    sprintf(property_name, "%s_%d.%d", xsname, g+1, g_prime+1);
    HASH_FIND_STR(material->property_datums, property_name, property_data);
    if (property_data != NULL) {
      *expr = &property_data->expr;
      return;
    }

    // probamos con underscore con flecha
    // no vale tener un menos ni un mayor en el nombre!
/*    
    sprintf(property_name, "%s_%d->%d", xsname, g+1, g_prime+1);
    HASH_FIND_STR(material->property_datums, property_name, property_data);
    if (property_data != NULL) {
      *expr = &property_data->expr;
      return;
    }
*/
  }

  return;


}

int milonga_assembly_objects(MatAssemblyType type) {

  // TODO: que orden es mejor?
  
  MatAssemblyBegin(milonga.R, type);
  if (milonga.has_fission) {
    MatAssemblyBegin(milonga.F, type);
  }
  if (milonga.has_sources) {
    VecAssemblyBegin(milonga.S);
  }
  
  MatAssemblyEnd(milonga.R, type);
  if (milonga.has_fission) {
    MatAssemblyEnd(milonga.F, type);
  }
  if (milonga.has_sources) {
    VecAssemblyEnd(milonga.S);
  }
  

  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "milonga_resolve_xs_values"
int milonga_resolve_xs_values(void) {
  
  int i;
  int g, g_prime;
  
  cell_t *cell;
  xs_t *material_xs;
  
  for (i = 0; i < milonga.spatial_unknowns; i++) {
    
    cell = &wasora_mesh.main_mesh->cell[i];
    material_xs = (xs_t *)(cell->element->physical_entity->material->ext);
    
    wasora_value(wasora_mesh.vars.x) = cell->x[0];
    wasora_value(wasora_mesh.vars.y) = cell->x[1];
    wasora_value(wasora_mesh.vars.z) = cell->x[2];
    
    for (g = 0; g < milonga.groups; g++) {
      
      // SigmaT_g
      if (material_xs->SigmaT[g]->n_tokens != 0) {
        material_xs->xs_values.SigmaT[g] = wasora_evaluate_expression(material_xs->SigmaT[g]);
      } else {
        material_xs->xs_values.SigmaT[g] = wasora_evaluate_expression(material_xs->SigmaA[g]);
        for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
          material_xs->xs_values.SigmaT[g] += wasora_evaluate_expression(material_xs->SigmaS0[g][g_prime]);
        }
      }
      
      // SigmaA_g
      if (material_xs->SigmaA[g]->n_tokens != 0) {
        material_xs->xs_values.SigmaA[g] = wasora_evaluate_expression(material_xs->SigmaA[g]);
      } 
      
      // nuSigmaF_g
      if (material_xs->nuSigmaF[g]->n_tokens != 0) {
        if ((material_xs->xs_values.nuSigmaF[g] = wasora_evaluate_expression(material_xs->nuSigmaF[g])) != 0) milonga.has_fission = 1;
      }
      
      // eSigmaF_g
      if (material_xs->eSigmaF[g]->n_tokens != 0) {
        material_xs->xs_values.eSigmaF[g] = wasora_evaluate_expression(material_xs->eSigmaF[g]);
      }
      
      // SigmaS0_g.g_prime
      for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
        if (material_xs->SigmaS0[g][g_prime]->n_tokens != 0) {
          material_xs->xs_values.SigmaS0[g][g_prime] = wasora_evaluate_expression(material_xs->SigmaS0[g][g_prime]);
        }
      }
      
      // S_g
      if (material_xs->S[g]->n_tokens != 0) {
        if ((material_xs->xs_values.S[g] = wasora_evaluate_expression(material_xs->S[g])) != 0) milonga.has_sources = 1;
      }
      
    }
  }
  
  return WASORA_RUNTIME_OK;
}
