/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's parsing routines
 *
 *  Copyright (C) 2010--2015 jeremy theler
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
#include "milonga.h"

#undef  __FUNCT__
#define __FUNCT__ "plugin_parse_line"
int plugin_parse_line(char *line) {
  
  char *token;

  if ((token = wasora_get_next_token(line)) != NULL) {

// ---------------------------------------------------------------------
///kw+MILONGA_SOLVER+usage MILONGA_SOLVER
///kw+MILONGA_SOLVER+desc Sets options related to the eigen-solver.
///kw+MILONGA_SOLVER+detail 
    if (strcasecmp(token, "MILONGA_SOLVER") == 0) {

      while ((token = wasora_get_next_token(NULL)) != NULL) {

///kw+MILONGA_SOLVER+usage [ ROUTINE <loadable_routine> ]
        if (strcasecmp(token, "ROUTINE") == 0) {
          if ((token = wasora_get_next_token(NULL)) == NULL) {
            wasora_push_error_message("expected solver name");
            return WASORA_PARSER_ERROR;
          }

          if ((milonga.user_provided_eigensolver = wasora_get_loadable_routine(token)) == NULL) {
            wasora_push_error_message("unknown routine '%s'", token);
            return WASORA_PARSER_ERROR;
          }

///kw+MILONGA_SOLVER+usage [ SPECTRUM { largest_eigenvalue | smallest_eigenvalue } ]
        } else if (strcasecmp(token, "SPECTRUM") == 0) {
          char *keywords[] = {"largest_eigenvalue", "smallest_eigenvalue", ""};
          int values[] = {spectrum_largest_eigenvalue, spectrum_smallest_eigenvalue, 0};
          wasora_call(wasora_parser_keywords_ints(keywords, values, (int *)&milonga.spectrum));
          
///kw+MILONGA_SOLVER+usage [ EPS_TYPE { krylovschur | gd | jd | power | arnoldi | subspace | ... } ]
///kw+MILONGA_SOLVER+detail List of `EPS_TYPE`s <http://www.grycap.upv.es/slepc/documentation/current/docs/manualpages/EPS/EPSType.html>
///kw+MILONGA_SOLVER+detail          
        } else if (strcasecmp(token, "EPS_TYPE") == 0) {
          wasora_call(wasora_parser_string(&milonga.eps_type));

///kw+MILONGA_SOLVER+usage [ ST_TYPE { sinvert | shift | cayley | precond } ]
///kw+MILONGA_SOLVER+detail List of `ST_TYPE`s <http://www.grycap.upv.es/slepc/documentation/current/docs/manualpages/ST/STType.html>
///kw+MILONGA_SOLVER+detail          
        } else if (strcasecmp(token, "ST_TYPE") == 0) {
          wasora_call(wasora_parser_string(&milonga.st_type));

///kw+MILONGA_SOLVER+usage [ KSP_TYPE { gmres | bcgs | bicg | richardson | chebyshev | ... } ]
///kw+MILONGA_SOLVER+detail List of `KSP_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html>
///kw+MILONGA_SOLVER+detail          
        } else if (strcasecmp(token, "KSP_TYPE") == 0) {
          wasora_call(wasora_parser_string(&milonga.ksp_type));

///kw+MILONGA_SOLVER+usage [ PC_TYPE { lu | none | sor | bjacobi | cholesky | ... } ]
///kw+MILONGA_SOLVER+detail List of `PC_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html>
///kw+MILONGA_SOLVER+detail          
        } else if (strcasecmp(token, "PC_TYPE") == 0) {
          wasora_call(wasora_parser_string(&milonga.pc_type));

///kw+MILONGA_SOLVER+usage [ SUBSPACE_DIM <expr> ]
        } else if (strcasecmp(token, "SUBSPACE_DIM") == 0) {
          wasora_call(wasora_parser_expression(&milonga.eps_ncv));

///kw+MILONGA_SOLVER+usage [ ST_SHIFT <expr> ]
        } else if (strcasecmp(token, "ST_SHIFT") == 0) {
          wasora_call(wasora_parser_expression(&milonga.st_shift));
          
///kw+MILONGA_SOLVER+usage [ ST_ANTI_SHIFT <expr> ]
        } else if (strcasecmp(token, "ST_ANTI_SHIFT") == 0) {
          wasora_call(wasora_parser_expression(&milonga.st_anti_shift));

        } else {
          wasora_push_error_message("undefined keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }
      }

      return WASORA_PARSER_OK;

// ---------------------------------------------------------------------
///kw+MILONGA_PROBLEM+usage MILONGA_PROBLEM
///kw+MILONGA_PROBLEM+desc Defines the number of spatial dimensions and groups of neutron energies.      
///kw+MILONGA_PROBLEM+desc It also selects the formulation of the neutronic problem to be solved (i.e. diffusion or tranport)
///kw+MILONGA_PROBLEM+desc and the spatial discretization scheme (i.e. finite volumes or finite elements).
///kw+MILONGA_PROBLEM+desc If several meshes are defined, it selects over which one it is that the neutronic problem is solved.      
    } else if (strcasecmp(token, "MILONGA_PROBLEM") == 0) {
      
      double xi;
      
      while ((token = wasora_get_next_token(NULL)) != NULL) {

///kw+MILONGA_PROBLEM+usage [ DIMENSIONS <expr> ]
        if (strcasecmp(token, "DIMENSIONS") == 0) {
          wasora_call(wasora_parser_expression_in_string(&xi));
          milonga.dimensions = (int)(xi);
          if (milonga.dimensions < 1 || milonga.dimensions > 3)  {
            wasora_push_error_message("either one, two or three dimensions should be selected instead of '%d'", milonga.dimensions);
            return WASORA_PARSER_ERROR;
          }
        
///kw+MILONGA_PROBLEM+usage [ GROUPS <expr> ]
        } else if (strcasecmp(token, "GROUPS") == 0) {
          wasora_call(wasora_parser_expression_in_string(&xi));
          milonga.groups = (int)(xi);
          if (milonga.groups < 1)  {
            wasora_push_error_message("a positive number of groups should be given instead of '%d'", milonga.groups);
            return WASORA_PARSER_ERROR;
          }

///kw+MILONGA_PROBLEM+usage [ MESH <identifier> ]
/*
        } else if (strcasecmp(token, "MESH") == 0) {
          char *mesh_name;
          
          wasora_call(wasora_parser_string(&mesh_name));
          if ((milonga.mesh = wasora_get_mesh_ptr(mesh_name)) == NULL) {
            wasora_push_error_message("unknown mesh '%s'", mesh_name);
            free(mesh_name);
            return WASORA_PARSER_ERROR;
          }
          free(mesh_name);
 */

///kw+MILONGA_PROBLEM+usage [ SCHEME { volumes | elements } ]
        } else if (strcasecmp(token, "SCHEME") == 0) {
          char *keywords[] = {"volumes", "elements", ""};
          int values[] = {scheme_volumes, scheme_elements, 0};
          wasora_call(wasora_parser_keywords_ints(keywords, values, (int *)&milonga.scheme));

///kw+MILONGA_PROBLEM+usage [ FORMULATION { diffusion | s2 | s4 | s6 | s8 | moc} ]
        } else if (strcasecmp(token, "FORMULATION") == 0) {
          
          char *formulation;
          wasora_call(wasora_parser_string(&formulation));
          
          // si milonga.SN == 0 entonces le creemos al input
          // si milonga.SN != 0 entonces le creemos a la linea de comandos
          if (milonga.SN == 0) {
            if (strcasecmp(formulation, "diffusion") == 0) {
              milonga.formulation = formulation_diffusion;
            } else if (strcasecmp(formulation, "s2") == 0) {
              milonga.formulation = formulation_sn;
              milonga.SN = 2;
            } else if (strcasecmp(formulation, "s4") == 0) {
              milonga.formulation = formulation_sn;
              milonga.SN = 4;
            } else if (strcasecmp(formulation, "s6") == 0) {
              milonga.formulation = formulation_sn;
              milonga.SN = 6;
            } else if (strcasecmp(formulation, "s8") == 0) {
              milonga.formulation = formulation_sn;
              milonga.SN = 8;
            } else if (strcasecmp(formulation, "moc") == 0) {
              milonga.formulation = formulation_moc;
            } else {
              wasora_push_error_message("unknown formulation '%s' (either diffusion, s2, s4, s6, s8 or moc)", formulation);
              free(formulation);
              return WASORA_PARSER_ERROR;
            }
          }
          
          free(formulation);
          
///kw+MILONGA_PROBLEM+usage [ TRACK <identifier> ]
        } else if (strcasecmp(token, "TRACK") == 0) {
          
          char *track_name;
          
          wasora_call(wasora_parser_string(&track_name));
          if ((milonga.moc_solver.tracks = milonga_get_ray_tracing_ptr(track_name)) == NULL) {
            wasora_push_error_message("unknown track '%s'", track_name);
            free(track_name);
            return WASORA_PARSER_ERROR;
          }
          free(track_name);
          
///kw+MILONGA_PROBLEM+usage [ POLAR_QUADRATURE <identifier> ]
        } else if (strcasecmp(token, "POLAR_QUADRATURE") == 0) {
          
          char *polar_quad_name;
          
          wasora_call(wasora_parser_string(&polar_quad_name));
          if ((milonga.moc_solver.polar_quadrature = milonga_get_polar_quadrature_ptr(polar_quad_name)) == NULL) {
            wasora_push_error_message("unknown polar quadrature '%s'", polar_quad_name);
            free(polar_quad_name);
            return WASORA_PARSER_ERROR;
          }
          free(polar_quad_name);
          
///kw+MILONGA_PROBLEM+usage [ DO_NOT_LINEARIZE_EXPONENTIALS ]
        } else if (strcasecmp(token, "DO_NOT_LINEARIZE_EXPONENTIALS") == 0) {
          milonga.moc_solver.do_not_linearize_exp = 1;
          
///kw+MILONGA_PROBLEM+usage [ VOLHOM ]
        } else if (strcasecmp(token, "VOLHOM") == 0) {
          milonga.volhom = 1;
          
        } else {
          wasora_push_error_message("undefined keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }
      }

      // si no nos dieron explicitamente la malla, ponemos la principal
      if (wasora_mesh.main_mesh == NULL) {
        wasora_push_error_message("unknown mesh for MILONGA_PROBLEM (no MESH keyword)", token);
        return WASORA_PARSER_ERROR;
      }
      
      // por si ya nos dieron mesh, usamos las dimensiones de la malla si no nos las dieron aca
      if (milonga.dimensions == 0 && wasora_mesh.main_mesh->bulk_dimensions != 0) {
        milonga.dimensions = wasora_mesh.main_mesh->bulk_dimensions;
      }
      // al reves, si ya nos la dieron, se las damos a todas las mallas
      // TOOD: esto no me termina de convencer
      if (wasora_mesh.main_mesh != NULL && wasora_mesh.main_mesh->bulk_dimensions == 0 && milonga.dimensions != 0) {
        mesh_t *mesh;
        for (mesh = wasora_mesh.main_mesh; mesh != NULL; mesh = mesh->hh.next) {
          if (mesh->bulk_dimensions == 0) {
            mesh->bulk_dimensions = milonga.dimensions;
          }
        }
      }
      
      // si estamos en formulacion moc
      if (milonga.formulation == formulation_moc) {
        
        // TODO: corregir las condiciones con los nuevos criterios
        
        // si no nos dieron los tracks
        if (milonga.moc_solver.tracks == NULL) {
          if (milonga_tracks.main_tracks != NULL) {
            milonga.moc_solver.tracks = milonga_tracks.main_tracks;
          } else {
            wasora_push_error_message("unknown TRACK for MILONGA_PROBLEM", token);
            return WASORA_PARSER_ERROR;
          }
        }
        
        // si no nos dieron la cuadratura polar, ponemos la principal, aunque si tampoco hay
        // ninguna principal no nos preocupamos porque mas adelante metemos una por defecto
        if (milonga.moc_solver.polar_quadrature == NULL) {
          if (milonga_polar_quadratures.main_polar_quadrature != NULL) {
            milonga.moc_solver.polar_quadrature = milonga_polar_quadratures.main_polar_quadrature;
          }
        }
      }
      
      wasora_call(milonga_define_result_functions());
      
      return WASORA_PARSER_OK;
      
// ---------------------------------------------------------------------
///kw+TRACK_MESH+usage TRACK_MESH
///kw+TRACK_MESH+desc Performs a ray tracing (tracking) over a mesh to later implement the Method of Characteristics.
///kw+TRACK_MESH+desc The number of azimuthal angles must be specified as a multiple of 4, as well as the track density or azimutal spacing.
///kw+TRACK_MESH+desc If several meshes are defined, it selects over which one it is that the tracking is performed.
    } else if (strcasecmp(token, "TRACK_MESH") == 0) {
      
      tracks_t *ray_tracing;
      char *name = NULL;
      mesh_t *mesh = NULL;
      expr_t *expr_n_azim = calloc(1, sizeof(expr_t));
      expr_t *expr_track_dens = calloc(1, sizeof(expr_t));
      expr_t *expr_track_spacing = calloc(1, sizeof(expr_t));
      expr_t *expr_tiny_step = calloc(1, sizeof(expr_t));
      int do_not_correct_volumes = 0;
      int debug = 0;
      
      while ((token = wasora_get_next_token(NULL)) != NULL) {
        
///kw+TRACK_MESH+usage [ NAME <expr> ]
        if (strcasecmp(token, "NAME") == 0) {
          
          wasora_call(wasora_parser_string(&name));
          
///kw+TRACK_MESH+usage { N_AZIM_ANGLES <expr> }
        } else if (strcasecmp(token, "N_AZIM_ANGLES") == 0 || strcasecmp(token, "N_AZIM") == 0){
          
          wasora_call(wasora_parser_expression(expr_n_azim));
          
///kw+TRACK_MESH+usage { TRACK_DENS <expr> |
        } else if (strcasecmp(token, "TRACK_DENS") == 0){
          
          wasora_call(wasora_parser_expression(expr_track_dens));
///kw+TRACK_MESH+usage TRACK_SPACING <expr> }
        } else if (strcasecmp(token, "TRACK_SPACING") == 0) {
          
          wasora_call(wasora_parser_expression(expr_track_spacing));
          
///kw+TRACK_MESH+usage [ TINY_STEP <expr> ]
        } else if (strcasecmp(token, "TINY_STEP") == 0) {
          
          wasora_call(wasora_parser_expression(expr_tiny_step));
          
///kw+TRACK_MESH+usage [ DO_NOT_CORRECT_VOLUMES ]
        } else if (strcasecmp(token, "DO_NOT_CORRECT_VOLUMES") == 0) {
          
          do_not_correct_volumes = 1;
          
///kw+TRACK_MESH+usage [ DEBUG ]
        } else if (strcasecmp(token, "DEBUG") == 0) {
          
          debug = 1;
          
///kw+TRACK_MESH+usage { MESH <identifier> }
        } else if (strcasecmp(token, "MESH") == 0) {
          char *mesh_name;
          
          wasora_call(wasora_parser_string(&mesh_name));
          if ((mesh = wasora_get_mesh_ptr(mesh_name)) == NULL) {
            wasora_push_error_message("unknown mesh '%s'", mesh_name);
            free(mesh_name);
            return WASORA_PARSER_ERROR;
          }
          free(mesh_name);
          
        } else {
          wasora_push_error_message("undefined keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }
      }
      
      // si no nos dieron el numero de angulos azimutales
      if (expr_n_azim->n_tokens == 0) {
        wasora_push_error_message("N_AZIM_ANGLES keyword should be specified");
        return WASORA_PARSER_ERROR;
      }
      
      //si nos dieron TRACK_SPACING y TRACK_DENS
      if (expr_track_spacing->n_tokens != 0 && expr_track_dens->n_tokens != 0) {
        wasora_push_error_message("both TRACK_SPACING and TRACK_DENS keywords specified");
        return WASORA_PARSER_ERROR;
        //o si no nos dieron ninguno
      } else if (expr_track_spacing->n_tokens == 0 && expr_track_dens->n_tokens == 0) {
        wasora_push_error_message("TRACK_SPACING or TRACK_DENS keyword should be specified");
        return WASORA_PARSER_ERROR;
      }  
      
      // si no nos dieron explicitamente la malla, ponemos la principal, pero si no hay, chau
      if (mesh == NULL && (mesh = wasora_mesh.main_mesh) == NULL) {
        wasora_push_error_message("unknown mesh for TRACK_MESH (no MESH keyword)", token);
        return WASORA_PARSER_ERROR;
      }
      
      if ((ray_tracing = milonga_define_ray_tracing(name, mesh, expr_n_azim, expr_track_dens, expr_track_spacing, expr_tiny_step, do_not_correct_volumes, debug)) == NULL) {
        return WASORA_PARSER_ERROR;
      }
      
      if (wasora_define_instruction(milonga_instruction_track, ray_tracing) == NULL) {
        return WASORA_PARSER_ERROR;
      }
      
      return WASORA_PARSER_OK;
      
// ---------------------------------------------------------------------
///kw+POLAR_QUADRATURE+usage POLAR_QUADRATURE
///kw+POLAR_QUADRATURE+desc Selects the number of polar angles and the polar quadrature type for the Method of Characteristics
    } else if (strcasecmp(token, "POLAR_QUADRATURE") == 0) {
      
      polar_quadrature_t *polar_quadrature;
      char *name = NULL;
      int polar_quad_type = 0;
      expr_t *expr_n_polar = calloc(1, sizeof(expr_t));
      
      while ((token = wasora_get_next_token(NULL)) != NULL) {
        
///kw+POLAR_QUADRATURE+usage [ NAME <expr> ]
        if (strcasecmp(token, "NAME") == 0) {
          
          wasora_call(wasora_parser_string(&name));
        
///kw+POLAR_QUADRATURE+usage { TYPE { equal_weight | equal_angle | gauss_legendre | leonard | tabuchi_yamamoto } }
        } else if (strcasecmp(token, "TYPE") == 0) {
          
          char *keywords[] = {"tabuchi_yamamoto", "equal_weight", "equal_angle", "gauss_legendre", "leonard", ""};
          int values[] = {tabuchi_yamamoto, equal_weight, equal_angle, gauss_legendre, leonard, 0};
          
          wasora_call(wasora_parser_keywords_ints(keywords, values, &polar_quad_type));
          
///kw+POLAR_QUADRATURE+usage { N_POLAR_ANGLES <expr> }
        } else if (strcasecmp(token, "N_POLAR_ANGLES") == 0 || strcasecmp(token, "N_POLAR") == 0){
          
          wasora_call(wasora_parser_expression(expr_n_polar));
          
        } else {
          wasora_push_error_message("undefined keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }
      }
      
      // si no nos dieron el numero de angulos polares
      if (expr_n_polar->n_tokens == 0) {
        wasora_push_error_message("N_POLAR_ANGLES keyword should be specified");
        return WASORA_PARSER_ERROR;
      }
      
      if ((polar_quadrature = milonga_define_polar_quadrature(name, expr_n_polar, polar_quad_type)) == NULL) {
        return WASORA_PARSER_ERROR;
      }
      
      if (wasora_define_instruction(milonga_instruction_polar_quadrature, polar_quadrature) == NULL) {
        return WASORA_PARSER_ERROR;
      }
      
      return WASORA_PARSER_OK;
// ---------------------------------------------------------------------
///kw+IMPLICIT_BC+usage IMPLICIT_BC
    } else if (strcasecmp(token, "IMPLICIT_BC") == 0) {


///kw+IMPLICIT_BC+usage { NONE | ALLOWED }
      char *keywords[] = {"NONE", "ALLOWED", ""};
      int values[] = {1, 0, 0};
      wasora_call(wasora_parser_keywords_ints(keywords, values, &milonga.implicit_bc_none));

      return WASORA_PARSER_OK;

// ---------------------------------------------------------------------
///kw+MILONGA_DEBUG+usage MILONGA_DEBUG
///kw+MILONGA_DEBUG+desc Generates debugging and benchmarking output and/or dumps the matrices into files or the screen.
    } else if ((strcasecmp(token, "MILONGA_DEBUG") == 0)) {
      
      debug_t *debug;
      debug = calloc(1, sizeof(debug_t));
      LL_APPEND(milonga.debugs, debug);

      while ((token = wasora_get_next_token(NULL)) != NULL) {
        
///kw+MILONGA_DEBUG+usage [ FILE <file_id> | 
        if (strcasecmp(token, "FILE") == 0) {
          wasora_call(wasora_parser_file(&debug->file));
          
///kw+MILONGA_DEBUG+usage [ FILE_PATH <file_path> ]
        } else if (strcasecmp(token, "FILE_PATH") == 0) {
            wasora_call(wasora_parser_file_path(&debug->file, "w"));
          
///kw+MILONGA_DEBUG+usage [ MATRICES_ASCII ]
        } else if (strcasecmp(token, "MATRICES_ASCII") == 0) {
          debug->matrices |= DEBUG_MATRICES_ASCII;
///kw+MILONGA_DEBUG+usage [ MATRICES_ASCII_STRUCTURE ]
        } else if (strcasecmp(token, "MATRICES_ASCII_STRUCTURE") == 0) {
          debug->matrices |= DEBUG_MATRICES_ASCII_STRUCT;
///kw+MILONGA_DEBUG+usage [ MATRICES_PETSC_BINARY ]
        } else if (strcasecmp(token, "MATRICES_PETSC_BINARY") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_BINARY;
///kw+MILONGA_DEBUG+usage [ MATRICES_PETSC_COMPRESSED_BINARY ]
        } else if (strcasecmp(token, "MATRICES_PETSC_COMPRESSED_BINARY") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_COMPRESSED_BINARY;
///kw+MILONGA_DEBUG+usage [ MATRICES_PETSC_ASCII ]
        } else if (strcasecmp(token, "MATRICES_PETSC_ASCII") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_ASCII;
///kw+MILONGA_DEBUG+usage [ MATRICES_PETSC_OCTAVE ]
        } else if (strcasecmp(token, "MATRICES_PETSC_OCTAVE") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_OCTAVE;
///kw+MILONGA_DEBUG+usage [ MATRICES_PETSC_DENSE ]
        } else if (strcasecmp(token, "MATRICES_PETSC_DENSE") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_DENSE;
///kw+MILONGA_DEBUG+usage [ MATRICES_X ]
        } else if (strcasecmp(token, "MATRICES_X") == 0) {
          debug->matrices |= DEBUG_MATRICES_X;
///kw+MILONGA_DEBUG+usage [ MATRICES_SNG ]
        } else if (strcasecmp(token, "MATRICES_SNG") == 0) {
          debug->matrices |= DEBUG_MATRICES_SNG;
///kw+MILONGA_DEBUG+usage [ MATRICES_SNG_STRUCT ]
        } else if (strcasecmp(token, "MATRICES_SNG_STRUCT") == 0) {
          debug->matrices |= DEBUG_MATRICES_SNG_STRUCT;
          
///kw+MILONGA_DEBUG+usage [ MATRICES_SIZE <expr> ]
        } else if (strcasecmp(token, "MATRICES_SIZE") == 0 || strcasecmp(token, "MATRICES_X_SIZE") == 0) {
          wasora_call(wasora_parser_expression(&debug->matrices_size));
          
///kw+MILONGA_DEBUG+usage [ MATRICES_STRIDE <expr> ]
        } else if (strcasecmp(token, "MATRICES_STRIDE") == 0) {
          wasora_call(wasora_parser_expression(&debug->matrices_stride));
          
///kw+MILONGA_DEBUG+usage [ INCLUDE_INPUT ]
        } else if (strcasecmp(token, "INCLUDE_INPUT") == 0) {
          debug->include_input = 1;
        } else {
          wasora_push_error_message("unknown keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }
      }

      // si pidieron DEBUG, le pedimos a petsc que loguee cosas
#if PETSC_VERSION_LT(3,7,0)
      PetscLogBegin();
#else
      PetscLogDefaultBegin();
#endif
      PetscMemorySetGetMaximumUsage();
      wasora_define_instruction(milonga_instruction_debug, debug);

      return WASORA_PARSER_OK;

// ---------------------------------------------------------------------
///kw+MILONGA_STEP+usage MILONGA_STEP
///kw+MILONGA_STEP+desc Solves the linear eigenvalue problem.
    } else if (strcasecmp(token, "MILONGA_STEP") == 0) {

      milonga_step_t *milonga_step = calloc(1, sizeof(milonga_step_t));
      
      if (wasora_mesh.main_mesh == NULL) {
        wasora_push_error_message("no mesh found! (MILONGA_STEP before MESH)");
        return WASORA_PARSER_ERROR;
      }

      // chequeo de dimensiones
      if (milonga.dimensions == 0 && wasora_mesh.main_mesh->bulk_dimensions == 0) {
        wasora_push_error_message("no spatial dimensions given neither in MILONGA_PROBLEM nor in MESH");
        return WASORA_PARSER_ERROR;
      }
      
      // si alguna es cero, la rellenamos con la otra
      // TODO: no hay que rellenar la dimension de la malla con la de milonga!
      if (milonga.dimensions == 0) {
        milonga.dimensions = wasora_mesh.main_mesh->bulk_dimensions;
      } else if (wasora_mesh.main_mesh->bulk_dimensions == 0) {
        wasora_mesh.main_mesh->bulk_dimensions = milonga.dimensions;
      }
      
      // si son diferentes nos quejamos
      if (milonga.dimensions != wasora_mesh.main_mesh->bulk_dimensions) {
        wasora_push_error_message("inconsistent dimensions (MILONGA_PROBLEM = %d, MESH = %d)", milonga.dimensions, wasora_mesh.main_mesh->bulk_dimensions);
        return WASORA_PARSER_ERROR;
      }
      
      if (milonga.functions.phi == NULL) {
        wasora_call(milonga_define_result_functions());
      }

      while ((token = wasora_get_next_token(NULL)) != NULL) {
///kw+MILONGA_STEP+usage [ JUST_BUILD |
        if (strcasecmp(token, "JUST_BUILD") == 0) {
          milonga_step->do_not_build = 0;
          milonga_step->do_not_solve = 1;
///kw+MILONGA_STEP+usage JUST_SOLVE ]
        } else if (strcasecmp(token, "JUST_SOLVE") == 0) {
          milonga_step->do_not_build = 1;
          milonga_step->do_not_solve = 0;
///kw+MILONGA_STEP+usage [ VERBOSE ]
        } else if (strcasecmp(token, "VERBOSE") == 0) {
          milonga_step->verbose = 1;
        } else {
          wasora_push_error_message("unknown keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }

      }
      
      if (milonga.formulation == formulation_moc) {
        wasora_define_instruction(milonga_instruction_moc_step, milonga_step);
      } else {
        wasora_define_instruction(milonga_instruction_step, milonga_step);
      }
      
      return WASORA_PARSER_OK;

// ---------------------------------------------------------------------
///kw+FLUX_POST+usage FLUX_POST
///kw+FLUX_POST+desc Writes a post-processing file with total and partial fluxes and optionally XS distributions.
    } else if (strcasecmp(token, "FLUX_POST") == 0) {
      
      int flux = 1;
      int xs = 0;
      mesh_post_t *mesh_post = calloc(1, sizeof(mesh_post_t));
      
      // con esto le pegamos la mayor parte de las veces
//      mesh_post->cell_centered = wasora_mesh.default_cell_centered;
      
      while ((token = wasora_get_next_token(NULL)) != NULL) {
      
///kw+FLUX_POST+usage { FILE <name> |
        if (strcasecmp(token, "FILE") == 0) {
          wasora_call(wasora_parser_file(&mesh_post->file));
///kw+FLUX_POST+usage FILE_PATH <file_path> }
        } else if (strcasecmp(token, "FILE_PATH") == 0) {
          char *file_path;
          wasora_call(wasora_parser_string(&file_path));
          if ((mesh_post->file = wasora_define_file(file_path, file_path, 0, NULL, "w", 0)) == NULL) {
            return WASORA_RUNTIME_ERROR;
          }
          free(file_path);

///kw+FLUX_POST+usage [ XS ]
        } else if (strcasecmp(token, "XS") == 0) {
          xs = 1;
          
///kw+FLUX_POST+usage [ NO_MESH ]
        } else if (strcasecmp(token, "NOMESH") == 0 || strcasecmp(token, "NO_MESH") == 0) {
          mesh_post->no_mesh = 1;
      
///kw+FLUX_POST+usage [ FORMAT { gmsh | vtk } ]
        } else if (strcasecmp(token, "FORMAT") == 0) {
          char *keywords[] = {"gmsh", "vtk", ""};
          int values[] = {post_format_gmsh, post_format_vtk, 0};
          wasora_call(wasora_parser_keywords_ints(keywords, values, (int *)&mesh_post->format));
        }
      }
    
    
      mesh_post->mesh = wasora_mesh.main_mesh;
      mesh_post->flags = POST_INCLUDE_FLUX*flux + POST_INCLUDE_XS*xs;
      if (mesh_post->file == NULL) {
        wasora_push_error_message("neither FILE not FILE_PATH given (use explicitly 'stdout' if you intend to)");
        return WASORA_PARSER_ERROR;
      }
      
      if (mesh_post->format == post_format_fromextension) {
        char *ext = mesh_post->file->format + strlen(mesh_post->file->format) - 4;
        
        if (strcasecmp(ext, ".pos") == 0 || strcasecmp(ext, ".msh") == 0) {
          mesh_post->format = post_format_gmsh;
        } else if (strcasecmp(ext, ".vtk") == 0) {
          mesh_post->format = post_format_vtk;
        } else {
          wasora_push_error_message("unknown extension '%s' and no FORMAT given", ext);
          return WASORA_PARSER_ERROR;
        }
      }
        
      switch (mesh_post->format) {
        case post_format_gmsh:
          mesh_post->write_header = mesh_gmsh_write_header;
          mesh_post->write_mesh = mesh_gmsh_write_mesh;
          mesh_post->write_scalar = mesh_gmsh_write_scalar;
          mesh_post->write_vector = mesh_gmsh_write_vector;
        break;
        case post_format_vtk:
          mesh_post->write_header = mesh_vtk_write_header;
          mesh_post->write_mesh = mesh_vtk_write_mesh;
          mesh_post->write_scalar = mesh_vtk_write_scalar;
          mesh_post->write_vector = mesh_vtk_write_vector;
        break;
        default:
          return WASORA_PARSER_ERROR;
        break;
      }
      
      LL_APPEND(wasora_mesh.posts, mesh_post);
      wasora_define_instruction(wasora_instruction_mesh_post, mesh_post);
      return WASORA_PARSER_OK;
    
// ---------------------------------------------------------------------
///kw+TRACK_POST+usage TRACK_POST
///kw+TRACK_POST+desc Writes a post-processing file with ray tracing information.
    } else if (strcasecmp(token, "TRACK_POST") == 0) {
      
      double xi;
      track_post_t *track_post = calloc(1, sizeof(track_post_t));
      
      while ((token = wasora_get_next_token(NULL)) != NULL) {
        
///kw+TRACK_POST+usage { FILE <name> |
        if (strcasecmp(token, "FILE") == 0) {
          wasora_call(wasora_parser_file(&track_post->file));
///kw+TRACK_POST+usage FILE_PATH <file_path> }
        } else if (strcasecmp(token, "FILE_PATH") == 0) {
          char *file_path;
          wasora_call(wasora_parser_string(&file_path));
          if ((track_post->file = wasora_define_file(file_path, file_path, 0, NULL, "w", 0)) == NULL) {
            return WASORA_RUNTIME_ERROR;
          }
          free(file_path);

///kw+TRACK_POST+usage [ NO_MESH ]
        } else if (strcasecmp(token, "NOMESH") == 0 || strcasecmp(token, "NO_MESH") == 0) {
          track_post->no_mesh = 1;
      
///kw+TRACK_POST+usage [ FORMAT { gnuplot } ]
        } else if (strcasecmp(token, "FORMAT") == 0) {
          char *keywords[] = {"gnuplot", ""};
          int values[] = {post_format_gnuplot, 0};
          wasora_call(wasora_parser_keywords_ints(keywords, values, (int *)&track_post->format));
          
///kw+TRACK_POST+usage { TRACK <identifier> }
        } else if (strcasecmp(token, "TRACK") == 0) {
          char *track_name;
          
          wasora_call(wasora_parser_string(&track_name));
          if ((track_post->tracks = milonga_get_ray_tracing_ptr(track_name)) == NULL) {
            wasora_push_error_message("unknown track '%s'", track_name);
            free(track_name);
            return WASORA_PARSER_ERROR;
          }
          free(track_name);
          
///kw+TRACK_POST+usage [ TRACK_ID ]
        } else if (strcasecmp(token, "TRACK_ID") == 0) {
          wasora_call(wasora_parser_expression_in_string(&xi));
          track_post->track_id = (int)(xi);
          if (track_post->track_id < 1)  {
            wasora_push_error_message("a positive number of track id should be given instead of '%d'", track_post->track_id);
            return WASORA_PARSER_ERROR;
          }
          
///kw+TRACK_POST+usage [ N_AZIM_ID ]
        } else if (strcasecmp(token, "AZIM_ID") == 0) {
          wasora_call(wasora_parser_expression_in_string(&xi));
          track_post->azim_id = (int)(xi);
          if (track_post->azim_id < 1)  {
            wasora_push_error_message("a positive number of azim id should be given instead of '%d'", track_post->azim_id);
            return WASORA_PARSER_ERROR;
          }
          
        } else {
          wasora_push_error_message("unknown keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }
      }
      
      // si no nos dieron un track, nos vamos
      if (track_post->tracks == NULL) {
        wasora_push_error_message("missing ray tracing structure to plot (TRACK keyword)");
        return WASORA_PARSER_ERROR;
      }
      
      if (track_post->file == NULL) {
        wasora_push_error_message("neither FILE not FILE_PATH given (use explicitly 'stdout' if you intend to)");
        return WASORA_PARSER_ERROR;
      }
      
      if (track_post->azim_id > 0 && track_post->track_id > 0) {
        wasora_push_error_message("both AZIM_ID and TRACK_ID keywords specified");
        return WASORA_PARSER_ERROR;
      }
      
      if (track_post->format == post_format_from_extension) {
        char *ext = track_post->file->format + strlen(track_post->file->format) - 4;
        
        if (strcasecmp(ext, ".plt") == 0 || strcasecmp(ext, ".gnu") == 0 || strcasecmp(ext, ".gpi") == 0) {
          track_post->format = post_format_gnuplot;
        } else {
          wasora_push_error_message("unknown extension '%s' and no FORMAT given", ext);
          return WASORA_PARSER_ERROR;
        }
      }
      
      // dejo todo preparado por si alguien quiere agregar un formato
      switch (track_post->format) {
        case post_format_gnuplot:
          track_post->write_header = track_gnuplot_write_header;
          track_post->write_mesh = track_gnuplot_write_mesh;
          track_post->write_tracks = track_gnuplot_write_tracks;
          break;
        default:
          wasora_push_error_message("unknown TRACK_POST format");
          return WASORA_PARSER_ERROR;
      }
      
      wasora_define_instruction(milonga_instruction_track_post, track_post);
      return WASORA_PARSER_OK;
    }
  }
  
  return WASORA_PARSER_UNHANDLED;
}


#undef  __FUNCT__
#define __FUNCT__ "milonga_define_result_functions"
int milonga_define_result_functions(void) {
  
  char name[32];
  int n, g;
  
  // las definimos solo si ya sabemos cuantas dimensiones tiene el problema
  if (milonga.dimensions == 0) {
    return WASORA_PARSER_OK;
  }

  if (milonga.formulation == formulation_sn) {
    if (milonga.SN == 0) {
      wasora_push_error_message("internal inconsistency (SN = 0)");
      return WASORA_RUNTIME_ERROR;
    }

    // ecuacion 6.19 de stammler
    switch(milonga.dimensions) {
      case 1:
        milonga.directions = milonga.SN;
        break;
      case 2:
        milonga.directions = 0.5*milonga.SN*(milonga.SN+2);
        break;
      case 3:
        milonga.directions = milonga.SN*(milonga.SN+2);
        break;
    }

    // las psi son solo para transporte
    milonga.functions.psi = malloc(milonga.directions * sizeof(function_t *));
    for (n = 0; n < milonga.directions; n++) {
      milonga.functions.psi[n] = malloc(milonga.groups * sizeof(function_t *));
      for (g = 0; g < milonga.groups; g++) {
        sprintf(name, "psi%d.%d", n+1, g+1);
        if ((milonga.functions.psi[n][g] = wasora_define_function(name, milonga.dimensions)) == NULL) {
          return WASORA_PARSER_ERROR;
        }
      }
    }
  } else if (milonga.formulation == formulation_diffusion) {
    // ponemos esto para difusion para que sea facil calcular los grados de libertad
    milonga.directions = 1;
  } else if (milonga.formulation == formulation_moc) {
    // TODO: luego hacer algo similar a sn con moc
    milonga.directions = 1;
  }

  // las phi para todos  
  milonga.functions.phi = malloc(milonga.groups * sizeof(function_t *));
  for (g = 0; g < milonga.groups; g++) {
    sprintf(name, "phi%d", g+1);
    if ((milonga.functions.phi[g] = wasora_define_function(name, milonga.dimensions)) == NULL) {
      return WASORA_PARSER_ERROR;
    }
  }

  if ((milonga.functions.pow = wasora_define_function("pow", milonga.dimensions)) == NULL) {
    return WASORA_PARSER_ERROR;
  }

  return WASORA_PARSER_OK;
}
