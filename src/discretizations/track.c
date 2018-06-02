/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's unstructured and structured ray tracing routines
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

int milonga_instruction_track(void *arg) {
  
  tracks_t *tracks = (tracks_t *)arg;
  
  if (tracks->initialized) {
    return WASORA_RUNTIME_OK;
  }
  
  if (tracks->mesh->bulk_dimensions != 2) {
    wasora_push_error_message("tracking is only supported for 2D geometries :( wanna code it?");
    return WASORA_RUNTIME_ERROR;
  }
  
  // leemos los parametros de entrada
  wasora_call(track_read_tracking_parameters(tracks));
  
  // inicializamos las cuadraturas
  wasora_call(track_init_quadratures_before_tracking(tracks));
  
  // computamos los tracks y rellenamos alguna informacion de cuadraturas
  wasora_call(track_compute_tracks(tracks));
  
  // completamos el calculo de las cuadraturas
  wasora_call(track_init_quadratures_after_tracking(tracks));
  
  // llevamos los tracks al origen real de la malla
  wasora_call(track_recalibrate_coords(tracks));
  
  // segmentamos los tracks segun las regiones que atraviesan
  wasora_call(track_segmentize_tracks(tracks));
  
  // si no se pide lo contrario, alteramos las longitudes de los tracks tal q $A_{i,m} = A_i$
  if (tracks->do_not_correct_volumes == 0) {
    wasora_call(track_correct_volumes(tracks));
  }
  
  tracks->initialized = 1;
  
  return WASORA_RUNTIME_OK;
}

int track_read_tracking_parameters(tracks_t *tracks) {
   
  // evaluamos los parametros de la cuadratura azimutal
  tracks->n_azim = (int) wasora_evaluate_expression(tracks->expr_n_azim);
  if (tracks->n_azim <= 0) {
    wasora_push_error_message("number of azimuthal angles '%d' cannot be negative or zero", tracks->n_azim);
    return WASORA_RUNTIME_ERROR;
  } else if (tracks->n_azim % 4 != 0) {
    wasora_push_error_message("number of azimuthal angles '%d' is not a multiple of 4", tracks->n_azim);
    return WASORA_RUNTIME_ERROR;
  } else {
    tracks->n_azim_2 = tracks->n_azim / 2;
    tracks->n_azim_4 = tracks->n_azim / 4;
  }
  
  // evaluamos los parametros del tracking
  if (tracks->expr_track_dens->n_tokens != 0) {
    tracks->track_dens = wasora_evaluate_expression(tracks->expr_track_dens);
    if (tracks->track_dens <= 0) {
      wasora_push_error_message("track density cannot be negative or zero for track '%s", tracks->name);
      return WASORA_RUNTIME_ERROR;
    }
    tracks->track_spacing = 1 / tracks->track_dens;
  } else {
    tracks->track_spacing = wasora_evaluate_expression(tracks->expr_track_spacing);
    if (tracks->track_spacing <= 0) {
      wasora_push_error_message("track spacing cannot be negative or zero for track '%s", tracks->name);
      return WASORA_RUNTIME_ERROR;
    }
    tracks->track_dens = 1 / tracks->track_spacing;
  }
  
  if (tracks->expr_tiny_step->n_tokens != 0) {
    tracks->tiny_step = wasora_evaluate_expression(tracks->expr_tiny_step);
  } else {
    // por defecto
    tracks->tiny_step = 1e-8;
  }
  
  return WASORA_RUNTIME_OK;
}

int track_init_quadratures_before_tracking(tracks_t *tracks) {
  
  // quadrature
  tracks->quadrature = calloc(1, sizeof(quadrature_t));
  
  // azimuthal quadrature
  tracks->quadrature->azimuthal = calloc(1, sizeof(azimuthal_quadrature_t));
  tracks->quadrature->azimuthal->n_azim = tracks->n_azim;
  tracks->quadrature->azimuthal->n_azim_2 = tracks->n_azim_2;
  tracks->quadrature->azimuthal->n_azim_4 = tracks->n_azim_4;
  tracks->quadrature->azimuthal->phi = calloc(tracks->quadrature->azimuthal->n_azim_2, sizeof(double));
  tracks->quadrature->azimuthal->w = calloc(tracks->quadrature->azimuthal->n_azim_2, sizeof(double));
  
  // effective spacing
  tracks->quadrature->effective_spacing = calloc(tracks->n_azim_2, sizeof(double));
  
  // polar quadrature
  //tracks->quadrature->polar = calloc(1, sizeof(polar_quadrature_t));
  
  return WASORA_RUNTIME_OK;
}

int track_compute_tracks(tracks_t * tracks) {
  
  int i, j;
  int id;
  
  double phi;
  
  double width_x;
  double width_y;
  
  double m;
  double exit_point[3];
  
  double *dx_eff;
  double *dy_eff;
  double *d_eff;
  
  track_t *track;
  
  tracks->n_tracks_x = calloc(tracks->n_azim_2, sizeof(int));
  tracks->n_tracks_y = calloc(tracks->n_azim_2, sizeof(int));
  tracks->n_tracks = calloc(tracks->n_azim_2, sizeof(int));
  
  // deltas entre los tracks en cada eje y perpendicular a los mismos
  dx_eff = calloc(tracks->n_azim_2, sizeof(double));
  dy_eff = calloc(tracks->n_azim_2, sizeof(double));
  d_eff = calloc(tracks->n_azim_2, sizeof(double));
  
  width_x = tracks->mesh->bounding_box_max.x[0] - tracks->mesh->bounding_box_min.x[0];
  width_y = tracks->mesh->bounding_box_max.x[1] - tracks->mesh->bounding_box_min.x[1];
  
  // barro unicamente el primer cuadrante
  for (i = 0; i < tracks->n_azim_4; i++) {
    
    phi = M_PI / tracks->n_azim_2 * (0.5 + i);
    
    tracks->n_tracks_x[i] = (int) (fabs(width_x / tracks->track_spacing * sin(phi))) + 1;
    tracks->n_tracks_y[i] = (int) (fabs(width_y / tracks->track_spacing * cos(phi))) + 1;
    tracks->n_tracks[i] = tracks->n_tracks_x[i] + tracks->n_tracks_y[i];
    
    // angulo efectivo, cerca del que deseamos
    phi = atan((width_y * tracks->n_tracks_x[i]) / (width_x * tracks->n_tracks_y[i]));
    tracks->quadrature->azimuthal->phi[i] = phi;
    
    // deltas efectivos
    dx_eff[i] = width_x / tracks->n_tracks_x[i];
    dy_eff[i] = width_y / tracks->n_tracks_y[i];
    d_eff[i] = dx_eff[i] * sin(phi);
    
    tracks->quadrature->effective_spacing[i] = d_eff[i];
    
    // y los angulos suplementarios
    tracks->n_tracks_x[tracks->n_azim_2-i-1] = tracks->n_tracks_x[i];
    tracks->n_tracks_y[tracks->n_azim_2-i-1] = tracks->n_tracks_y[i];
    tracks->n_tracks[tracks->n_azim_2-i-1] = tracks->n_tracks[i];
    
    tracks->quadrature->azimuthal->phi[tracks->n_azim_2-i-1] = M_PI - phi;
    
    dx_eff[tracks->n_azim_2-i-1] = dx_eff[i];
    dy_eff[tracks->n_azim_2-i-1] = dy_eff[i];
    d_eff[tracks->n_azim_2-i-1] = d_eff[i];
    
    tracks->quadrature->effective_spacing[tracks->n_azim_2-i-1] = d_eff[i];
  }
  
  // empiezo a allocar la matriz de tracks
  tracks->track = calloc(tracks->n_azim_2, sizeof(track_t *));
  
  // puntos de inicio y final de cada track (referenciados al (0,0) ficticio en el borde inferior izquierdo)
  // barro los dos primeros cuadrantes
  for (i = 0; i < tracks->n_azim_2; i++) {
    
    // recupero el angulo
    phi = tracks->quadrature->azimuthal->phi[i];

    // alloco segun la cantidad de tracks dado el angulo i
    tracks->track[i] = calloc(tracks->n_tracks[i], sizeof(track_t));

    // puntos de inicio para cada track sobre el eje x
    for (j = 0; j < tracks->n_tracks_x[i]; j++) {
      
      track = &tracks->track[i][j];
      
      // si el track apunta a la derecha
      if (i < tracks->n_azim_4) {
        wasora_set_point_coords(track->p_i, dx_eff[i] * (tracks->n_tracks_x[i] - j - 0.5), 0, tracks->z_value);
        
        // o si el track apunta a la izquierda
      } else {
        wasora_set_point_coords(track->p_i, dx_eff[i] * (0.5 + j), 0, tracks->z_value);
      }
    }

    // puntos de inicio sobre el eje y
    for (j = 0; j < tracks->n_tracks_y[i]; j++) {

      track = &tracks->track[i][tracks->n_tracks_x[i]+j];
      
      // tracks que inician en x = 0 (apuntan a la derecha)
      if (i < tracks->n_azim_4) {
        wasora_set_point_coords(track->p_i, 0, dy_eff[i] * (0.5 + j), tracks->z_value);
        
        // tracks que inician en x = width_x (apuntan a la izquierda)
      } else {
        wasora_set_point_coords(track->p_i, width_x, dy_eff[i] * (0.5 + j), tracks->z_value);
      }
    }
    
    // recorremos todos los tracks para una dada direccion azimutal y buscamos el punto de salida
    for (j = 0; j < tracks->n_tracks[i]; j++) {
      
      track = &tracks->track[i][j];
      
      // seteamos indice azimutal y angulo de cada track
      track->azim_index = i;
      track->phi = phi;
      
      m = tan(track->phi);
      
      // cualquier track, indistintamente a donde apunte, puede salir por y = width_y
      wasora_set_point_coords(exit_point, track->p_i[0] - (track->p_i[1] - width_y) / m, width_y, tracks->z_value);
      
      // nos fijamos si encontramos el punto de salida
      if (exit_point[0] > 0.0 && exit_point[0] < width_x) {
        
        wasora_set_point_coords(track->p_o, exit_point[0], exit_point[1], exit_point[2]);
        
        // sino, seguimos buscando el punto de salida
      } else {
        
        // aquellos tracks que apuntan a la derecha
        if (track->azim_index < tracks->n_azim_4) {
          
          // pueden salir tambien por x = width_x
          wasora_set_point_coords(exit_point, width_x, track->p_i[1] + m * (width_x - track->p_i[0]), tracks->z_value);
          
          // aquellos tracks que apuntan a la izquierda
        } else {
          
          // pueden salir tambien por x = 0
          wasora_set_point_coords(exit_point, 0.0, track->p_i[1] - m * track->p_i[0], tracks->z_value);
        }
        
        // nos fijamos si encontramos el punto de salida
        if (exit_point[1] > 0.0 && exit_point[1] < width_y) {
          wasora_set_point_coords(track->p_o, exit_point[0], exit_point[1], exit_point[2]);
        } else {
          wasora_push_error_message("track %d exit point was not found", track->id);
          return WASORA_RUNTIME_ERROR;
        }
      }
      
      // una vez que encontramos el punto de salida, calculamos la longitud del track
      track->length = mesh_subtract_module(track->p_o, track->p_i);
    }
  }
  
  free(dx_eff);
  free(dy_eff);
  free(d_eff);
  
  // una vez que tengo los tracks, los cuento y ...
  for (tracks->n_total_tracks = 0, i = 0; i < tracks->n_azim_2; i++) 
    tracks->n_total_tracks += tracks->n_tracks[i];
  
  // ... armo un vector de tracks en funcion de un unico indice para hacer easier loops
  tracks->track_by_id = calloc(tracks->n_total_tracks, sizeof(tracks_t *));
  for (id = 0, i = 0; i < tracks->n_azim_2; i++) {
    for (j = 0; j < tracks->n_tracks[i]; j++) {
      tracks->track_by_id[id] = &tracks->track[i][j];
      tracks->track_by_id[id]->id = id + 1;
      id++;
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int track_init_quadratures_after_tracking(tracks_t *tracks) {
  
  // calculamos la azimuthal quadrature
  wasora_call(track_init_azimuthal_weights(tracks->quadrature->azimuthal));
  
  
  // CREO QUE ESTA MAL CONCEPTUALMENTE LLAMAR ACA ESTO, YA QUE EL RAY TRACING NO
  // DEPENDE REALMENTE DE LO POLAR (ES 2D), ENTONCES BASTA CON LLAMARLA EN LA 
  // INSTRUCCION DE POLAR QUADRATURE O EN MILONGA PROBLEM
  // ESTO INCLUSO ME PERMITIRIA PASAR EL RAY TRACING A WASORA, DONDE DEBE ESTAR!
  // calculamos una polar quadrature por defecto
  // wasora_call(track_init_polar_quadrature(tracks->quadrature->polar));
  
  return WASORA_RUNTIME_OK;
}

int track_init_azimuthal_weights(azimuthal_quadrature_t *azimuthal) {
  
  int i;
  
  for (i = 0; i < azimuthal->n_azim_4; i++) {
    if (i == 0) {
      azimuthal->w[i] = 0.5 * (azimuthal->phi[i+1] + azimuthal->phi[i]);
    } else if (i == (azimuthal->n_azim_4-1)) {
      azimuthal->w[i] = M_PI_2 - 0.5 * (azimuthal->phi[i] + azimuthal->phi[i-1]);
    } else {
      azimuthal->w[i] = 0.5 * (azimuthal->phi[i+1] - azimuthal->phi[i-1]);
    }
    // normalizo
    azimuthal->w[i] /= (2.0 * M_PI);
    // e igualo los pesos que correspondan
    azimuthal->w[azimuthal->n_azim_2-i-1] = azimuthal->w[i];
  }
  
  return WASORA_RUNTIME_OK;  
}

int track_init_polar_quadrature(polar_quadrature_t *polar) {
  
  // si no se dio una expresion
  if (polar->expr_n_polar == NULL) {
  //if (polar == NULL) {
    
    //polar = calloc(1, sizeof(polar_quadrature_t));
    
    // por defecto, tabuchi-yamamoto con 6 polar angles
    polar->quad_type = tabuchi_yamamoto;
    polar->n_polar = 6;
    polar->n_polar_2 = 6 / 2;
    
  } else {
    
    // no se si hace falta esto cuado este cambiando todo
    // if (polar->expr_n_polar != NULL)
    
    polar->n_polar = (int) wasora_evaluate_expression(polar->expr_n_polar);
    
    if (polar->n_polar <= 0) {
      wasora_push_error_message("number of polar angles '%d' cannot be negative or zero", polar->n_polar);
      return WASORA_RUNTIME_ERROR;
    } else if (polar->n_polar % 2 != 0) {
      wasora_push_error_message("number of polar angles '%d' is not a multiple of 2", polar->n_polar);
      return WASORA_RUNTIME_ERROR;
    } else {
      polar->n_polar_2 = polar->n_polar / 2;
    }
  }
  
  wasora_call(track_init_polar_quadrature_values(polar));
  
  return WASORA_RUNTIME_OK;
}

int track_init_polar_quadrature_values(polar_quadrature_t *polar_quad) {
  
  int i;
  
  polar_quad->theta = calloc(polar_quad->n_polar, sizeof(double *));
  polar_quad->sin_theta = calloc(polar_quad->n_polar, sizeof(double *));
  polar_quad->w = calloc(polar_quad->n_polar, sizeof(double *));
  
  switch (polar_quad->quad_type) {
    case tabuchi_yamamoto:
      // del Handbook de Ingenieria Nuclear pg. 1154
      if (polar_quad->n_polar == 2) {
        
        polar_quad->sin_theta[0] = 0.798184;
        polar_quad->w[0] =   1.0 / 2.0;
        
      } else if (polar_quad->n_polar == 4) {
        
        polar_quad->sin_theta[0] = 0.363900;
        polar_quad->w[0] =   0.212854 / 2.0;
        
        polar_quad->sin_theta[1] = 0.899900;
        polar_quad->w[1] =   0.787146 / 2.0;
        
      } else if (polar_quad->n_polar == 6) {
        
        polar_quad->sin_theta[0] = 0.166648;
        polar_quad->w[0] =   0.046233 / 2.0;
        
        polar_quad->sin_theta[1] = 0.537707;
        polar_quad->w[1] =   0.283619 / 2.0;
        
        polar_quad->sin_theta[2] = 0.932954;
        polar_quad->w[2] =   0.670148 / 2.0;
        
      } else {
        wasora_push_error_message("'tabuchi_yamamoto' quadrature not suported for N_POLAR_ANGLES superior than 6");
        return WASORA_RUNTIME_ERROR;
      }
      break;
    case equal_weight:
      wasora_push_error_message("'equal_weight' polar quadrature not yet suported");
      return WASORA_RUNTIME_ERROR;
    case equal_angle:
      wasora_push_error_message("'equal_angle' polar quadrature not yet suported");
      return WASORA_RUNTIME_ERROR;
    case gauss_legendre:
      // del Handbook de Ingenieria Nuclear pg. 1153
      if (polar_quad->n_polar == 2) {
        
        polar_quad->theta[0]   = acos(0.5773502691);
        polar_quad->w[0] = 1.0 / 2.0;
        
      } else if (polar_quad->n_polar == 4) {
        
        polar_quad->theta[0]   = acos(0.3399810435);
        polar_quad->w[0] = 0.6521451549 / 2.0;
        
        polar_quad->theta[1]   = acos(0.8611363115);
        polar_quad->w[1] = 0.3478548451 / 2.0;
        
      } else if (polar_quad->n_polar == 6) {
        
        polar_quad->theta[0]   = acos(0.2386191860);
        polar_quad->w[0] = 0.4679139346 / 2.0;
        
        polar_quad->theta[1]   = acos(0.6612093864);
        polar_quad->w[1] = 0.3607615730 / 2.0;
        
        polar_quad->theta[2]   = acos(0.9324695142);
        polar_quad->w[2] = 0.1713244924 / 2.0;
        
      } else if (polar_quad->n_polar == 8) {
        
        polar_quad->theta[0]   = acos(0.1834346424);
        polar_quad->w[0] = 0.3626837834 / 2.0;
        
        polar_quad->theta[1]   = acos(0.5255324099);
        polar_quad->w[1] = 0.3137066459 / 2.0;
        
        polar_quad->theta[2]   = acos(0.7966664774);
        polar_quad->w[2] = 0.2223810344 / 2.0;
        
        polar_quad->theta[3]   = acos(0.9602898564);
        polar_quad->w[3] = 0.1012285363 / 2.0;
        
      } else if (polar_quad->n_polar == 10) {
        
        polar_quad->theta[0]   = acos(0.1488743387);
        polar_quad->w[0] = 0.2955242247 / 2.0;
        
        polar_quad->theta[1]   = acos(0.4333953941);
        polar_quad->w[1] = 0.2692667193 / 2.0;
        
        polar_quad->theta[2]   = acos(0.6794095682);
        polar_quad->w[2] = 0.2190863625 / 2.0;
        
        polar_quad->theta[3]   = acos(0.8650633666);
        polar_quad->w[3] = 0.1494513492 / 2.0;
        
        polar_quad->theta[4]   = acos(0.9739065285);
        polar_quad->w[4] = 0.0666713443 / 2.0;
        
      } else if (polar_quad->n_polar == 12) {
        
        polar_quad->theta[0]   = acos(0.1252334085);
        polar_quad->w[0] = 0.2491470458 / 2.0;
        
        polar_quad->theta[1]   = acos(0.3678314989);
        polar_quad->w[1] = 0.2334925365 / 2.0;
        
        polar_quad->theta[2]   = acos(0.5873179542);
        polar_quad->w[2] = 0.2031674267 / 2.0;
        
        polar_quad->theta[3]   = acos(0.7699026741);
        polar_quad->w[3] = 0.1600783286 / 2.0;
        
        polar_quad->theta[4]   = acos(0.9041172563);
        polar_quad->w[4] = 0.1069393260 / 2.0;
        
        polar_quad->theta[5]   = acos(0.9815606342);
        polar_quad->w[5] = 0.0471753364 / 2.0;
        
      } else {
        wasora_push_error_message("'gauss_legendre' quadrature not suported for N_POLAR_ANGLES superior than 12");
        return WASORA_RUNTIME_ERROR;
      }
      break;
    case leonard:
      // ver de que referencia sacar esto (paper tabuchi es diferente a openmoc y handbook tiene otras optimizadas por Herbert)
      // estas se corresponden con el paper de tabuchi
      if (polar_quad->n_polar == 2) {
        
        polar_quad->sin_theta[0] = 0.752244;
        polar_quad->w[0] =   1.0 / 2.0;
        
      } else if (polar_quad->n_polar == 4) {
        
        polar_quad->sin_theta[0] = 0.273658;
        polar_quad->w[0] =   0.139473 / 2.0;
        
        polar_quad->sin_theta[1] = 0.865714;
        polar_quad->w[1] =   0.860527 / 2.0;
        
      } else if (polar_quad->n_polar == 6) {
        
        polar_quad->sin_theta[0] = 0.103840;
        polar_quad->w[0] =   0.020530 / 2.0;
        
        polar_quad->sin_theta[1] = 0.430723;
        polar_quad->w[1] =   0.219161 / 2.0;
        
        polar_quad->sin_theta[2] = 0.905435;
        polar_quad->w[2] =   0.760309 / 2.0;
        
      } else {
        wasora_push_error_message("'leonard' quadrature not suported for N_POLAR_ANGLES superior than 6");
        return WASORA_RUNTIME_ERROR;
      }
      break;
  }
  
  // si nos dieron los senos
  if (polar_quad->quad_type == tabuchi_yamamoto || polar_quad->quad_type == leonard) {
    // seteamos las componentes que faltan: rellenamos los sin_theta, los pesos y calculamos los angulos 
    for (i = 0; i < polar_quad->n_polar_2; i++) {
      polar_quad->sin_theta[polar_quad->n_polar-i-1] = polar_quad->sin_theta[i];
      polar_quad->w[polar_quad->n_polar-i-1] = polar_quad->w[i];
      
      polar_quad->theta[i] = asin(polar_quad->sin_theta[i]);
      polar_quad->theta[polar_quad->n_polar-i-1] = M_PI - polar_quad->theta[i];
    }
  }
  
  // si nos dieron los angulos
  if (polar_quad->quad_type == gauss_legendre) {
    // seteamos las componentes que faltan: rellenamos los theta, los pesos y calculamos los sin_theta 
    for (i = 0; i < polar_quad->n_polar_2; i++) {
      polar_quad->theta[polar_quad->n_polar-i-1] = M_PI - polar_quad->theta[i];
      polar_quad->w[polar_quad->n_polar-i-1] = polar_quad->w[i];
      
      polar_quad->sin_theta[i] = sin(polar_quad->theta[i]);
      polar_quad->sin_theta[polar_quad->n_polar-i-1] = sin(polar_quad->theta[polar_quad->n_polar-i-1]);
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int track_recalibrate_coords(tracks_t *tracks) {
  
  int i;
  
  track_t *track;
  
  // llevamos al origen posta todos los tracks
  for (i = 0; i < tracks->n_total_tracks; i++) {
    
    track = tracks->track_by_id[i];
    
    track->p_i[0] += tracks->mesh->bounding_box_min.x[0];
    track->p_i[1] += tracks->mesh->bounding_box_min.x[1];
    track->p_i[2] += 0.0;
    
    track->p_o[0] += tracks->mesh->bounding_box_min.x[0];
    track->p_o[1] += tracks->mesh->bounding_box_min.x[1];
    track->p_o[2] += 0.0;
  }
  
  return WASORA_RUNTIME_OK;
}

int track_segmentize_tracks(tracks_t *tracks) {
  
  int i;
  
  double x0[3];
  
  double **segment_points;
  double p_i[3];
  double p_o[3];
  double eps;
  double length;
  
  element_t *element;
  element_t *prev_element;
  segment_t *segment;
  track_t *track;
  segment_list_item_t *associated_segment;
  
  segment_points = calloc(2, sizeof(double *));
  for (i = 0; i < 3; i++) {
    segment_points[i] = calloc(3, sizeof(double));
  }
  
  // el ray tracing necesita que point_in_element dentro de mesh_find_element 
  // tenga una tolerancia nula para no retornar elementos que no contienen al
  // punto
  eps = wasora_var(wasora_mesh.vars.eps);
  wasora_var(wasora_mesh.vars.eps) = 0.0;
  
  // recorremos cada track a segmentar
  for (i = 0; i < tracks->n_total_tracks; i++) {
    
    track = tracks->track_by_id[i];
    
    // la ecuacion general del track
    track_compute_general_equation(track->ABC, track->p_i, track->p_o);
    
    // perturbamos el punto inicial del track para caer dentro de un elemento
    wasora_set_point_coords(x0, track->p_i[0], track->p_i[1], track->p_i[2]);
    wasora_adjust_point_coords(x0, tracks->tiny_step, M_PI_2, track->phi);
    
    prev_element = NULL;
    while ((element = mesh_find_element(tracks->mesh, x0)) != NULL) {
      
      // si estamos afuera del dominio y aun asi wasora nos mando un elemento,
      // nos vamos: tener presente que el segmento definido en el paso anterior
      // posee como punto final la superficie del dominio.
      if (x0[0] >= tracks->mesh->bounding_box_max.x[0] || x0[0] <= tracks->mesh->bounding_box_min.x[0] || 
          x0[1] >= tracks->mesh->bounding_box_max.x[1] || x0[1] <= tracks->mesh->bounding_box_min.x[1])
        break;
      
      // si, por ejemplo, el tiny step fue cortito y wasora retorno nuevamente
      // el anterior elemento, nos movemos un paso mas
      if (prev_element != NULL && element == prev_element) {
        wasora_adjust_point_coords(x0, tracks->tiny_step, M_PI_2, track->phi);
        continue;
      }
      
      // puede ocurrir que wasora retorne un elemento que no contenga al punto.
      // esto puede ocurrir por dos motivos:
      //   1. point_in_element retorna true cuando en realidad deberia ser false;
      //   2. la malla esta "deformada" y el radio de busqueda de nodos no es lo
      //      suficientemente grande
      // 1. se soluciona mediante el uso de eps = 0.0
      // 2. se soluciona incrementando el radio de busqueda.
      // como aun no hemos aplicado la solucion de 2. y ademas, por completitud, 
      // realizamos un checkeo aqui
      if (!(element->type->point_in_element(element, x0))) {
        wasora_push_error_message("point in element did not return the correct element for point (%.4e,%.4e), please increase search radius", x0[0], x0[1]);
        return WASORA_RUNTIME_ERROR;
      }
      
      // recursivamente encontramos los puntos de definicion del segmento
      track_compute_element_intersections(element, track, tracks->tiny_step, 1e-16, &segment_points);
      
      // se determina si el point es de entrada o salida
      if (track->phi < M_PI_2) {
        if (segment_points[0][0] < segment_points[1][0]) {
          wasora_set_point_coords(p_i, segment_points[0][0], segment_points[0][1], segment_points[0][2]);
          wasora_set_point_coords(p_o, segment_points[1][0], segment_points[1][1], segment_points[1][2]);
        } else {
          wasora_set_point_coords(p_i, segment_points[1][0], segment_points[1][1], segment_points[1][2]);
          wasora_set_point_coords(p_o, segment_points[0][0], segment_points[0][1], segment_points[0][2]);
        }
      } else {
        if (segment_points[0][0] > segment_points[1][0]) {
          wasora_set_point_coords(p_i, segment_points[0][0], segment_points[0][1], segment_points[0][2]);
          wasora_set_point_coords(p_o, segment_points[1][0], segment_points[1][1], segment_points[1][2]);
        } else {
          wasora_set_point_coords(p_i, segment_points[1][0], segment_points[1][1], segment_points[1][2]);
          wasora_set_point_coords(p_o, segment_points[0][0], segment_points[0][1], segment_points[0][2]);
        }
      }
      
      // estamos en condiciones de crear un segmento
      track_create_segment(&segment, p_i, p_o, element);
      track_append_segment_to_list(&track->associated_segments, segment);
      
      // y actualizamos los valores
      wasora_set_point_coords(x0, p_o[0], p_o[1], p_o[2]);
      wasora_adjust_point_coords(x0, tracks->tiny_step, M_PI_2, track->phi);
      
      prev_element = element;
    }
    
    // computamos la longitud del track a partir de la sumatoria de las longitudes de los segmentos
    length = 0.0;
    LL_FOREACH(track->associated_segments, associated_segment) {
      segment = associated_segment->segment;
      length += segment->length;
    }
    
    // y luego verificamos si se cumple la igualdad con la longitud del track (para ver si nos salteamos algun elemento)
    if (!(gsl_fcmp(length, track->length, tracks->tiny_step) == 0)) {
      
      if (tracks->debug) {
        
        track_post_t *track_post = calloc(1, sizeof(track_post_t));
        
        char id[10];
        sprintf(id, "%d", track->id);
        
        char *file_path = malloc(6 + strlen(tracks->name) + 7 + strlen(id) + 3 + 1);
        sprintf(file_path, "debug-%s-track-%d.gp", tracks->name, track->id);
        
        if ((track_post->file = wasora_define_file(file_path, file_path, 0, NULL, "w", 0)) == NULL) {
          return WASORA_RUNTIME_ERROR;
        }
        
        track_post->tracks = tracks;
        track_post->track_id = track->id;
        track_post->write_header = track_gnuplot_write_header;
        track_post->write_mesh = track_gnuplot_write_mesh;
        track_post->write_tracks = track_gnuplot_write_tracks;
        
        milonga_instruction_track_post(track_post);
        
        wasora_push_error_message("Check out the file '%s' with gnuplot to inspect the problem.", file_path);
        
        free(file_path);
        free(track_post);
      } else {
        wasora_push_error_message("Use DEBUG option when doing the raytracing to get a gnuplot file where the problem can be inspected.");
      }
      
      wasora_push_error_message("track %d has been poorly segmented in '%s': real length is %.7e while approximate length is %.7e (try using a different TINY_STEP).", track->id, tracks->name, track->length, length);
      
      return WASORA_RUNTIME_ERROR;
    }
  }
  
  // restauramos el valor
  wasora_var(wasora_mesh.vars.eps) = eps;
  
  for (i = 0; i < 3; i++) {
    free(segment_points[i]);
  }
  free(segment_points);
  
  return WASORA_RUNTIME_OK;
}

int track_compute_element_intersections(element_t *element, track_t *track, double tiny_step, double eps, double ***segment_points) {
  
  int i, j;
  int n_int;
  int parallel_side;
  int parallel_found;
  
  double x_int[3] = { 0 };
  double int_points[4][3];
  double ABC[3];
  
  double length;
  
  double *x1;
  double *x2;
  
  // primero chekeamos no habernos ido de tema
  if (eps > tiny_step) {
    wasora_push_error_message("recursive ray tracing algorithm has failed the segmentation for track '%d' since search tolerance %.4e is bigger than tiny step %.4e", track->id, eps, tiny_step);
    return WASORA_RUNTIME_ERROR;
  }
  
  // contador de intersecciones
  n_int = 0;
  
  // asumimos que el track no es paralelo a ninguna cara del elemento
  parallel_found = 0;
  
  // nos movemos sobre las caras del elemento
  for (i = 0; i < element->type->faces; i++) {
    
    x1 = element->node[i]->x;
    x2 = element->node[(i == (element->type->faces-1)) ? 0 : i+1]->x;
    
    // la ecuacion general a partir de los nodos
    track_compute_general_equation(ABC, x1, x2);
    
    // la interseccion entre las rectas
    track_compute_intersection(x_int, track->ABC, ABC, &parallel_side);
    
    // si no hubo interseccion (rectas paralelas), nos vamos
    if (parallel_side == 1) {
      
      parallel_found = 1;
      continue;
      
      // si la interseccion no se da en el segmento, no cuenta
    } else if (!(wasora_point_in_segment(x1, x2, x_int, eps))) {
      
      continue;
      
      // sino, guardo el punto de interseccion
    } else {
      
      wasora_set_point_coords(int_points[n_int], x_int[0], x_int[1], x_int[2]);
      
      // incrementamos el numero de intersecciones
      n_int++;
    }
  }
  
  // si hay 3 o 4 puntos de interseccion, elegimos los dos mas alejados
  if (n_int == 3 || n_int == 4) {
    
    for (length = 0.0, i = 1; i < n_int; i++) {
      for (j = i; j < n_int; j++) {
        
        x1 = int_points[i-1];
        x2 = int_points[j];
        
        if (fabs(mesh_subtract_module(x1, x2)) > length) {
          wasora_set_point_coords((*segment_points)[0], x1[0], x1[1], x1[2]);
          wasora_set_point_coords((*segment_points)[1], x2[0], x2[1], x2[2]);
          length = mesh_subtract_module(x1, x2);
        }
      }
    }
    
    return WASORA_RUNTIME_OK;
    
    // si hay dos puntos de interseccion y ademas un segmento era paralelo al track
  } else if (n_int == 2 && parallel_found == 1) {
    
    wasora_set_point_coords((*segment_points)[0], int_points[0][0], int_points[0][1], int_points[0][2]);
    wasora_set_point_coords((*segment_points)[1], int_points[1][0], int_points[1][1], int_points[1][2]);
    
    return WASORA_RUNTIME_OK;
    
  } else if (n_int == 2 && parallel_found == 0) {
    
    // si las intersecciones no estan separadas al menos un tiny step, volvemos a calcular
    if (mesh_subtract_module(int_points[0], int_points[1]) < tiny_step) {
      
      track_compute_element_intersections(element, track, tiny_step, 10*eps, segment_points);
      
    // si no, esta todo piola
    } else {
      
      wasora_set_point_coords((*segment_points)[0], int_points[0][0], int_points[0][1], int_points[0][2]);
      wasora_set_point_coords((*segment_points)[1], int_points[1][0], int_points[1][1], int_points[1][2]);
      
      return WASORA_RUNTIME_OK;
    }
    
    // pero tambien puede pasar esta desgracia y hay que volver a calcular
  } else if (n_int == 1 || n_int == 0) {
    
    track_compute_element_intersections(element, track, tiny_step, 10*eps, segment_points);
  }
  
  return WASORA_RUNTIME_OK;
}

int track_correct_volumes(tracks_t *tracks) {
  
  int i, j, k;
  
  double alpha_im;
  
  track_t *track;
  segment_list_item_t *associated_segment;
  segment_t *segment;
  
  tracks->volumes = calloc(tracks->mesh->n_cells, sizeof(double));
  
  for (i = 0; i < tracks->n_azim_2; i++) {
    for (j = 0; j < tracks->n_tracks[i]; j++) {
      
      track = &tracks->track[i][j];
      
      LL_FOREACH(track->associated_segments, associated_segment) {
        segment = associated_segment->segment;
        tracks->volumes[segment->element->cell->id - 1] += tracks->quadrature->effective_spacing[track->azim_index] * segment->length;
      }
    }
    
    // primero un check
    for (k = 0; k < tracks->mesh->n_cells; k++) {
      
      // si no tengo el volumen real, lo computo
      if (tracks->mesh->cell[k].volume == 0) {
        tracks->mesh->cell[k].volume = tracks->mesh->cell[k].element->type->element_volume(tracks->mesh->cell[k].element);
      }
      
      if (tracks->volumes[k] == 0) {
        wasora_push_error_message("cell %d has been untracked for number %d azimutal direction with the defined ray tracing parameters in '%s'", tracks->mesh->cell[k].id, i, tracks->name);
        return WASORA_RUNTIME_ERROR;
      }
    }
      
    // recorro todos los tracks corrigiendo las longitudes de los segmentos para esta direccion azimutal i (corrijo tal que Ai,m = Ai)
    for (j = 0; j < tracks->n_tracks[i]; j++) {
      
      track = &tracks->track[i][j];
      
      LL_FOREACH(track->associated_segments, associated_segment) {
        
        segment = associated_segment->segment;
        alpha_im = tracks->mesh->cell[segment->element->cell->id - 1].volume / tracks->volumes[segment->element->cell->id - 1];
        
        // por ahi este checkeo lo vuelo a la mierda, y solo corrijo sin checkear si corregir es una buena aproximacion
        // checkeo estar dentro del 15% de variacion respecto al volumen real
/*
        if (!(gsl_fcmp(alpha_im, 1.0, 0.15 / 2) == 0)) {
          wasora_push_error_message("cell %d has been poorly tracked with the defined tracking parameters in '%s':\n \
real volume \t %e \n \
appr volume \t %e \n \
alpha_im    \t %e", tracks->mesh->cell[segment->element->cell->id - 1].id, tracks->name, tracks->mesh->cell[segment->element->cell->id - 1].volume , tracks->volumes[segment->element->cell->id - 1], alpha_im);
          return WASORA_RUNTIME_ERROR;
        }
*/
        
        // actualizo la longitud
        segment->length *= alpha_im;
      }
    }
    
    // y reinicio el volumen
    for (k = 0; k < tracks->mesh->n_cells; k++) {
      tracks->volumes[k] = 0;
    }
  }
  
  // computamos los volumenes con los tracks corregidos, y los guardamos
  for (i = 0; i < tracks->n_total_tracks; i++) {
    
    track = tracks->track_by_id[i];
    
    LL_FOREACH(track->associated_segments, associated_segment) {
      segment = associated_segment->segment;
      // el peso azimutal por dos debido a que no barremos en la direccion bwd
      tracks->volumes[segment->element->cell->id - 1] += 2 * tracks->quadrature->azimuthal->w[track->azim_index] * tracks->quadrature->effective_spacing[track->azim_index] * segment->length;
    }
  }
  
  // checkeamos que efectivamente se cumpla la igualdad
  for (i = 0; i < tracks->mesh->n_cells; i++) {
    
    // should I use tiny step as eps here?
    if (!(gsl_fcmp(tracks->volumes[i], tracks->mesh->cell[i].volume, 1e-10) == 0)) {
      wasora_push_error_message("trivialmente no deberias entrar aca");
      return WASORA_RUNTIME_ERROR;
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int tracks_set_tracks_boundary_conditions(tracks_t *tracks) {
  
  int i, j;
  int parallels;
  
  double ABC[3];
  double x_int[3];
  double eps;
  double d1, d2;
  double *x1;
  double *x2;
  
  track_t *track;
  element_t *element;
  mesh_t *mesh = tracks->mesh;
  
  // esto no es tan util...
  x_int[2] = tracks->z_value;
  
  // para cada track
  for (i = 0; i < tracks->n_total_tracks; i++) {
    
    track = tracks->track_by_id[i];
    
    // recorro todos los elementos
    for (j = 0; j < mesh->n_elements; j++) {
      
      // nos vamos si ya seteamos las condiciones de contorno de entrada y salida del track
      if (track->boundary_start != NULL && track->boundary_end != NULL) {
        break;
      }
      
      element = &mesh->element[j];
      
      // pero solo tomo los elementos de superficie
      if (element->type->dim == (tracks->mesh->bulk_dimensions-1)) {
        
        x1 = element->node[0]->x;
        x2 = element->node[1]->x;
        
        // computo su ecuacion general
        track_compute_general_equation(ABC, x1, x2);
        
        // primero asumimos que el track no es paralelo al elemento de superficie
        // tener presente que un track nunca es paralelo a un elemento de superficie (siempre hay interseccion)
        parallels = 0;
        
        // calculo la interseccion entre el track y el elemento de superficie
        // obs: siempre hay interseccion, solo resta saber si esta cae dentro del elemento de sup.
        track_compute_intersection(x_int, track->ABC, ABC, &parallels);
        
        // vemos que nunca se verifique que un track es paralelo a alguna superficie
        if (parallels) {
          wasora_push_error_message("track '%d' is parallel to surface element '%d'", track->id, element->id);
          return WASORA_RUNTIME_ERROR;
        }
        
        // chekeamos si la interseccion se da sobre el dominio del elemento:
        // no importa si exactamente no es este el elemento que contiene a la entrada\salida del track
        // ya que todos los elementos sobre una superficie de la bounding box poseen la misma boundary
        eps = 1e-10;
        if (!(wasora_point_in_segment(x1, x2, x_int, eps))) continue;
        
        // se determina si la condicion de contorno es de entrada o salida
        d1 = mesh_subtract_squared_module2d(x_int, track->p_i);
        d2 = mesh_subtract_squared_module2d(x_int, track->p_o);
        if (d1 < d2) {
          track->boundary_start = element->physical_entity;
        } else {
          track->boundary_end = element->physical_entity;
        }
      }
    }
  }
  
  for (i = 0; i < tracks->n_total_tracks; i++) {
    
    track = tracks->track_by_id[i];
    
    // check de que se hayan puesto todas las boundary
    if (track->boundary_start == NULL) {
      wasora_push_error_message("forward boundary condition unassigned to track %d", track->id);
      return WASORA_RUNTIME_ERROR;
    } else if (track->boundary_end == NULL) {
      wasora_push_error_message("backward boundary condition unassigned to track %d", track->id);
      return WASORA_RUNTIME_ERROR;
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int tracks_set_tracks_next_tracks(tracks_t *tracks) {
  
  int i, j, k;
  
  track_t *track;
  
  
  for (i = 0; i < tracks->n_azim_2; i++) {
    
    // el indice complementario
    k = tracks->n_azim_2 - i - 1;
    
    for (j = 0; j < tracks->n_tracks[i]; j++) {
      
      track = &tracks->track[i][j];
      
      // dado un track, busco su next track en direccion fwd.
      // miramos los que llegan al eje "y" que apunta en su sentido fwd
      if (j < tracks->n_tracks_y[i]) {
        track->dir_next_track_fwd = 0;
        if(track->boundary_end->bc_type_phys == BC_PERIODIC) {
          track->next_track_fwd = &tracks->track[i][j + tracks->n_tracks_x[i]];
          // y si la bc es reflectiva o vacio
        } else {
          track->next_track_fwd = &tracks->track[k][j + tracks->n_tracks_x[i]];
        }
        
        // miramos los que llegan al eje "x" superior (siguiendo su sentido fwd)
      } else {
        if(track->boundary_end->bc_type_phys == BC_PERIODIC) {
          track->dir_next_track_fwd = 0;
          track->next_track_fwd = &tracks->track[i][j - tracks->n_tracks_y[i]];
          // y si la bc es reflectiva o vacio
        } else {
          track->dir_next_track_fwd = 1;
          track->next_track_fwd = &tracks->track[k][tracks->n_tracks[i] + tracks->n_tracks_y[i] - j - 1];
        }
      }
      
      // dado un track, busco su next track en direccion bwd.
      // miramos los que llegan al eje "x" inferior (siguiendo su sentido bwd)
      if (j < tracks->n_tracks_x[i]) {
        if(track->boundary_end->bc_type_phys == BC_PERIODIC) {
          track->dir_next_track_bwd = 1;
          track->next_track_bwd = &tracks->track[i][j + tracks->n_tracks_y[i]];
          // y si la bc es reflectiva o vacio
        } else {
          track->dir_next_track_bwd = 0;
          track->next_track_bwd = &tracks->track[k][tracks->n_tracks_x[i] - j - 1];
        }
        
        // miramos los que llegan al eje "y" que apunta en su sentido bwd
      } else {
        track->dir_next_track_bwd = 1;
        if(track->boundary_end->bc_type_phys == BC_PERIODIC) {
          track->next_track_bwd = &tracks->track[i][j - tracks->n_tracks_x[i]];
          // y si la bc es reflectiva o vacio
        } else {
          track->next_track_bwd = &tracks->track[k][j - tracks->n_tracks_x[i]];
        }
      }
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int track_compute_segments_tau(tracks_t *tracks) {
  
  int i, g, g_prime;
  
  double sigma_t;
  
  track_t *track;
  segment_list_item_t *associated_segment;
  segment_t *segment;
  cell_t *cell;
  xs_t *material_xs;
  
  tracks->max_tau = 0;
  for (i = 0; i < tracks->n_total_tracks; i++) {
    
    track = tracks->track_by_id[i];
    
    LL_FOREACH(track->associated_segments, associated_segment) {
      
      segment = associated_segment->segment;
      cell = segment->element->cell;
      material_xs = (xs_t *) (segment->element->physical_entity->material->ext);
      
      // por si tenemos expresiones dependientes de x,y
      wasora_value(wasora_mesh.vars.x) = cell->x[0];
      wasora_value(wasora_mesh.vars.y) = cell->x[1];
      wasora_value(wasora_mesh.vars.z) = cell->x[2];
      
      // si es la primera vez que entra el segmento aqui, hay que inicializar
      if (!segment->initialized_tau) segment->tau = calloc(milonga.groups, sizeof(double));
      
      for (g = 0; g < milonga.groups; g++) {
        
        //computo sigma total
        if (material_xs->SigmaT[g]->n_tokens != 0) {
          sigma_t = wasora_evaluate_expression(material_xs->SigmaT[g]);
        } else {
          sigma_t = wasora_evaluate_expression(material_xs->SigmaA[g]);
          for (g_prime = 0; g_prime < milonga.groups; g_prime++) {
            sigma_t += wasora_evaluate_expression(material_xs->SigmaS0[g][g_prime]);
          }
        }
        
        segment->tau[g] = sigma_t * segment->length;
        segment->initialized_tau = 1;
        
        // y guardo el maximo tau
        if (segment->tau[g] > tracks->max_tau) tracks->max_tau = segment->tau[g];
      }
    }
  }
  
  return WASORA_RUNTIME_OK;
}

int track_segmentize_segments(tracks_t *tracks, double max_tau) {
  
  int i, j;
  int g;
  int n_sub_segments;
  
  double max_segment_tau;
  double new_length;
  double dummy;
  
  double new_segment_start_point[3];
  double new_segment_end_point[3];
  
  track_t *track;
  segment_list_item_t *associated_segment;
  segment_list_item_t *temporary_segment;
  segment_t *segment;
  segment_t *new_segment;
  
  for (i = 0; i < tracks->n_total_tracks; i++) {
    
    track = tracks->track_by_id[i];
    
    LL_FOREACH_SAFE(track->associated_segments, associated_segment, temporary_segment) {
      
      segment = associated_segment->segment;
      
      max_segment_tau = 0;
      for (g = 0; g < milonga.groups; g++) {
        // nos quedamos con el tau maximo del segmento (para cualquier g)
        if (segment->tau[g] > max_segment_tau) max_segment_tau = segment->tau[g]; 
      }
      
      // si este camino optico es mayor al permitido, hay que hacer un split del segmento
      if (max_segment_tau > max_tau) {
        n_sub_segments = ceil(max_segment_tau / max_tau);
      } else {
        // sino, no tenemos que segmentar este segmento
        n_sub_segments = 1;
        // y nos vamos a otro segmento
        continue;
      }
      
      // sino dividimos al segmento en j subsegmentos
      new_length = segment->length / n_sub_segments;
      dummy =  mesh_subtract_module(segment->p_o, segment->p_i) / n_sub_segments;
      for (j = 0; j < n_sub_segments; j++) {
        
        // computamos los puntos de inicio y final del nuevo segmento
        wasora_set_point_coords(new_segment_start_point, segment->p_i[0], segment->p_i[1], segment->p_i[2]);
        wasora_adjust_point_coords(new_segment_start_point, j * dummy, M_PI_2, track->phi);
                
        wasora_set_point_coords(new_segment_end_point, new_segment_start_point[0], new_segment_start_point[1], new_segment_start_point[2]);
        wasora_adjust_point_coords(new_segment_end_point, dummy, M_PI_2, track->phi);
        
        // creamos un nuevo segmento y lo insertamos antes del que estamos subdividiendo
        track_create_segment(&new_segment, new_segment_start_point, new_segment_end_point, segment->element);
        // pero le sobreescribimos una longitud ficticia acomodada segun la correccion 
        new_segment->length = new_length;
        track_prepend_elem_segment_to_list(&track->associated_segments, &associated_segment, new_segment);
      }
      
      // eliminamos el segmento segmentizado
      LL_DELETE(track->associated_segments, associated_segment);
    }  
  }
  
  return WASORA_RUNTIME_OK;
}

int track_compute_general_equation(double ABC[], const double *x1, const double *x2) {
  
  double norm;
  
  // ecuacion general del track
  wasora_set_point_coords(ABC, 
         x1[1] - x2[1],
         x2[0] - x1[0],
         x1[0] * x2[1] - x2[0] * x1[1]);
  
  // normalizo para evitar problemas numericos
  norm = gsl_hypot(ABC[0], ABC[1]);
  
  wasora_set_point_coords(ABC,
         ABC[0] / norm, 
         ABC[1] / norm, 
         ABC[2] / norm);
  
  return WASORA_RUNTIME_OK;
}

int track_compute_intersection(double x_int[], const double ABC_1[], const double ABC_2[], int *parallels) {
  
  // calculamos el determinante
  double det = ABC_1[1] * ABC_2[0] - ABC_2[1] * ABC_1[0];
  
  // si existe (no son paralelas ni coincidentes), computo la interseccion
  if (gsl_finite(1.0 / det)) {
    x_int[0] = (ABC_1[2] * ABC_2[1] - ABC_2[2] * ABC_1[1]) / det;
    x_int[1] = (ABC_1[0] * ABC_2[2] - ABC_2[0] * ABC_1[2]) / det;
    x_int[2] = 0;
    *parallels = 0;
  } else {
    *parallels = 1;
  }
  
  return WASORA_RUNTIME_OK;
}

int track_create_segment(segment_t **segment, const double *start_point, const double *end_point, element_t *element) {
  
  *segment = calloc(1, sizeof(segment_t));
  
  wasora_set_point_coords((*segment)->p_i, start_point[0], start_point[1], start_point[2]);
  wasora_set_point_coords((*segment)->p_o, end_point[0], end_point[1], end_point[2]);
  
  (*segment)->length = mesh_subtract_module((*segment)->p_o, (*segment)->p_i);
  
  (*segment)->element = element;
  
  return WASORA_RUNTIME_OK;
}

int track_append_segment_to_list(segment_list_item_t **list, segment_t *segment) {
  segment_list_item_t *item = calloc(1, sizeof(segment_list_item_t));
  item->segment = segment;
  LL_APPEND(*list, item);
  
  return WASORA_RUNTIME_OK;
}

int track_prepend_elem_segment_to_list(segment_list_item_t **list, segment_list_item_t **element, segment_t *segment) {
  segment_list_item_t *item = calloc(1, sizeof(segment_list_item_t));
  item->segment = segment;
  LL_PREPEND_ELEM(*list, *element, item);
  
  return WASORA_RUNTIME_OK;
}

int milonga_instruction_track_post(void *arg) {
  
  track_post_t *track_post = (track_post_t *)arg;
  
  if (track_post->file->pointer == NULL) {
    if (wasora_instruction_open_file(track_post->file) != WASORA_RUNTIME_OK) {
      return WASORA_RUNTIME_ERROR;
    }
  }
  
  // solo hacemos algo si el archivo es nuevo
  if (ftell(track_post->file->pointer) == 0) {
    wasora_call(track_post->write_header(track_post));
    if (track_post->no_mesh == 0) {
      wasora_call(track_post->write_mesh(track_post));
    }
    wasora_call(track_post->write_tracks(track_post));
  }
  
  return WASORA_RUNTIME_OK;
}


int track_gnuplot_write_header (track_post_t *track_post) {
  
  fprintf(track_post->file->pointer, "# milonga ray tracing file for gnuplot\n\n");
  
  if (track_post->no_mesh == 0) {
    fprintf(track_post->file->pointer, "plot '-' u 1:2 w l lw 1 lt rgb 'black' ti 'mesh %s', '-' u 1:2:3 w lp palette ti 'tracks segments'\n", track_post->tracks->mesh->name);
  } else {
    fprintf(track_post->file->pointer, "plot '-' u 1:2:3 w lp palette ti 'tracks segments'\n");
  }
  
  fprintf(track_post->file->pointer, "\n");
  
  return WASORA_RUNTIME_OK;
}

int track_gnuplot_write_mesh(track_post_t *track_post) {
  
  int i, j;
  
  double *x1;
  double *x2;
  
  element_t *element;
  
  // imprimo los elementos
  // por ej: tendremos lineas de superficie y triangulos. imprimo sus caras (por lo que se repetiran segmentos)
  for (i = 0; i < track_post->tracks->mesh->n_elements; i++) {
    
    element = &track_post->tracks->mesh->element[i];
    
    for (j = 0; j < element->type->faces; j++) {
      
      x1 = element->node[j]->x;
      x2 = element->node[(j == (element->type->faces-1)) ? 0 : j+1]->x;
      
      fprintf(track_post->file->pointer, "%e\t%e\n", x1[0], x1[1]);
      fprintf(track_post->file->pointer, "%e\t%e\n", x2[0], x2[1]);
      fprintf(track_post->file->pointer, "\n");
    }
  }
  
  fprintf(track_post->file->pointer, "end\n\n");
  
  return WASORA_RUNTIME_OK;
}

int track_gnuplot_write_tracks(track_post_t *track_post) {
  
  int i, j;
  
  track_t *track;
  segment_list_item_t *associated_segment;
  segment_t * segment;
  
  // si no se pide ni un track ni un angulo azimutal particular, imprimimos todos
  if (track_post->track_id == 0 && track_post->azim_id == 0) {
    
    for (i = 0; i < track_post->tracks->n_azim_2; i++) {
      for (j = 0; j < track_post->tracks->n_tracks[i]; j++) {
        
        track = &track_post->tracks->track[i][j];
        
/*
        fprintf(track_post->file->pointer, "%e\t%e\t%d\n", track->p_i[0], track->p_i[1], track->id);
        fprintf(track_post->file->pointer, "%e\t%e\t%d\n", track->p_o[0], track->p_o[1], track->id);
        fprintf(track_post->file->pointer, "\n");
*/
        
        LL_FOREACH(track->associated_segments, associated_segment) {
          segment = associated_segment->segment;
          fprintf(track_post->file->pointer, "%e\t%e\t%d\n", segment->p_i[0], segment->p_i[1], segment->element->id);
          fprintf(track_post->file->pointer, "%e\t%e\t%d\n", segment->p_o[0], segment->p_o[1], segment->element->id);
          fprintf(track_post->file->pointer, "\n");
        }
      }
    }
    
    // si se pide un track_id especifico, imprimimos solo ese track
  } else if (track_post->track_id != 0 && track_post->azim_id == 0) {
    
    if (track_post->track_id > track_post->tracks->n_total_tracks)  {
      wasora_push_error_message("the given track id '%d' exceeds the amount of total tracks '%d'", track_post->track_id, track_post->tracks->n_total_tracks);
      return WASORA_RUNTIME_ERROR;
    }
    
    track = track_post->tracks->track_by_id[track_post->track_id - 1];
    LL_FOREACH(track->associated_segments, associated_segment) {
      segment = associated_segment->segment;
      fprintf(track_post->file->pointer, "%e\t%e\t%d\n", segment->p_i[0], segment->p_i[1], segment->element->id);
      fprintf(track_post->file->pointer, "%e\t%e\t%d\n", segment->p_o[0], segment->p_o[1], segment->element->id);
      fprintf(track_post->file->pointer, "\n");
    }
    
    // pero si se pide un azim id especifico, imprimimos esos tracks
  } else if (track_post->track_id == 0 && track_post->azim_id != 0) {
    
    if (track_post->azim_id > track_post->tracks->n_azim)  {
      wasora_push_error_message("the given azimuthal id '%d' exceeds the amount of total azimuthal angles '%d'", track_post->azim_id, track_post->tracks->n_azim);
      return WASORA_RUNTIME_ERROR;
    }
    
    // corregimos el angulo azimutal si estamos mirando el 3er o 4to cuadrante
    if (track_post->azim_id > track_post->tracks->n_azim_2) {
      track_post->azim_id /= 2;
    }
    
    // imprimimos los tracks con este angulo azimutal
    for (j = 0; j < track_post->tracks->n_tracks[track_post->azim_id - 1]; j++) {
      
      track = &track_post->tracks->track[track_post->azim_id - 1][j];
      LL_FOREACH(track->associated_segments, associated_segment) {
        segment = associated_segment->segment;
        fprintf(track_post->file->pointer, "%e\t%e\t%d\n", segment->p_i[0], segment->p_i[1], segment->element->id);
        fprintf(track_post->file->pointer, "%e\t%e\t%d\n", segment->p_o[0], segment->p_o[1], segment->element->id);
        fprintf(track_post->file->pointer, "\n");
      }
    }
  }
  
  fprintf(track_post->file->pointer, "end\n");
  fprintf(track_post->file->pointer, "pause -1");
  
  return WASORA_RUNTIME_OK;
}

int milonga_instruction_polar_quadrature(void *arg) {
  
  polar_quadrature_t *polar_quad = (polar_quadrature_t *)arg;
  
  // en este caso la cuadratura no es por defecto
  wasora_call(track_init_polar_quadrature(polar_quad));
  
  return WASORA_RUNTIME_OK;
}

int wasora_set_point_coords(double *xyz, const double x, const double y, const double z) {
  
  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = z;
  
  return WASORA_RUNTIME_OK;
}

int wasora_adjust_point_coords(double *xyz, const double delta, const double theta, const double phi) {
  
  if (theta == M_PI_2) {
    xyz[0] += delta * cos(phi);
    xyz[1] += delta * sin(phi);
  } else {
    xyz[0] += delta * sin(theta) * cos(phi);
    xyz[1] += delta * sin(theta) * sin(phi);
    xyz[2] += delta * cos(theta);
  }
  
  return WASORA_RUNTIME_OK;
}

int wasora_point_in_segment(const double *x1, const double *x2, const double *xp, double eps) {
  
  return (((gsl_fcmp(x1[0], xp[0], eps) <= 0) && (gsl_fcmp(xp[0], x2[0], eps) <= 0)  || 
           (gsl_fcmp(x2[0], xp[0], eps) <= 0) && (gsl_fcmp(xp[0], x1[0], eps) <= 0)) &&
          ((gsl_fcmp(x1[1], xp[1], eps) <= 0) && (gsl_fcmp(xp[1], x2[1], eps) <= 0)  || 
           (gsl_fcmp(x2[1], xp[1], eps) <= 0) && (gsl_fcmp(xp[1], x1[1], eps) <= 0)));
    
}
