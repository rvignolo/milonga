/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's polar quadrature routines
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

#include "milonga.h"

polar_quadrature_t *milonga_define_polar_quadrature(char *name, expr_t *expr_n_polar, int polar_quad_type) {
  
  polar_quadrature_t *polar_quadrature;
  
  if (name == NULL) {
    name = strdup("polar_quadrature");
  }
  
  if (milonga_get_polar_quadrature_ptr(name) != NULL) {
    wasora_push_error_message("there already exists a polar quadrature named '%s'", name);
    return NULL;
  }
  
  polar_quadrature = calloc(1, sizeof(polar_quadrature_t));
  polar_quadrature->name = strdup(name);
  polar_quadrature->expr_n_polar = expr_n_polar;
  polar_quadrature->quad_type = polar_quad_type;
  
  // la agregamos al hash
  HASH_ADD_KEYPTR(hh, milonga_polar_quadratures.polar_quadratures, polar_quadrature->name, strlen(polar_quadrature->name), polar_quadrature);
  // y seteamos la principal como esta
  milonga_polar_quadratures.main_polar_quadrature = polar_quadrature;
  
  return polar_quadrature;
}

polar_quadrature_t *milonga_get_polar_quadrature_ptr(const char *name) {
  polar_quadrature_t *polar_quadrature;
  HASH_FIND_STR(milonga_polar_quadratures.polar_quadratures, name, polar_quadrature);
  return polar_quadrature;
}

int milonga_instruction_polar_quadrature(void *arg) {
  
  polar_quadrature_t *polar_quadrature = (polar_quadrature_t *)arg;
  
  // en este caso la cuadratura no es por defecto
  wasora_call(milonga_init_polar_quadrature(polar_quadrature));
  
  return WASORA_RUNTIME_OK;
}

int milonga_init_polar_quadrature(polar_quadrature_t *polar_quadrature) {
  
  // si tenemos que armar una por defecto
  if (polar_quadrature->expr_n_polar == NULL) {
    
    // por defecto, tabuchi-yamamoto con 6 polar angles
    polar_quadrature->quad_type = tabuchi_yamamoto;
    polar_quadrature->n_polar = 6;
    polar_quadrature->n_polar_2 = 6 / 2;
    
  } else {
    
    // no se si hace falta esto cuado este cambiando todo
    // if (polar->expr_n_polar != NULL)
    
    polar_quadrature->n_polar = (int) wasora_evaluate_expression(polar_quadrature->expr_n_polar);
    
    if (polar_quadrature->n_polar <= 0) {
      wasora_push_error_message("number of polar angles '%d' cannot be negative or zero", polar_quadrature->n_polar);
      return WASORA_RUNTIME_ERROR;
    } else if (polar_quadrature->n_polar % 2 != 0) {
      wasora_push_error_message("number of polar angles '%d' is not a multiple of 2", polar_quadrature->n_polar);
      return WASORA_RUNTIME_ERROR;
    } else {
      polar_quadrature->n_polar_2 = polar_quadrature->n_polar / 2;
    }
  }
  
  wasora_call(milonga_init_polar_quadrature_values(polar_quadrature));
  
  return WASORA_RUNTIME_OK;
}

int milonga_init_polar_quadrature_values(polar_quadrature_t *polar_quadrature) {
  
  int i;
  
  polar_quadrature->theta = calloc(polar_quadrature->n_polar, sizeof(double *));
  polar_quadrature->sin_theta = calloc(polar_quadrature->n_polar, sizeof(double *));
  polar_quadrature->w = calloc(polar_quadrature->n_polar, sizeof(double *));
  
  switch (polar_quadrature->quad_type) {
    case tabuchi_yamamoto:
      // del Handbook de Ingenieria Nuclear pg. 1154
      if (polar_quadrature->n_polar == 2) {
        
        polar_quadrature->sin_theta[0] = 0.798184;
        polar_quadrature->w[0] =   1.0 / 2.0;
        
      } else if (polar_quadrature->n_polar == 4) {
        
        polar_quadrature->sin_theta[0] = 0.363900;
        polar_quadrature->w[0] =   0.212854 / 2.0;
        
        polar_quadrature->sin_theta[1] = 0.899900;
        polar_quadrature->w[1] =   0.787146 / 2.0;
        
      } else if (polar_quadrature->n_polar == 6) {
        
        polar_quadrature->sin_theta[0] = 0.166648;
        polar_quadrature->w[0] =   0.046233 / 2.0;
        
        polar_quadrature->sin_theta[1] = 0.537707;
        polar_quadrature->w[1] =   0.283619 / 2.0;
        
        polar_quadrature->sin_theta[2] = 0.932954;
        polar_quadrature->w[2] =   0.670148 / 2.0;
        
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
      if (polar_quadrature->n_polar == 2) {
        
        polar_quadrature->theta[0]   = acos(0.5773502691);
        polar_quadrature->w[0] = 1.0 / 2.0;
        
      } else if (polar_quadrature->n_polar == 4) {
        
        polar_quadrature->theta[0]   = acos(0.3399810435);
        polar_quadrature->w[0] = 0.6521451549 / 2.0;
        
        polar_quadrature->theta[1]   = acos(0.8611363115);
        polar_quadrature->w[1] = 0.3478548451 / 2.0;
        
      } else if (polar_quadrature->n_polar == 6) {
        
        polar_quadrature->theta[0]   = acos(0.2386191860);
        polar_quadrature->w[0] = 0.4679139346 / 2.0;
        
        polar_quadrature->theta[1]   = acos(0.6612093864);
        polar_quadrature->w[1] = 0.3607615730 / 2.0;
        
        polar_quadrature->theta[2]   = acos(0.9324695142);
        polar_quadrature->w[2] = 0.1713244924 / 2.0;
        
      } else if (polar_quadrature->n_polar == 8) {
        
        polar_quadrature->theta[0]   = acos(0.1834346424);
        polar_quadrature->w[0] = 0.3626837834 / 2.0;
        
        polar_quadrature->theta[1]   = acos(0.5255324099);
        polar_quadrature->w[1] = 0.3137066459 / 2.0;
        
        polar_quadrature->theta[2]   = acos(0.7966664774);
        polar_quadrature->w[2] = 0.2223810344 / 2.0;
        
        polar_quadrature->theta[3]   = acos(0.9602898564);
        polar_quadrature->w[3] = 0.1012285363 / 2.0;
        
      } else if (polar_quadrature->n_polar == 10) {
        
        polar_quadrature->theta[0]   = acos(0.1488743387);
        polar_quadrature->w[0] = 0.2955242247 / 2.0;
        
        polar_quadrature->theta[1]   = acos(0.4333953941);
        polar_quadrature->w[1] = 0.2692667193 / 2.0;
        
        polar_quadrature->theta[2]   = acos(0.6794095682);
        polar_quadrature->w[2] = 0.2190863625 / 2.0;
        
        polar_quadrature->theta[3]   = acos(0.8650633666);
        polar_quadrature->w[3] = 0.1494513492 / 2.0;
        
        polar_quadrature->theta[4]   = acos(0.9739065285);
        polar_quadrature->w[4] = 0.0666713443 / 2.0;
        
      } else if (polar_quadrature->n_polar == 12) {
        
        polar_quadrature->theta[0]   = acos(0.1252334085);
        polar_quadrature->w[0] = 0.2491470458 / 2.0;
        
        polar_quadrature->theta[1]   = acos(0.3678314989);
        polar_quadrature->w[1] = 0.2334925365 / 2.0;
        
        polar_quadrature->theta[2]   = acos(0.5873179542);
        polar_quadrature->w[2] = 0.2031674267 / 2.0;
        
        polar_quadrature->theta[3]   = acos(0.7699026741);
        polar_quadrature->w[3] = 0.1600783286 / 2.0;
        
        polar_quadrature->theta[4]   = acos(0.9041172563);
        polar_quadrature->w[4] = 0.1069393260 / 2.0;
        
        polar_quadrature->theta[5]   = acos(0.9815606342);
        polar_quadrature->w[5] = 0.0471753364 / 2.0;
        
      } else {
        wasora_push_error_message("'gauss_legendre' quadrature not suported for N_POLAR_ANGLES superior than 12");
        return WASORA_RUNTIME_ERROR;
      }
      break;
    case leonard:
      // ver de que referencia sacar esto (paper tabuchi es diferente a openmoc y handbook tiene otras optimizadas por Herbert)
      // estas se corresponden con el paper de tabuchi
      if (polar_quadrature->n_polar == 2) {
        
        polar_quadrature->sin_theta[0] = 0.752244;
        polar_quadrature->w[0] =   1.0 / 2.0;
        
      } else if (polar_quadrature->n_polar == 4) {
        
        polar_quadrature->sin_theta[0] = 0.273658;
        polar_quadrature->w[0] =   0.139473 / 2.0;
        
        polar_quadrature->sin_theta[1] = 0.865714;
        polar_quadrature->w[1] =   0.860527 / 2.0;
        
      } else if (polar_quadrature->n_polar == 6) {
        
        polar_quadrature->sin_theta[0] = 0.103840;
        polar_quadrature->w[0] =   0.020530 / 2.0;
        
        polar_quadrature->sin_theta[1] = 0.430723;
        polar_quadrature->w[1] =   0.219161 / 2.0;
        
        polar_quadrature->sin_theta[2] = 0.905435;
        polar_quadrature->w[2] =   0.760309 / 2.0;
        
      } else {
        wasora_push_error_message("'leonard' quadrature not suported for N_POLAR_ANGLES superior than 6");
        return WASORA_RUNTIME_ERROR;
      }
      break;
  }
  
  // si nos dieron los senos
  if (polar_quadrature->quad_type == tabuchi_yamamoto || polar_quadrature->quad_type == leonard) {
    // seteamos las componentes que faltan: rellenamos los sin_theta, los pesos y calculamos los angulos 
    for (i = 0; i < polar_quadrature->n_polar_2; i++) {
      polar_quadrature->sin_theta[polar_quadrature->n_polar-i-1] = polar_quadrature->sin_theta[i];
      polar_quadrature->w[polar_quadrature->n_polar-i-1] = polar_quadrature->w[i];
      
      polar_quadrature->theta[i] = asin(polar_quadrature->sin_theta[i]);
      polar_quadrature->theta[polar_quadrature->n_polar-i-1] = M_PI - polar_quadrature->theta[i];
    }
  }
  
  // si nos dieron los angulos
  if (polar_quadrature->quad_type == gauss_legendre) {
    // seteamos las componentes que faltan: rellenamos los theta, los pesos y calculamos los sin_theta 
    for (i = 0; i < polar_quadrature->n_polar_2; i++) {
      polar_quadrature->theta[polar_quadrature->n_polar-i-1] = M_PI - polar_quadrature->theta[i];
      polar_quadrature->w[polar_quadrature->n_polar-i-1] = polar_quadrature->w[i];
      
      polar_quadrature->sin_theta[i] = sin(polar_quadrature->theta[i]);
      polar_quadrature->sin_theta[polar_quadrature->n_polar-i-1] = sin(polar_quadrature->theta[polar_quadrature->n_polar-i-1]);
    }
  }
  
  return WASORA_RUNTIME_OK;
}