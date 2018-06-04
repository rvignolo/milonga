/*------------ -------------- -------- --- ----- ---   --       -            -
 * milonga's polar quadrature protoypes
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

typedef struct polar_quadrature_t polar_quadrature_t;

struct polar_quadrature_t {
  char *name;
  
  expr_t *expr_n_polar;
  int n_polar;                // el numero de angulos polares en (0, 1*pi)
  int n_polar_2;              // el numero de angulos polares en (0, pi/2)
  
  double *theta;              // arreglo de angulos polares
  double *sin_theta;          // arreglo de senos de angulos polares
  double *w;                  // peso polar
  
  enum {                      // diferentes tipos de cuadraturas polares en la literatura
    tabuchi_yamamoto,
    equal_weight,
    equal_angle,
    gauss_legendre,
    leonard
  } quad_type;
  
  UT_hash_handle hh;
};

struct {
  polar_quadrature_t *main_polar_quadrature;
  
  polar_quadrature_t *polar_quadratures;       // tabla de hash de polar quadratures
} milonga_polar_quadratures;

extern int milonga_instruction_polar_quadrature(void *);

extern int milonga_init_polar_quadrature(polar_quadrature_t *);
extern int milonga_init_polar_quadrature_values(polar_quadrature_t *);

extern polar_quadrature_t *milonga_define_polar_quadrature(char *, expr_t *, int);
extern polar_quadrature_t *milonga_get_polar_quadrature_ptr(const char *);