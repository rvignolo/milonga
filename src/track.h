/*------------ -------------- -------- --- ----- ---   --       -            -
 * milonga's unstructured tracking protoypes
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

#define MAX_TAU         10.0
#define EXP_PRECISION   1e-5

typedef struct exp_evaluator_t exp_evaluator_t;
typedef struct segment_t segment_t;
typedef struct segment_list_item_t segment_list_item_t;
typedef struct azimuthal_quadrature_t azimuthal_quadrature_t;
typedef struct polar_quadrature_t polar_quadrature_t;
typedef struct quadrature_t quadrature_t;
typedef struct track_t track_t;
typedef struct tracks_t tracks_t;
typedef struct track_post_t track_post_t;

struct exp_evaluator_t {
  int size;                 // tamanyo de la tabla
  
  double max_tau_allowed;   // maximo camino optico permitido
  double max_err_allowed;   // maximo error para la aproximacion
  
  double delta;             // longitud del intervalo de interpolacion
  
  double *table;            // tabla para interpolar
  
  double (*compute_exponential)(double, int);   // apuntador a donde corresponda
};

struct segment_t {
  double p_i[3];           // punto inicial del segmento
  double p_o[3];           // punto final del segmento
  
  double length;           // longitud del segmento
  
  int initialized_tau;     // flag para indiciar que el vector tau ha sido inicializado (se agregan nuevos segmentos en segmentize_segments)
  double *tau;             // camino optico del segmento para cada grupo
  
  element_t *element;      // asociando al elemento tengo informacion de materiales, ids, boundaries, etc
};

struct segment_list_item_t {
  segment_t *segment;
  segment_list_item_t *next;
};

struct azimuthal_quadrature_t {
  int n_azim;                 // el numero de angulos azimutales en (0, 2*pi)
  int n_azim_2;               // el numero de angulos azimutales en (0, 1*pi)
  int n_azim_4;               // el numero de angulos azimutales en (0, pi/2)
  
  double *phi;                // arreglo de angulos azimutales
  double *w;                  // peso azimutal
};

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

struct quadrature_t {
  
  azimuthal_quadrature_t *azimuthal;  // azimuthal quadrature 
  
  polar_quadrature_t *polar;          // polar quadrature
  
  double *effective_spacing;          // espacio entre tracks (cm) para un dado angulo azimutal
  
  double **w_total;                   // para llamar comodamente a los pesos (es funcion del indice azimutal y polar)
};

struct track_t {
  int id;                     // id universal del track
  
  int azim_index;             // indice azimutal
  
  double ABC[3];              // ecuacion general del track
  double phi;                 // angulo azimutal del track
  
  double p_i[3];              // punto de partida del track
  double p_o[3];              // punto de salida del track
  double length;              // longitud del track
  
  segment_list_item_t *associated_segments;  // linked list de segmentos asociados al track
  
  physical_entity_t *boundary_start;         // boundary en el p_i del track
  physical_entity_t *boundary_end;           // boundary en el p_o del track
  
  track_t *next_track_fwd;    // proximo track al recorrer el track actual en direccion forward  y dependiendo de la bc
  track_t *next_track_bwd;    // proximo track al recorrer el track actual en direccion backward y dependiendo de la bc
  
  int dir_next_track_fwd;     // si el proximo track en direccion fwd es recorrido en su sentido fwd (0) o bwd (1) por como elegi angular_index
  int dir_next_track_bwd;     // si el proximo track en direccion bwd es recorrido en su sentido fwd (0) o bwd (1) por como elegi angular_index
};

struct tracks_t {
  char *name;
  int initialized;
  
  mesh_t *mesh;
  
  expr_t *expr_n_azim;
  int n_azim;                   // el numero de angulos azimutales en (0, 2*pi)
  int n_azim_2;                 // el numero de angulos azimutales en (0, 1*pi)
  int n_azim_4;                 // el numero de angulos azimutales en (0, pi/2)
  
  expr_t *expr_track_dens;
  double track_dens;            // densidad de lineas (1/cm)
  expr_t *expr_track_spacing;
  double track_spacing;         // espacio entre lineas (cm)
  
  expr_t *expr_tiny_step;
  double tiny_step;             // distancia (cm) recorrida hacia dentro de un elemento en la direccion del track
  
  int *n_tracks_x;              // cantidad de tracks que arrancan en el eje x para un dado angulo azimutal
  int *n_tracks_y;              // cantidad de tracks que arrancan en el eje y para un dado angulo azimutal
  int *n_tracks;                // cantidad de tracks para un dado angulo azimutal
  int n_total_tracks;           // cantidad total de tracks
  
  track_t **track;              // depende del indice del angulo azimutal y del numero de track
  track_t **track_by_id;        // para buscar cada track pero con un universal id (hacemos loops de una)
  
  int do_not_correct_volumes;   // flag que inidica si checkear el volumen real frente al aprox por el ray tracing para cada celda
  double *volumes;              // volumenes de cada celda aproximados con los tracks
  
  int debug;                    // si se detecta un error en el ray tracing, esta opcion permite dibujar el track con dicha inconsistencia
  
  double max_tau;               // camino optico maximo del ray tracing antes y despues de segmentar segmentos. \
                                   al principio nos sirve para ver si segmentamos y luego se actualiza para saber \
                                   cuantos intervalos de interpolacion necesitamos
  
  double z_value;               // coordenada axial (esto es mas dummy... quizas lo vuele)
  
  quadrature_t *quadrature;     // cuadratura para computar integrales
  
  UT_hash_handle hh;
};

struct track_post_t {
  file_t *file;             // archivo donde se guarda la informacion
  
  tracks_t *tracks;         // ray tracing a plotear (tambien contiene la malla del ray tracing)
  
  int no_mesh;              // flag que indica si ploteamos o no la malla
  int track_id;             // si se pide un unico track, este es su id
  int azim_id;              // si se pide un id azimutal especifico
  
  enum  {
    post_format_from_extension,
    post_format_gnuplot,
  } format;
  
  int (*write_header)(track_post_t *);   // escribe las instrucciones para plotear
  int (*write_mesh)(track_post_t *);     // escribe los segmentos de la malla
  int (*write_tracks)(track_post_t *);   // escribe los segmentos de los tracks
};

struct {
  tracks_t *main_ray_tracing;
  polar_quadrature_t *main_polar_quadrature;
  
  tracks_t *ray_tracings;                      // tabla de hash de trackings
  polar_quadrature_t *polar_quadratures;       // tabla de hash de polar quadratures
} tracking;

extern int milonga_instruction_track(void *);

extern int track_read_tracking_parameters(tracks_t *);
extern int track_init_quadratures_before_tracking(tracks_t *);
extern int track_compute_tracks(tracks_t *);
extern int track_init_quadratures_after_tracking(tracks_t *);
extern int track_init_azimuthal_weights(azimuthal_quadrature_t *);
extern int track_init_polar_quadrature(polar_quadrature_t *);
extern int track_init_polar_quadrature_values(polar_quadrature_t *);
extern int track_recalibrate_coords(tracks_t *);
extern int track_segmentize_tracks(tracks_t *);
extern int track_compute_element_intersections(element_t *, track_t *, double, double, double ***);
extern int track_correct_volumes(tracks_t *);
extern int tracks_set_tracks_boundary_conditions(tracks_t *);
extern int tracks_set_tracks_next_tracks(tracks_t *);
extern int track_compute_segments_tau(tracks_t *);
extern int track_segmentize_segments(tracks_t *, double);
extern int track_compute_general_equation(double [], const double *, const double *);
extern int track_compute_intersection(double [], const double [], const double [], int *);
extern int track_create_segment(segment_t **, const double *, const double *, element_t *);
extern int track_append_segment_to_list(segment_list_item_t **, segment_t *);
extern int track_prepend_elem_segment_to_list(segment_list_item_t **, segment_list_item_t **, segment_t *);


extern int milonga_instruction_track_post(void *);

extern int track_gnuplot_write_header(track_post_t *);
extern int track_gnuplot_write_mesh(track_post_t *);
extern int track_gnuplot_write_tracks(track_post_t *);


extern int milonga_instruction_polar_quadrature(void *);


extern tracks_t *milonga_define_ray_tracing(char *, mesh_t *, expr_t *, expr_t *, expr_t *, expr_t *, int, int);
extern tracks_t *milonga_get_ray_tracing_ptr(const char *);
extern polar_quadrature_t *milonga_define_polar_quadrature(char *, expr_t *, int);
extern polar_quadrature_t *milonga_get_polar_quadrature_ptr(const char *);


extern int wasora_set_point_coords(double *, const double, const double, const double);
extern int wasora_adjust_point_coords(double *, const double, const double, const double);
extern int wasora_point_in_segment(const double *, const double *, const double *, double);
