#include <stdio.h>
#include <math.h>
#include "idl_export.h"
#ifdef _OPENMP
	#include <omp.h>
#endif

IDL_VPTR build_beam_worker(int argc, IDL_VPTR argv[]) {
/* Declare variables from IDL */
IDL_VPTR drift_dec, drift_err, obs, map;		// For arrays
IDL_VPTR pnt_dec_raw, sigma_raw;			// For scalars
double pnt_dec, sigma;
int n_drifts, n_map, n_rec;
double *drift_dec_d, *drift_err_d, *obs_d, *map_d;	// Pointers to hold/access data

/* Moving arrays to C variables */
drift_dec    = argv[0];
drift_err    = argv[1];
pnt_dec_raw  = argv[2];
obs          = argv[3];
sigma_raw    = argv[4];
map          = argv[5];

/* Data type checks */
/** Arrays **/
IDL_ENSURE_SIMPLE(drift_dec);
IDL_ENSURE_ARRAY(drift_dec);
IDL_ENSURE_SIMPLE(drift_err);
IDL_ENSURE_ARRAY(drift_err);
IDL_ENSURE_SIMPLE(obs);
IDL_ENSURE_ARRAY(obs);
IDL_ENSURE_SIMPLE(map);
IDL_ENSURE_ARRAY(map);
/** Scalars **/
IDL_ENSURE_SIMPLE(pnt_dec_raw);
IDL_ENSURE_SCALAR(pnt_dec_raw);
IDL_ENSURE_SIMPLE(sigma_raw);
IDL_ENSURE_SCALAR(sigma_raw);

/* Move scalars to type and then C variables */
if( pnt_dec_raw->type != IDL_TYP_DOUBLE ) {
	pnt_dec_raw = IDL_CvtDbl(1, &pnt_dec_raw);
}
pnt_dec = (double) pnt_dec_raw->value.d;

if( sigma_raw->type != IDL_TYP_DOUBLE ) {
	sigma_raw = IDL_CvtDbl(1, &sigma_raw);
}
sigma = (double) sigma_raw->value.d;

/* Gather together array dimensions */
n_drifts = (int) (drift_dec->value.arr->n_elts / 7);
n_map = (int) sqrt( (int) (map->value.arr->n_elts) );
n_rec = (int) ((obs->value.arr->n_elts) / (drift_dec->value.arr->n_elts));

/* Define local variables */
int i,j,l,d,b,rp;
long obs_count, dec_count;
double dec_offset, loc_y, r, cosd, x2, y2, d2, w, map_pnt;
double cnt_dec, cnt_err, cnt_obs, map_tot=0.0;

/* Define constants */
double degrad=M_PI / 180.0;

/* Load in array data into pointers */
drift_dec_d = (double *) drift_dec->value.arr->data;
drift_err_d = (double *) drift_err->value.arr->data;
obs_d = (double *) obs->value.arr->data;
map_d = (double *) map->value.arr->data;

/* Create output array */
//map_d = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE, map->value.arr->n_dim, map->value.arr->dim, IDL_ARR_INI_NOP, &map);

/* New for-loop order here.  (i,j)->i and i->l, j->m 
 * Possibility of running in parallel with OpenMP if this module is compiled that way.  
 * It doesn't work on GCC <  4.3 because of how the OpenMP library is compiled (cannot 
 * dlopen it) but it does work with ICC 10.1. For ICC I use the options:
 *	-O3 -openmp -I/usr/local/itt/idl/external/include/ -m64 -shared -fPIC -lm -vec_report2
 * to get OpenMP working and everything set so that IDL can use the library.
 */
#ifdef _OPENMP
	// These variables need to be private to each thread if we are compiled to use OpenMP
	#pragma omp parallel default(shared) private(i,j,d,b,rp,map_pnt,obs_count,dec_count,cnt_dec,cnt_err,cnt_obs,dec_offset,loc_y,r,cosd,x2,y2,w)
#endif
{
	/* Our parallel for-loop.  NOTE:  The number of threads to use is set by the
	 * shell enviroment variable OMP_NUM_THREADS. */
	#ifdef _OPENMP
		#pragma omp for schedule(dynamic, 200)
	#endif
	for(l=0; l<(n_map*n_map); l++) {
		i = l%n_map;	// old i value
		j = l/n_map;	// old j value
		
		map_pnt = 0.0;
		obs_count = 0;
		dec_count = 0;
		for(d=0; d<n_drifts; d++) {
			for(b=0; b<7; b++) {
				cnt_dec = *(drift_dec_d+dec_count);
				cnt_err = *(drift_err_d+dec_count++);
	
				dec_offset = (cnt_dec - pnt_dec) * 60.0;
				loc_y = dec_offset - 0.5*cnt_err*0.05 + (double) (n_map / 2);
				
				for(rp=0; rp<n_rec; rp++) {
					cnt_obs = *(obs_d+obs_count++);

					r = ((double) rp - 1.5) / 4.0;
					cosd = cos(cnt_dec*degrad);
					x2 = pow((r-i)*cosd,2.0);
					y2 = pow(loc_y-j,2.0);

					w = exp(-(x2+y2)/(2.0*pow(sigma,2.0)));
					if(b==0) {
						w *= 1.27;
					}
					if(x2 > 25.0 && y2 > 25.0) {
						w = 0.0;
					}

					map_pnt = map_pnt + cnt_obs*w;
				}
			}
		}

		*(map_d+l) = map_pnt;
	}
}

/* Normalize by first totaling the beam map and then dividing every element
 * by the beam total.  This way, total(beam) = 1. Again, we can do this in a
 * parallel fashion. */
# ifdef _OPENMP
	#pragma omp parallel for shared(map_tot)
#endif
for(l=0; l<(n_map*n_map); l++) {
	map_tot = map_tot + *(map_d+l);
}

#ifdef _OPENMP
	#pragma omp parallel for schedule(static)
#endif
for(l=0; l<(n_map*n_map); l++) {
	*(map_d+l) = *(map_d+l) / map_tot;
}

return(map);

}
