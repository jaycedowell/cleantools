#include <stdio.h>
#include <math.h>
#include "idl_export.h"
#ifdef _OPENMP
	#include <omp.h>
#endif

long whereMax(double* data, long offset, long elements) {
	long i;
	double max = *(data+offset);
	long index = 0;

	for(i=1; i<elements; i++) {
		if( *(data+offset + i) > max ) {
			max = *(data+offset + i);
			index = i;
		}
	}

	return index;
}

IDL_VPTR alfa_clean_worker(int argc, IDL_VPTR argv[]) {
/* Declare variables from IDL */
IDL_VPTR clean, work, beams, exit_status, beam_norm;			// For arrays
IDL_VPTR cstart_raw, cstop_raw, flimit_raw, gain_raw;			// For scalars - deconvolution control
IDL_VPTR ilimit_raw, beam_size_raw, clean_size_raw, n_chan_raw;		// For scalars - deconvolution control/array shape

double flimit, gain;
int c, cstart, cstop, i, ilimit, beam_size, clean_size, n_chan;
int *exit_status_d;
double *clean_d, *work_d, *beams_d, *beam_norm_d;	// Pointers to hold/access data in arrays

/* Moving arrays to C variables */
clean      = argv[0];
work       = argv[1];
beams      = argv[2];
cstart_raw = argv[3];
cstop_raw  = argv[4];
flimit_raw = argv[5];
ilimit_raw = argv[6];
gain_raw   = argv[7];
beam_norm  = argv[8];
beam_size_raw  = argv[9];
clean_size_raw = argv[10];
n_chan_raw     = argv[11];
exit_status    = argv[12];

/* Data type checks */
/** Arrays **/
IDL_ENSURE_SIMPLE(clean);
IDL_ENSURE_ARRAY(clean);
IDL_ENSURE_SIMPLE(work);
IDL_ENSURE_ARRAY(work);
IDL_ENSURE_SIMPLE(beams);
IDL_ENSURE_ARRAY(beams);
IDL_ENSURE_SIMPLE(exit_status);
IDL_ENSURE_ARRAY(exit_status);
IDL_ENSURE_SIMPLE(beam_norm);
IDL_ENSURE_ARRAY(beam_norm);
/** Scalars **/
IDL_ENSURE_SIMPLE(cstart_raw);
IDL_ENSURE_SCALAR(cstart_raw);
IDL_ENSURE_SIMPLE(cstop_raw);
IDL_ENSURE_SCALAR(cstop_raw);
IDL_ENSURE_SIMPLE(flimit_raw);
IDL_ENSURE_SCALAR(flimit_raw);
IDL_ENSURE_SIMPLE(ilimit_raw);
IDL_ENSURE_SCALAR(ilimit_raw);
IDL_ENSURE_SIMPLE(gain_raw);
IDL_ENSURE_SCALAR(gain_raw);
IDL_ENSURE_SIMPLE(beam_size_raw);
IDL_ENSURE_SCALAR(beam_size_raw);
IDL_ENSURE_SIMPLE(clean_size_raw);
IDL_ENSURE_SCALAR(clean_size_raw);
IDL_ENSURE_SIMPLE(n_chan_raw);
IDL_ENSURE_SCALAR(n_chan_raw);

/* Move scalars to type and then C variables */
if( cstart_raw->type != IDL_TYP_LONG ) {
	cstart_raw = IDL_CvtLng(1, &cstart_raw);
}
cstart = cstart_raw->value.l;

if( cstop_raw->type != IDL_TYP_LONG ) {
	cstop_raw = IDL_CvtLng(1, &cstop_raw);
}
cstop = cstop_raw->value.l;

if( flimit_raw->type != IDL_TYP_DOUBLE ) {
	flimit_raw = IDL_CvtDbl(1, &flimit_raw);
}
flimit = flimit_raw->value.d;

if( ilimit_raw->type != IDL_TYP_LONG ) {
	ilimit_raw = IDL_CvtLng(1, &ilimit_raw);
}
ilimit = ilimit_raw->value.l;

if( gain_raw->type != IDL_TYP_DOUBLE ) {
	gain_raw = IDL_CvtDbl(1, &gain_raw);
}
gain = gain_raw->value.d;

if( beam_size_raw->type != IDL_TYP_LONG ) {
	beam_size_raw = IDL_CvtLng(1, &beam_size_raw);
}
beam_size = beam_size_raw->value.l;

if( clean_size_raw->type != IDL_TYP_LONG ) {
	clean_size_raw = IDL_CvtLng(1, &clean_size_raw);
}
clean_size= clean_size_raw->value.l;

if( n_chan_raw->type != IDL_TYP_LONG ) {
	n_chan_raw = IDL_CvtLng(1, &n_chan_raw);
}
n_chan = n_chan_raw->value.l;

/* Define local variables */
long peak, n_work, chan_offset;
int clean_size2;
int peak_x, peak_y, x, y, n_clean, n_beam, channel_exit, l, m;
int work_x_lo, work_x_hi, work_y_lo, work_y_hi, beam_x_lo, beam_y_lo;
double peak_value;

/* Load in array dimensions */
n_work = (long) (work->value.arr->n_elts) / n_chan;
clean_size2 = (int) (clean->value.arr->n_elts) / n_chan / clean_size;

/* Load in array data into pointers */
clean_d = (double *) clean->value.arr->data;
work_d  = (double *) work->value.arr->data;
beams_d = (double *) beams->value.arr->data;
beam_norm_d = (double *) beam_norm->value.arr->data;
exit_status_d = (int *) exit_status->value.arr->data;

/* Loop over channels
 * Possibility of running in parallel with OpenMP if this module is compiled that way.  
 * It doesn't work on GCC <  4.3 because of how the OpenMP library is compiled (cannot 
 * dlopen it) but it does work with ICC 10.1. For ICC I use the options:
 *	-O3 -openmp -I/usr/local/itt/idl/external/include/ -m64 -shared -fPIC -lm -vec_report2
 * to get OpenMP working and everything set so that IDL can use the library. */
#ifdef _OPENMP
	// These variables need to be private to each thread if we are compiled to use OpenMP
	#pragma omp parallel default(shared) private(chan_offset,channel_exit, i,peak,peak_value,peak_x,peak_y, work_x_lo,work_x_hi,work_y_lo,work_y_hi,beam_x_lo,beam_y_lo, l,m)
#endif
{
	/* Our parallel for-loop.  By using shedule(static), the work is dividely ~evenly 
	 * across the different threads.  NOTE:  The number of threads to use is set by the
	 * shell enviroment variable OMP_NUM_THREADS. */
	#ifdef _OPENMP
		#pragma omp for schedule(static)
	#endif
	for(c=cstart; c<=cstop; c++) {
		chan_offset = (long) c * clean_size*clean_size2;
		channel_exit = -1;
	
		/* Loop over iterations */
		for(i=0; i<ilimit; i++) {
			// Idenfity the location and value of the global maximum
			peak = whereMax( work_d, chan_offset, n_work);
			peak_value = *(work_d+chan_offset + peak);
		
			// Check to see if we have reached the flux limit for cleaning
			if( peak_value < flimit ) {
				channel_exit = i+1;
				goto IterDone;
			}
	
			/* If we have made it this far, get the variables set up to define the region we
			 * need to clean.  We also need to make sure that we don't wrap around in the 
			 * array.  To do this, we need to set a work x and y range as well as a beam x 
			 * and y range. */
			peak_x = peak % clean_size;
			peak_y = peak / clean_size;

			/* work_x_lo -> left x coordinate where we need to subtract the map
			 * work_x_hi -> right x coordinate where we need to subtract the map
			 * work_y_lo -> lower y coordinate where we need to subtract the map
			 * work_y_hi -> upper y coordinate where we need to subtract the map
			 * beam_x_lo -> left x coordinate of the beam we need
			 * beam_y_lo -> lower y coordinate of the beam we need */
			work_x_lo = IDL_MAX((peak_x - beam_size/2), 0);
			work_x_hi = IDL_MIN((peak_x + beam_size/2), (clean_size-1));
			work_y_lo = IDL_MAX((peak_y - beam_size/2), 0);
			work_y_hi = IDL_MIN((peak_y + beam_size/2), (clean_size2-1));
			beam_x_lo = beam_size/2 - (peak_x - work_x_lo);
			beam_y_lo = beam_size/2 - (peak_y - work_y_lo);
	
			// Subtract out the scaled beam from the `work` map
			for(m=work_y_lo; m<=work_y_hi; m++) {
				for(l=work_x_lo; l<=work_x_hi; l++) {
					*(work_d+chan_offset + l+m*clean_size) = *(work_d+chan_offset + l+m*clean_size) - *(beams_d+peak + ((beam_x_lo+l-work_x_lo)+(beam_y_lo+m-work_y_lo)*beam_size)*n_work) * gain*peak_value;
				}
			}
	
			/* Add that component to `clean` map.  Relative to the subtraction above, this is 
			 * easy because the dimensions of *(clean_d+chan_offset) and *(beam_norm_d) are the
			 * same. */
			*(clean_d+chan_offset + peak) = *(clean_d+chan_offset + peak) + *(beam_norm_d + peak) * gain*peak_value;
		}
		IterDone:
	
		// Save the channel exit status (-1 or 1)
		*(exit_status_d+c) = channel_exit;
	}
}

return(clean);

}
