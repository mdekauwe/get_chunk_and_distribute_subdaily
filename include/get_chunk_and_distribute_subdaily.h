#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>


#define STRING_LENGTH 2000
#define TRUE 1
#define FALSE 0

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))



typedef struct  {
    int    row_start;
    int    row_end;
    int    col_start;
    int    col_end;
    int    nrows_in_slice;
    int    ncols_in_slice;
    int    ncols;
    int    nrows;
    int    num_land_pixels;
    char   land_mask_fn[STRING_LENGTH];
    int    root_processor;
    int    rank;
    int    size;
    int    nsize;
    int    remainder;
    float  land_id;
    float  cellsize;
    float  xllcorner;
    float  yllcorner;
    int    start_yr;
    int    end_yr;
    int    start_yr_forcing;
    int    end_yr_forcing;
    int    start_yr_rad;
    int    end_yr_rad;
    char   fdir[STRING_LENGTH];
} control;

typedef struct  {
    int   *tmax_dates;
    int   *tmin_dates;
    int   *rain_dates;
    int   *rad_dates;
    int   *vph09_dates;
    int   *vph15_dates;
    int   tmax_size;
    int   tmin_size;
    int   rain_size;
    int   rad_size;
    int   vph09_size;
    int   vph15_size;
    int   tmax_ndays;
    int   tmin_ndays;
    int   rain_ndays;
    int   rad_ndays;
    int   vph09_ndays;
    int   vph15_ndays;
    float *tmax_slice;
    float *tmin_slice;
    float *rain_slice;
    float *rad_slice;
    float *vph09_slice;
    float *vph15_slice;
 } met;


void   clparser(int, char **, control *);
void   initialise_stuff(control *);
void   mask_ij(control *, float *, int *);
void   read_met_data_slice(control *, met *, int *);
void   get_data(control *, char *, int, float **, int **, int *);
int    distribute_ij(control *, int *, int **);
int    distribute(control *, int *, float *, float **, int, int);
void   build_radiation_clim(control *, int *, float *, float **, float **);
void write_spinup_file(int, int, control *, met *, float *, float *,
                       float *, float *, float *, float *, float *);
void write_forcing_file(int, int, control *, met *, float *, float *,
                        float *, float *, float *, float *, float *, float *);



float  calc_day_length(int, int, float);
void   calc_tam_tpm(float *, float *, float, float, float, float);
float  calc_vpd(float, float);
int    is_leap_year(int);
