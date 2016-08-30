/*--------------------------------------------------------------------
 *	$Id: grdfilter.c,v 1.72 2009/08/11 12:03:00 guru Exp $
 *
 *	Copyright (c) 1991-2009 by P. Wessel and W. H. F. Smith
 *	See COPYING file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; version 2 of the License.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/*
 * grdfilter.c  reads a grid file and creates filtered grid file
 *
 * user selects boxcar, gaussian, cosine arch, custom, median, or mode filter
 * user selects distance scaling as appropriate for map work, etc.
 *
 * Author:  W.H.F. Smith
 * Date: 	16 Aug 88
 *
 * Modified:	27 Sep 88 by W.H.F. Smith, to use new median routine.
 *
 * Updated:	PW: 13-Jun-1991 to v2.0
 *		PW: 13-Jul-1992 to actually do v2 i/o.
 *		PW: 15-Jul-1992 to handle arbitrary new -R -I
 *		PW: 03-Jan-1995 to offer mode-filter (from filter1d)
 *		WS: 21-May-1998 to make warnings silent unless -V on.
 *		PW: 03-Jun-1998 upgrade to GMT 3.1
 *		PW: 02-Jun-1999 upgrade to GMT 3.3 + added SINCOS option
 *		PW: 18-Oct-1999 Use sincos directly
 *		PW: 18-JUN-2000 3.3.5
 * Version:	4
 *		PW: 27-MAY-2004 Added extreme values filter options l, L, u, U
 *		PW: 16-NOV-2005 Added spherical filtering support & highpass option
 *		PW: 27-APR-2007 Added -D5 to handle Mercator (img) grids
 *              RR (Roelof Rietbroek): 31-Aug-2009 Added the use of a custom kernel from ascii file (1D/2D)
/*  		or gmt_grid file (2D). */
/* 		Fixed 3 bugs: 1) fixed F.nx and F.ny such that they are always odd */
/* 		2): Changed inner loop over filter y direction, such that the loop continues ( instead of breaking) when encountering a pole which has already been used. The previous version caused the loop to ignore all remaining (and possibly valid) filterpoints when encountering an input point which was exactly on the poles */
/* 		3): weighted pole values with pole_weight*weight[ij_wt] instead of with pole_weight only */ 

#define GMT_WITH_NO_PS
#include "gmt.h"

struct GRDFILTER_CTRL {
	struct D {	/* -D<distflag> */
		BOOLEAN active;
		GMT_LONG mode;
	} D;
	struct F {	/* <type>[-]<filter_width>[<mode>] */
	  BOOLEAN active;
	  BOOLEAN highpass;
	  char filter;	/* Character codes for the filter */
 	  double width,width_x,width_y;
	  GMT_LONG mode;
	} F;
	struct G {	/* -G<file> */
		BOOLEAN active;
		char *file;
	} G;
	struct I {	/* -Idx[/dy] */
		BOOLEAN active;
		double xinc, yinc;
	} I;
	struct N {	/* -Np|i|r */
		BOOLEAN active;
		GMT_LONG mode;	/* 0 is default (i), 1 is replace (r), 2 is preserve (p) */
	} N;
	struct T {	/* -T */
		BOOLEAN active;
	} T;
};

#define IMG2LAT(img) (2.0*atand(exp((img)*D2R))-90.0)
#define LAT2IMG(lat) (R2D*log(tand(0.5*((lat)+90.0))))

#define GRDFILTER_WIDTH		0
#define GRDFILTER_HALF_WIDTH	1
#define GRDFILTER_X_SCALE	2
#define GRDFILTER_Y_SCALE	3
#define GRDFILTER_INV_R_SCALE	4

#define GRDFILTER_N_FILTERS	10 /* incremented RR*/

#define NAN_IGNORE	0
#define NAN_REPLACE	1
#define NAN_PRESERVE	2

struct FILTER_INFO {
  GMT_LONG nx;			/* The max number of filter weights in x-direction */
  GMT_LONG ny;			/* The max number of filter weights in y-direction */
  GMT_LONG x_half_width;	/* Number of filter nodes to either side needed at this latitude */
  GMT_LONG y_half_width;	/* Number of filter nodes above/below this point (ny_f/2) */
  GMT_LONG d_flag;
  double dx, dy;		/* Grid spacing in original units */
  double *x, *y;		/* Distances in original units along x and y to distance nodes */
  /*	double par[5];	removed by RR: unused  */
  double x_off, y_off;	/* Offsets relative to original grid */
  PFD weight_func, radius_func, azimuth_func; /*added azimuth function for 2D kernels*/
  double *Kern_x,*Kern_y; /* added RR* array holding custom kernel data*/
  float *Kernel;
  struct GRD_HEADER Kern_h; /* Add a grid header for the 2D filter (needed for using gmt_bcr interpolation routines */
  struct GMT_BCR Kern_bcr; /* Add a interpolation structure for 2D interpolation */
  struct GMT_EDGEINFO Kern_edge; /* Edge info structure */
  char * Kernfile;/*  Filename containing the external filter */
  BOOLEAN ISO; /*Added RR* Boolean whether kernel is 1 or 2D*/
  BOOLEAN gmtgrd;/*  Parameter determining whether the kernel file is a gmt grid */
};

/* Moved by RR from main to global such that it is accessible by all methods:*/
struct	FILTER_INFO F;


int main (int argc, char **argv)
{
	GMT_LONG	nx_out, ny_out, n_in_median, n_nan = 0, i_west, i_east, j_s_pole, j_n_pole;
	GMT_LONG	j_origin, i_out, j_out, half_nx, i_orig = 0, n_bad = 0, nm;
	GMT_LONG	i_in, j_in, ii, jj, i, j, ij_in, ij_out, ij_wt, effort_level;
	GMT_LONG	filter_type, one_or_zero = 1, GMT_n_multiples = 0;
	
	GMT_LONG *i_origin;


	BOOLEAN	error, new_range, fast_way, slow = FALSE, same_grid = FALSE;
	BOOLEAN wrap_case_x = FALSE, wrap_case_y = FALSE, pole_check = FALSE;
	BOOLEAN full_360, full_180, duplicate_check, i_west_used, i_east_used;

	float	*input, *output;

	double	west_new, east_new, south_new, north_new, merc_range, lat_out;
	double	x_scale = 1.0, y_scale = 1.0, x_width, y_width, y, pole_weight;
	double  par[5];/*[0] is filter width, [1] is 0.5*filter_width, [2] is xscale, [3] is yscale, [4] is 1/r_half for filter */
	double	x_out, y_out, wt_sum, value, last_median, this_median, xincnew2, yincnew2;
	double	xincold2, yincold2, y_shift = 0.0, x_fix = 0.0, y_fix = 0.0, max_lat;
	double	*weight, *work_array = VNULL, *x_shift = VNULL;

	char	*fin = CNULL, c, filter_code[GRDFILTER_N_FILTERS] = {'b', 'c', 'g','e', 'm', 'p', 'l', 'L', 'u', 'U'};
	char	*filter_name[GRDFILTER_N_FILTERS] = {"Boxcar", "Cosine Arch", "Gaussian","Custom", "Median", "Mode", "Lower", "Lower+", "Upper", "Upper-"};

	struct	GRD_HEADER h, test_h;
	struct	GRDFILTER_CTRL *Ctrl;

	void   set_weight_matrix (double *weight, double output_lat, double par[], double x_off, double y_off);
	double UnitWeight (double r, double par[]);
	double CosBellWeight (double r, double par[]);
	double GaussianWeight (double r, double par[]);
	double CartRadius (double x0, double y0, double x1, double y1, double par[]);
	double CartScaledRadius (double x0, double y0, double x1, double y1, double par[]);
	double FlatEarthRadius (double x0, double y0, double x1, double y1, double par[]);
	double SphericalRadius (double x0, double y0, double x1, double y1, double par[]);
	/* new functions by RR*/
	double CustomWeight_1D(double r, double par[]);
	double CustomWeight_2D(double r, double az, double par[]);
	double GMT_az_backaz_sphere(double lonE, double latE, double lonS, double latS, BOOLEAN baz);
	double GMT_az_backaz_flatearth(double lonE, double latE, double lonS, double latS, BOOLEAN baz);
	double GMT_az_backaz_cartesian(double lonE, double latE, double lonS, double latS, BOOLEAN baz);
	void Load_Customfilter(char *file, double * Fwidth_x,double * Fwidth_y);

	char *loc_colon;
	char *loc_slash;
	void *New_Grdfilter_Ctrl (), Free_Grdfilter_Ctrl (struct GRDFILTER_CTRL *C);
	int nfilt,w_nan,npole;
	BOOLEAN SP_used,NP_used;
	GMT_LONG padding[4]={0,0,0,0};
	argc = GMT_begin (argc, argv);

	Ctrl = (struct GRDFILTER_CTRL *)New_Grdfilter_Ctrl ();	/* Allocate and initialize a new control structure */

/* 	Some defaults */
	error = new_range = FALSE;
	fin = NULL;
	west_new = east_new = 0.0;
	filter_type = -1;
	F.ISO=TRUE; /* Default uses isotropic filter 22-11-2009 fixed bug (was set to false)*/
	F.gmtgrd=FALSE; /* Assume ascii input first */
	F.Kern_h.y_inc=F.Kern_h.y_inc=0;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'R':
				case 'V':
				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west_new, &east_new, &south_new, &north_new);
					break;

				case 'D':
					Ctrl->D.active = TRUE;
					Ctrl->D.mode = atoi(&argv[i][2]);
					break;
				case 'F':
					if (strchr ("bcgeEmpLlUu", argv[i][2])) {	/* OK filter code */
						Ctrl->F.active = TRUE;
						Ctrl->F.filter = argv[i][2];
						if (Ctrl->F.filter == 'p') {	/* Check for some futher info in case of mode filtering */
							c = argv[i][strlen(argv[i]-1)];
							if (c == '-') Ctrl->F.mode = -1;
							if (c == '+') Ctrl->F.mode = +1;
						}
						if(Ctrl->F.filter =='E'){
						  F.gmtgrd=TRUE;/*  Input filter is a gmt grid */
						  Ctrl->F.filter ='e';/*  Reset filter to plain external */
						}
/* 						Check for a colon: (and possibly a /) */
						loc_colon=strchr(argv[i],':');
						loc_slash=strchr(argv[i],'/');
/* 						Replace locations with a terminating null character */
						if(loc_slash != NULL){
						  *loc_slash='\0';
						  F.Kern_h.y_inc=atof(loc_slash+1);
						}
						if(loc_colon != NULL){
						  *loc_colon='\0';
						  F.Kern_h.x_inc=atof(loc_colon+1);
						}


						if (Ctrl->F.filter =='e'){ /*ADDED RR: Get filename from external kernel*/
						  if(argv[i][3]=='-')
						    {
						      Ctrl->F.highpass = TRUE;
						      F.Kernfile = strdup(&argv[i][4]);
						      
						    }
						  else
						    {
						      F.Kernfile = strdup(&argv[i][3]);
						    }

						  /* Read in filter file and determine filterwidth*/
						  Load_Customfilter(F.Kernfile,&Ctrl->F.width_x,&Ctrl->F.width_y);	
						  /*Get minimum of the filter widths (this will only be used in the filterwidth check)*/
						  Ctrl->F.width= (Ctrl->F.width_x < Ctrl->F.width_y) ? Ctrl->F.width_x :Ctrl->F.width_y;

						  	

						}else{ /*  No filtername but filterwidth */
						  Ctrl->F.width = atof (&argv[i][3]);
						  if (Ctrl->F.width < 0.0) Ctrl->F.highpass = TRUE;
						  Ctrl->F.width = fabs (Ctrl->F.width);
/* 						  Set Filter length in X and Y equal */
						  Ctrl->F.width_x=Ctrl->F.width_y=Ctrl->F.width;

						}
					}
					else {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -F option.  Correct syntax: -FX[<width>|<filterfile>], X one of bcgemplLuU\n", GMT_program);
						error++;
					}
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					break;
				case 'N':
					if (!argv[i][2]) {	/* Pixel registration OBSOLETE but BACKWARD COMPATIBLE for now */
						one_or_zero = 0;
					}
					else {	/* Treatment of NaNs */
						Ctrl->N.active = TRUE;
						switch (argv[i][2]) {
							case 'i':
								Ctrl->N.mode = NAN_IGNORE;	/* Default */
								break;
							case 'r':
								Ctrl->N.mode = NAN_REPLACE;	/* Replace */
								break;
							case 'p':
								Ctrl->N.mode = NAN_PRESERVE;	/* Preserve */
								break;
							default:
								fprintf (stderr, "%s: GMT SYNTAX ERROR -N option.  Correct syntax: -Ni|p|r\n", GMT_program);
								break;
						}
					}
					break;
				case 'T':	/* Toggle registration */
					Ctrl->T.active = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			fin = argv[i];
	}


	if (argc == 1 || GMT_give_synopsis_and_exit) {
	  fprintf (stderr, "%s %s - Filter a 2-D grid file in the space (or time) domain\n\n",GMT_program, GMT_VERSION);
	  fprintf(stderr,"usage: %s input_file -D<distance_flag> -F<type>[-]<filter_width>|<filter_file>[:<xinc>[/<yinc>]][<mode>]\n",GMT_program);
		fprintf(stderr,"\t-G<output_file> [%s] [-Ni|p|r] [%s] [-T] [-V] [%s]\n", GMT_I_OPT, GMT_Rgeo_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf(stderr,"\tDistance flag determines how grid (x,y) maps into distance units of filter width as follows:\n");
		fprintf(stderr,"\t   -D0 grid x,y same units as <filter_width>, cartesian Distances.\n");
		fprintf(stderr,"\t   -D1 grid x,y in degrees, <filter_width> in km, cartesian Distances.\n");
		fprintf(stderr,"\t   -D2 grid x,y in degrees, <filter_width> in km, x_scaled by cos(middle y), cartesian Distances.\n");
		fprintf(stderr,"\t   These first three options are faster; they allow weight matrix to be computed only once.\n");
		fprintf(stderr,"\t   Next three options are slower; weights must be recomputed for each scan line.\n");
		fprintf(stderr,"\t   -D3 grid x,y in degrees, <filter_width> in km, x_scale varies as cos(y), cartesian Distances.\n");
		fprintf(stderr,"\t   -D4 grid x,y in degrees, <filter_width> in km, spherical Distances.\n");
		fprintf(stderr,"\t   -D5 grid x,y in Mercator units (-Jm1), <filter_width> in km, spherical Distances.\n");
		fprintf(stderr,"\t-F sets the low-pass filter type and full diameter (6 sigma) filter_width.  Choose between\n");
		fprintf(stderr,"\t   convolution-type filters which differ in how weights are assigned and geospatial\n");
		fprintf(stderr,"\t   filters that seek to return a representative value.\n");
		fprintf(stderr,"\t   Give negative filter width to select highpass filtering [lowpass].\n");
		fprintf(stderr,"\t   Convolution filters:\n");
		fprintf(stderr,"\t     b: Boxcar : a simple averaging of all points inside filter radius.\n");
		fprintf(stderr,"\t     c: Cosine arch : a weighted averaging with cosine arc weights\n");
		fprintf(stderr,"\t     g: Gaussian : weighted averaging with Gaussian weights.\n");
		fprintf(stderr,"\t     e|E: External filter from file: Append the name of the file containing the equidistant filter\n");
		fprintf(stderr,"\t        For 'e' the file contains either 2 columns (isotropic kernel with radius and weight)\n");
		fprintf(stderr,"\t        Or three columns with xyz triplets ( 2D filter). Unless <xinc> and optionally <yinc> are\n");
		fprintf(stderr,"\t        given the files must contain a square filter and may not contain data gaps.\n");
		fprintf(stderr,"\t        When 'E' is specified the filename is assumed to be a GMT_grid file.\n");
		fprintf(stderr,"\t        The filter_width is derived from the data and should NOT be given for the e|E option.\n");
		fprintf(stderr,"\t        The filtervalues are derived by cubic (bi) spline interpolation from the original (1D/2D)filter.\n");
		fprintf(stderr,"\t        So they may be defined on a different grid. The filters must be centered on the origin.\n");
		fprintf(stderr,"\t   Geospatial filters:\n");
		fprintf(stderr,"\t     l: Lower : return minimum of all points.\n");
		fprintf(stderr,"\t     L: Lower+ : return minimum of all +ve points.\n");
		fprintf(stderr,"\t     m: Median : return the median value of all points.\n");
		fprintf(stderr,"\t     p: Maximum likelihood probability estimator : return mode of all points.\n");
		fprintf(stderr,"\t        By default, we return the average if more than one mode is found.\n");
		fprintf(stderr,"\t        Append - or + to the width to instead return the smallest or largest mode.\n");
		fprintf(stderr,"\t     u: Upper : return maximum of all points.\n");
		fprintf(stderr,"\t     U: Upper- : return maximum of all -ve points.\n");
		fprintf(stderr,"\t-G sets output filename for filtered grid\n");
		fprintf(stderr, "\n\tOPTIONS:\n");
		GMT_inc_syntax ('I', 0);
		fprintf(stderr,"\t   The new xinc and yinc should be divisible by the old ones (new lattice is subset of old).\n");
		fprintf(stderr, "\t-N specifies how NaNs in the input grid should be treated.  There are three options:\n");
		fprintf(stderr, "\t   -Ni skips all NaN values and returns a filtered value unless all are NaN [Default]\n");
		fprintf(stderr, "\t   -Np sets filtered output to NaN is any NaNs are found inside filter circle.\n");
		fprintf(stderr, "\t   -Nr sets filtered output to NaN if the corresponding input node was NaN.\n");
		fprintf(stderr, "\t      (only possible if the input and output grids are coregistered).\n");
		fprintf(stderr, "\t-T Toggles between grid and pixel registration for output grid [Default is same as input registration]\n");
		fprintf(stderr, "\t-R for new Range of output grid; enter <WESN> (xmin, xmax, ymin, ymax) separated by slashes.\n");
		GMT_explain_option ('V');
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, NULL, &Ctrl->I.active);

	if (!Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option:  Must specify output file\n", GMT_program);
		error++;
	}
	if (!fin) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
		error++;
	}
	if (Ctrl->D.mode < 0 || Ctrl->D.mode > 5) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -D option:  Choose from the range 0-5\n", GMT_program);
		error++;
	}
	if (!Ctrl->F.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: -F option is required:\n", GMT_program);
		error++;
	}
	if (Ctrl->F.active && Ctrl->F.width == 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -F option:  filter fullwidth must be nonzero:\n", GMT_program);
		error++;
	}
	if (Ctrl->I.active && (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error++;
	}
	if (Ctrl->T.active && one_or_zero == 0) {	/* Both -N and -T set, not good */
		fprintf (stderr, "%s: GMT SYNTAX ERROR -T option:  Not allowed with obsolete -N option\n", GMT_program);
		error++;
	}
	if (project_info.region_supplied && Ctrl->I.active && Ctrl->F.highpass) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -F option:  Highpass filtering requires original -R -I\n", GMT_program);
		error++;
	}
	
	if (error) exit (EXIT_FAILURE);

	/* Assign filter_type number */
	
	for (filter_type = 0; filter_type < GRDFILTER_N_FILTERS && filter_code[filter_type] != Ctrl->F.filter; filter_type++);
	

	if (project_info.region_supplied) new_range = TRUE;

	GMT_err_fail (GMT_read_grd_info (fin, &h), fin);
	GMT_grd_init (&h, argc, argv, TRUE);	/* Update command history only */

	if (Ctrl->T.active) {	/* Make output grid of the opposite registration */
		one_or_zero = (h.node_offset) ? 1 : 0;
	}
	else
		one_or_zero = (h.node_offset) ? 0 : 1;

	/* Read the input grid file and close it  */

	nm = GMT_get_nm (h.nx, h.ny);
	input = (float *) GMT_memory (VNULL, (size_t)nm, sizeof(float), GMT_program);

	GMT_err_fail (GMT_read_grd (fin, &h, input, 0.0, 0.0, 0.0, 0.0, padding, FALSE), fin);

	full_360 = (Ctrl->D.mode && GMT_360_RANGE (h.x_max, h.x_min));	/* Periodic geographic grid */
	full_180 = (Ctrl->D.mode && GMT_180_RANGE (h.y_min, h.y_max));	/* Full latitude range for geographic grid */

	last_median = 0.5 * (h.z_min + h.z_max);

	/* Check range of output area and set i,j offsets, etc.  */

	if (!new_range) {
		west_new = h.x_min;
		east_new = h.x_max;
		south_new = h.y_min;
		north_new = h.y_max;
	}
	if (!Ctrl->I.active) {
		Ctrl->I.xinc = h.x_inc;
		Ctrl->I.yinc = h.y_inc;
	}

	if (!full_360) {
		if (west_new < h.x_min) error = TRUE;
		if (east_new > h.x_max) error = TRUE;
	}
	if (south_new < h.y_min) error = TRUE;
	if (north_new > h.y_max) error = TRUE;
	if (Ctrl->I.xinc <= 0.0) error = TRUE;
	if (Ctrl->I.yinc <= 0.0) error = TRUE;

	if (error) {
		fprintf(stderr,"%s: New WESN incompatible with old.\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	/* Make sure output grid is kosher */

	test_h.x_min = west_new;	test_h.x_max = east_new;	test_h.x_inc = Ctrl->I.xinc;
	test_h.y_min = south_new;	test_h.y_max = north_new;	test_h.y_inc = Ctrl->I.yinc;
	GMT_RI_prepare (&test_h);	/* Ensure -R -I consistency and set nx, ny */
	GMT_err_fail (GMT_grd_RI_verify (&test_h, 1), Ctrl->G.file);
	/* Copy back in case test_h were changed */
	west_new  = test_h.x_min;	east_new  = test_h.x_max;	Ctrl->I.xinc = test_h.x_inc;
	south_new = test_h.y_min;	north_new = test_h.y_max;	Ctrl->I.yinc = test_h.y_inc;

	/* We can save time by computing a weight matrix once [or once pr scanline] only
	   if new grid spacing is a multiple of old spacing */

	fast_way = (fabs (fmod (Ctrl->I.xinc / h.x_inc, 1.0)) < GMT_SMALL && fabs (fmod (Ctrl->I.yinc / h.y_inc, 1.0)) < GMT_SMALL);
	same_grid = !(new_range  || Ctrl->I.active || h.node_offset == one_or_zero);
	if (!fast_way && gmtdefs.verbose) {
		fprintf (stderr, "%s: Warning - Your output grid spacing is such that filter-weights must\n", GMT_program);
		fprintf (stderr, "be recomputed for every output node, so expect this run to be slow.  Calculations\n");
		fprintf (stderr, "can be speeded up significantly if output grid spacing is chosen to be a multiple\n");
		fprintf (stderr, "of the input grid spacing.  If the odd output grid is necessary, consider using\n");
		fprintf (stderr, "a \'fast\' grid for filtering and then resample onto your desired grid with grdsample.\n");
	}
	if (Ctrl->N.mode == NAN_REPLACE && !same_grid) {
		fprintf (stderr, "%s: Warning: -Nr requires co-registered input/output grids, option is ignored\n", GMT_program);
		Ctrl->N.mode = NAN_IGNORE;
	}
	nx_out = one_or_zero + irint ( (east_new - west_new) / Ctrl->I.xinc);
	ny_out = one_or_zero + irint ( (north_new - south_new) / Ctrl->I.yinc);

	nm = GMT_get_nm (nx_out, ny_out);
	output = (float *) GMT_memory (VNULL, (size_t)nm, sizeof(float), GMT_program);
	i_origin = (GMT_LONG *) GMT_memory (VNULL, (size_t)nx_out, sizeof(GMT_LONG), GMT_program);
	if (!fast_way) x_shift = (double *) GMT_memory (VNULL, (size_t)nx_out, sizeof(double), GMT_program);

	xincnew2 = (one_or_zero) ? 0.0 : 0.5 * Ctrl->I.xinc;
	yincnew2 = (one_or_zero) ? 0.0 : 0.5 * Ctrl->I.yinc;
	xincold2 = (h.node_offset) ? 0.5 * h.x_inc : 0.0;
	yincold2 = (h.node_offset) ? 0.5 * h.y_inc : 0.0;

/* 	Moved and changed by RR downward */
/* 	if (fast_way && h.node_offset == one_or_zero) {	/\* multiple grid but one is pix, other is grid *\/ */
/* 		x_fix = 0.5 * h.x_inc; */
/* 		y_fix = 0.5 * h.y_inc; */
/* 	} */

	/* Set up the distance scalings for lon and lat, and assign pointer to distance function  */

	switch (Ctrl->D.mode) {
		case 0:	/* Plain, unscaled isotropic Cartesian distances */
			x_scale = y_scale = 1.0;
			F.radius_func = CartRadius;
			F.azimuth_func = GMT_az_backaz_cartesian;
			break;
		case 1:	/* Plain, scaled (degree to km) isotropic Cartesian distances */
			x_scale = y_scale = project_info.KM_PR_DEG;
			F.radius_func = CartScaledRadius;
			F.azimuth_func = GMT_az_backaz_cartesian;
			break;
		case 2:	/* Flat Earth Cartesian distances, xscale fixed at mid latitude */
			x_scale = project_info.KM_PR_DEG * cosd (0.5 * (north_new + south_new));
			y_scale = project_info.KM_PR_DEG;
			F.radius_func = FlatEarthRadius;
			F.azimuth_func = GMT_az_backaz_flatearth;
			break;
		case 3:	/* Flat Earth Cartesian distances, xscale reset for each latitude */
			x_scale = project_info.KM_PR_DEG * ((fabs (south_new) > north_new) ? cosd (south_new) : cosd (north_new));
			y_scale = project_info.KM_PR_DEG;
			F.radius_func = FlatEarthRadius;
			F.azimuth_func = GMT_az_backaz_flatearth;
			break;
		case 4:	/* Great circle distances */
			x_scale = 0.0;
			y_scale = project_info.KM_PR_DEG;
			F.radius_func = SphericalRadius;
			F.azimuth_func = GMT_az_backaz_sphere;
			if (full_360) wrap_case_x = TRUE;		/* For spherical filtering */
			if (full_180) wrap_case_y = wrap_case_x;	/* For spherical filtering */
			break;
		case 5:	/* Great circle distances with Mercator coordinates */
			/* Get the max |lat| extent of the grid */
			max_lat = IMG2LAT (MAX (fabs (h.y_min), fabs (h.y_max)));
			merc_range = LAT2IMG (max_lat + (0.5 * Ctrl->F.width_y / project_info.KM_PR_DEG)) - LAT2IMG (max_lat);
			x_scale = y_scale = 0.5 * Ctrl->F.width_y / merc_range;
			F.radius_func = SphericalRadius;
			F.azimuth_func = GMT_az_backaz_sphere;
			if (full_360) wrap_case_x = TRUE;		/* For spherical filtering */
			break;
	}

	switch (filter_type) {
		case 1:	/*  Cosine-bell filter weights */
			par[GRDFILTER_INV_R_SCALE] = 2.0 / Ctrl->F.width;
			F.weight_func = CosBellWeight;
			break;
		case 2:	/*  Gaussian filter weights */
			par[GRDFILTER_INV_R_SCALE] = -18.0 / (Ctrl->F.width * Ctrl->F.width);
			F.weight_func = GaussianWeight;
			break;
		case 3:/*Custom kernel from file*/
		  if(F.ISO){
		    F.weight_func = CustomWeight_1D;
		  }else{
		    F.weight_func = CustomWeight_2D;
		  }
		  break;
		default:	/* Everything else uses unit weights */
		  F.weight_func = UnitWeight;
		  break;
	}

	/* Set up miscellaneous filter parameters needed when computing the weights */
	
	par[GRDFILTER_WIDTH] = Ctrl->F.width;
	par[GRDFILTER_HALF_WIDTH] = 0.5 * Ctrl->F.width;
	par[GRDFILTER_X_SCALE] = x_scale;
	par[GRDFILTER_Y_SCALE] = (Ctrl->D.mode == 5) ? project_info.KM_PR_DEG : y_scale;
	F.d_flag = Ctrl->D.mode;
	F.dx = h.x_inc;
	F.dy = h.y_inc;
	x_width = Ctrl->F.width_x / (h.x_inc * x_scale);
	y_width = Ctrl->F.width_y / (h.y_inc * y_scale);
	F.y_half_width = (GMT_LONG) (ceil(y_width) / 2.0);
	F.x_half_width = (GMT_LONG) (ceil(x_width) / 2.0);
	F.nx = 2 * F.x_half_width + 1;
	F.ny = 2 * F.y_half_width + 1;
	if (x_scale == 0.0 || F.nx < 0 || F.nx > h.nx) {	/* Safety valve when x_scale -> 0.0 */
	  F.nx = 2*(h.nx/2) + 1; /* RR: removed bug: ensure that F.nx and F.ny are ALWAYS odd */
	  F.x_half_width = F.nx / 2;
	}
	if (F.ny < 0 || F.ny > h.ny) {
	  F.ny = 2*(h.ny/2)+1;
		F.y_half_width = F.ny / 2;
	}
	F.x = (double *) GMT_memory (VNULL, (size_t)(F.x_half_width+1), sizeof (double), GMT_program);
	F.y = (double *) GMT_memory (VNULL, (size_t)(F.y_half_width+1), sizeof (double), GMT_program);
	for (i = 0; i <= F.x_half_width; i++) F.x[i] = i * F.dx;
	for (j = 0; j <= F.y_half_width; j++) F.y[j] = j * F.dy;
	
	weight = (double *) GMT_memory (VNULL, (size_t)(F.nx*F.ny), sizeof(double), GMT_program);

	if (filter_type >= 4) {	/* These filters are not convolutions; they require sorting or comparisons */
		slow = TRUE;
		work_array = (double *) GMT_memory (VNULL, (size_t)(F.nx*F.ny), sizeof(double), GMT_program);
	}

	if (gmtdefs.verbose) {
		fprintf(stderr,"%s: Input nx,ny = (%d %d), output nx,ny = (%ld %ld), filter nx,ny = (%ld %ld)\n", GMT_program, h.nx, h.ny, nx_out, ny_out, F.nx, F.ny);
		fprintf(stderr,"%s: Filter type is %s.\n", GMT_program, filter_name[filter_type]);
	}

	/* Compute nearest xoutput i-indices and shifts once */

	for (i_out = 0; i_out < nx_out; i_out++) {
		x_out = west_new + i_out * Ctrl->I.xinc + xincnew2;
		i_origin[i_out] = GMT_x_to_i (x_out, h.x_min, h.x_inc, h.xy_off, h.nx);
		if (!fast_way) x_shift[i_out] = x_out - (h.x_min + i_origin[i_out] * h.x_inc + xincold2);
	}
	
/* 	Calculate constant gridpoint offset between input and output grid (RR: statement moved and removed bug) */
	if( fast_way ) {
	  /* Compare first point from the left of the output grid  */
	  x_out = west_new + xincnew2; 
	  x_fix=x_out - (h.x_min + i_origin[0] * h.x_inc + xincold2);

/* 	  Same for the most northern point of the output grid */
	  y_out = north_new  - yincnew2;
	  j_origin = GMT_y_to_j (y_out, h.y_min, h.y_inc, h.xy_off, h.ny);
	  y_fix = y_out - (h.y_max - j_origin * h.y_inc - yincold2);
	}
	


	/* Determine how much effort to compute weights:
		1 = Compute weights once for entire grid
		2 = Compute weights once per scanline
		3 = Compute weights for every output point [slow]
	*/

	if (fast_way && Ctrl->D.mode <= 2)
		effort_level = 1;
	else if (fast_way && Ctrl->D.mode > 2)
		effort_level = 2;
	else 
		effort_level = 3;

	if (effort_level == 1) set_weight_matrix ( weight, 0.0, par, x_fix, y_fix);
	half_nx = (h.node_offset) ? h.nx / 2 : (h.nx - 1) / 2;
	pole_weight = 0.25 * M_PI * (h.y_inc / h.x_inc);	/* This is the weight of the pole point for gridregistered grids */
	pole_check = (Ctrl->D.mode && !h.node_offset);		/* Must check to make sure we only use the N and S pole value once */
	duplicate_check = (full_360 && !h.node_offset);		/* Must avoid using the duplicated Greenwich node twice */
	j_s_pole = (h.y_min == -90.0) ? h.ny - 1 : INT_MIN;	/* row number of S pole (if in the grid)*/
	j_n_pole = (h.y_max == 90.0) ? 0 : INT_MIN;		/* row number of N pole (if in the grid)*/
	i_west = 0;		/* column number of western edge */
	i_east = h.nx - 1;	/* column number of eastern edge */

	for (j_out = 0; j_out < ny_out; j_out++) {

		if (gmtdefs.verbose) fprintf (stderr, "%s: Processing output line %ld\r", GMT_program, j_out);

		y_out = north_new - j_out * Ctrl->I.yinc - yincnew2;
		lat_out = (Ctrl->D.mode == 5) ? IMG2LAT (y_out) : y_out;
		j_origin = GMT_y_to_j (y_out, h.y_min, h.y_inc, h.xy_off, h.ny);
		if (Ctrl->D.mode == 3) par[GRDFILTER_X_SCALE] = project_info.KM_PR_DEG * cosd (lat_out);	/* Update flat-earth longitude scale */

		if (Ctrl->D.mode > 2) {	/* Update max filterweight nodes to deal with at this latitude */
			y = fabs (lat_out);

			if (Ctrl->D.mode == 4) y += (Ctrl->F.width_y /(2.0 * par[GRDFILTER_Y_SCALE]));
			F.x_half_width = (y < 90.0) ? MIN (F.nx / 2, irint (Ctrl->F.width_x / (2.0 * F.dx * par[GRDFILTER_Y_SCALE] * cosd (y)))) : F.nx / 2;
		}

		if (effort_level == 2) set_weight_matrix ( weight, y_out, par, x_fix, y_fix);
		
		if (!fast_way) y_shift = y_out - (h.y_max - j_origin * h.y_inc - yincold2);

		for (i_out = 0; i_out < nx_out; i_out++) {
		  
			if (effort_level == 3) set_weight_matrix ( weight, y_out, par, x_shift[i_out], y_shift);
			wt_sum = value = 0.0;
			n_in_median = n_bad = 0;
			ij_out = j_out * nx_out + i_out;
			if (Ctrl->N.mode == NAN_REPLACE && GMT_is_fnan (input[ij_out])) continue;	/* Since output will be NaN we bypass the filter loop */

			/* Now loop over the filter domain and collect those points that should be considered by the filter operation */

/*RR: 			Tag to track whether ( North and South pole has been used for a given output point) */
			SP_used=FALSE;
			NP_used=FALSE;
			nfilt=0;
			w_nan=0;
			npole=0;
			i_west_used = i_east_used = FALSE;
			for (ii = -F.x_half_width; ii <= F.x_half_width; ii++) {
				i_in = i_origin[i_out] + ii;

				if (wrap_case_x) {	/* Just wrap around the globe */
					if (i_in < 0) i_in += h.nx;
					if (i_in >= h.nx) i_in -= h.nx;
					i_orig = i_in;
				}
				else if ( (i_in < 0) || (i_in >= h.nx)){continue;}

				if (i_in == i_west) i_west_used = TRUE;
				if (i_in == i_east) i_east_used = TRUE;
				if (duplicate_check && ((i_in == i_east && i_west_used) || (i_in == i_west && i_east_used))) continue;	/* Do not use the node at 360 if the one at 0 was used, or vice versa */
				
				for (jj = -F.y_half_width; jj <= F.y_half_width; jj++) {
/* 				  j_in = j_origin + jj; removed bug should be -jj since grid starts in the north while weight starts in the south*/
				  j_in = j_origin - jj;
				  
				  if (pole_check && (j_in == j_n_pole || j_in == j_s_pole)) {
				    /*RR removed bug: jj loop should CONTINUE on encountering a pole, which is already used, NOT break the jj loop */
				    if(SP_used || NP_used)continue;/*  Continue loop when Point is on the poles and poles are already used*/

				    
				    ij_in = j_in*h.nx;	/* All points on pole rows are the same; pick the first one */
				    if (GMT_is_fnan (input[ij_in])) {	/* Whaddaya know; the pole is NaN */
				      n_bad++;
				      continue;
				    }
				    /* Get here when pole value is usable  */
				    if (slow) {
				      work_array[n_in_median] = input[ij_in];
				      n_in_median++;
				    }
				    else {
				      
				      ij_wt = (jj + F.y_half_width) * F.nx + ii + F.x_half_width;

				      if (GMT_is_dnan(weight[ij_wt])){w_nan++;continue;}/*  Adapted by RR: also allow negative weights but no NaN*/
				      value += input[ij_in] * pole_weight*weight[ij_wt];
				      wt_sum += weight[ij_wt] * pole_weight;
				      nfilt++;
				      npole+=F.nx-1;
/* 				      value += input[ij_in] * pole_weight; */
/* 				      wt_sum += pole_weight; */
				    }
				    /* Tag N or S pole as used */
				    if(j_in == j_n_pole){
				      NP_used=TRUE;
				    }else{
				      SP_used=TRUE;
				    }

				    continue;	/* Done with this filter latitude */
				  }
				  
				  if (wrap_case_y) {	/* Just wrap over the pole */
				    i_in = i_orig;
				    if (j_in < 0) j_in = -j_in, i_in += half_nx;
				    if (j_in >= h.ny) j_in = 2*h.ny - j_in - 2, i_in += half_nx;
				    if (i_in < 0) i_in += h.nx;
				    if (i_in >= h.nx) i_in -= h.nx;
				  }
				  else if ( (j_in < 0) || (j_in >= h.ny) ) continue;

				  ij_wt = (jj + F.y_half_width) * F.nx + ii + F.x_half_width;
				  if (GMT_is_fnan(weight[ij_wt])){w_nan++;continue;}/*  Adapted by RR: also allow negative weights but no NaN*/
				  
				  ij_in = j_in*h.nx + i_in;
				  if (GMT_is_fnan (input[ij_in])) {
				    n_bad++;
				    continue;
				  }

				  
				  /* Get here when point is usable  */
				  if (slow) {
				    work_array[n_in_median] = input[ij_in];
				    n_in_median++;
				  }
				  else {
				 
				    value += input[ij_in] * weight[ij_wt];
				    wt_sum += weight[ij_wt];
/* 				    fprintf(stderr,"%f %f %f %f \n",ii*0.5,jj*0.5,input[ij_in],weight[ij_wt]); */
				    nfilt++;
				  }
				}
			}
/* 			if (gmtdefs.verbose) fprintf (stderr, "value %f normweight %f bad points %ld filtered points %d\n", value, wt_sum,n_bad,nfilt); */
			/* Now we have done the convolution and we can get the value  */
/* 			fprintf(stderr,"valid, rejected and pole_duplicates %d %d %d total; %d badpoints %ld\n",nfilt,w_nan,npole,npole+nfilt+w_nan,n_bad); */
			if (Ctrl->N.mode == NAN_REPLACE && GMT_is_fnan (input[ij_out])) {
				output[ij_out] = GMT_f_NaN;
				n_nan++;
			}
			else if (Ctrl->N.mode == NAN_PRESERVE && n_bad) {	/* -Np in effect and there were NaNs inside circle */
				output[ij_out] = GMT_f_NaN;
				n_nan++;
			}
			else if (slow) {
				if (n_in_median) {
					switch (filter_type) {
						case 4:	/* Median */
							GMT_median (work_array, n_in_median, h.z_min, h.z_max, last_median, &this_median);
							last_median = this_median;
							break;
						case 5:	/* Mode */
							GMT_mode (work_array, n_in_median, n_in_median/2, TRUE, Ctrl->F.mode, &GMT_n_multiples, &this_median);
							break;
						case 6:	/* Lowest of all */
							this_median = GMT_extreme (work_array, n_in_median, DBL_MAX, 0, -1);
							break;
						case 7:	/* Lowest of positive values */
							this_median = GMT_extreme (work_array, n_in_median, 0.0, +1, -1);
							break;
						case 8:	/* Upper of all values */
							this_median = GMT_extreme (work_array, n_in_median, -DBL_MAX, 0, +1);
							break;
						case 9:	/* Upper of negative values */
							this_median = GMT_extreme (work_array, n_in_median, 0.0, -1, +1);
							break;
					}
					output[ij_out] = (float)this_median;
				}
				else {
					output[ij_out] = GMT_f_NaN;
					n_nan++;
				}
			}
			else {
				if (wt_sum == 0.0) {	/* Assign value = GMT_f_NaN */
					output[ij_out] = GMT_f_NaN;
					n_nan++;
				}
				else {
				  output[ij_out] = (float)(value / wt_sum);
				}
			}
		}
	}

		
		
	
	if (gmtdefs.verbose) fprintf (stderr, "\n");

	if (Ctrl->F.highpass) {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Subtracting lowpass-filtered data from grid to obtain high-pass filtered data\n", GMT_program);
		for (ij_out = 0; ij_out < nx_out * ny_out; ij_out++) output[ij_out] = input[ij_out] - output[ij_out];
	}
	
	/* At last, that's it!  Output: */

	if (n_nan && gmtdefs.verbose) fprintf (stderr, "%s: Unable to estimate value at %ld nodes, set to NaN\n", GMT_program, n_nan);
	if (GMT_n_multiples > 0 && gmtdefs.verbose) fprintf (stderr, "%s: WARNING: %ld multiple modes found\n", GMT_program, GMT_n_multiples);

	h.nx = nx_out;
	h.ny = ny_out;
	h.x_min = west_new;
	h.x_max = east_new;
	h.y_min = south_new;
	h.y_max = north_new;
	h.x_inc = Ctrl->I.xinc;
	h.y_inc = Ctrl->I.yinc;
	h.node_offset = !one_or_zero;
	GMT_err_fail (GMT_write_grd (Ctrl->G.file, &h, output, 0.0, 0.0, 0.0, 0.0, padding, FALSE), Ctrl->G.file);

	GMT_free ((void *) input);
	GMT_free ((void *) output);
	GMT_free ((void *) weight);
	GMT_free ((void *) i_origin);
	GMT_free ((void *) F.x);
	GMT_free ((void *) F.y);
	if (slow) GMT_free ((void *) work_array);
	if (!fast_way) GMT_free ((void *) x_shift);

	Free_Grdfilter_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void set_weight_matrix (double *weight, double output_lat, double par[], double x_off, double y_off)
{
	/* x_off and y_off give offset between output node and 'origin' input node for this window (0,0 for integral grids).
	 * fast is TRUE when input/output grids are offset by integer values in dx/dy.
	 * Here, par[0] = filter_width, par[1] = filter_width / 2, par[2] = x_scale, part[3] = y_scale, and
	 * par[4] is the normalization distance needed for the Cosine-bell or Gaussian weight function.
	 */

	GMT_LONG	i, j, ij;
	double	x, y, yc, y0, x_jump = 0.0, normalize = 1.0, r, dy = 1.0, dy_half = 0.5, dy_shrink = 1.0;
	double az;
	yc = y0 = output_lat - y_off;
	if (F.d_flag == 5) {
		yc = IMG2LAT (yc);	/* Recover actual latitudes */
		dy = F.y[1] - F.y[0];
		dy_half = 0.5 * dy;
	}
	for (j = -F.y_half_width; j <= F.y_half_width; j++) {
		y = y0 + ((j < 0) ? -F.y[-j] : F.y[j]);
		if (F.d_flag) {	/* Must normalize based on shrinking area representation and check for going over the pole */
			if (F.d_flag == 5) {
				dy_shrink = (IMG2LAT (y + dy_half) - IMG2LAT (y - dy_half)) / dy;		/* Recover actual latitudes */
				y = IMG2LAT (y);	/* Recover actual latitudes */
			}
			if (fabs (y) > 90.0) {	/* Must find the point across the pole */
				y = copysign (180.0 - fabs (y), y);
				x_jump = 180.0;
			}
			else
				x_jump = 0.0;
			normalize = dy_shrink * cosd (y);
		}
		for (i = -F.x_half_width; i <= F.x_half_width; i++) {
			x = (i < 0) ? -F.x[-i] : F.x[i];
			ij = (j + F.y_half_width) * F.nx + i + F.x_half_width;
			r = F.radius_func (x_off, yc, x + x_jump, y, par);
			if(F.ISO){ /* Don't need azimuth calculation for Isotropic filters*/
/* 			  fprintf(stderr,"ISO FILTER  \n");					 */
			  weight[ij] = (r > par[GRDFILTER_HALF_WIDTH]) ? GMT_f_NaN : F.weight_func (r, par);

			}else{/* Do need azimuth calculation for 2D filter */
/* 			  fprintf(stderr,"NON-ISO FILTER  \n");					 */
			  az=F.azimuth_func (x_off, yc, x + x_jump, y, TRUE);
/* 			  fprintf(stderr,"x0,y0,x1,y1,r,az %f %f %f %f %f %f\n",x_off, yc, x + x_jump, y,r,az); */
			  weight[ij] = F.weight_func (r,az,par);
			}
			if (F.d_flag) weight[ij] *= normalize;	/* Adjust for variation in area with latitude */
		}
	}

}


double UnitWeight (double r, double par[])
{
	/* Return unit weight since we know r is inside radius */
	return (1.0);
}

double CosBellWeight (double r, double par[])
{
	/* Return the cosine-bell filter weight for given r.
	 * The parameter r_f_half is the 5th parameter passed.
	 */

	return (1.0 + cos (M_PI * r * par[GRDFILTER_INV_R_SCALE]));
}

double GaussianWeight (double r, double par[])
{
	/* Return the Gaussian filter weight for given r.
	 * The parameter sig_2 is the 5th parameter passed.
	 */

	return (exp (r * r * par[GRDFILTER_INV_R_SCALE]));
}

void *New_Grdfilter_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDFILTER_CTRL *C;
	
	C = (struct GRDFILTER_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDFILTER_CTRL), "New_Grdfilter_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->D.mode = -1;	
	return ((void *)C);
}

void Free_Grdfilter_Ctrl (struct GRDFILTER_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}


  /* Return the interpolated value ( Cubic spline) from the loaded 1D Kernel contained in FILTER_INFO structure F*/
double CustomWeight_1D(double r, double par[]){
  GMT_LONG err;
  double out;

  if(GMT_intpol (F.Kern_x, F.Kern_y, F.Kern_h.nx, (GMT_LONG) 1, &r, &out, (GMT_LONG) 2) != GMT_NOERROR)
    {
      fprintf(stderr,"Error in linear interpolation of Custom Kernel, exiting\n");
      exit(1);
    }
  return (out);
}


  /* Return the interpolated value ( Cubic spline) from the loaded 2D Kernel contained in FILTER_INFO structure F*/
double CustomWeight_2D(double r, double az, double par[]){
  GMT_LONG err;
  double out,xx,yy;
  
  if(r < GMT_SMALL && GMT_is_dnan(az)){/*  Then it is the origin of the grid */
    xx=0;
    yy=0;
  } else{
  xx=r*sind(az);
  yy=r*cosd(az);
  }
/*   Interpolate from grid suplied in F.Kernel (returns NaN when outside grid)*/
  out= GMT_get_bcr_z (&F.Kern_h, xx, yy, F.Kernel, &F.Kern_edge, &F.Kern_bcr);
  return (out);
}

/* NOTE by RR: in the XXXXRadius routines below changed x0-x1 and y0-y1 to x1-x0 and y1-y0 for better readability ( positive distances)) */

double CartRadius (double x0, double y0, double x1, double y1, double par[])
{	/* Plain Cartesian distance */
	return (hypot (x1 - x0, y1 - y0));
}

double CartScaledRadius (double x0, double y0, double x1, double y1, double par[])
{	/* Plain scaled Cartesian distance (xscale = yscale) */
	return (par[GRDFILTER_X_SCALE] * hypot (x1 - x0, y1 - y0));
}

double FlatEarthRadius (double x0, double y0, double x1, double y1, double par[])
{	/* Cartesian radius with different scales */
	return (hypot (par[GRDFILTER_X_SCALE] * (x1 - x0), par[GRDFILTER_Y_SCALE] * (y1 - y0)));
}

double SphericalRadius (double x0, double y0, double x1, double y1, double par[])
{	/* Great circle distance with polar wrap-around test on 2nd point */
	if (fabs (y1) > 90.0) {	/* Must find the point across the pole */
		y1 = copysign (180.0 - fabs (y1), y1);
		x1 += 180.0;
	}
/* 	fprintf(stderr,"spherical inmput %f %f %f %f\n",x0, y0, x1, y1); */
	return (GMT_great_circle_dist_km (x0, y0, x1, y1));
}


/* Function to read in data from file and set up some filter parameters*/
void Load_Customfilter(char *file, double * Fwidth_x,double * Fwidth_y) {
  FILE *fid;
  int err,ncol,nalloc,ij;
  int nchunk=1200; /*  (dividible by 2 and 3) */
  int ndat,row,col,i;
  double *dat;
  double min1,min2,max1,max2,range,inc,test_dblx,test_dbly;
  double tol;
  char line[120];
  
/*   Some initializations  */
  GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;	/* Add padding for 2D interpolation */
  tol=0.05; /* Use a 5 percent tolerance on determining whether a lat,lon value is close enough to a grid node*/

  if(!F.gmtgrd){ /* Read from ascii file */
    /*   initialize data array with zeros (use initially ndat  */
    dat=(double*) calloc((size_t)nchunk,sizeof(double));
    
  
    /*open file with read only*/
    fid=fopen(file,"r");
    if(fid==NULL){
      fprintf(stderr,"Load_Customfilter: ERROR opening file, exiting\n");
      exit(1);
    }
    
    /*Read first line for checking the amount of columns*/
    if(fgets(line,120,fid)==NULL){
      fprintf(stderr,"Load_Customfilter: ERROR EMPTY File?\n");
      exit(1);
    }
    
    /*   Determine the amount of columns and read in first data line ( try reading 4 values) */
    ncol=sscanf(line,"%lf %lf %lf %lf",&dat[0],&dat[1],&dat[2],&dat[3]);
    
    if(ncol ==2){ /* Isotropic filter */
      F.ISO=TRUE;
    }else if(ncol==3){ /* Anisotropic filter */
      F.ISO=FALSE; 
    }else{
      fprintf(stderr,"Load_Customfilter: ERROR: file has %d columns??\n",ncol);
      exit(1);
    } 

    /*   Read remaining data */
    ndat=ncol;
    nalloc=nchunk;
    while(TRUE){
      /*     Reallocate a chunk of memory if needed */
      if(nalloc <= ndat){
	/*       fprintf(stderr,"reallocating %ld\n",(nalloc+nchunk)*sizeof(double)); */
	dat=(double *)realloc(dat,sizeof(double)*(nalloc+nchunk));
	nalloc+=nchunk;
      }
      
      
      /*     Read additional value */
      err=fscanf(fid,"%lf",&dat[ndat]);
      /*    fprintf(stderr,"read %d %f",ndat,dat[ndat]); */
      if(err==EOF)break;/*  Exit loop on end of file */
      ndat++;
    }

  
    /*close file*/
    if(fclose(fid)!=0){
      fprintf(stderr,"Load_Customfilter: ERROR closing file, exiting\n");
      exit(1);
    }
  


/*   Now find the extend of the filter to determine the filterwidth */
/*   Intialize values to extremes */
    min1=9999999;
    min2=min1;
    
    max1=-9999999;
    max2=max1;
    
    /* Find filter extend and resort 2Dkernel data */
    if(F.ISO){
      for( i=0;i<ndat;i+=ncol){
	/*       R(adius) values only */
	min1= (min1>dat[i]) ?dat[i] : min1;
	max1= (max1<dat[i]) ?dat[i] : max1;
      }
      /*     Test for nonsense */
      if(min1!=0 || max1<0){
	fprintf(stderr,"Load_Customfilter:ERROR minimum and/or maximumfilter distances must be [0,->] for isotropic case , exiting\n");
	exit(1);
      }
      /*     Set filterwidth ( equally long in x and y direction) */
      *Fwidth_x=*Fwidth_y=2*(max1-min1);

      /* Put data in Kernel variable (possibly resort)*/
      range=max1-min1;

      /*     Calculate increment from range if needed */
      if(F.Kern_h.x_inc==0){
	F.Kern_h.nx=ndat/ncol;
	F.Kern_h.x_inc=range/(F.Kern_h.nx-1);
      
      }else{ /* calculate nx from provided increment (allowing for data gaps (set to zero)) */
	F.Kern_h.nx=floor(range/F.Kern_h.x_inc+0.5)+1;
      }

      /*     Initialize Kernel array to zero */
      F.Kern_x=(double *)calloc(F.Kern_h.nx,sizeof(double));
      F.Kern_y=(double *)calloc(F.Kern_h.nx,sizeof(double));
      
      for( i=0;i<ndat;i+=ncol){
	/*       Try to catch errors for nonequidistant data */
	test_dblx=(dat[i]-min1)/F.Kern_h.x_inc;
	if(fabs(test_dblx-floor(test_dblx+0.5)) > tol){
	  fprintf(stderr,"Load_Customfilter:ERROR Provided filter is not equidistant?, exiting\n");
	  exit(1);
	}
	
	row=(int)((dat[i]-min1)/F.Kern_h.x_inc);
	F.Kern_y[row]=dat[i+1];
	F.Kern_x[row]=dat[i];
      }
    }else{
      for( i=0;i<ndat;i+=ncol){
	/*       X values */
	min1= (min1>dat[i]) ?dat[i] : min1;
	max1= (max1<dat[i]) ?dat[i] : max1;
	/*       Y values */
	min2= (min2>dat[i+1]) ?dat[i+1] : min2;
	max2= (max2<dat[i+1]) ?dat[i+1] : max2;
      }
      
      /*     Test for nonsense */
      /* ( origin must be contained and in the center The filter must also be square equidistant) */
      if((min1+max1)/2!=0 || (min2+max2)/2!=0 ){
	fprintf(stderr,"Load_Customfilter:ERROR 2D kernel grid must be centered on the origin, exiting\n");
	exit(1);
      }
      
/*       Retrieve maximum extends of the filter */
      *Fwidth_x = max1-min1;
      *Fwidth_y = max2-min2;
      
      range = max1-min1;
      /*     Calculate X increment if not provided */
      if(F.Kern_h.x_inc == 0){ /* Assume filter is square */
	F.Kern_h.nx = sqrt(ndat/ncol);
	F.Kern_h.x_inc = range/(F.Kern_h.nx-1);
      }else{/*  calculate number of xdata from increment allowing data gaps*/
	F.Kern_h.nx = floor(range/F.Kern_h.x_inc+0.5)+1;
	
      }

/*       Same stuff for Y */
      range=max2-min2;
      /*     Calculate Y increment if not provided */
      if(F.Kern_h.y_inc == 0){ /* Assume filter is square */
	F.Kern_h.ny = F.Kern_h.nx;
	F.Kern_h.y_inc=range/(F.Kern_h.ny-1);
      }else{/*  calculate number of xdata from increment( allowing data gaps)*/
	F.Kern_h.ny=floor(range/F.Kern_h.y_inc+0.5)+1;
      }
      


      /* Set up remaining values in the grd header structure (needed for interpolation routines)*/
      F.Kern_h.node_offset=0;
      F.Kern_h.x_min=min1;
      F.Kern_h.y_min=min2;    
      F.Kern_h.x_max=max1;
      F.Kern_h.y_max=max2;
      F.Kern_h.xy_off=0.0;
      F.Kern_h.z_scale_factor=1;

   
      /*     Initialize Kernel array to zero */
      F.Kernel=(float *)calloc((size_t)((F.Kern_h.nx+GMT_pad[0]+GMT_pad[1])*(F.Kern_h.ny+GMT_pad[2]+GMT_pad[3])),sizeof(float));
      F.Kern_x=(double *)calloc((size_t)(F.Kern_h.nx),sizeof(double));
      F.Kern_y=(double *)calloc((size_t)(F.Kern_h.ny),sizeof(double));
      
      
      for( i=0; i < ndat; i+=ncol){ /* Loop over all data rows */
	/*       Try to catch errors for nonequidistant data */
	test_dblx=(dat[i]-min1)/F.Kern_h.x_inc;
	test_dbly=(dat[i+1]-min2)/F.Kern_h.y_inc;
	if(fabs(test_dblx-floor(test_dblx+0.5)) > tol){
	  fprintf(stderr,"Load_Customfilter:ERROR Provided filter is not equidistant in x?,entries %f %f\n",test_dblx,floor(test_dblx+0.5));
	  exit(1);
	}
	else if(fabs(test_dbly-floor(test_dbly+0.5))  > tol){
	  fprintf(stderr,"Load_Customfilter:ERROR Provided filter is not equidistant in y?, entries %f %f\n",test_dbly,floor(test_dbly+0.5));
	  exit(1);
	}

	col=GMT_x_to_i(dat[i],F.Kern_h.x_min,F.Kern_h.x_inc,F.Kern_h.xy_off,F.Kern_h.nx);
	row=GMT_y_to_j(dat[i+1],F.Kern_h.y_min,F.Kern_h.y_inc,F.Kern_h.xy_off,F.Kern_h.ny);
	/*       Store data in same scheme as in gmt grids */
	
	ij=GMT_IJ(row+GMT_pad[2],col+GMT_pad[0],F.Kern_h.nx+GMT_pad[0]+GMT_pad[1]);
 	F.Kernel[ij]=(float)dat[i+2]; 

	F.Kern_x[col]=dat[i];
	F.Kern_y[row]=dat[i+1];
      }
      
    
    }
    
    /* Free memory used by dat */
    free(dat);
  
  }else{/*  Read filter from gridfile */
    /*       Header info  */
    GMT_err_fail (GMT_read_grd_info (file, &F.Kern_h), file);    
    F.ISO=FALSE; 
/*     Check if filter is centered around the origin */
    if((F.Kern_h.x_min+F.Kern_h.x_max)/2!=0 || (F.Kern_h.y_min+F.Kern_h.y_max)/2!=0 ){
      fprintf(stderr,"Load_Customfilter:ERROR 2D kernel grid must be centered on the origin, exiting\n");
      exit(1);
    }


/*       Retrieve maximum extends of the filter */
      *Fwidth_x = F.Kern_h.x_max-F.Kern_h.x_min;
      *Fwidth_y = F.Kern_h.y_max-F.Kern_h.y_min;


    /*       Allocate space */
      F.Kernel=(float *)calloc((size_t)((F.Kern_h.nx+GMT_pad[0]+GMT_pad[1])*(F.Kern_h.ny+GMT_pad[2]+GMT_pad[3])),sizeof(float));
/*       read in grid data */
    GMT_err_fail (GMT_read_grd (file, &F.Kern_h, F.Kernel, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), file);

  }
  

    /*  Inititialize interpolation structure (use interpolation as a start */

  if(!F.ISO){
    /* Setup edge info structure */
    GMT_boundcond_init(&F.Kern_edge);
    GMT_boundcond_param_prep (&F.Kern_h, &F.Kern_edge);
/*       Use appropriate padding */
    GMT_boundcond_set (&F.Kern_h, &F.Kern_edge, GMT_pad, F.Kernel);
/*     Prepare fro interpolation */
    GMT_bcr_init (&F.Kern_h, GMT_pad, BCR_BSPLINE, 1.0, &F.Kern_bcr);
  }
}
