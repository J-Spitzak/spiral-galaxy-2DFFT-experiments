//
// P2MAP.CPP - This program will generate the polar mapped projection for FITS
//             galaxy images.  The output files will have the prefix M_ added
//             to avoid collisions.  It will also produce a text table showing
//             the mapping of the polar coordinates to cartesian X, Y values.
//
//
//  Version 1.2: 03-May-2019
//
//
//  2DFFT (Original) Author: Dr. Ivanio Puerari
//                Instituto Nacional de Astrofisca, Optica y Electronica
//                Santa Maria Tonantzintla, Puebla, Mexico
//
//  2DFFT (Revised) Lead Author: Dr. Marc Seigar
//                University of Minnesota, Duluth
//                Duluth, MN USA
//
//  2DFFT (Progenitor of P2DFFT) Lead Author: Dr. Benjamin Davis
//                Swinburne University of Technology,
//                Centre for Astrophysics and Supercomputing,
//                Melbourne, Victoria, Australia
//                http://www.d.umn.edu/~msseigar/2DFFT/2DFFT.tar.gz
//
//  P2DFFT By: Ian Hewitt & Dr. Patrick Treuthardt
//                NC Museum of Natural Sciences
//                Astronomy & Astrophysics Lab,
//                Raleigh, NC USA
//                https://github.com/treuthardt/P2DFFT
//
//  LICENSE
//
//  P2DFFT Spiral Galaxy Arm Pitch Angle Analysis Suite
//  Copyright (c) 2016-2019  Ian B. Hewitt & Dr. Patrick Treuthardt
//
//  The program is free software:  you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the Free
//  Software Foundation, version 3.
//
//  This program is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//  more details.
//
//  You should have received a copy of the GNU General Public License along with
//  this program.  If not, see < https://www.gnu.org/licenses >.
//
//  The authors can be contacted at:
//
//      North Carolina Museum of Natural Sciences
//      Astronomy & Astrophysics Laboratory
//      11 West Jones Street
//      Raleigh, NC, 27601  USA
//      +1.919.707.9800
//
//      -- or --
//
//      patrick.treuthardt@naturalsciences.org
//
//
//  Usage: p2map [-i|--input <file>] [-v|--verbose] [<args>]
// 
//         There are several non-mandatory options:
//              -i|--input  : Will read file names, results file, and radius
//                            from the file specified with this option instead
//                            of standard input.
//              -v|--verbose: Prints status messages during the
//                            processing (good for those who like to see
//                            things during a run).
//
//
//  Version History:
//
//      1.2  03-May-2019 - Correct header comments
//      1.1  13-Nov-2018 - Update error messages to be more consistent
//                       - Correct error handling bug
//                       - Make the matrix file have a M_ prefix
//                       - Add ability to print a text table of values with a
//                         T_ prefix
//      1.0  28-May-2018 - Initial version
//


//
// GLOBAL CONSTANTS
//

#include    "globals.h"


//
// INCLUDE FILES
//

#include    <math.h>
#include    <stdio.h>
#include    <stdlib.h>
#include    <string.h>
#include    <errno.h>
#include    <sys/types.h>
#include    <sys/stat.h>
#include    <unistd.h>
#include    <getopt.h>

//
// Include the Astro Functions Class
//

#include    "astro_class.h"

//
// Version number
//

#define     VERSION     "1.2/20190503"

//
// Set this flag to #define to get a data matrix debugging information.  This
//   is independent of the DEBUG flag in globals.h because it produces a lot
//   of output.  Highly recommend redirecting stdout to a file when using this
//   option.
//

#undef      DEBUG_DAT

//
// Set this flag to #define to get a input matrix debugging information.  This
//   is independent of the DEBUG flag in globals.h because it produces a lot
//   of output.  Highly recommend redirecting stdout to a file when using this
//   option.
//

#undef       DEBUG_MAT

//
//  VARIABLES
//

int     c;                 /* Return value for command line options parser   */
int     fn;                /* Argument index                                 */
int     lim;               /* Initialization size for fft_sum                */
int     num;               /* Number of threads on the host machine          */
int     a, b;              /* Cartesian coordinates of ln(r)/theta in image  */
int     pflag;             /* Table printing flag                            */
int     msize;             /* Binary FITS file data size                     */
int     jm, im;            /* Local index variables                          */
int     radius;            /* Radius for current file                        */
int     status;            /* Return value for scanx() and system() calls    */
int     items=0;           /* Number of files processed                      */
int     i, j, k;           /* Index variables                                */
int     counter;           /* Counter for output array                       */
int     core_val;          /* Core brightness value                          */
int     x_0, y_0;          /* Carteian coordinates for the image center      */
int     verbose=0;         /* Flag for printing of status messages           */
int     count_theta;       /* Counter for theta steps in degrees             */
int     proc_error=0;      /* Input file error count                         */
int     input_file=0;      /* Flag to indicate if input file is used         */
int     x_dim, y_dim;      /* The cartesian dimensions of the input file     */
int     count_radians;     /* Counter for theta steps in radians             */

char    *tmp;              /* Pointer to string for integer conversion       */
char    cmd[128];          /* Buffer for system(2) commands                  */
char    infile[80];        /* Input filename for -i                          */
char    fname[128];        /* FITS filename                                  */

FILE    *table;            /* Output mapping table file pointer              */
FILE    *fp_inp;           /* Input file pointer                             */
    
float   lnr;               /* Natural log of radius for a certain point      */
float   x, y;              /* Relative cartesian coordinates of ln(r)/theta  */
float   **mat;             /* 2D cartesian image data                        */
float   *data;             /* Polar mapped image data matrix                 */
float   **polar;           /* 2D cartesian output image data                 */
float   log_tmp;           /* The natural log of the current radius value    */
float   log_rad;           /* The natural log of the current radius value    */
float   log_itrad;         /* The natural log of the maximum radius value    */
float   freq_counter;      /* Frequency counter value                        */
float   theta_degrees;     /* Current theta (polar angle) in degrees         */
float   theta_radians;     /* Current theta (polar angle) in radians         */

const   float   radstep=2.0*PI/STEP_P/DIM_RAD;    /*                         */
const   float   theta_step=2.0*PI/GR_RAD/DIM_THT; /*                         */

astro   ast;               /* Instantiation of astro_class functions         */
        

//
// MAIN() CODE BLOCK
//

int main(int argc, char **argv)
    {
    printf("p2map version: %s\n", VERSION);
    ast.version();

//
// Parse the command line options, if any, and set the flags associated
//   with the options
//

    static struct option long_options[] =
        {
        {"verbose", no_argument,     0, 'v'},
        /* These options require an argument. */
        {"input", optional_argument, 0, 'i'},
        {0, 0, 0, 0}
        };

    int option_index = 0;

    while ((c = getopt_long (argc, argv, "vi:", long_options, &option_index)
) != -1)
        {
        switch (c)
            {
            case 'v':
                {
                verbose = 1;
                break;
                }
            case 'i':
                {
                input_file = 1;
                if (!ast.file_exists(optarg))
                    {
                    printf("ERROR: Input File %s Not Found...Exiting\n",optarg);
                    exit(-1);
                    }
                strcpy(infile, optarg);
                break;
                }
            default:
                {
                fprintf(stderr, "Usage: p2map [-v|--verbose] [-i <file>] [<n>[,<n>...]]\n");
                exit(-1);
                break;
                }
            }
        }

//
// Allocate the Cartesian data arrays.  Also, zero out the first cell of mat 
//   because FITS image indices start at 1.  Please Note:  ArrayAlloc()
//   allocates a contiguous C-style array.
//

    if (verbose) printf("Allocating Cartesian mat[] Array...\n");

    if (!(mat=ast.ArrayAlloc(MAX_DIM, MAX_DIM)))
        {
        if (mat != NULL) free(mat);
        printf("ERROR: Memory allocation failed while allocating for mat[]/n");
        exit(-1);
        }

    if (verbose) printf("Allocating Cartesian polar[] Array...\n");

    if (!(polar=ast.ArrayAlloc(DIM_RAD, DIM_THT)))
        {
        if (polar != NULL) free(polar);
        printf("ERROR: Memory allocation failed while allocating for polar[]/n");
        exit(-1);
        }

//
// Read the input parameters for the analysis.  The input parameters will 
//   include:
//
//     * All filenames to be processed
//     * Any keywords associated with those files (optional)
//     * Any radius values for the files (optional)
//
// Input can come from one of the following:
//
//     * Input file specified with -i
//     * Command line arguments
//

    if (input_file)
        {
        if (1)
            {
            std::cout << "ERROR: Can't Read File Name: " << infile << std::endl;
            exit(-1);
            }
        if ((items==0))
            {
            std::cout << "ERROR: No Valid Items in Input File: " << infile << std::endl;
            exit(-1);
            }
        }
    else
        {
//
// Check if there are command line arguments
//

        if (DEBUG) printf("optind=%d, argc=%d\n",optind,argc);

        if (optind >= argc)
            {
//
// No command line arguments for files, so read from stdin.  Assume old
//   style 2DFFT input format here.
//

            std::cout << "ERROR: No valud arguments...Exiting" << std::endl;
            exit(-1);
            }
//
// Get the command line arguments and put them in vector of items
//

        for (fn=optind; fn < argc; fn++)
            {
            items++;
            if (DEBUG) printf("argv[%d]=%s\n",fn,argv[fn]);

            if (!(ast.file_exists(argv[fn])))
                {
                std::cout << "WARNING: " << argv[fn] << " Does Not Exist...Skipping" << std::endl;
                proc_error++;
                continue;
                }

            if (ast.file_type(argv[fn]) != ASTRO_BIN_FILE)
                {
                std::cout << "WARNING: Can't Get File Type: " << argv[fn] << " Skipping..." << std::endl;
                proc_error++;
                continue;
                }

            x_dim=0;
            y_dim=0;

            std::cout << "Processing Entry - Name: " << argv[fn] << std::endl;

// 
// Read the data from the image and determine the radius
//

            
            if (!(data=ast.fits_read(argv[fn], &msize)))
                {
//
// Read Failure
//

                std::cout << "WARNING: Can't Read FITS Binary File: " << argv[fn] << " Skipping..." << std::endl;
                proc_error++;
                continue;
                }

//
// Get the radius
//

            if (ast.fits_dims(argv[fn],&x_dim, &y_dim))
                {
//
// Failure to get size from Header
//

                std::cout << "ERROR: Can't Read FITS Dimensions for " << argv[fn] << " Skipping..." << std::endl;
                proc_error++;
                continue;
                }

            printf("FITS DIMS: X_DIM=%d, Y_DIM=%d\n",x_dim, y_dim);

//
// Find radius
//

            if ( x_dim < y_dim )
                {
                radius=(x_dim-1)/2;
                }
            else
                {
                radius=(y_dim-1)/2;
                }

//
// Zero output array
//

            for (i=0; i < DIM_THT; i++)
                {
                for (j=0; j < DIM_RAD; j++)
                    {
                    polar[j][i]=0.0;
                    }
                }

//
// Copy the FITS data into the mat 2D Cartesian array.  Need to take this from
//   FITS ordering (column major) to C ordering (row major).   In addition,
//   note that the X-axis is also reversed.
//

#ifdef DEBUG_DAT
            for(i=0;i<(DIM_RAD*DIM_THT*2+2); i++)
                {
                printf("DEBUG: data[%d]=%f\n",i,data[i]);
                }
#endif

            counter=0;
            for(j=1;j<=y_dim;j++) 
                {
                for(i=1;i<=x_dim;i++)
                    {
                    mat[i][j]=data[counter++];

#ifdef DEBUG_MAT
                    printf("DEBUG: mat[%d][%d]=%f\n",i,j,mat[i][j]);
#endif

                    }
                }

//
// Print out the intermediate matrix.  NOTE: since the row/column has been
//   changed, it will have a rotation in the image.
//

            sprintf(fname,"!M_%s.fits",argv[fn]);

            if (ast.fits_write(fname, mat[0], MAX_DIM, MAX_DIM, 1, "p2map/",VERSION))
                {
                printf("ERROR: fits_write(%s) Failed\n",fname);
                proc_error++;
                }

            sprintf(fname,"T_%s.txt",argv[fn]);

            if ((table=fopen(fname,"w"))==NULL)
                {
                printf("ERROR: Could Not Write T_%s.txt\n",fname);
                exit(1);
                }

            fprintf(table,"File Mapping: %s\n", fname);

            if (verbose) std::cout << "Processing Entry - Name: " << argv[fn] << " Radius: " << radius << std::endl;

//
// Use (dim-1)/2 for each dimension.  This makes it work for both odd and even
//   sized images.  Need to calculate both because image may not be rectangular.
//

            x_0=((x_dim-1)/2)+1;
            y_0=((y_dim-1)/2)+1;
            fprintf(table,"X_0=%d, Y_0=%d\n",x_0,y_0);
            fprintf(table,"Radius\tln(R)\tX\tY\tRel X\tRel Y\n");
            fprintf(table,"------\t-----\t-\t-\t-----\t-----\n");

//
// log() functions are computationally expensive, so calculate the logs
//   outside of the loop.
//

            log_rad=log((double)radius);
            log_tmp=log((double)76);

//
// Step around theta angles (360 degrees in 0.35 steps)
//

            core_val=mat[x_0][y_0];
            pflag=1;
            count_theta=0;
            for(theta_degrees=0.0;count_theta<DIM_THT;theta_degrees+=theta_step) 
                {

//
// Convert the degrees to radians
//

                theta_radians=theta_degrees*GR_RAD;	
                count_radians=0;

                for(lnr=0.0;count_radians<DIM_RAD;lnr+=radstep) 
                    {
                    if (lnr>log_rad)
                        {
                        count_radians++;
                        continue;
                        }
                    x=expf(lnr)*cosf(theta_radians);
                    y=expf(lnr)*sinf(theta_radians);

                    a=(int)x+x_0;
                    b=(int)y+y_0;

                    if (mat[a][b] < core_val - 3)
                        polar[count_radians][count_theta]=(float) mat[a][b];
                    
                    if (pflag) fprintf(table,"%f\t%f\t%d\t%d\t%d\t%d\n",expf((float)radius),lnr,a,b,a-x_0,b-y_0);
                    count_radians++;
                    }
                pflag=0;
                count_theta++;
                }

//
// Do a reverse mapping Step around theta angles (360 degrees in 0.35 steps)
//

            count_theta=0;
            for(theta_degrees=0.0;count_theta<DIM_THT;theta_degrees+=theta_step) 
                {

//
// Convert the degrees to radians
//

                theta_radians=theta_degrees*GR_RAD;	
                count_radians=0;

                for(lnr=0.0;count_radians<DIM_RAD;lnr+=radstep) 
                    {
                    if (lnr>log_rad)
                        {
                        count_radians++;
                        continue;
                        }
                    x=lnr*cosf(theta_radians);
                    y=lnr*sinf(theta_radians);

                    a=(int)x+x_0;
                    b=(int)y+y_0;

                    mat[a][b]=polar[count_radians][count_theta];
                    
                    count_radians++;
                    }
                count_theta++;
                }
//
// Now create a FITS file using the astro_class libraries.  Please note that 
//   when the new flag is set to 1 it will overwrite any existing file with
//   with the same name.
//
  
            if (verbose) printf("  --- Write P_%s.fits File\n",argv[fn]);

            sprintf(fname,"!P_%s.fits",argv[fn]);

            ast.set_warn(1);
            if (ast.fits_write(fname, polar[0], 1024, 2048, 1, "p2map/",VERSION))
                {
                printf("ERROR: fits_write(%s) Failed\n",fname);
                proc_error++;
                }
            sprintf(fname,"!R_%s.fits",argv[fn]);

            ast.set_warn(1);
            if (ast.fits_write(fname, mat[0], MAX_DIM, MAX_DIM, 1, "p2map/",VERSION))
                {
                printf("ERROR: fits_write(%s) Failed\n",fname);
                proc_error++;
                }
            fclose(table);
            }
        }

    printf("-------------------------------\n");
    printf("Successfuly Processed        %d\n",items-proc_error);
    printf("Errors                       %u\n",proc_error);
    }
