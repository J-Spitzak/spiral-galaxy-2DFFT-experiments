//
// P2DFFT.CPP - This program will calculate the 2D FFT for determining pitch
//              angles.
//
//  Version 5.9: 20-Jun-2019
//
//           The program was derived from 2DFFT at the NC Museum of Natural
//           Sciences Astronmy & Astrophysics Research lab.  It calculates the
//           pitch angles of spiral galaxy arms using the Fourier tranform
//           method, but uses a different output/interpretation method for the
//           FFT output data.  It runs the calculations in parallel, which
//           makes for much faster execution
//           on multi-threaded machines.  In addition, there have been some 
//           changes to the input to the programs to increase usability. 
//           Finally, the program has changed the FFT algorithm to the 
//           Harvard "Fastest Fourier Transform in the West" (FFTW3) library.
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
//  Usage: p2dfft [-i|--input <file>] [-v|--verbose] [-w|--warn]  [-r|--reverse]
//                [-f|--fixed <size>] [-p|--polar] [-z|--zero] [-m|--mask 0,1]
//                [-h|--highpass] [<args>]
// 
//         p2dfft will process a list of files.  These files can come from 
//         standard input, the command line, or an input file.  The files can
//         be binary FITS files or ASCII text FITS files. By default, no
//         radius or output file prefix needs to be provided as these will
//         determined by the program.
//
//         There are several non-mandatory options:
//              -i|--input  : Will read file names, results file, and radius
//                            from the file specified with this option instead
//                            of standard input.
//              -v|--verbose: Prints status messages during the
//                            processing (good for those who like to see
//                            things during a run).
//              -w|--warn   : Causes warnings to be printed when individual
//                            values are not rational or minor processing 
//                            problems are discovered.  This will slow 
//                            performance and is mainly useful if unexpected
//                            results are seen.
//              -r|--reverse: By default the pitch angle is calculated using
//                            different annuli with increasing inner radius,
//                            this option performs the calculations with
//                            annuli with decreasing outer radius.
//              -f|--fixed  : Use annuli with a fixed radius (specified by
//                            the argument to fixed.  The inner radius of
//                            the annuli will start from 1 and increase but
//                            the end will alway be that value plus the 
//                            argument to fixed.
//              -p|--polar  : Generate a FITS file with the logarithmic polar
//                            mapping used by P2DFFT. 
//              -z|--zero   : Generate a zero filled padding around the polar
//                            projection to simulate an FFT window
//              -m|--mask   : Option to maks bright values (core and bluge).
//                            A value of 0 will mask out any values >= to the
//                            value at the center.  A mask value of 1 will
//                            calculate the highest radius of continuous
//                            values that have a value >= the center value
//                            and will eliniate (zero) all values in the polar
//                            plot below that radius.  This eliminates 
//                            structural artefacts that could affect the
//                            data when a mask value of 0 is used.
//              -h|--highpass: Apply a high pass filter to the results after
//                             the low pass filter is applied (experimental)
//
//
//  Input formats:
//
//         Command Line Arguments - The file names of binary FITS and ASCII
//         text files are specified on the command line.  The result file for
//         each file will be all characters of the file name before the first
//         occurrence of a "."  The radius will be the maximum radius based on
//         the size of the image.
//
//         Input File - The format of the input file is:
//
//           image_file_1,result_file_1,outer_radius_1
//           image_file_2,result_file_2,outer_radius_2
//           image_file_3,result_file_3,outer_radius_3
//           image_file_etc,result_file_etc,outer_radius_etc
//       
//        The result file name and radius are optional and so are the commas
//        for those fields (e.g. you can just have a line like "image_file_1"
//        or "image_file_1,result_file_1".  If they are not specified, the
//        program will calculate them.  Previous scripter input files can be
//        used as inputs to this version of p2dfft.  Blank lines are allowed
//        in the file, but are no longer needed and will be ignored.
//
//        Standard Input - The standard input format is for compatibility
//        with older scripts generated by scripter.  The expected format is
//        (one per line):
//
//           image_file_name
//           result_file
//           keyword
//           radius
//
//        Only one file can be read from standard input.  The keyword is kept
//        for compatiability reasons, but not used in P2DFFT v5.0+
//
//  Version History:
//
//      5.9  20-Jun-2019 - Remove incorrect FFT row sense fix
//      5.8  02-Jun-2019 - Change 2D FFT row sense to correct bug left over
//                         from FFTW migrations and to improve accuracy
//                       - Remove deprecated leftover -s argument code
//      5.7  21-May-2019 - Fix incorrect original author attribution with
//                         apologies to those previously left out of the list
//      5.6  08-May-2019 - Eliminate the -s option so the program will always
//                         create the FFT data files
//                       - Add additional debugging to show .rip file mapping
//      5.5  18-Dec-2018 - Fix bug where save files are not created in the 
//                         current working directory if fits file is not
//                         located in the current working directory
//                       - Fix bug where failure to write the the .dat or
//                         .rip files would cause a crash
//      5.4  20-Oct-2018 - Change the version printing to only happen when -v
//                         is specified
//                       - Add -m|--mask option to mask bright core values
//                       - Add -z|--zero option to add artificial FFT window
//                       - add -p|--polar option to generate polar FITS image
//                       - Change the mapping from the fits file to account
//                         for the FITS/C++ row ordering difference.
//                       - Change version printing to only occur with -v
//                       - Bug fix for failure when using non local arguments
//                       - Update comments and remove unused code blocks
//      5.3  19-Jun-2018 - Fix two bugs with sum output to address potential
//                         race conditions which caused indeterminate output
//                       - Fix bug in summation when multiple files specified
//                       - Fix bug in reverse mode where annuli changing in
//                         wrong direction
//      5.2  27-May-2018 - Add -r|--reverse and -f|--fixed options
//                       - Simplify center location logic
//                       - Remove unused variables
//                       - Fix bug that resulted in center being off by 1
//                         for certain sizes
//                       - Update comments/usage info
//                       - Remove unused variables
//                       - Move fftw plan creation so it's only done once
//                       - instead of for every file
//                       - Add summed FFT data files to output
//                       - Remove code to limit threads because it was too
//                         complex for marginal usefulness
//                       - Remove fftw3 multithreaded code and replace with
//                         OpenMP code for better utilization of threads
//                       - Replace fixed frequency values with FREQ constants
//      5.1  15-Mar-2018 - Add support for sn() and fwhm() functions and adding
//                         their results to the output
//                       - Clean up error handling for pitch_class to return
//                         less extraneous output
//                       - Add -w|--warning option to allow/supress warning
//                         messages being printed
//                       - Fix small bug in getopt_long() call
//                       - Preface each printed message with aa class (warning,
//                         debug, or error)
//                       - Fix bug where an error would not cause an orderly
//                         exit
//                       - Print the total items and how many were successfully
//                         processed at the end of the program
//                       - Significant comment updates
//      5.0  05-Feb-2018 - Change name and comments from 2dfft to p2dfft
//                       - Significant changes to file I/O scheme.  Will no
//                         longer write/use temporary files (.rip & .data) 
//                         if saverip is not enabled.
//                       - saverip data files are no longer created in the
//                         current working directory, but are created in the
//                         <result_file> subdirectory.
//                       - Remove external give_maximum_pitch_phase executable
//                         and replace it with a pitch_class function.  This
//                         eliminates the need for FORTRAN and a multitude
//                         of interim working files.
//                       - Remove the code to automatically remove working
//                         files in the current directory and change the code
//                         which creates the <result_file> output files
//                       - Add comments for variables and change some variable
//                         names to de-obsfucate the code :-)
//                       - Clean up version print code and add class versions
//                       - Move global constants to a separate include file
//                       - Remove requirement that images files are square
//                       - Add capability to override default number of threads
//                       - Change DEBUG to be a variable, not a comp. cond.
//                       - Add comments for array mapping information
//                       - Add error handling for output files, including 
//                         deselecting saverip if data files can't be opened
//      4.0  28-Aug-2017 - Rewrite program to use FFTW3 libraries
//                       - Make some changes to mapping algorithm to improve
//                         performance
//                       - Eliminate creates of the gal.* file
//                       - Updated header comments
//                       - Add code for more detailed debugging output
//                       - Updated/fixed comments
//                       - Fixed bug in reading from input file where no files
//                         would be found, even in a valid input file
//                       - Minor changes to eliminate compiler warnings on
//                         some Linux distributions
//      3.4  27-Jul-2017 - Print limited status information even when -v is
//                         not specified
//      3.3  09-Jun-2017 - Correct an error in the mapping of ASCII fits files
//                       - Update comments
//                       - Fix some error messages
//                       - Remove deprecated -f/--force_read logic
//                       - Add code to print out version when run
//                       - Correct help printout error
//      3.2  07-Jun-2017 - Print version number to stdout
//      3.1  21-Feb-2017 - Remove unused variable opt_err
//      3.0  20-Feb-2017 - Significant rewrite to add the following changes:
//                          * Move setup blocks of code and rearrange for
//                            readability.
//                          * Move from C to C++ (to get classes)
//                          * Program always builds for for multiple threads
//                            (will still work on single threaded machine)
//                          * Make other compile time options into command line
//                            options (STATUS and SAVE_RIP)
//                          * Add astro class functions for working with
//                            FITS files (see astro_class.cpp)
//                          * Eliminate the need for scripter.  2dfft can
//                            now read the input files originally meant for
//                            scripter.  This is done with the -i option.
//                          * If the input file does not have the result 
//                            file (used to be keyword) or radius, it will
//                            be calculated automatically.
//                          * 2dfft can now take image file names via the
//                            command line.
//                          * 2dfft now reads both ASCI FITS image and binary
//                            FITS images (of any sort).
//      2.5  03-Jan-2017 - Add check during reading of input file.  If the file
//                         has a FITS header, it will loop forever.  Now it
//                         will abort with an error message and move on to the
//                         next file.
//      2.4  30-Dec-2016 - Add check for EOF when looking for input file.  This
//                         prevents 2DFFT from getting in an infinite loop when
//                         the input file does not exist.
//      2.3  27-Nov-2016 - Change the program to now look for the
//                         give_maximum_pitch_phase program based on the
//                         BIN_DIR entry in the makefile.  The program will:
//                          * Look in the current directory and use that one
//                          * Look in the BIN_DIR and use that one
//                          * If not found, search some common locations
//                         The first executable version found will be used.
//      2.2  27-Sep-2016 - Some changes to address:
//                          * In parallel mode, add reading of the final output
//                            filename so we can make the *_m[0-6] files 
//                            correctly (so it works with scripter.cpp v2.2)
//                          * Reduced loop to calculate one less iteration in
//                            parallel mode to give exact match for old results.
//			    * Fixed bug in non-parallel mode output by setting
//                            radius = r_min.  This allowed a number of 
//                            conditional compile statements to be removed, 
//                            which simplified the code.
//			    * Fix path name for give_maximum_pitch_phase which
//                            was hard-coded. Now it just relies on it being in
//                            the command path. <This change is deprecated>
//      2.1  26-Sep-2016 - Fixed two subtle race conditions in the parallel
//                         loop.  Also cleaned up STATUS messages to be
//                         cleaner and print less (for better performance)
//      2.0  12-Sep-2016 - Update for openmp parallel threads and to add
//                         compile time options (See Makefile and Readme.txt)
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
#include    <omp.h>
#include    <fftw3.h>
#include <libgen.h>
//
// GLOBAL CONSTANTS
//

#include    "globals.h"

//
// Include the Astro Functions Class from the NCNMS/NRC
//

#include    "astro_class.h"
#include    "pitch_class.h"

//
// Version number definition
//

#define     VERSION     "5.9/20190620"

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

#undef      DEBUG_MAT

//
//  VARIABLES
//

int     c;                 /* Return value for command line options parser   */
int     lim;               /* Initialization size for fft_sum                */
int     num;               /* Number of threads on the host machine          */
int     msize;             /* Binary FITS file data size                     */
int     mask=0;            /* Flag for masking only high values              */
int     zero=0;            /* Flag ot insert zero padding in FFT data        */
int     warn=0;            /* Flag to indicate if warnings are printed       */
int     dindex;            /* Counter for debug statement counting           */
int     status;            /* Return value for scanx() and system() calls    */
int     fixed=0;           /* Flag for fixed annuli for calculations         */
int     polar=0;           /* Flag to control if polar proj image created    */
int     i, j;              /* Index variables                                */
int     counter;           /* Index variable for the polar data array        */
int     x_0, y_0;          /* Carteian coordinates for the image center      */
int     offset=0;          /* Index for start of image data in input array   */
int     reverse=0;         /* Flag to control if inner or outer radis varies */
int     verbose=0;         /* Flag for printing of status messages           */
int     proc_error;        /* Input file error count                         */
int     high_pass=0;       /* Flag for applying high pass filter             */
int     mask_line=0;       /* Flag for masking on an even line               */
int     input_file=0;      /* Flag to indicate if input file is used         */
int     x_dim, y_dim;      /* The cartesian dimensions of the input file     */

unsigned    int     it;    /* Files vector index variable                    */

char    *tmp;              /* Pointer to string for integer conversion       */
char    *base;             /* Base name for FFT data output files            */
char    *fname;            /* FITS filename                                  */
char    cmd[128];          /* Buffer for system(2) commands                  */
char    pfile[80];         /* Input filename for -i                          */
char    infile[80];        /* Input filename for -i                          */
char    keyword[80];       /* String for intermediate data file prefix       */
char    outfile[80];       /* String for intermediate file name              */
char    tmpofile[80];      /* Intermediate data file file name               */
char    resultfile[80];    /* Summary file (*_m[0-6]) file name              */

FILE    *fp_inp;           /* Input file pointer                             */
FILE    *sum_out;          /* Output file pointer for per mode summed data   */
FILE    *mode_out;         /* Output file pointer for per mode peak data     */
    
float   **mat;             /* 2D cartesian image data                        */
float   *data;             /* Polar mapped image data matrix                 */
float   *proj;             /* Polar mapped image data matrix                 */
float   ctr_val;           /* Core brightness for masking                    */
float   log_rad;           /* The natural log of the current radius value    */
float   log_bar;           /* The natural log of the bar radius value        */
float   log_itrad;         /* The natural log of the maximum radius value    */
float   freq_counter;      /* Frequency counter value                        */

const   float   radstep=2.0*PI/STEP_P/DIM_RAD;    /*                         */
const   float   theta_step=2.0*PI/GR_RAD/DIM_THT; /*                         */

astro   ast;               /* Instantiation of astro_class functions         */
pitch   pit;               /* Instantiation of pitch_class functions         */
        
fftw_plan   plan;          /* FFTW execution plan variable                   */

std::vector  <file_rec>    items; /* Vector of input files                   */

struct  result_pa   mode_data[M_FIN+1][(MAX_DIM/2)+1];   /* FFT analysis data*/

//
// SUBROUTINES
//


//
// REMOVE_EXTENSION() - Takes a filename and strips the extension and returns
//                      that value.  If not extension is found, the original
//                      argument was returned.
//
// Arguments:
//      filename - string with path/filname
//
// Return Value:
//      String with filename without extension, or if no extension was present,
//      the original filename is returned
//

std::string remove_extension(std::string& filename)

    {
//
// Locate potential extension
//

    size_t lastdot = filename.find_last_of(".");

//
// If not dot, we are done
//

    if (lastdot == std::string::npos) return filename;

//
// Need to check that this dot is not part of the path
//

    size_t slashpos = filename.find_last_of("/");

    if (slashpos != std::string::npos)
        {
        if (slashpos > lastdot) return filename;
        }

//
// Strip extension
//

    return filename.substr(0, lastdot); 
    }


//
// READ_STD_INPUT() - Reads parameters for processing/analysis from std input.
//                    NOTE: In this version, it is assumed that the format is
//                    is the same that was used in the original 2DFFT.
//
// Arguments:
//      NONE
//
// Return Value:
//      NONE    - Just sets global parameters.  Exits program if error.
//

void    read_std_input()
    {
    int     r_max;             /* Radius value from standard input           */

    if (verbose) puts("--- reading filename");

//
// Get input filename
//

    do 
        {
        if (verbose) puts("--- waiting...");

//
// Check for EOF here.  This prevents p2dfft running forever if the input 
//   file cannot be found.
//

        if (scanf("%s",infile)==EOF)
            {
            printf("ERROR: Could Not Read Input File: %s\n",infile);
            exit(-1);
            }
        } while(!ast.file_exists(infile));

//
// Get the result file name, keyword for the output file, and the radius
//

    if (scanf("%s",resultfile)==EOF)
        {
        printf("ERROR: Unexpected End of std input while getting resultfile\n");
        exit(-1);
        }

    if (scanf("%s",keyword)==EOF)
        {
        printf("ERROR: Unexpected End of std input stream while getting keyword\n");
        exit(-1);
        }

    if (scanf("%d",&r_max)==EOF)
        {
        printf("ERROR: Unexpected End of std input stream\n");
        exit(-1);
        }

//
// Create a C++ vector for each file entry
//

    file_rec    ff;

    ff.name=std::string(infile);
    if ((ff.binary=ast.file_type(ff.name)) == -1) exit(-1);

    ff.result=std::string(resultfile);
    ff.keyword=std::string(keyword);
    ff.radius=r_max;
    ff.valid=1;
    items.push_back(ff);
    }


//
// FIND_BAR() - This routine is used with the mask option.  It will start by
//              in the center of the image and search for the LARGEST radius
//              which has a pixel value greater than the limit provided.  That
//              will be assume ot be the bar radius.
//
// Arguments:
//      rad -     Outer radius of image
//      x_org -   X coordinate of center point
//      y_org -   Y coordinate of center point
//      lim_val - Limit value for masking
//
// Return Value:
//      Radius of esimated bar
//

float   find_bar(int rad, int x_org, int y_org, float lim_val)
    {
    int     skip;          /* Loop variable set to 1 after low value found   */
    int     aa, bb;        /* Cartesian coordinates of ln(r)/theta in image  */
    int     cnt_rad;       /* Counter for theta steos in radians             */
    int     cnt_tht=1;     /* Counter for theta steps in degrees             */

    float   r;             /* Natural log of radius for a certain point      */
    float   lb;            /* Largest bar radius value                       */
    float   curr;          /* Pixel value at current radius (for slope)      */
    float   prev1;         /* Pixel value at current radius -1 (for slope)   */
    float   prev2;         /* Pixel value at current radius -1 (for slope)   */
    float   prev3;         /* Pixel value at current radius -1 (for slope)   */
    float   xx, yy;        /* Relative cartesian coordinates of ln(r)/theta  */
    float   log_edge;      /* Natural log of current value of radius         */
    float   tht_deg;       /* Current theta (polar angle) in degrees         */
    float   tht_rad;       /* Current theta (polar angle) in radians         */

    printf("Rad=%d, X_org=%d, Y_org=%d, Lim_val=%f\n",rad,x_org,y_org,lim_val);

//
// log() functions are computationally expensive, so calculate the logs
//   outside of the loop.
//

    log_edge=log((double) rad);
    printf("Log_edge=%g\n",log_edge);
    lb=0.0;

//
// Step around theta angles (360 degrees in 0.35 steps)
//

    for(tht_deg=0.0;cnt_tht<=DIM_THT;tht_deg+=theta_step)
        {
        cnt_tht++;
        tht_rad=tht_deg*GR_RAD;
        cnt_rad=1;
        skip=0;
        curr=-65535.0;
        prev1=-65535.0;
        prev2=-65535.0;
        prev3=-65535.0;

//
// Step out along radius steps
//

        for(r=0.0;cnt_rad<=DIM_RAD;r+=radstep)
            {
            cnt_rad++;

            if (skip) continue;
            if (r > log_edge) continue;

            xx=expf(r)*cosf(tht_rad);
            yy=expf(r)*sinf(tht_rad);

            aa=(int)xx+x_org;
            bb=(int)yy+y_org;

//
// Set the slope values for calculation - TBA function
//

            prev3=prev2;
            prev2=prev1;
            prev1=curr;
            curr=mat[aa][bb];
            

            if (DEBUG) printf("R=%f, Mat[%d][%d]=%f\n",r, aa,bb,mat[aa][bb]);
            if (mat[aa][bb] >= lim_val)
                {
                if (r > lb) lb=r;
                }
            else
                {
                skip=1;
                }
            }
        }

    printf("--- bar length: %d (%f)\n",(int) expf(lb),lb);
    return(lb);
    }


//
// MAIN() CODE BLOCK
//

int main(int argc, char **argv)
    {
//
// Parse the command line options, if any, and set the flags associated
//   with the options
//

    static struct option long_options[] =
        {
        {"data", no_argument,        0, 'd'},
        {"polar", no_argument,       0, 'p'},
        {"zero",  no_argument,       0, 'z'},
        {"warning", no_argument,     0, 'w'},
        {"verbose", no_argument,     0, 'v'},
        {"reverse", no_argument,     0, 'r'},
        {"highpass", no_argument,    0, 'h'},
        /* These options require an argument. */
        {"mask",  optional_argument, 0, 'm'},
        {"fixed", optional_argument, 0, 'f'},
        {"input", optional_argument, 0, 'i'},
        {0, 0, 0, 0}
        };

    int option_index = 0;

    while ((c = getopt_long (argc, argv, "pzwvrhm:f:i:", long_options, &option_index)
) != -1)
        {
        switch (c)
            {
            case 'p':
                {
                polar = 1;
                break;
                }
            case 'z':
                {
                zero = 1;
                break;
                }
            case 'v':
                {
                verbose = 1;
                break;
                }
            case 'r':
                {
                reverse = 1;
                break;
                }
            case 'h':
                {
                high_pass = 1;
                break;
                }
            case 'w':
                {
                warn = 1;
                pit.set_warn(1);
                ast.set_warn(1);
                break;
                }
            case 'm':
                {
                if (atoi(optarg) != 0)
                    {
                    mask_line=1;
                    }
                else
                    {
                    mask=1;
                    }
                break;
                }
            case 'f':
                {
                fixed=atoi(optarg);
                if ((fixed > MAX_WINDOW) || (fixed < MIN_WINDOW))
                    {
                    printf("ERROR: Window Size Must Be Between %d and %d...Exiting\n",MIN_WINDOW,MAX_WINDOW);
                    exit(-1);
                    }
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
                fprintf(stderr, "Usage: p2dfft [-i|--input <file>] [-v|--verbose] [-w|--warn]  [-r|--reverse] [-f|--fixed <size>] [-p|--polar] [-z|--zero] [-m|--mask 0|1] [<args>]\n");
                exit(-1);
                break;
                }
            }
        }

//
// Output version numbers
//

    if (verbose)
        {
        printf("p2dfft version: %s\n", VERSION);
        ast.version();
        pit.version();
        }

//
// Check for conflicting arguments
//

    if (fixed && reverse)
        {
        printf("ERROR: Cannot specify -r|-reverse and -f|--fixed...Exiting\n");
        exit(-1);
        }

//
// Allocate the Cartesian data array.  Also, zero out the first cell of mat 
//   because FITS image indices start at 1.  Please Note:  ArrayAlloc()
//   allocates a C-style array, not a FFTW3 style array.
//

    if (verbose) printf("Allocating Cartesian mat[] Array...\n");

    if (!(mat=ast.ArrayAlloc(MAX_DIM, MAX_DIM)))
        {
        if (mat != NULL) free(mat);
        printf("ERROR: Memory allocation failed while allocating for mat[]/n");
        exit(-1);
        }

    mat[0][0]=0.0;

//
// Get number of threads for this machine.  By default this should return
//   a value = #cores * threads per core.
//

    num=omp_get_max_threads();

//
// Allocate structure array for sum of FFT outputs and initialize the
//   frequencies
//

    i=M_FIN-M_INI+1;
    lim=(int)((fabs(FREQ_START)+fabs(FREQ_END))*4.0)+1.0;

    struct  fft_out     fft_sum[i][lim];  /* Array for sum of FFT outputs */

    freq_counter=FREQ_START;

    for (i=0; i<lim; i++)
        {
        for(j=M_INI; j<=M_FIN; j++)
            {
            fft_sum[j][i].freq=freq_counter;
            }
        freq_counter+=STEP_P;
        }

//
// Allocate the FFT arrays.  These need to be allocated with fftw_ functions
//   since they are not C-style 2D arrays and the fact they need to be aligned
//   on 16 byte boundaries if the target machine has SIMD support.
//
// Note we need to allocate one set per thread, since this program is managing
//   the threads and not the FFTW library.
//

    if (polar) proj = (float *) malloc((DIM_RAD*DIM_THT+1) * sizeof(float));
    struct  fft_out     fft_data[num][DIM_RAD+2];  /* FFT output data array */

    fftw_complex    *in_data[num];
    fftw_complex    *out_data[num];

    for ( i=0; i < num; i++ )
        {
        in_data[i] = (fftw_complex *) fftw_malloc((DIM_RAD*DIM_THT+1) * sizeof(fftw_complex));
    if(NULL == in_data[i])
            {
            printf("ERROR: FFTW Memory allocation failed for in_data[%d]/n",i);
            exit(-1);
            }

        out_data[i] = (fftw_complex *) fftw_malloc((DIM_RAD*DIM_THT+1) * sizeof(fftw_complex));
        if(NULL == out_data[i])
            {
            printf("ERROR: FFTW Memory allocation failed for out_data[%d]/n",i);
            exit(-1);
            }
        }
        
//
// Read the input parameters for the analysis.  The input parameters will 
//   include:
//
//     * All filenames to be processed
//     * Any keywords associated with those files (optional)
//     * Any radius values for the files (optional)
//
// Input can come from one of the following sources (in 
//  priority order):
//
//     * Input file specified with -i
//     * Command line arguments
//     * Std input
//

    if (input_file)
        {
        if (ast.read_lines(std::string(infile), &items))
            {
            std::cout << "ERROR: Can't Read File Name: " << infile << std::endl;
            exit(-1);
            }
        if ((items.size()==0))
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

            read_std_input();
            }
        else
            {
//
// Get the command line arguments and put them in vector of items
//

            for (i=optind; i < argc; i++)
                {
                if (DEBUG) printf("argv[%d]=%s\n",i,argv[i]);

                if (ast.file_exists(argv[i]))
                    {
                    file_rec    f;
                   
                    f.name=std::string(argv[i]);
                    f.result=remove_extension(f.name);
                    f.keyword="outi";
                    f.radius=-1;
                    f.valid=0;
//
// Up to this point, don't assume it's a rational file, just populate the
//   defaults.  This next bit tests if it's a binary file and if that is 
//   true, we can put it on the list.
//
                    if ((f.binary=ast.file_type(f.name))!= -1) items.push_back(f);
                    }
                }
            }
        }

//
// Final check to make sure we have items.   No Reason to Fail Here, but......
//

    if (items.size() == 0)
        {
        printf("ERROR: No Valid Files to Process (Empty work list)\n");
        exit(-1);
        }
    else
        {
        printf("Total files to Process:    %u\n",(unsigned int)items.size());
        }

    proc_error=0;

//
// Build the plan for the FFT transform
//

    if (verbose) printf("Building plan for FFTW...");
    plan=fftw_plan_dft_2d( (int) DIM_THT, (int) DIM_RAD, in_data[0], out_data[0], FFTW_FORWARD, FFTW_MEASURE);
    if ( plan == NULL )
        {
        printf("ERROR: FFTW Plan (%d) Build Failed\n",i);
        exit(1);
        }
    if (verbose) printf("Done\n");


//
// MAIN PROCESSING LOOP
//

//
// Now we have the list of files.  Loop through all fo them and process.
//

    for ( it = 0; it < items.size();  it++)
        {
//
// Zero out x_dim and y_dim.  This is important for the logic to 
//   determine the radius
//

        x_dim=0;
        y_dim=0;

        std::cout << "Processing Entry - Name: " << items[it].name << std::endl;
        if (DEBUG) std::cout << " Result: " << items[it].result << " Keyword: " << items[it].keyword << " Radius: " << items[it].radius << " Binary: " << items[it].binary << " Valid: " << items[it].valid << std::endl;

//
// Zero out the fft_sum array abs values
//

        lim=(int)((fabs(FREQ_START)+fabs(FREQ_END))*4.0)+1.0;

        for (i=0; i<lim; i++)
            {
            for(j=M_INI; j<=M_FIN; j++)
                {
                fft_sum[j][i].abs=0.0;
                }
            }

// 
// Read the data from the image.  P2DFFT can read either a FITS ASCII .txt
//   file or binary FITS file.  Also determine the radius, if needed.
//

        if (items[it].binary)
            {
//
// It's a binary FITS file - Data will start at location 0 from fits_read()
//
            
            offset=0;
            fname=(char *) items[it].name.c_str();
            if (!(data=ast.fits_read(fname, &msize)))
                {
//
// Read Failure
//

                std::cout << "WARNING: Can't Read Binary File: " << items[it].name << " Skipping..." << std::endl;
                proc_error++;
                continue;
                }
//
// Get the radius
//

            if (ast.fits_dims(items[it].name.c_str(),&x_dim, &y_dim))
                {
//
// Failure to get size from Header
//

                std::cout << "ERROR: Can't Read Binary File Size: " << items[it].name << " Skipping..." << std::endl;
                proc_error++;
                continue;
                }

//
// Find radius.  Images are no longer required to be square so find the 
//   shortest dimension for the radius.
//

            if (!items[it].valid)
                {
                if ( x_dim < y_dim )
                    {
                    items[it].radius=(x_dim-1)/2;
                    items[it].valid=1;
                    }
                else
                    {
                    items[it].radius=(y_dim-1)/2;
                    items[it].valid=1;
                    }
                }
            }
        else
            {
//
// It's a ASCII FITS file -- IMPORTANT NOTE: These type of files must have
//   two bytes for size information.  The bytes can be zero, but must be
//   be there or the first two bytes of the data will be ignored and the
//   alignment of the other data incorrect, which will lead to changes
//   in the output values.
//

//
// Read ASCII File into data
//

            if (verbose) puts("--- reading image");

            if((fp_inp=fopen(items[it].name.c_str(),"r"))==NULL)
                {
                std::cout << "WARNING: Problem Reading ASCII FITS File: " << items[it].name << std::endl;
                proc_error++;
                continue;
                }

            data=(float *) malloc((DIM_RAD*DIM_THT*2+2) * sizeof(float *));

            if ( data == NULL )
                {
                printf("ERROR: malloc() failed for data array, Skipping...\n");
                exit(1);
                }

            i=0;
            status=fscanf(fp_inp,"%f",&data[i]);
            do
                {
                i++;
                if (i >= (MAX_DIM * MAX_DIM))
                    {
                    std::cout << "ERROR: File Exceeded Maximum Size " << items[it].name << std::endl;
                    exit(1);
                    }
                } while((fscanf(fp_inp,"%f",&data[i])) != EOF);
        
            i--;

//
// Try to read the size from the first two bytes
//

            if(data[0]==data[1] && data[0]>0.0 && data[1]>0.0)
                {
                x_dim=data[0];
                y_dim=data[1];
                if (verbose) printf("--- dimensions (read) : xdim=%d : ydim=%d\n",x_dim,y_dim);
                }

//
//  If there were problems reading the size from the file, or the force read, 
//    calculate the size
//

            if ((x_dim == 0) || (y_dim == 0))
                {
                x_dim=sqrt(i-1);
                y_dim=sqrt(i-1);
                if (verbose) printf("--- dimensions (not read) : xdim=%d : ydim=%d\n",x_dim,y_dim);
                }

            items[it].radius=(x_dim-1)/2;
            items[it].valid=1;
            offset=2;
            }

//
// Copy the FITS data into the mat 2D Cartesian array
//

#ifdef DEBUG_DAT
//        for(i=0;i<(DIM_RAD*DIM_THT*2+2); i++)
        for(i=0;i<msize; i++)
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

        if (verbose) std::cout << "Processing Entry - Name: " << items[it].name << " Result: " << items[it].result << " Keyword: " << items[it].keyword << " Radius: " << items[it].radius << " Binary: " << items[it].binary << " Valid: " << items[it].valid << std::endl;

        if (verbose) puts("--- transforming X x Y -> Theta x ln r");

//
// Use (dim-1)/2 for each dimension.  This makes it work for both odd and even
//   sized images.
//

        x_0=((x_dim-1)/2)+1;
        y_0=((y_dim-1)/2)+1;

//
// Determine the masking value by determining the core brightness
//

        ctr_val=mat[x_0][y_0];
        if (mask_line)
            {
            if (verbose) printf("Center Value %f\n",ctr_val);
            log_bar=find_bar(items[it].radius,x_0,y_0,ctr_val);
            printf("Bar is %f\n",expf(log_bar));
            }
        else
            {
            log_bar=0.0;
            }

//
// log() functions are computationally expensive, so calculate the logs
//   outside of the loop.
//

        log_itrad=log((double)items[it].radius);

//
// Create the directory for the FFT output data
//

        base=basename((char *)items[it].result.c_str());
        sprintf(cmd,"mkdir -p %s\n",base);
        status=system(cmd);

//
//  This is the parallel version of the code.  All the inner radius values for
//    each annuli will be caculated in groups of parallel threads starting here.
//    Even if there is is only one thread, this code will still work.
//

#pragma omp parallel for

        for (int radius = 1; radius < items[it].radius; radius++)
            {
//
// VERY IMPORTANT - current is unique to each thread, so it must be defined 
//   inside this block of code.  Since these threads run in parallel, any
//   variable must be unique per thread.  We do this by making them arrays
//   and using the thread number as the index.  This is current and it will
//   appear often below, except for any variable defined here (which by virtue
//   of it's location will be unique per thread.
//

int     current=omp_get_thread_num(); /* Current index for arrays for thread */

//
// Other definitions that are unique instances per thread.
//

int     a, b;              /* Cartesian coordinates of ln(r)/theta in image  */
int    	mode;              /* Mode index value                               */
int     jm, im;            /* Local index variables                          */
int     cont_p;            /* Index for remapping output data in fft_data    */
int     status;            /* Pitch_class return value                       */
int     sum_ptr;           /* Index for FFT summed data strcuture            */
int     counter=0;         /* FFT array index value                          */
int     count_theta=1;     /* Counter for theta steps in degrees             */
int     count_radians;     /* Counter for theta steos in radians             */

char    outfile1[80];      /* Intermediate .rip file name string             */
char    outfile2[80];      /* Intermediate .dat file name string             */

FILE    *fp_out1;          /* Intermediate .rip file pointer                 */
FILE    *fp_out2;          /* Intermediate .dat file pointer                 */

float   lnr;               /* Natural log of radius for a certain point      */
float   x, y;              /* Relative cartesian coordinates of ln(r)/theta  */
float   log_lo;            /* Natural log of inside of fixed annuli          */
float   log_hi;            /* Natural log of outside of fixed annuli         */
float   log_rad;           /* Natural log of current value of radius         */
float   norma=0.0;         /* Normalization value (sum of number of values)  */
float   freq_save;         /* Current frequency calculation value            */
float   theta_degrees;     /* Current theta (polar angle) in degrees         */
float   theta_radians;     /* Current theta (polar angle) in radians         */


            if (reverse)
                {
                log_rad=log((double)(items[it].radius-radius+1));
                }
            else
                {
                log_rad=log((double)radius);
                }

            if (fixed && ((radius <= (fixed/2)) || (radius >= items[it].radius-(fixed/2)))) continue; 

//
// Zero out the arrays.  This is really important to getting the correct
//   results.
//

            for (im=0; im < DIM_RAD*DIM_THT+1; im++)
                {
                in_data[current][im][0]=0.0;
                in_data[current][im][1]=0.0;
                out_data[current][im][0]=0.0;
                out_data[current][im][1]=0.0;
                }
        
//
// Step around theta angles (360 degrees in 0.35 steps)
//

            for(theta_degrees=0.0;count_theta<=DIM_THT;theta_degrees+=theta_step) 
                {
                count_theta++;

//
// Convert the degrees to radians
//

                theta_radians=theta_degrees*GR_RAD;	
                count_radians=1;

                for(lnr=0.0;count_radians<=DIM_RAD;lnr+=radstep) 
                    {
                    count_radians++;

                    if ((zero) && (count_theta < 4 || count_theta > 1021))
                        {
                        in_data[current][counter][0]=0.0;
                        in_data[current][counter++][1]=0.0;
                        continue;
                        }

                    if ((mask_line) && (lnr <= log_bar))
                        {
                        in_data[current][counter][0]=0.0;
                        in_data[current][counter++][1]=0.0;
                        continue;
                        }
                       
//
// Here's the bit that controls what get mapped and what is set to zero.
//   These tests depend on the value of reverse and fixed.
//
    
                    if (reverse && (lnr>log_rad || lnr>log_itrad))
                        {
                        in_data[current][counter][0]=0.0;
                        in_data[current][counter++][1]=0.0;
                        continue;
                        }
                       
                    if (fixed && (lnr>log_hi || lnr<log_lo)) 
                        {
                        in_data[current][counter][0]=0.0;
                        in_data[current][counter++][1]=0.0;
                        continue;
                        }
    
                    if (!reverse && !fixed && (lnr>log_itrad || lnr<log_rad))
                        {

                        in_data[current][counter][0]=0.0;
                        in_data[current][counter++][1]=0.0;
                        continue;
                        }

                    x=expf(lnr)*cosf(theta_radians);
                    y=expf(lnr)*sinf(theta_radians);

                    a=(int)x+x_0;
                    b=(int)y+y_0;

                    if ((mask) && (mat[a][b] >= ctr_val))
                        {
                        in_data[current][counter][0]=0.0;
                        }
                    else
                        {
                        in_data[current][counter][0]=(double) mat[a][b];
                        norma+=in_data[current][counter][0];
                        }
                    in_data[current][counter++][1]=0.0;
                    }
                }

            if (verbose) printf("--- calculating 2DFFT: %d/%d\n",radius, items[it].radius);

#ifdef DEBUG_DAT
            if (radius<5)
                {
                printf("RADIUS: %d\n",radius);
                for(im=0;im<=counter;im++) 
                    {
                    printf("DEBUG: In Data[%d][0]=%f\n",im,in_data[current][im][0]);
                    printf("DEBUG: In Data[%d][1]=%f\n",im,in_data[current][im][1]);
                    }
                }
#endif

//
// Save the polar mapped image if the -p option was specified
//

            if ((polar) && (radius==1))
                {
                counter=0;
                for (jm=0; jm < DIM_RAD; jm++)
                    {
                    for (im=0; im < DIM_THT; im++)
                        {
                        proj[counter++]=(float) in_data[current][(im*2048)+jm+1][0];
                        }
                    }
                fname=(char *) items[it].name.c_str();
                if (verbose) printf("  --- Write P_%s File\n",fname);

                sprintf(pfile,"!P_%s",fname);
                if (ast.fits_write(pfile, proj, DIM_THT, DIM_RAD, 1, "p2dfft/",VERSION))
                    {
                    printf("WARNING: fits_write(%s) Failed\n",pfile);
                    }
                }

//
// Perform the FFT using the plan
//

            fftw_execute_dft(plan,in_data[current],out_data[current]);

//
// Normalize the output data
//

            for(im=0;im<=counter;im++) 
                {
#ifdef DEBUG_DAT
                printf("DEBUG: Out Data[%d][0]=%f\n",im,out_data[current][im][0]);
                printf("DEBUG: Out Data[%d][1]=%f\n",im,out_data[current][im][1]);
#endif
                out_data[current][im][0]=out_data[current][im][0]/(double)norma;
                out_data[current][im][1]=out_data[current][im][1]/(double)norma;
                }

//
// Loop for each mode
//

            for(mode=M_INI;mode<=M_FIN;mode++) 
                {
                counter=mode*DIM_RAD;

//
// If data files are being generated, open them and write the initial data
//

                base=basename((char *)items[it].result.c_str());
                sprintf(outfile1,"%s/%s%d_m%1d.rip",base,items[it].keyword.c_str(),radius,mode);
                sprintf(outfile2,"%s/%s%d_m%1d.dat",base,items[it].keyword.c_str(),radius,mode);
                if ((fp_out1=fopen(outfile1,"w"))==NULL)
                    {
                    if (warn) printf("WARNING: Could Not Write %s\n",outfile1);
                    }
                else
                    {
                    fprintf(fp_out1,"%d\n",x_dim/2);
                    fprintf(fp_out1,"%e\n",norma);
                    }

                if ((fp_out2=fopen(outfile2,"w"))==NULL)
                    {
                    if (warn) printf("WARNING: Could Not Write %s\n",outfile2);
                    }

//
// Extract the FFT output components for -50 to +50 Hz and populate them in
//   the fft_data structure.  P2DFFT uses a different order than FFTW uses
//   for it's output.   The mapping is:
//
//     Data Array   Description               fft_data Index   rip File
//     ----------   ---------------------     --------------   --------
//        0,0       Real Mid Freq (0/DC)           1025          403
//        0,1       Imag Mid Freq                  1025          404
//        1,0       Real Min Pos. Freq Start       1026          405
//        1,1       Imag Min Pos. Freq Start       1026          406
//       200,0      Real Pos. Freq                 1225          803
//       200,1      Imag Pos. Freq                 1225          804
//       201,0      Real Pos. Freq                 1226          N/S
//       201,1      Imag Pos. Freq                 1226          N/S
//      1024,0      Real Max Pos. Freq End         2049          N/S
//      1024,1      Imag Max Pos. Freq End         2049          N/S
//      1025,0      Real Min Neg. Freq                2          N/S
//      1025,1      Imag Max Neg. Freq                2          N/S
//      1847,0      Real Neg. Freq                  824          N/S
//      1847,1      Imag Neg. Freq                  824          N/S
//      1848,0      Real Neg. Freq                  825            3
//      1848,1      Imag Neg. Freq                  825            4
//      2047,0      Real Max Neg. Freq End         1024          401
//      2047,1      Imag Max Neg. Freq End         1024          402
//
//   This mapping is also for m0 only, other modes would have the same fft_data
//   and rip values, but the data indices would have mode*2048 added to them. 
//
//   Finally, this table also assumes a frequency mapping of -50 to +50, if
//   thats different, then the start and end will be different.
//
//   Also note that we multiply the imaginary component by -1.0 because FFTW3
//   returns a sign reversed value compared to the previous algorithm.
//

                for(cont_p=0;cont_p<DIM_RAD/2;cont_p++) 
                    {
                    fft_data[current][cont_p+(DIM_RAD/2)+1].real=out_data[current][counter][0];
                    fft_data[current][cont_p+(DIM_RAD/2)+1].imag=-1.0*out_data[current][counter][1];
                    fft_data[current][cont_p+DIM_RAD/2+1].abs=sqrt(pow(out_data[current][counter][0],2.0)+pow(out_data[current][counter][1],2.0));

                    if (DEBUG && radius==1) printf("DEBUG: Map out_data[%d][1] to fft_data[current][%d].real/imag/abs\n",counter,cont_p+(DIM_RAD/2)+1);

                    ++counter;
                    }

                fft_data[current][DIM_RAD+1].real=out_data[current][counter][0];
                fft_data[current][DIM_RAD+1].imag=-1.0*out_data[current][counter][1];
                fft_data[current][DIM_RAD+1].abs=sqrt(pow(out_data[current][counter][0],2.0)+pow(out_data[current][counter][1],2.0));

//
// This was in the original code.  Not sure if it is still needed.
//

                fft_data[current][1].abs=sqrt(pow(out_data[current][counter][0],2.0)+pow(out_data[current][counter][1],2.0));

                if (DEBUG && radius==1) printf("DEBUG: Map out_data[%d][1] to fft_data[%d].real/imag/abs\n",counter,cont_p+(DIM_RAD/2)+1);

                ++counter;

                for(cont_p=(-1)*(DIM_RAD/2)+1;cont_p<=-1;cont_p++) 
                    {
                    fft_data[current][cont_p+(DIM_RAD/2)+1].real=out_data[current][counter][0];
                    fft_data[current][cont_p+(DIM_RAD/2)+1].imag=-1.0*out_data[current][counter][1];
                    fft_data[current][cont_p+(DIM_RAD/2)+1].abs=sqrt(pow(out_data[current][counter][0],2.0)+pow(out_data[current][counter][1],2.0));

                    if (DEBUG && radius==1) printf("DEBUG: Map out_data[%d][1] to fft_data[%d].real/imag/abs\n",counter,cont_p+(DIM_RAD/2)+1);

                    ++counter;
                    }

//
// Add frequency values to the fft_data array, the summed data array, and
//   optionally the intermediate FFT output files
//

                sum_ptr=0;
		dindex=2;
                for(jm=1;jm<=DIM_RAD+1;jm++) 
                    {
                    freq_save=(-1)*STEP_P*DIM_RAD/2+(jm-1)*STEP_P;
                    if(freq_save>=FREQ_START && freq_save<=FREQ_END) 
                        {
                        if (fft_data[current][jm].abs == fft_data[current][jm].abs)
                            {
#pragma omp critical
                            {
                            fft_sum[mode][sum_ptr].abs+=fft_data[current][jm].abs;
                            }
                            }
                        sum_ptr++;
                        fft_data[current][jm].freq=freq_save;
                	    if (DEBUG && radius==1) printf("DEBUG: Map fft_data[%d][%d] to RIP Index=%d\n",current,jm,dindex);
			dindex++;
                        if (high_pass && (freq_save < ((float)mode*0.25)) && (freq_save > ((float)mode*-0.25)))
                            {
                            fft_data[current][jm].abs=0.0;
                            fft_data[current][jm].real=0.0;
                            fft_data[current][jm].imag=0.0;
                            }
#pragma omp critical
                            {
                            fprintf(fp_out2,"%f %e\n",freq_save,fft_data[current][jm].abs);
                            fprintf(fp_out1,"%e\n",fft_data[current][jm].real);
                            fprintf(fp_out1,"%e\n",fft_data[current][jm].imag);
                            }
                        }
                    }

//
// This set of output files are complete, so close them
//

                fclose(fp_out1);
                fclose(fp_out2);

//
// Call the pitch class function to determine the dominant pitch angle for this
//   radius.
//

                status=pit.pitch_phase(fft_data[current],mode,&mode_data[mode][radius]);

//
// The pitch_phase routine should have populated the mode_data structure with
//   all the data for this radius.  However, if the FFT returned NaN's (due
//   to low signal or monochromatic space), need to handle that before calling
//   the higher order analysis functions.   NOTE:  NaN's are not an error
//   (this can happen if the image is small in the frame).
//

                if ((status == PITCH_RET_ERR) || (status == PITCH_RET_NAN))
                    {
//
//  This is some error in the pitch angle calculation.  Should not happen with
//    rational input data.  However, we can calculate other meaningful
//    parameters.
//

                    if (warn) printf("WARNING: pitch_phase() failed (%d) for radius %d and mode %d\n",pit.get_err(),radius,cont_p);
                    mode_data[mode][radius].index=0;
                    mode_data[mode][radius].freq=NAN;
                    mode_data[mode][radius].amp=NAN;
                    mode_data[mode][radius].avg_amp=NAN;
                    mode_data[mode][radius].pa=NAN;
                    mode_data[mode][radius].phase=NAN;
                    mode_data[mode][radius].snr=NAN;
                    mode_data[mode][radius].fwhm=NAN;
                    }
                else
                    {
                    status=pit.snr(fft_data[current],&mode_data[mode][radius]);
                    if (status==PITCH_RET_ERR)
                        {
                        if (warn) printf("WARNING: snr() failed (%d) for radius %d and mode %d\n",pit.get_err(),radius,cont_p);
                        mode_data[mode][radius].avg_amp=NAN;
                        mode_data[mode][radius].snr=NAN;
                        mode_data[mode][radius].fwhm=NAN;
                        }
                    else
                        {
                        status=pit.fwhm(fft_data[current],&mode_data[mode][radius]);
                        if (status==PITCH_RET_ERR)
                            {
                            if (warn) printf("WARNING: fwhm() failed (%d) for radius %d and mode %d\n",pit.get_err(),radius,cont_p);
                            mode_data[mode][radius].fwhm=NAN;
                            }
                        }
                    }
                if (DEBUG) printf("DEBUG: Pitch Phase Angle=%f, SNR=%f, FWHM=%f\n",mode_data[mode][radius].pa,mode_data[mode][radius].snr,mode_data[mode][radius].fwhm);
                }
            }

// **** END OF PARALLEL THREAD FOR LOOP

//
// Now that all radii are complete, write the per mode and summed output files
//

        for (i = M_INI; i <= M_FIN; i++)
            {
            sprintf(outfile,"%s_m%1d",items[it].result.c_str(),i);
            if ((mode_out=fopen(outfile,"w"))==NULL)
                {
                printf("ERROR: Could Not Write %s\n",outfile);
                exit(1);
                }

            for (j = 1; j <= items[it].radius; j++)
                {
                sprintf(tmpofile,"%s%d_m%1d",items[it].keyword.c_str(),j,i);
                fprintf(mode_out,"%6d%11s%8.2f%12.3f%9.2f%11.3f%11.3f%11.3f\n",i,tmpofile,mode_data[i][j].freq,mode_data[i][j].amp,mode_data[i][j].pa,mode_data[i][j].phase,mode_data[i][j].snr,mode_data[i][j].fwhm);
                }
            fclose(mode_out);

            sprintf(outfile,"%s_sum_m%1d",items[it].result.c_str(),i);
            if ((sum_out=fopen(outfile,"w"))==NULL)
                {
                printf("ERROR: Could Not Write %s\n",outfile);
                exit(1);
                }

            for (j = 0; j < lim; j++)
                {
                fprintf(sum_out,"%6.2f     %f\n",fft_sum[i][j].freq,fft_sum[i][j].abs);
                }
            fclose(sum_out);
            }
        }
    printf("-------------------------------\n");
    it=(unsigned int)items.size()-(unsigned int)proc_error;
    printf("Successfuly Processed        %d\n",it);
    printf("Errors                       %u\n",proc_error);
    }
