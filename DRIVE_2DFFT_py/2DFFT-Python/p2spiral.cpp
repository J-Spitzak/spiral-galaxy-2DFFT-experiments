//
// P2SPIRAL.CPP - Generate Spiral Galaxy Models for Testing
//
//            This program will generate FITS files of model galaxies.  The
//            files can be either binary floating point images of ASCII FITS
//            files.  Multiple features of the image can be specified through
//            the options (see below).
//
//
// Version 4.1: 13-Dec-2018
//
//
// Authors:  Ian Hewitt & Dr. Patrick Treuthardt,
//           NC Museum of Natural Sciences,
//           Astronomy & Astrophysics Lab,
//           Raleigh, NC USA.
//           http://github.com/treuthardt/P2DFFT
//
//
// LICENSE
//
// P2DFFT Spiral Galaxy Arm Pitch Angle Analysis Suite
// Copyright (c) 2016-2019  Ian B. Hewitt & Dr. Patrick Treuthardt
//
// The program is free software:  you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, version 3.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see < https://www.gnu.org/licenses >.
//
// The authors can be contacted at:
//
//     North Carolina Museum of Natural Sciences
//     Astronomy & Astrophysics Laboratory
//     11 West Jones Street
//     Raleigh, NC, 27601  USA
//     +1.919.707.9800
//
//     -- or --
//
//     patrick.treuthardt@naturalsciences.org
//
//
// Algorithms:
//
//            The basic formula for the spirals is:
//
//
//                                        <theta> tan( <phi> )
//                 Polar form:  r = r0 * e
//
//                      Where r is the radius and thera is the central angle
//                            r0 is the initial radius at theta 0/-90
//                            phi is the pitch angle (can be variable)
//
//            the luminosity chage (decrease/increase) can be either linear or
//            logarithmic.  The logarithmic change follows the form of 
//            exponential growth:
//
//                            -kr
//                  f(r) = C e
//
//            Where r is the radius, and both C and k are constants.
//
//            ---------------------------------------------------------------
//
//            The cores are drawn as circles with formula:
//
//                   2    2    2
//                  x  + y  = r
//
//            Where x and y are the coordinates (0,0 origin) and r is the
//            radius.
//
//            ---------------------------------------------------------------
//
//            The bars are generated using the ellipse formula:
//
//            __                       __2  __                       __2
//            |(x*cos(phi))+(y*sin(phi))|   |(y*cos(phi))-(x*sin(phi))|
//            |-------------------------| + |-------------------------| = 1.0
//            |            a            |   |            b            |
//
//            Where x and y are the coordinates (0,0 origin), a is the semi-
//            major aixs, b is the semi-minor axis and phi is the angle of
//            rotation.
//
//
// References: 
//           "Measurement of Galactic Logarithmic Spiral Arm Pitch Angle Using
//            Two-Dimensional Fast Fourier Transform Decomposition."
//
//           Davis, Benjamin L., Berrier, Joel C., Shields, Douglas W., 
//           Kennefick, Julia, Kennefick, Daniel, Seigar, Marc S., Lacy, 
//           Claud H.S., & Puerari, Ivanio, 2012, ApJS, 199, 33
//
//           - or -
//
//           http://arxiv.org/abs/1202.4870
//
//
// Usage: p2spiral [-i|--input <file>] [-v|--verbose] [-t|--text] [-p|--print]
// 
//        If p2spiral is called with no -i argument, the user will be
//        prompted for all the parameters (see below) pitch angle,     
//        number of arms, and total file size (which must be odd).  Up to 64 
//        files can be specified (this can be changed by altering MAX_FILES).
//        Entering <ctrl>-D will stop the input and start the generation of
//        the files.  For each entry, up to two files will be generated: 
//        <base name>.fits will be a 32-bit floating point FITS file and 
//        (optionally) <base name>.txt will be an ASCII text FITS image file
//        suitable for input into 2DFFT.
//
//        Optionally, an input file can be specified which lists the parameters
//        for the files to be created (see below for data file description).  
//        The complete list of options is:
//
//              -i|--input : Will read file names and keywords from the file
//                           specified with this option instead of standard 
//                           input.
//              -v|--verbose : Causes status messages to be printed.
//              -t|--text  : Generate ASCII FITS files as well as binary
//              -p|--print : Print a listing of pitch angle by radius to stdout
//
//
// INPUT FILE FORMAT:
//
//   The input is in the form of a comma, space or tab delimited text file with
//   the following input fields, all on a single line (delimiters not shown for
//   clarity):
//
//   <name> <pa> <arms> <hsize> <vsize> <feath> <sweep> <rot> <r0> <core>
//        <bara> <barb> <mar> <fg> <bg> <delta> <lum> <log> <arm_lum> <noise>
//
//      <name> -    The prefix of the file name to be generated (not including
//                  the .fits suffix)
//      <pa> -      Pitch angle of the arms.  Can be a positive or negative
//                  floating point number with -75.0 < pa < 75.0
//      <arms> -    The integer number of spiral arms in the image (from 1
//                  to 6)
//      <hsize> -   The integer height of the image
//      <vsize> -   The integer height of the image
//      <feath> -   Value from 0 to 10 setting the added thickness of arm
//                  points in pixels.
//      <sweep> -   Value from 0.0 to 720.0 which indicates how far in polar 
//                  coordinates to map the arms in degrees, e.g. 720.0 is
//                  two complete circles.
//      <rot> -     Rotation of the arm start (-90 to +90).  By default (0.0)
//                  the first arm will start at noon and the others even 
//                  spaced around that.
//      <r0> -      Initial integer radius of the core/start of arms from 0
//                  to <size>-1
//      <core> -    Value of 2 to fill in the central core with values twice
//                  the default arm value.  A value of 1 fills core the
//                  starting value of arms.   A value of 0 means the core i
//                  not filled in at all.
//      <bara> -    Semi-major aixs of bar (0=no bar), must be an integer. 
//      <barb> -    Semi-minor aixs of bar (integer). Ignored if bara=0.
//      <mar> -     Outside margin in pixels.  Value must be an integer. 
//                  Nothing will be mapped beyond this, no matter what 
//                  other settings are used.
//      <fg> -      Foreground value for arms
//      <bg> -      Background value in the image.  This is equivalent to the
//                  bias, and any noise will be added to this value.
//      <delta> -   Change in pitch angle (negative or positive) in degrees of
//                  the arms over their length.  This value is in degrees
//                  from -60.0 and 60.0.
//      <lum> -     Change in luminosity (negative or positive) of an arm over
//                  it's length.  This is a percentage value between -0.99
//                  and 0.99
//      <log> -     Change in luminosity iis logaritmic (1) or linear (0)
//      <arm_lum> - Apply change in luminosity over the width of an arm as well
//                  as it's length. A value of 0 means constant width
//                  luminosity and 1 means decreasing luminosity.
//      <noise> -   This is a floating point value controlling the amount of
//                  Poisson (shot) noise in the image.  This value sets the 
//                  largest possible value of the noise in any pixel and the
//                  actual noise values used in any given pixel will be
//                  random up to this amount. Please note that noise will not
//                  be added to the core or arms and those values should be
//                  > that the noise value.
//
//      A an example is:
//
//      Sample,25.0,3,935,935,3,180.0,0,25,2,0,0,20,128,10,0,0,1,0,50
//
//
// Revision History:
//      4.1  13-Dec-2018: - Fix bug in feathering code to make more consistent
//                          arm widths
//      4.0  10-Jun-2018: - Add parameter to add/specify a bar
//                        - Change to C++ (only minor changes needed, so 
//                          keep revision numbering scheme)
//                        - Change FITS writing code to use astro_class
//                        - Change array declarations to match FITS standard
//                        - Update comments/usage information
//                        - Fix bug with non-comma delimiters not working
//                        - Swich array allocation to use the one in 
//                          astro_class
//                        - Move MAX_FILES to globals.h
//                        - Allow different vertical and horizontal dimensions
//                        - Increase maximum dimensions
//      3.0  16-Mar-2018: - Correct bug in command line option strcuture
//                        - Add input field for type of brightness change, 
//                          arm width brightness change, and one to control
//                          core brightness value.
//                        - Add simulated Poisson noise to images
//                        - Signifcantly update comments
//                        - Make arm sweep part of input files rather than 
//                          command line option
//                        - Add support for comments and blank lines in input
//                          file to improve readability
//                        - Redorder parameters into logical order (group
//                          similar ones together)
//                        - Remove optional parameters, all parameters need
//                          to be specified for all images (this is clearer)
//                        - Add subroutines for input parameters to simplify
//                          the code
//                        - Increased most default values to allow more varied
//                          output files
//                        - Fix bug in command line option processing for -p
//                        - Fix bug in linear brightness change algorithm with
//                          increasing brightness
//                        - Fix bug in variable pitch angle calculation
//                        - Fix bug where specifying too many input files would
//                          cause unpredictable behavior instead of giving error
//      2.1  21-Jan-2018: - Add option to specify change in luminosity with 
//                          radius
//                        - Add command line option for arm sweep
//      2.0  04-Jan-2018: - Add input option for rotation of spiral (-90 to 90
//                          degrees) for each file
//                        - Add input option to specify initial radius (r0) for
//                          each file
//                        - Add input option specify if the core area should be
//                          filled in with a value twice the value of the arm
//                          foreground pixel value
//                        - Add input option to specify feather value
//                        - Add input options to specify default foreground and
//                          background pixel values
//                        - Add input option to specify the change in pitch
//                          angle from beginning to end
//                        - Add input option to specify outer margin (i.e. the
//                          number of pixels at the edge where nothing is
//                          plotted)
//                        - Remove restriction on file size being odd
//                        - Change input scheme to include defaults and not
//                          fail, but re-ask if invalid value is given
//                        - Change theta value from a function of size to 360
//                          degrees
//                        - Make ASCII FITS file generation optional by
//                          specifying a command line flag
//                        - Add option to list pitch angle by radius
//                        - Fix FITS generation code to use CFITSIO function
//                          built-in code to remove existing file
//                        - Update/correct comments
//      1.4  28-Aug-2017: - Fix bug in file close statement
//                        - Fix four small errors that led to compiler warnings
//                          on some Linux distributions (like Ubuntu)
//      1.3  18-Jan-2017: - Added the optional parameters to control values to
//                          adjust signal to noise ratio in the images
//                        - Changed delimiter in input file from tab to comma
//                        - Fixed error in comments/help test
//      1.2  13-Dec-2016: - Fixed bug where chirality was reversed
//      1.1  07-Dec-2016: - Increased maximum number of files from 64 to 256
//      1.0  06-Dec-2016:   Initial version
//
//

#define     VERSION "4.1/20181213"

//
// HEADER FILES
//

#include    <math.h>
#include    <stdio.h>
#include    <stdlib.h>
#include    <string.h>
#include    <getopt.h>
#include    <sys/stat.h>
#include    <sys/types.h>

#include    "globals.h"
#include    "astro_class.h"

//
// CONSTANTS - These are the default values for the parameters
//

#define     DEF_PA    20.0 /* Default pitch angle (degrees)                  */
#define     DEF_ARMS     2 /* Default number of arms                         */
#define     DEF_SIZE   255 /* Default size of image                          */
#define     DEF_FTHR     5 /* Default fethering (fuzziness) value            */
#define     DEF_SWEEP  360.0 /* Degrees used for polar mapping of arms       */
#define     DEF_ROT    0.0 /* Rotation (in degress) of 0 degree point        */
#define     DEF_R0    20.0 /* Initial Radius in Pixels at 0 degrees          */
#define     DEF_CORE     1 /* Default value for filling in the core          */
#define     DEF_BARA   0.0 /* Default value for bar semi-major axis          */
#define     DEF_BARB   0.0 /* Default value for bar semi-minor axis          */
#define     DEF_MAR     20 /* Default outer margin value                     */
#define     DEF_FG   255.0 /* Default value used for spiral components       */
#define     DEF_BG     0.0 /* Default background (bias) value                */
#define     DEF_DELTA  0.0 /* Change in Pitch Angle                          */
#define     DEF_LUM    0.0 /* Change in luminosity                           */
#define     DEF_LOG    1.0 /* Change in luminosity algorithm                 */
#define     DEF_ARM_LUM 1.0 /* Change in luminosity across arm width         */
#define     DEF_NOISE  0.0 /* Max value of Poisson (Shot) Noise added to bg  */

//
// CONSTANTS - Can change these to alter the program limits
//

#undef      DEBUG

#define     STR_SIZE    64 /* Maximum characters in input line               */
#define     MAX_FILES 1024 /* The limit on the number of files specified     */

#define     MIN_PA   -75.0 /* Minimum pitch angle (floating point value)     */
#define     MAX_PA    75.0 /* Maximum pitch angle (floating point value)     */
#define     MIN_ARM      1 /* Minimum number of spiral arms                  */
#define     MAX_ARM      6 /* Maximum number of spiral arms                  */
#define     MIN_SIZE    50 /* Minimum file size (should be at least 50)      */
#define     MAX_SIZE  2048 /* Maximum files size ( < 2049 for 2DFFT)         */
#define     MIN_FTHR     0 /* Minimum feather setting (integer value)        */
#define     MAX_FTHR    15 /* Maximum feather setting (integer value)        */
#define     MIN_SWEEP 90.0 /* Minimum arm sweep (< 45 is not meaningful)     */
#define     MAX_SWEEP 720.0 /* Maximum arm mapping angle                     */
#define     MIN_ROT  -90.0 /* Minimum rotation angle (floating point value)  */
#define     MAX_ROT   90.0 /* Maximum rotation angle (floating point value)  */
#define     MIN_R0     1.0 /* Minimum initial radius (floating point value)  */
#define     MAX_R0  1000.0 /* Maximum initial radius (floating point value)  */
#define     MIN_CORE     0 /* Lowest core mapping value                      */
#define     MAX_CORE     2 /* Highest core mapping value                     */
#define     MIN_BARA     0 /* Minimum bar semi-major axis                    */
#define     MAX_BARA 1024.0 /* Maximum bar semi-major axis                   */
#define     MIN_BARB     0 /* Minimum bar semi-minor axis                    */
#define     MAX_BARB 1024.0 /* Maximum bar semi-minor axis                   */
#define     MIN_MAR      0 /* Minimum outer margin of blank space            */
#define     MAX_MAR    200 /* Maximum outer margin of blank space            */
#define     MIN_PIXEL -1024.0 /* The minimum value for fg and bg variables   */
#define     MAX_PIXEL  1024.0 /* The maximum value for fg and bg variables   */
#define     MIN_DELTA -60.0 /* Minimum pitch angle change (floating point)   */
#define     MAX_DELTA  60.0 /* Maximum pitch angle change (floating point)   */
#define     MIN_LUM   -0.99 /* Minimum brightness change (floating point)    */
#define     MAX_LUM    0.99 /* Maximum brightness change (floating point)    */
#define     MIN_LOG      0 /* Lowest brightness option                       */
#define     MAX_LOG      1 /* Highest brightness option                      */
#define     MIN_ARM_LUM  0 /* Lowest arm width luminosity change value       */
#define     MAX_ARM_LUM  1 /* Highest arm width luminosity change value      */
#define     MIN_NOISE -512.0 /* The minimum value for fg and bg variables    */
#define     MAX_NOISE  512.0 /* The maximum value for fg and bg variables    */

//
// VARIABLES
//

int     c;                 /* Getopt_long return value                       */
int     r2;                /* Initial core radius squared                    */
int     i, j;              /* Index variables                                */
int     s, t;              /* Index variables                                */
int     x, y;              /* Array index variables                          */
int     ctr;               /* Counter for number of entries/line in txt file */
int     mode;              /* Index used in the formula to process each arm  */
int     txt=0;             /* Flag for creating ASCII FITS files             */
int     outer;             /* Estimate of arm length                         */
int     longr;             /* Estimate of longest radius                     */
int     list=0;            /* Flag for listing the pitch angles by radius    */
int     errcnt;            /* Total number of input errors encountered       */
int     starti;            /* Starting radius for arms (integer)             */
int     status;            /* Result from system() command                   */
int     centerx;           /* Center coordinate for x (horizontal) axis      */
int     centery;           /* Center coordinate for x (vertical) axis        */
int     counter;           /* Main processing loop variable                  */
int     verbose;           /* Flag for verbose mode (1=true)                 */
int     starti_s;          /* Starting radius for semi-minor axis            */
int     num_files;         /* Total number of files to process - 1           */
int     mar[MAX_FILES];    /* Outer margin value                             */
int     arm[MAX_FILES];    /* Number of arms for each file                   */
int     hsize[MAX_FILES];  /* Horizontal size in pixels for each file        */
int     vsize[MAX_FILES];  /* Vertical size in pixels for each file          */
int     core[MAX_FILES];   /* Flag for filling in the core area              */
int     feath[MAX_FILES];  /* Feathering/fuzziness value for each file       */
int     linear[MAX_FILES]; /* Brightness change algorithm flag               */
int     arm_lum[MAX_FILES];  /* Flag for changing brightness over arm width  */

char    *item;             /* Token parsed from input line                   */
char    cmd[STR_SIZE];     /* String for command to remove old file versions */
char    key[256];          /* String for FITS keywords                       */
char    line[256];         /* String for reading input file lines            */
char    keys[5][32];       /* FITS header key names                          */
char    items[5][80];      /* FITS header key values                         */
char    fname[STR_SIZE];   /* Input file name string                         */
char    fname2[STR_SIZE];  /* Input file name string                         */
char    entry[STR_SIZE];   /* Input file name string                         */

char    base[MAX_FILES][64]; /* Base names for generated files               */

float   r;                 /* R value for the polar coordinates              */
float   rn;                /* calculated noise value                         */
float   mod;               /* Modifier for the chirality (direction)         */
float   brt;               /* Current pixel value for foreground             */
float   dist;              /* Cartesian distance from center                 */
float   mb,ma;             /* Bar (ellipse) mapping coordinates              */
float   si,co;             /* Sine and cosine of ellipse mapping coordinates */
float   theta;             /* Loop variable for theta angles                 */
float   pitch;             /* Loop variable for changing pitch angle         */
float   **mat;             /* Pointer for the Cartesion mapping of image     */
float   change;            /* Rate of change for varying pitch angles        */
float   startf;            /* Starting radius for arms (float)               */
float   lum_rate;          /* Luminosity change per radius step              */
float   newpitch;          /* Updated pitch value                            */
float   avg_pitch;         /* Average pitch value                            */
float   max_pitch;         /* Maximum pitch value                            */
float   min_pitch;         /* Minimum pitch value                            */
float   num_pitch;         /* Number of pitch values used                    */
float   separation;        /* Angular separation arms (polar coordinates)    */
        
float   r0[MAX_FILES];     /* Initial radius at 0 degrees for each file      */
float   fg[MAX_FILES];     /* Foreground FITS pixel value                    */
float   bg[MAX_FILES];     /* Background FITS pixel value                    */
float   pa[MAX_FILES];     /* Pitch angles of files                          */
float   lum[MAX_FILES];    /* Brightness change over radius                  */
float   rot[MAX_FILES];    /* Rotation (degrees) value for each file         */
float   bara[MAX_FILES];   /* Pitch angle change value for each file         */
float   barb[MAX_FILES];   /* Pitch angle change value for each file         */
float   delta[MAX_FILES];  /* Pitch angle change value for each file         */
float   sweep[MAX_FILES];  /* Pitch angle range for mapping                  */
float   noise[MAX_FILES];  /* Maximum noise ceiling value                    */

long    naxis=2;           /* CFITSIO number of axes - always 2              */
long    naxes[2];          /* Size of array give to CFITSIO                  */

FILE    *ofile;            /* File stream for output .txt file               */
FILE    *file_list;        /* File stream for input file                     */

astro   ast;               /* Instantiation of astro_class                   */

//
// SUBROUTINES
//


// 
// READ_TOKEN() - Runs a loop to get the input value for a parameter from a
//                file line.  Parses the next token in the line and compares
//                it for validity.
//
// Prerequisites:  strtok() must have alread been run on the input line prior
//                 to first execution of this routine.
//
// Globals:
//      errcnt  - Count of lines with errors
//
// Arguments:
//      fname   - Character string with base name of the file
//      name    - Character string with name of the field
//      min     - Minimum value for any valid input
//      max     - Minimum value for any valid input
//
// Return Value:
//      ret     - Floating point value or -2048.0 if error
//
            
float   read_token(char *fname, const char *name,float min,float max)

    {

    float   ret;

//
// Get next token in string (line from file)
//

    if ((item=strtok(NULL,",\t "))!=NULL)
        {
        ret=(float)atof(item);

//
// Check it
//

        if ( (ret < min) || (ret > max) )
            {
            printf("WARNING: Invalid %s %f for File %s\n", name, ret, fname);
            errcnt++;
            ret= -2048.0;
            }
        }
    else
        {
        printf("ERROR: No %s for File %s\n",name, fname);
        errcnt++;
        ret= -2048.0;
        }

    return(ret);
    }

// 
// GET_INPUT() - Runs a loop to get the input value for a parameter from stdin.
//               Will keep looping until a <ctrl-d>, <cr>, or valud entry is
//               received and will return the appropriate code.  All inputs
//               and outputs are floating point.
//
// Arguments:
//      prompt  - Character string with prompt for user.  Should have one %f
//                  field to show the default value.
//      min     - Minimum value for any valid input
//      max     - Minimum value for any valid input
//      def     - Default value assumed if only <cr> is read
//
// Return Value:
//      ret     - Value to be used from input.  If <ctrl-d> is received, it
//                will be -2048.0
//

float   get_input(const char *prompt, float min, float max, float def)

    {
    float   ret;

//
// Loop until we get a valid entry
//

    while (1)
        {
        printf(prompt,def);

//
// If a <ctrl-d> set return to -2048.0
//

        if (fgets(entry,STR_SIZE,stdin) == NULL)
            {
            printf("\n");
            return(-2048.00);
            }

//
// If a <cr> set return to the default
//

        entry[strcspn(entry,"\n")]=0;
        if (strlen(entry)==0)
            {
            ret=def;
            }
        else
            {
            ret=atof(entry);
            }

//
// Check the values against min/max
//

        if ((ret < (float) min) || (ret > (float) max))
            {
            printf("WARNING: Bad Value %s\n",entry);
            continue;
            }
        else
            {
            break;
            }
        }
    return(ret);
    }


//
// MAIN ROUTINE
//

int main(int argc, char **argv)
    {

//
// Get and process the command line arguments.
//

//
// Define the command line options, see getopt_long(3) for details
//

    static struct option long_options[] =
        {
        {"verbose", no_argument,         0, 'v'},
        {"text", no_argument,            0, 't'},
        {"print", no_argument,           0, 'p'},
        /* These options require an argument. */
        {"input", optional_argument,     0, 'i'},
        {0, 0, 0, 0}
        };
      
    int option_index = 0;

    while ((c = getopt_long (argc, argv, "vtpi:", long_options, &option_index)) 
!= -1)
        {
        switch (c)
            {
            case 'v':
                {
                verbose = 1;
                break;
                }
            case 't':
                {
                txt = 1;
                break;
                }
            case 'p':
                {
                list = 1;
                break;
                }
            case 'i':
                {
                strcpy(fname,optarg);
                break;
                }
            default:
                {
                fprintf(stderr, "Usage: p2spiral [-i|--input <file>] [-v|--verbose] [-t|--text] [-p|--print]\n");
                exit(1);
                break;
                }
            }
        }

//
// Open the input file name if given
//

    if ((fname[0] != '\0')&&(fname[0] != '\1'))
        {
        file_list=fopen(fname,"r");
        if ( file_list == NULL )
            {
            printf("ERROR: Cannot open input file - %s\n",fname);
            exit(1);
            }

//
// Open the input file and read the values into the arrays.  This code seems a
//   bit tedious because of all the error checking.
//

        errcnt=0;
        num_files=0;
        while (fgets(line,256,file_list)!=NULL)
            {
//
// Check for comment or blank line
//

            if ((line[0]=='#') || (strlen(line) < 2)) continue;

//
// If we have reached the maximum file number, we need to flag an error.
//

            if (num_files == MAX_FILES)
                {
                printf("ERROR: Too many input lines!\n");
                exit(1);
                }

//
// Start breaking the line into individual fields.  Need to run strtok() at
//   least once before calling read_token().
//

            if ((item=strtok(line,", \t"))!=NULL)
                {
                strcpy(base[num_files],item);
                }
            else
                {
                printf("WARNING: Invalid Keyword on Line %d\n",num_files);
                errcnt++;
                continue;
                }

            if ((pa[num_files]=read_token(base[num_files],"Pitch Angle",MIN_PA,MAX_PA)) < -2000.0) continue;
            
            if ((arm[num_files]=(int)read_token(base[num_files],"Arm Number",MIN_ARM,MAX_ARM)) < -2000.0) continue;
            
            if ((hsize[num_files]=(int)read_token(base[num_files],"Horizontal File Size",MIN_SIZE,MAX_SIZE)) < -2000.0) continue;
            
            if ((vsize[num_files]=(int)read_token(base[num_files],"Vertical File Size",MIN_SIZE,MAX_SIZE)) < -2000.0) continue;
            
            if ((feath[num_files]=(int)read_token(base[num_files],"Feather",MIN_FTHR,MAX_FTHR)) < -2000.0) continue;
            
            if ((sweep[num_files]=read_token(base[num_files],"Sweep Angle",MIN_SWEEP,MAX_SWEEP)) < -2000.0) continue;
            
            if ((rot[num_files]=read_token(base[num_files],"Rotation Angle",MIN_ROT,MAX_ROT)) < -2000.0) continue;
            
            if ((r0[num_files]=read_token(base[num_files],"Initial Radius",MIN_R0,MAX_R0)) < -2000.0) continue;
           
            if ((core[num_files]=(int)read_token(base[num_files],"Core Setting",MIN_CORE,MAX_CORE)) < -2000.0) continue;
            
            if ((bara[num_files]=read_token(base[num_files],"Bar Semi-Major Axis",MIN_BARA,MAX_BARA)) < -2000.0) continue;
           
            if ((barb[num_files]=read_token(base[num_files],"Bar Semi-Minor Axis",MIN_BARB,MAX_BARB)) < -2000.0) continue;

            if (bara[counter] && (barb[num_files] < 1.0))
                {
                printf("WARNING: Semi-Minor Axis Must Be At Least 1.0...Ignoring\n");
                bara[num_files]=0.0;
                barb[num_files]=0.0;
                }
 
            if (barb[num_files] > bara[num_files])
                {
                printf("WARNING: Semi-Major Axis Must Be >= Than Semi-Minor Axis...Skipping\n");
                continue;
                }
 
            if (bara[num_files] && (r0[num_files] >= bara[num_files]))
                {
                printf("WARNING: Semi-Major Axis Must Be > Than Initial Radius...Ingoring Bar Values\n");
                bara[num_files]=0.0;
                barb[num_files]=0.0;
                }
 
            if ((mar[num_files]=(int)read_token(base[num_files],"Outer Margin",MIN_MAR,MAX_MAR)) < -2000.0) continue;
 
            if ((fg[num_files]=read_token(base[num_files],"Foreground",MIN_PIXEL,MAX_PIXEL)) < -2000.0) continue;
           
            if ((bg[num_files]=read_token(base[num_files],"Background (Bias)",MIN_PIXEL,MAX_PIXEL)) < -2000.0) continue;
           
            if ((delta[num_files]=read_token(base[num_files],"Pitch Angle Change",MIN_DELTA,MAX_DELTA)) < -2000.0) continue;
           
            if ((lum[num_files]=read_token(base[num_files],"Luminosity Change",MIN_LUM,MAX_LUM)) < -2000.0) continue;
           
            if ((linear[num_files]=(int)read_token(base[num_files],"Brightness Algorithm",MIN_LOG,MAX_LOG)) < -2000.0) continue;
            
            if ((arm_lum[num_files]=(int)read_token(base[num_files],"Arm Width Luminosity Change",MIN_ARM_LUM,MAX_ARM_LUM)) < -2000.0) continue;
            
            if ((noise[num_files]=read_token(base[num_files],"Noise (Shot)",MIN_NOISE,MAX_NOISE)) < -2000.0) continue;
           
//
// Success!  Carry on with next item.
//

            num_files++;
            continue;
            }
        }
    else
//
// No input file specified, so read the parameters interactively from stdin
//

        {
        num_files=0;

//
// Forever loop, will break out when user is done and enters <CTRL>-D
//   No default for the base name.
//

        while (1)
            {
            printf("\nBase File Name: ");

//
// Here's where we bail on the loop if <CTRL>-D
//

            if ( fgets(entry,STR_SIZE,stdin) == NULL ) break;
            entry[strcspn(entry,"\n")]=0; 
            if (strlen(entry)==0)
                {
                printf("WARNING: Invalid Keyword %s\n",base[num_files]);
                continue;
                }
            strcpy(base[num_files],entry);

//
// For the rest of the values, use the subroutine
//

            if ((pa[num_files]=get_input("Pitch Angle [%f]: ", MIN_PA, MAX_PA, DEF_PA)) < -2000.0) break;

            if ((arm[num_files]=(int)get_input("Arms [%f]: ", MIN_ARM, MAX_ARM, DEF_ARMS)) < -2000.0) break;

            if ((hsize[num_files]=(int)get_input("Horizontal Size [%f]: ", MIN_SIZE, MAX_SIZE, DEF_SIZE)) < -2000.0) break;

            if ((vsize[num_files]=(int)get_input("Vertical Size [%f]: ", MIN_SIZE, MAX_SIZE, DEF_SIZE)) < -2000.0) break;

            if ((feath[num_files]=(int)get_input("Feather [%f]: ", MIN_FTHR, MAX_FTHR, DEF_FTHR)) < -2000.0) break;

            if ((sweep[num_files]=get_input("Sweep Angle[%f]: ", MIN_SWEEP, MAX_SWEEP, DEF_SWEEP)) < -2000.0) break;

            if ((rot[num_files]=get_input("Rotation Angle[%f]: ", MIN_ROT, MAX_ROT, DEF_ROT)) < -2000.0) break;

            if ((r0[num_files]=get_input("Initial Radius [%f]: ", MIN_R0, MAX_R0, DEF_R0)) < -2000.0) break;

            if ((core[num_files]=(int)get_input("Core Setting [%f]: ", MIN_CORE, MAX_CORE, DEF_CORE)) < -2000.0) break;

            if ((bara[num_files]=get_input("Initial Radius [%f]: ", MIN_BARA, MAX_BARA, DEF_BARA)) < -2000.0) break;

            if (bara[num_files])
                {
                if ((barb[num_files]=get_input("Initial Radius [%f]: ", MIN_BARB, MAX_BARB, DEF_BARB)) < -2000.0) break;
                }
            else
                {
                if ((barb[num_files]=get_input("Initial Radius [%f]: ", MIN_BARB, MAX_BARB, DEF_BARB+1.0)) < -2000.0) break;
                }

            if ((mar[num_files]=(int)get_input("Outer Margin [%f]: ", MIN_MAR, MAX_MAR, DEF_MAR)) < -2000.0) break;

            if ((fg[num_files]=get_input("Foreground [%f]: ", MIN_PIXEL, MAX_PIXEL, DEF_FG)) < -2000.0) break;

            if ((bg[num_files]=get_input("Background (Bias) [%f]: ", MIN_PIXEL, MAX_PIXEL, DEF_BG)) < -2000.0) break;

            if ((delta[num_files]=get_input("Pitch Angle Change[%f]: ", MIN_DELTA, MAX_DELTA, DEF_DELTA)) < -2000.0) break;

            if ((lum[num_files]=get_input("Luminosity Change [%f]: ", MIN_LUM, MAX_LUM, DEF_LUM)) < -2000.0) break;

            if ((linear[num_files]=(int)get_input("Brightness Change Algorithm [%f]: ", MIN_LOG, MAX_LOG, DEF_LOG)) < -2000.0) break;

            if ((arm_lum[num_files]=(int)get_input("Arm Width Luminosity Change Setting [%f]: ", MIN_ARM_LUM, MAX_ARM_LUM, DEF_ARM_LUM)) < -2000.0) break;

            if ((noise[num_files]=get_input("Noise (Shot) [%f]: ", MIN_NOISE, MAX_NOISE, DEF_NOISE)) < -2000.0) break;

            num_files++;
            }
        }

//
//  Need at least one entry or there is no valid input
//

    if (num_files < 1)
        {
        printf("No files to generate (%d)\n", num_files);
        exit (1);
        }

//
// Now create the requested FITS files
//

    for ( counter=0; counter < num_files; counter++ )
        {

//
// If the verbose flag is set, print out all input information/assumptions
//

        if (verbose)
            {
            printf("Processing File %d: Name=%s, Pitch Angle=%f\n",counter+1,base[counter],pa[counter]);

            printf("    Arms=%d, Hor. Size=%d, Ver. Size=%d, Feather=%d\n",arm[counter],hsize[counter],vsize[counter],feath[counter]);

            printf("    Sweep=%f, Rotation=%f, r0=%f, Core=%d, Bar Semi-Major=%f, Bar Semi-Minor=%f\n",sweep[counter],rot[counter],r0[counter],core[counter],bara[counter],barb[counter]);

            printf("    Margin=%d, Fg=%f, Bg=%f, Delta=%f, Lum=%f\n",mar[counter],fg[counter],bg[counter],delta[counter],lum[counter]);

            printf("    Log=%d, Arm_lum=%d, Noise=%f\n",linear[counter],arm_lum[counter],noise[counter]);
            }

//
// Initialize array values and populate it with random noise values, if needed.
//   PLEASE NOTE: The FITS array ordering is different than C so columns are
//   major and rows are minor (fastest varying).   The astro_class libraries
//   expect this as well.
//

        if (verbose) printf("  --- Generating Arrays\n");

        mat=ast.ArrayAlloc(vsize[counter], hsize[counter]);

        for (x=0; x < hsize[counter]; x++)
            {
            for (y=0; y < vsize[counter]; y++)
                {
                if (noise[counter] != 0.0)
                    {
                    rn=(((float)rand())/((float) RAND_MAX + 1.0))*noise[counter];
                    mat[y][x]=bg[counter]+rn;
                    }
                else
                    {
                    mat[y][x]=bg[counter];
                    }
                }
            }

//
// Set the chirality (direction) value for the equation based on neg/pos p.a.
//

        if (verbose) printf("  --- Set Chirality\n");

        if ( pa[counter] > 0 )
            {
            mod=-1.0;
            }
        else
            {
            mod=1.0;
            }

//
// Determine the arm separation.  If arm > 1 set to the angle of separation
//   between the arms.  If arm=1, then set to zero, which will cause the
//   separation term in the formula to go to 0.
//

        if (verbose) printf("  --- Set Arm Separation\n");

        if (arm[counter] > 1) 
            {
            separation=360.0/(float)arm[counter];
            }
        else
            {
            separation=0.0;
            }

//
// Determine the radius where the arms will start based on core size and bar
//   parameters
//

        if (bara[counter] > r0[counter])
            {
            startf=bara[counter];
            starti=(int) bara[counter];
            }
        else
            {
            startf=r0[counter];
            starti=(int) r0[counter];
            }

//
// Calculate the pitch angle rate of change (if any)
//

        longr=1;
        for (theta=0.0;theta<=sweep[counter];theta+=1.0)
            {
            r = startf*expf(tan(fabs(pa[counter]+delta[counter])*(M_PI/180.0))*theta*(M_PI/180.0));
            x=(hsize[counter]/2)+(int)(r*cos(mod*theta*(M_PI/180.0)));
            y=(vsize[counter]/2)+(int)(r*sin(mod*theta*(M_PI/180.0)));

            if ( (x > mar[counter]) && (x < (hsize[counter]-mar[counter])) && (y > mar[counter]) && (y < (vsize[counter]-mar[counter])) )
                {
                longr=r;
                }
            }

        if (verbose) printf("Longest r=%d,  ",longr);

//
// Check for sanity of parameters
//

        if (hsize[counter] < vsize[counter])
            {
            outer=hsize[counter]/2-mar[counter]-starti-feath[counter]-1;
            }
        else
            {
            outer=vsize[counter]/2-mar[counter]-starti-feath[counter]-1;
            }
        if (verbose) printf("Outer Arm=%d,  ",outer);
        if ((outer < 2)||(outer>(hsize[counter]/2))||(outer>(vsize[counter]/2)))
            {
            printf("ERROR: Input parameters inconsistent - arm length is %d\n",outer);
            errcnt++;
            continue;
            }

//
// Calculate the pitch angle change
//

        change=delta[counter]/((float)longr-startf);
        if (verbose) printf("Pitch Angle Incremental Change=%f\n",change);

//
// Calculate the brightness change for the arms
//

        if (lum[counter] == 0.0)
            {
            lum_rate=0.0;
            }
        else
            {
            if (linear[counter]==0) lum_rate=(fg[counter]-fabs(fg[counter]*lum[counter]))/((float)longr-startf);
                else lum_rate=-1.0*logf(fg[counter]/(fg[counter]+(fg[counter]*lum[counter])))/((float)(longr-1)-startf);
            if ((lum[counter]< 0)&&(linear[counter]==0))
                {
                lum_rate= -1.0*lum_rate;
                }
            if (verbose) printf("Brightness Incremental Change=%f\n",lum_rate);
            }

//
// Start mapping the polar coordinates.  Loop through theta from 1 to sweep.
//

        if (verbose) printf("  --- Map Coordinates\n");

        pitch=pa[counter];
        min_pitch=pitch;
        max_pitch=pitch;
        avg_pitch=0;
        num_pitch=0;

        for (theta=0.0;theta<=sweep[counter];theta+=1.0)
            {
//
// This loop is needed for multiple arms
//

            for(mode=0; mode < arm[counter]; mode++)
                {

//
// Calculate the true r value.  This is where using consecutively larger
//   values of theta makes the calculation easier.  Please note tan(3) on 
//   Linux expects values in radians and not degrees.
//

                r = startf*expf(tan(fabs(pitch)*(M_PI/180.0))*theta*(M_PI/180.0));

//
// Now map r and theta to the Cartesian values.  Note that we map
//   from the center of the image outward.  Also note, that we vary the actual
//   theta angle for each arm by the separation.
//

                x=(hsize[counter]/2)+(int)(r*cos(mod*(theta+rot[counter]+((float)mode*separation))*(M_PI/180.0)));
                y=(vsize[counter]/2)+(int)(r*sin(mod*(theta+rot[counter]+((float)mode*separation))*(M_PI/180.0)));

//
// Depending on the pitch angle value, X and/or Y can go outside of the array
//   bounds, so don't plot them if that's the case.  Please note, since the arm
//   lines are feathered, we need to include the padding.
//

                if ((x >= (mar[counter]+feath[counter])) && (x < (hsize[counter]-mar[counter]-feath[counter])) && (y >= (mar[counter]+feath[counter])) && (y < (vsize[counter]-mar[counter]-feath[counter])))
                    {
                    if (linear[counter]==0)
                        {
                        brt=fg[counter]+((r-1.0-startf)*lum_rate);
                        }
                    else
                        {
                        brt=fg[counter]*expf(lum_rate*(r-1.0-startf));
                        }
                    mat[y][x]=brt;
                    avg_pitch+=pitch;
                    num_pitch+=1.0;
                    if (list) printf("Radius: %f\t Pitch: %f Luminosity: %f\n",r,pitch,brt);
                    newpitch=pa[counter]+((int)(r-startf)*change);

//
// The pitch angle formula in unstable in the outer regions for variable
//   pitch angles.  This logic attempts to maintain the growth of the curve.
//

                    if (pa[counter] > 0.0)
                        {
                        if (delta[counter] > 0.0)
                            {
                            if (newpitch > pitch)
                                {
                                pitch=newpitch;
                                if (pitch > max_pitch) max_pitch=pitch;
                                if (pitch < min_pitch) min_pitch=pitch;
                                }
                            }
                        if (delta[counter] < 0.0)
                            {
                            if (newpitch < pitch) 
                                {
                                pitch=newpitch;
                                if (pitch > max_pitch) max_pitch=pitch;
                                if (pitch < min_pitch) min_pitch=pitch;
                                }
                            }
                        }

                    if (pa[counter] < 0.0)
                        {
                        if (delta[counter] > 0.0)
                            {
                            if (newpitch > pitch) 
                                {
                                pitch=newpitch;
                                if (pitch > max_pitch) max_pitch=pitch;
                                if (pitch < min_pitch) min_pitch=pitch;
                                }
                            }
                        if (delta[counter] < 0.0)
                            {
                            if (newpitch < pitch) 
                                {
                                pitch=newpitch;
                                if (pitch > max_pitch) max_pitch=pitch;
                                if (pitch < min_pitch) min_pitch=pitch;
                                }
                            }
                        }

//
// Make the lines thicker with a 2D feathering
//

                    if (feath[counter] > 0)
                        {
                        for ( t=1; t <= feath[counter]; t++)
                            {
                            for ( s=1; s <= feath[counter]; s++)
                                {

//
// Check the X & Y values again
//

                                if ((x > mar[counter]) && (x < (hsize[counter]-mar[counter])) && (y > mar[counter]) && (y < (vsize[counter]-mar[counter])))
                                    {
                                    mat[y-t][x]=brt;
                                    mat[y][x-s]=brt;
                                    mat[y+t][x-s]=brt;
                                    mat[y-t][x-s]=brt;
                                    mat[y+t][x]=brt;
                                    mat[y][x+s]=brt;
                                    mat[y+t][x+s]=brt;
                                    mat[y+t][x-s]=brt;
                                    }
                                }
                            }
                        }

                    }
                }
            }

//
// Calculate average pitch angle and other constants to complete drawing
//

        avg_pitch=avg_pitch/num_pitch;
        centerx=hsize[counter]/2;
        centery=vsize[counter]/2;
        r2=(int) r0[counter]*r0[counter];
        si=sin(rot[counter]*(M_PI/180.0));
        co=cos(rot[counter]*(M_PI/180.0));
        starti_s=(int) barb[counter];

//
// Fill in bar ellipse
//

        if (bara[counter])
            {
            brt=fg[counter];
            for (x=centerx-starti; x <= centerx+starti; x++)
                {
                for (y=centery-starti_s; y <= centery+starti_s; y++)
                    {
                    ma=(float)(x-centerx)*co+(float)(y-centery)*si;
                    mb=(float)(y-centery)*co-(float)(x-centerx)*si;
                    if ((pow(ma/bara[counter],2.0)+pow(mb/barb[counter],2.0)) <= 1.0)
                        {
                        mat[y][x]=brt;
                        }
                    }
                }
            }

//
// Fill in the core.  Cannot use the same polar to cartesian mapping as in
//   p2dfft because larger cores will have gaps and patterns develop that
//   can be interpreted as structure.
//

        if (core[counter])
            {
            brt=fg[counter]*(float)core[counter];
            for (x=centerx-r0[counter]; x <= centerx+r0[counter]; x++)
                {
                for (y=centery-r0[counter]; y <= centery+r0[counter]; y++)
                    {
                    if ((x-centerx)*(x-centerx)+(y-centery)*(y-centery) <= r2)
                        {
                        mat[y][x]=brt;
                        }
                    }
                }
            }

#ifdef DEBUG

    for (i=0;i<vsize[counter];i++)
        {
        for (j=0;j<hsize[counter];j++)
            {
            printf("DEBUG: mat[%d][%d]=%f\n",i,j,mat[i][j]);
            }
        }

#endif

//
// Now that we have a Cartesian matrix, create a FITS .txt file to write out.
//   Remove any old version first.
//

        if (txt)
            {
            if (verbose) printf("  --- Write %s.txt File\n",base[counter]);

            sprintf(fname,"%s.txt",base[counter]);
            sprintf(cmd,"rm -f %s",fname);
            status=system(cmd);

            ofile=fopen(fname,"w");

//
// Loop through file and write it in .txt format.  2DFFT expect 80 character
//   lines, which means five values per line.  ctr keeps track of how many
//   entries on the current line and when that gets to five, write a return
//   and start a new line.
//

            fprintf(ofile,"%14f%14f",(double)vsize[counter],(double)hsize[counter]);
            ctr=0; 
            for(i=0;i<vsize[counter];i++)
                {
                for(j=0;j<hsize[counter];j++) 
                    {
                    fprintf(ofile,"%14f",mat[i][j]);
                    ctr++;
                    if (ctr==5)
                        {
                        fprintf(ofile,"\n");
                        ctr=0;
                        }
                    }
                }

            fclose(ofile);
            }

//
// Now create a FITS file using the astro_class libraries.  Please note that 
//   when the new flag is set to 1 it will overwrite any existing file with
//   with the same name.
//
  
        if (verbose) printf("  --- Write %s.fits File\n",base[counter]);

        sprintf(fname,"%s.fits",base[counter]);
        sprintf(fname2,"!%s.fits",base[counter]);
        ast.set_warn(1);
        if (ast.fits_write(fname2, mat[0], hsize[counter], vsize[counter], 1, "p2spiral/",VERSION))
            {
            printf("ERROR: fits_write() Failed\n");
            free(mat);
            continue;
            }

//
// Add some extra key values to the FITS header.  The COLORSPC key is used by
//   some other analysis programs like SPaRcFiRe.  The other keys are for
//   convienience for testing.
//

        strcpy(keys[0],"COLORSPC");
        strcpy(items[0],"Grayscale");
        strcpy(keys[1],"ARMS");
        sprintf(items[1],"%d", arm[counter]);
        strcpy(keys[2],"AVGPITCH");
        sprintf(items[2],"%f", avg_pitch);
        strcpy(keys[3],"MINPITCH");
        sprintf(items[3],"%f", min_pitch);
        strcpy(keys[4],"MAX_PITCH");
        sprintf(items[4],"%f", max_pitch);

        if (ast.fits_header_write(fname, keys, items, 5))
            {
            printf("WARNING: fits_header_write() Failed\n");
            }

//
// Deallocate the array for this file
//

        free(mat);
        }

    printf("Total Files Processed: %d\n",num_files);
    printf("Total Errors: %d\n",errcnt);
    }
