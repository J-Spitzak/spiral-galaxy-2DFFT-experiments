//
// P2IFFT.CPP - This program will take the FFT output data from P2DFFT and 
//              generate a FITS image based on that data.  By default, the
//              sum of all data is used, but there are options to use only
//              some .rip files for a subset of the data.
//
//
// Version 3.4: 20-Jun-2019
//
//
// 2DFFT (original) Author: Dr. Ivanio Puerari
//                          Instituto Nacional de Astrofisica,
//                          Optica y Electronica,
//                          Santa Maria Tonantzintla,
//                          Puebla, Mexico
//
// 2DFFT (revised) Lead Author: Dr. Marc Seigar
//                              University of Minnesota Duluth,
//                              Duluth, MN USA
//
// 2DFFT (progenitor of P2DFFT) Lead Author: Dr. Benjamin Davis
//                                           Swinburne University of Technology.
//                                           Centre for Astrophysics and
//                                           Supercomputing
//                                           Melbourne, Victoria, Australia
//                                       http://d.umn.edu/~msseigar/2DFFT.html
//
// P2DFFT By: Ian Hewitt & Dr. Patrick Treuthardt,
//            NC Museum of Natural Sciences,
//            Astronomy & Astrophysics Lab,
//            Raleigh, NC USA.
//            http://github.com/treuthardt/P2DFFT
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
// Usage: p2ifft [-i|--input <file>] [-v|--verbose] [-m|--mode <n>[,<n>...]] 
//               [-s|--start <arg>] [-e|--end <arg>] [<file>[,<file>...]]
// 
//        If there is an input file specified with -i, that will be used for
//            the list of file names to be processed (one per line).  If no
//            -i option is specified and the command line arguements will
//            be used.
//
//        The file names can be the prefix of the files or can have the .fits
//        extension.  Either will work (e.g. <galaxy> and <galaxy>.fits
//        will both work.  The output files will have the form of 
//        I_<keyword>.fits or I<mode>_<keyword>.fits (if the -m option is used)
//
//        There are several command line options:
//
//              -i|--input : Will read file names and keywords from the file
//                           specified instead of the command line
//              -v|--verbose: Prints status messages during the processing.
//              -m|--mode  : Specifies the modes to be used in the output IFFT.
//                           By default, all modes (0-6) are used.  The -m 
//                           option specifies comma separated list of modes
//                           to be used in the plot, e.g. -m 1,3,5,6
//              -s|--start : Specify a starting inner radius (default is 1)
//              -e|--end   : Specify an ending inner radius (default is file
//                           size - 10%)
//
// Algorithm Notes:
//
//      The inverse image is created by reading the complex FFT frequency space
//      data for each annuli starting with an inner radius of 1 (or the start
//      value) to the ending inner radius annuli specified by end (or 10% less
//      than the maximum radius, per Davis et. al. 2012).  The FFT data is then
//      summed and a reverse two dimensional transform is done.   The resultant
//      data is then mapped from polar to cartesian space and a FITS image is
//      created.
//
// Revision History:
//      3.4  20-Jun-2019: - Fix small bug in ifft image generation
//                        - Correct/rework some DEBUG information printing
//                        - Fix bounds checking as isnan() did not detect -nan
//                        - Fix error in usage information
//                        - Clarify author/licensing information
//      3.3  03-Jun-2019: - Remove debug information printing
//                        - Change row/column sense for inverse FFT to match
//                          change in p2dfft.cpp
//                        - Fix bug in centering alogrithm for output image
//                        - Remove leftover -a option code
//      3.2  31-Jan-2019: - Fix bug introduced by compiler change in 18.04
//      3.1  18-Dec-2018: - Update comments and formatting
//                        - Fix version display bug
//                        - Update FITS mapping code to account for different
//                          ordering of CFITSIO vs C/C++
//                        - Remove ASCII FITS file output option
//                        - Fix one instance where error counter was not updated
//                        - Remove averaging for summed data
//                        - Exclude mode 0 data by default
//                        - Add checks for nan in FFT data and exclude them
//                        - Fix bug where default mode values not correct
//                        - Remove normalization as not needed because FFTW3
//                          scales the array which is fine for IFFT
//                        - Fix bug in input file handling that could cause
//                          a infinite loop
//                        - Fix small bug in verbose print statements
//                        - Fix bug in file naming for specfied modes
//                        - Fix bug where using an input file could result in
//                          no modes being selected
//      3.0  15-Sep-2018: - Remove multi-threaded code
//                        - Move directory check earlier to avoid potential
//                          crash
//                        - Standardize output message format
//                        - Move file processing message to require verbose
//                        - Add -s and -e options and update input file to
//                          support starting and ending radii for each file
//                        - Change algorithm to use sum file or create sum
//                          of FFT data to create the inverse image
//                        - Remove entering files from stdin
//      2.1  05-Feb-2018: - Remove unused variables and add global constants
//                        - Update comments/data table and check spelling
//                        - Remove unused variables
//      2.0  28-Aug-2017: - Update header comments/usage information
//                        - Change from C to C++
//                        - Major rewrite to switch from numerical recipes FFT
//                          algorithm to the libfftw3 routines.  This required
//                          major program flow changes around how parallelism
//                          was implemented as well as changing the input/output
//                          data type to be complex.
//                        - Made text file creation optional vi the -t option
//                        - Removed some old/deprecated/unused variables
//                        - Removed the -f (frequency) and (-p) phase mappings,
//                          as they were not really useful the way they were
//                          implemented and avg.py produces a more useful
//                          mapping.
//                        - Minor changes to eliminate compiler warnings on
//                          some distributions
//      1.6  29-Jul-2017: - Fix bug where mode 1 was not set as a default
//                        - Correct bug where results matrix not initialized
//                        - Add filter to prevent NAN values from IFFT from
//                          being mapped as they display improperly
//                        - Add some status messages to non-verbose mode
//                        - Add feature where -m creates mode specific filenames
//                        - Add feature to allow user to specify .fits
//                          filenames, instead of directories for ease of use.
//      1.5  23-Jan-2017: - Fix default mode to not include mode 0
//      1.4  17-Jan-2017: - Fix bug in mode selection option and made the mode
//                          bounds checking reference a constant.
//      1.3  05-Jan-2017: - Complete frequency mapping file output
//                        - Add phase angle mapping file output (-p)
//      1.2  03-Jan-2017: - Add capability to make  data file for frequency
//                          mapping (-f) -- PARTIAL
//                        - Update help message to include -v|--verbose
//                        - Add I_ prefix to output file to avoid collisions
//      1.1  06-Dec-2016: - Add averaging for data from all radii
//      1.0  29-Nov-2016: - Initial version from 2DFFT code for the first radius
//                          data.
//


//
// INCLUDE FILES
//

#include    <math.h>
#include    <stdio.h>
#include    <stdlib.h>
#include    <string.h>
#include    <getopt.h>
#include    <sys/stat.h>
#include    <sys/types.h>
#include    <omp.h>
#include    <fftw3.h>
#include    "fitsio.h"

//
// Include the Astro Functions Class from the NCNMS/NRC
//

#include    "astro_class.h"

#include    "globals.h"

//
// CONSTANTS
//

#define VERSION "3.4/20190620"

//
// Number of total frequency steps
//

#define FREQ_STEPS  200

//
// GLOBAL VARIABLES
//

int     m;        /* Mode data being processed                           */
int     t;        /* Calculated index value                              */
int     st;       /* User argument for starting radius                   */
int     en;       /* User argument for ending radius                     */
int     i,j;      /* Index variables for loops                           */
int     x,y;      /* Index variables for loops and matrices              */
int     ctr;      /* Counter for number of entries/line in .txt file out.*/
int     ind;      /* Index for extension in input file name (prefix)     */
int     val;      /* Mode value                                          */
int     dim;      /* Dimension of FITS data matrix for the current file  */
int     rmap;     /* Index used in mapping rip[] to out_data[][]         */
int     begin;    /* Beginning inner annuli value                        */
int     dummy;    /* Radius value place holder used in reading rip files */
int     finish;   /* Ending inner radius vlue                            */
int     dummy1;   /* Place holder for scanf call unused value            */
int     status;   /* CFITSIO return status                               */
int     maxrad;   /* Value of outer radius from the data files           */
int     radius;   /* Loop variable for the inner radius                  */
int     looper;   /* Loop variable for the processing loop               */
int     sindex;   /* Index for parsing mode string from input file       */
int     err_cnt;  /* Count of number of files where processing fails     */
int     counter;  /* Index value for mapping data_out[][] to mat[][]     */
int     arg_ptr;  /* Index for command line arguments                    */
int     verbose;  /* Flag to control whether status messages are output  */
int     mode_opt=0; /* Flag to indicate if mode option was specified     */
int     cmd_line; /* Flag to indicate if command line arguments are used */
int     maxrad90; /* Outer radius value - 10%                            */
int     num_files; /* Number of file to eb processed                     */
int     inp_mode=0; /* Glaf to indicate if an input file was used        */
int     count_theta;   /* Step counter for radial degrees                */
int     count_radians; /* Step counter for radial radians                */
int     option_index=0; /* Used for argument processing                  */

int     end[MAX_FILES];   /* Array for user specified starting radii     */
int     start[MAX_FILES]; /* Array for user specified ending radii       */
int     mode[MAX_FILES][M_FIN+2]; /* Flags for modes to use              */
        
long    naxis=2;       /* CFITSIO number of axes                         */
long    naxes[2];      /* CFITSIO axes size                              */

char    c;             /* Value from getopt_long(3)                      */
char    cval;          /* Character holder for mode                      */
char    *item;         /* Pointer to string for input file parsing       */
char    str[2];        /* String for stripping extension from input file */
char    tmp[64];       /* Temporary string for mode specific outfile     */
char    cstr[2];       /* Character string holder for mode               */
char    cmd[128];      /* String for command line to rm file             */
char    buff[256];     /* Line read from input file                      */
char    *last_line;    /* Pointer for find end radius in data file       */
char    fname[128];    /* Input file name from command line              */
char    infile[64];    /* Current file being read                        */
char    outfile[64];   /* Base name for output files (.txt and .fits)    */
char    modeinp[32];   /* Mode input from command line                   */
char    *last_newline; /* Ptr for end of next to last line in data file  */
char    files[MAX_ARGS]; /* Input file name from command line            */
char    base[MAX_FILES][128]; /* Input file name prefixes                */
char    mode_str[MAX_FILES][32]; /* Input mode value strings             */

float   lnr;           /* LN(R) value for polar-->Cartesian mapping      */
float   norma;         /* Normalization value from FFT                   */
float   fx, fy;        /* Floating point x,y Cartesian values            */
float   radstep;       /* Radian step increment                          */
float   **result;      /* Pointer for matrix which has final result      */
float   rip[805];      /* Array holding rip file contents                */
float   theta_step;    /* Theta angle increment                          */
float   theta_degrees; /* Value of polar mapping theta angle in degrees  */
float   theta_radians; /* Value of polar mapping theta angle in radians  */
float   mat[MAX_DIM][MAX_DIM]; /* Individual radius loop result matrix   */
float   vals[MAX_DIM][MAX_DIM]; /* Number of entries in a given x,y      */

size_t	num_read;

double  log_rad;       /* Log(2) of current radius                       */
double  log_maxrad;    /* Log(2) of maximum radius                       */

FILE    *tmp_file;    /* File input for first file to get sizes/norma    */
FILE    *ofileptr;    /* Text file output stream                         */
FILE    *file_list;   /* Directory listing from command line             */
FILE    *rip_ptr;     /* Rip file pointer                                */

astro   ast;          /* Class object for NCNMS astro_class library      */

struct  stat    sb;   /* Structure for stat command to check files/dir   */

fitsfile *fptr;       /* CFITSIO file pointer                            */


//
// MAIN ROUTINE
//

int main(int argc, char **argv)
    {
//
// Get and process the command line arguments
//

    static struct option long_options[] =
        {
        {"verbose", no_argument,     0, 'v'},
        /* These options require an argument. */
        {"start",  optional_argument, 0, 's'},
        {"end",  optional_argument, 0, 'e'},
        {"mode",  optional_argument, 0, 'm'},
        {"input", optional_argument, 0, 'i'},
        {0, 0, 0, 0}
        };
      
    while ((c = getopt_long (argc, argv, "vfs:e:i:m:", long_options, &option_index)) != -1)
        {
        switch (c)
            {
            case 'v':
                {
                verbose = 1;
                break;
                }
            case 's':
                {
                st=atoi(optarg);
                break;
                }
            case 'e':
                {
                en=atoi(optarg);
                break;
                }
            case 'm':
                {
                mode_opt=1;
                strcpy(modeinp,optarg);
                break;
                }
            case 'i':
                {
                strcpy(fname,optarg);
                break;
                }
            default:
                {
                fprintf(stderr, "Usage: p2ifft [-i|--input <file>] [-v|--verbose] [-s|--start <arg>] [-e|--end <arg>] [-m|--mode <n>[,<n>...]]\n");
                exit(1);
                break;
                }
            }
        }

    if (verbose) printf("p2ifft - Version: %s\n",VERSION);

//
// Check the start and end values, if given.  Will adjust the default endpoint
//   later, if it is too high.
//

    if ((st != 0) || (en != 0))
        {
        if (en < st)
            {
            printf("ERROR: Radius range %d to %d is invalid...Exiting\n",st,en);
            exit(1);
            }
        if ((st < 1) || (st > MAX_DIM))
            {
            printf("ERROR: Start value %d is invalid...Exiting\n",st);
            exit(1);
            }
        if ((en < 1) || (en > MAX_DIM))
            {
            printf("ERROR: End value %d is invalid...Exiting\n",en);
            exit(1);
            }
        }

//
// Check the input file name, if given, and read the parameters.  If no input
//   file, check command line arguments.  Either way, get number of files to
//   process fill arrrays with file information/parameters
//

    num_files=0;
    err_cnt=0;

    if ((fname[0] != '\0')&&(fname[0] != '\1'))
        {
        file_list=fopen(fname,"r");
        if ( file_list == NULL )
            {
            printf("ERROR: Cannot open input file - %s\n",fname);
            exit(1);
            }
        while (fgets(buff,256,file_list)!=NULL)
            {
//
// Check for comment or blank line
//

            if ((buff[0]=='#') || (strlen(buff) < 2)) continue;

//
// If we have reached the maximum file number, we need to flag an error.
//

            if (num_files == MAX_FILES)
                {
                printf("ERROR: Too many input lines!\n");
                exit(1);
                }
//
// Break the line into individual fields.
//

            if ((item=strtok(buff,",\t\n "))!=NULL)
                {
                strcpy(base[num_files],item);
                if ((item=strtok(NULL,",\t\n "))!=NULL)
                    {
                    strcpy(mode_str[num_files],item);
                    inp_mode=1;
                    if ((item=strtok(NULL,",\t\n "))!=NULL)
                        {
                        start[num_files]=atoi(item);
                        if (start[num_files] < 1)
                            {
                            printf("WARNING: Invalid Start on Line %d\n",num_files+1);
                            err_cnt++;
                            continue;
                            }
                        if ((item=strtok(NULL,",\t\n "))!=NULL)
                            {
                            end[num_files]=atoi(item);
                            if (end[num_files] < 1)
                                {
                                printf("WARNING: Invalid End on Line %d\n",num_files+1);
                                err_cnt++;
                                continue;
                                }
                            }
                        }
                    }
                else
                    {
                    sprintf(mode_str[num_files],"123456");
                    }
                }
            else
                {
                printf("WARNING: Invalid Filename on Line %d\n",num_files+1);
                err_cnt++;
                continue;
                }
            num_files++;
            }
        }
    else
        {
//
// No file, so go through command line arguments and build the file list
//

        if (optind < argc)
            {
            for (i=optind; i < argc; i++)
                {
                if ((item=strtok(argv[i],",\t "))!=NULL)
                    {
                    strcpy(base[num_files],item);
                    num_files++;
                    }
                }
            }
        else
            {
            printf("ERROR: No files specified\n");
            exit(1);
            }
        }

//
// If there was an input file, assign the modes here
//

    if (inp_mode == 1)
        {
        for (i=0; i <= num_files; i++)
            {
            if (mode_str[i] != NULL)
                {
                sindex=0;
                while ( mode_str[i][sindex] != '\000' )
                    {
                    cval=mode_str[i][sindex];
                    sprintf(cstr,"%c",cval);
                    val=atoi(cstr);
                    if ((val >= M_INI) || (val <= M_FIN))
                        {
                        mode[i][val]=1;
                        }
                    else
                        {
                        printf("WARNING: Unknown mode %c on line %d\n",cval,i);
                        }
                    sindex++;
                    }
                }
            else
                {
                for (j=M_INI+1; j <= M_FIN; j++) mode[i][j]=1;
                }
            }
        }


//
// Check the mode values from the command line or the input file
//
   
    if ((modeinp[0] != '\0')&&(modeinp[1] != '\1'))
        {
        item=strtok(modeinp,",");
        while ( item != NULL )
            {
            val=atoi(item);
            if ((val >= M_INI) || (val <= M_FIN))
                {
                for (i=0; i <= num_files; i++) mode[i][val]=1;
                }
            else
                {
                printf("ERROR: Unknown mode %s\n",item);
                exit(1);
                }
            item = strtok(NULL,",");
            }
        }
    else
        {
        if (inp_mode==0)
            {
            for (i=0; i <= num_files; i++)
                {
                for (j=1; j <= M_FIN; j++) mode[i][j]=1;
                }
            }
        }

    fftw_plan   plan;

//
// Allocate the FFT arrays.  These need to be allocated with fftw_ functions
//   since they are not C-style 2D arrays.
//

    fftw_complex *in_data;
    fftw_complex *out_data;

    if (verbose) printf("Allocating FFT Arrays...");

    in_data = (fftw_complex *) fftw_malloc((DIM_RAD*DIM_THT+1) * sizeof(fftw_complex));
    if(NULL == in_data)
        {
        printf("ERROR: FFTW Memory allocation failed for in_data[]/n");
        exit(-1);
        }

    out_data = (fftw_complex *) fftw_malloc((DIM_RAD*DIM_THT+1) * sizeof(fftw_complex));
    if(NULL == out_data)
        {
        printf("ERROR: FFTW Memory allocation failed for out_data[]/n");
        exit(-1);
        }

//
// Build a plan for the FFT transform
//

    if (verbose) printf("Building plan for FFT...");
    plan=fftw_plan_dft_2d( (int) DIM_THT, (int) DIM_RAD, in_data, out_data, FFTW_BACKWARD, FFTW_MEASURE);
    if ( plan == NULL )
        {
        printf("ERROR: FFTW plan failed, Exiting\n");
        exit(1);
        }
    else
        {
        if (verbose) printf("Done\n");
        }

// 
// MAIN LOOP through the list of input files in order to process them
//

    for(looper=0; looper < num_files; looper++)
        {

//
// Need to set start and end values appropriately
//

        begin=finish=-1;
        if ((st != 0)&&(en !=0))
            {
            begin=st;
            finish=en;
            }
        else
            {
            if (start[looper]>0) begin=start[looper];
            if (end[looper]>0) finish=end[looper];
            }

        if (begin < 0) begin=1;

//
// If the filename ends in .fits (for user convenience), strip it
//

        if (strlen(base[looper]) > 5)
            {
            ind=strlen(base[looper])-5;
            if (strstr(base[looper],".fits") == &base[looper][ind])
                {
                base[looper][ind]='\000';
                }
            }

//
// Read the first file to get the maximum radius.  If we fail to find this,
//   skip it and move on to the next file.
//

        if (verbose) printf("  --> Processing Files for %s\n",base[looper]);

        sprintf(infile,"%s_m1",base[looper]);

//
// If failure, just move on to next file
//

        if ((tmp_file=fopen(infile,"r")) == NULL)
            {
            printf("WARNING: Cannot Get Radius From %s\n...Skipping this directory\n",infile);
            err_cnt++;
            continue;
            }

//
// Read the outer radius value from the end of the file.  Skip to next file if
//   it's not a valid number.  Do this by seeking the end and then doig a 
//   reverse read until the newline.
//

        fseek(tmp_file, -81, SEEK_END); 
        num_read=fread(buff, 80, 1, tmp_file);
        fclose(tmp_file);

        buff[27] = '\0';
        last_newline = strrchr(buff, '\n');
        last_line = last_newline+1;

        status=sscanf(last_line,"%d outi%d_m%d",&dummy, &maxrad, &dummy1);

//
// The last 10% of the radius calculations tend to return invalid results (like
//   -nan values, so reduce radius by 10%
//

        maxrad90=(int)((float)maxrad*0.9);

//
// Do some checking on the radius to see if it's logical
//

        if ((maxrad90 < 1) || (maxrad90 > MAX_DIM/2))
            {
            printf("WARNING: Invalid radius %d in file %s...Skipping\n",maxrad90,infile);
            err_cnt++;
            continue;
            }

        if ((maxrad < 1) || (maxrad > MAX_DIM/2))
            {
            printf("WARNING: Invalid radius %d in file %s...Skipping\n",maxrad,infile);
            err_cnt++;
            continue;
            }

        if (verbose) printf("%s: Radius=%d:%d ",base[looper],maxrad,maxrad90);

//
// New set finish based on the radius
//

        if (finish<0)
            {
            finish=maxrad90;
            }
        else
            {
            if (finish > maxrad90)
                {
                finish=maxrad90;
                printf("WARNING: End radius beyond 90 percent for file %s...Trimming to %d\n",base[looper],finish);
                }
            }

//
// Check if the data files or directory exists for the file.  If not, then skip
//   it and move to the next.
//

        if ((stat(base[looper], &sb) != 0) || !S_ISDIR(sb.st_mode))
            {
            printf("WARNING: Directory %s does not exist -- Skipping...\n",base[looper]);
            continue;
            err_cnt++;
            }

//
// Zero out the mat and val arrays completely, even though we only use a
//   subset of them.  Needs to be done once for each new file/directory, not
//   each radius.
//
        for (i=0; i < MAX_DIM; i++)
            {
            for (j=0; j < MAX_DIM; j++)
                {
                mat[i][j]=0.0;
                vals[i][j]=0.0;
                }
            }

//
// Set up the dimension variables as these will be the same for each radius.
//

        dim=(maxrad*2)+1;

        radstep=2.*PI/STEP_P/DIM_RAD;
        theta_step=2.*PI/GR_RAD/DIM_THT;

//
// Must completely zero the data arrays each time or you get incorrect results
//

        for(x=0;x<DIM_RAD*DIM_THT+1;x++)
            {
            in_data[x][0]=0.0;
            in_data[x][1]=0.0;
            out_data[x][0]=0.0;
            out_data[x][1]=0.0;
            }

//
// Loop for each radius.  Do the logarithm of the max radius outside the loop
//   for better performance.
//

        log_maxrad=log((double)finish);

        for (radius = begin; radius <= finish; radius++)
            {
//
// Create filename base.  Use 7 for the mode now and will change that as needed.
//

            sprintf(infile,"%s/outi%d_m%d.rip",base[looper],radius,7);

//            
// Read in files by mode.  
//

            for (m=M_INI; m <= M_FIN; m++)
                { 
//
// Check if this mode is to be plotted
//

                if (!mode[looper][m])
                    {
                    continue;
                    }
//
// Update filename and try to open
//

                sprintf(str,"%d",m);
                infile[strlen(infile)-5]=str[0];
                if ((rip_ptr=fopen(infile,"r")) == NULL)
                    {
                    printf("WARNING: Cannot open %s\n...Skipping\n",infile);
                    continue;
                    }
            
                if (verbose) printf("--- Adding %s from File\n",infile);

//    
// Read header information.  We already had radius, so use a dummy variable, but
//   the normalization value is needed.
//

                status=fscanf(rip_ptr,"%d",&dummy);
                status=fscanf(rip_ptr,"%e",&norma);

                if (verbose) printf("Norma=%f\n",norma);

//
// Read data into the part of the array corresponding to the mode.  There
//  should be 802 entries in the rip file (not including the radius value
//  and the normalization number.  Since it's not a straight mapping of the
//  rip file values to the data file entries, it's read into the rip array.
//  Also, print a warning if the numbers don't add to the right amount.
//

                for (i=0; i< 805; i++) rip[i]=0.0;

                counter=1;
                while((fscanf(rip_ptr,"%e",&rip[counter-1]))!=EOF) counter++;
                if (verbose && (counter != (3+(FREQ_STEPS*4)))) printf("WARNING: Count for File %s was not %d, but %d.. Continuing Anyway....\n",infile,(FREQ_STEPS*2)+3,counter);
                fclose(rip_ptr);

                if (DEBUG) printf("Counter=%d\n",counter);

//
// Map the rip entries to the data array. P2DFFT writes these in frequency
//   order going from lowest to highest, but that's not the same order in 
//   the data file.  This comment is long, but it's important if you want
//   to make code changes (like adding a different frequency cutoff). Of
//   course, it is vital that this mapping match the P2DFFT program, or the
//   output will be meaningless.
//
//     Data Array   Description               Fft_data     Rip File   rip[]
//     ----------   ---------------------     --------     --------   ----
//        0,0       Real Mid Freq (0/DC)        1025        403       401
//        0,1       Imag Mid Freq               1025        404       402
//        1,0       Real Min Pos. Freq Start    1026        405       403
//        1,1       Imag Min Pos. Freq Start    1026        406       404
//       200,0      Real Pos. Freq              1225        803       801
//       200,1      Imag Pos. Freq              1225        804       802
//       201,0      Real Pos. Freq              1226        N/S       N/S
//       201,1      Imag Pos. Freq              1226        N/S       N/S
//      1024,0      Real Max Pos. Freq End      2049        N/S       N/S
//      1024,1      Imag Max Pos. Freq End      2049        N/S       N/S
//      1025,0      Real Min Neg. Freq             2        N/S       N/S
//      1025,1      Imag Max Neg. Freq             2        N/S       N/S
//      1847,0      Real Neg. Freq               824        N/S       N/S
//      1847,1      Imag Neg. Freq               824        N/S       N/S
//      1848,0      Real Neg. Freq               825          3         1
//      1848,1      Imag Neg. Freq               825          4         2
//      2047,0      Real Max Neg. Freq End      1024        401       399
//      2047,1      Imag Max Neg. Freq End      1024        402       400
//
//   A few notes on the table.  The fft_data column is for the indices of the
//   array used in the P2DFFT program.  These are used for an intermediate
//   transform (but not used in this program).
//
//   This mapping is also for m0 only, other modes would have the same fft_data
//   and rip values, but the data indices would have mode*2048 added to them. 
//
//   Finally, this table also assumes a fequency mapping of -50 to +50, if thats
//   different, then the start and end will be different.
//

                rmap=0;
                for (x = 0; x < FREQ_STEPS ; x++)
                    {
                    if (!(rip[rmap]!=rip[rmap])&& (!isinf(rip[rmap]))) in_data[x+1848+(m*2048)][0]+=rip[rmap];
                    if (!(rip[rmap+400]!=rip[rmap+400])&& (!isinf(rip[rmap+400]))) in_data[x+(m*2048)][0]+=rip[rmap+400];
                    rmap++;
                    if (!(rip[rmap]!=rip[rmap])&& (!isinf(rip[rmap]))) in_data[x+1848+(m*2048)][1]+=(-1.0)*rip[rmap];
                    if (!(rip[rmap+400]!=rip[rmap+400])&& (!isinf(rip[rmap+400]))) in_data[x+(m*2048)][1]+=(-1.0)*rip[rmap+400];
                    rmap++;

                    if (DEBUG && radius==1)
                        {
                        printf("Map rip[%d]=%e to in_data[%d][0]\n",rmap,rip[rmap],x+1848+(m*2048));
                        printf("Map rip[%d]=%e to in_data[%d][0]\n",rmap+400,rip[rmap+400],x+(m*2048));
                        printf("Map rip[%d]=%e to in_data[%d][1]\n",rmap,rip[rmap],x+1848+(m*2048));
                        printf("Map rip[%d]=%e to in_data[%d][1]\n",rmap+400,rip[rmap+400],x+(m*2048));
                        }
                    }
                if ((!(rip[800]!=rip[800])) && (!isinf(rip[800]))) in_data[200+(m*2048)][0]+=rip[800];
                if ((!(rip[801]!=rip[801])) && (!isinf(rip[801]))) in_data[200+(m*2048)][1]+=(-1.0)*rip[801];
                }

//
// This if statement is handy for debugging if you make changes
//

            if (DEBUG && radius == 1)
                {
                for (x=0;x< DIM_RAD*DIM_THT; x++)
                    {
                    printf ("In Data[%d][0]=%f\n",x,in_data[x][0]);
                    printf ("In Data[%d][1]=%f\n",x,in_data[x][1]);
                    }
                }
            }
                        
//
// This if statement is handy for debugging if you make changes
//

        if (DEBUG)
            {
            for (x=0;x< DIM_RAD*DIM_THT; x++)
                {
                printf ("All In Data[%d][0]=%f\n",x,in_data[x][0]);
                printf ("All In Data[%d][1]=%f\n",x,in_data[x][1]);
                }
            }
        fftw_execute(plan);

//
// Need to normalize the data from number of signal points
//

        for (x=0; x < DIM_RAD*DIM_THT; x++) out_data[x][0]/=(DIM_RAD*DIM_THT);

//
// This if statement is handy for debugging if you make changes
//

        if (DEBUG)
            {
            for (x=0;x< DIM_RAD*DIM_THT; x++)
                {
                printf ("Out Data[%d][0]=%f\n",x,out_data[x][0]);
                printf ("Out Data[%d][1]=%f\n",x,out_data[x][1]);
                }
            }

//
// Map the polar coordinates in out_data[][] to Cartesian. The original mapping 
//   in P2DFFT makes radial slices in small steps of theta, so this just
//   reverses the mapping.  Please note that some values are duplicated in the
//   polar version, so we need to account for this when mapping back to cart.
//

        if (verbose) printf("Transform data lnr theta ---> X,Y\n");
        counter=0;
        count_theta=1;
        log_rad=log((double)finish);

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
                if(lnr>(double)log_rad) 
                    {
                    ++counter;
                    continue;
                    }

                fx=exp(lnr)*cos(theta_radians);
                fy=exp(lnr)*sin(theta_radians);

                x=(int)fx+maxrad+1;
                y=(int)fy+maxrad+1;

//
// If the data is valid, add the result to the master matrix and increment the
//   number of values used.  This will be used to normalize the total value
//   later.  Invalid values tend to occur at outer radii, but can happen in
//   other places.  Having a NAN result in the image file is allowed, but
//   causes programs to ds9 to display the data in a less useful way.
//

                if (!(out_data[counter][0]!=out_data[counter][0]))
                    {
                     mat[x][y]+=out_data[counter][0];
    
                     vals[x][y]+=1.0;

                    if (DEBUG) printf("Radius: %d: Assign Mat[%d][%d]=%f,radius,vals[%d][%d]=%f, Index=%d\n",radius,x,y,out_data[counter][0],x,y,vals[x][y],counter);
                    }
                ++counter;
                }
            }

//
// Need to create a matrix with the exact size and average data from all runs.
//   This is needed because the CFITSIO libraries can return incorrect results
//   if part of a larger array is used.  In addition, the matrix must have a
//   contiguous memory space, so you the astro_class function to allocate.
//

        if (verbose) printf("Creating Output File...\n");

        result = ast.ArrayAlloc(dim,dim);

        for(i=0;i<dim;i++)
            {
            for(j=0;j<dim;j++)
                {
                result[i][j]=0.0;
                }
            }

//
// Normalize the values
//

        for(i=0;i<dim;i++)
            {
            for(j=0;j<dim;j++)
                {
                if (vals[j][i] != 0.0)
                    {
                    result[j][i] = mat[i][j] / vals[i][j];
                    if (DEBUG) printf("Result[%d][%d]=%f mat=%f vals=%f\n",j,i,result[j][i],mat[i][j],vals[i][j]);
                    }
                }
            }

// 
// Create binary FITS file
//

        naxes[0] = (long) dim;
        naxes[1] = (long) dim;
        
        if ((mode_opt == 1) || (inp_mode == 1))
            {
            sprintf(outfile,"%s", "I_");
            for (i = 1; i <= M_FIN; i++)
                {
                if (mode[looper][i]==1)
                    {
                    sprintf(tmp,"%d",i);
                    strcat(outfile,tmp);
                    }
                }
            strcat(outfile,"_");
            strcat(outfile,base[looper]);
            strcat(outfile,".fits");
            }
        else
            {
            sprintf(outfile,"I_%s.fits",base[looper]);
            }

//
// The CFITSIO libraries don't work if the file already exists, so remove it
//

        sprintf(cmd,"rm -f %s",outfile);
        status=system(cmd);

        fits_create_file(&fptr,outfile, &status);

        fits_create_img(fptr,FLOAT_IMG,naxis,naxes, &status);

        fits_write_img(fptr,TFLOAT,1,dim*dim,result[0],&status);

        fits_close_file(fptr, &status);

        fits_report_error(stderr,status);
    
        free(result);
        }

//
// Clean up open files and memory allocations
//

    fftw_destroy_plan(plan);
    fftw_free(in_data);
    fftw_free(out_data);

    if (file_list) fclose(file_list);
    if (verbose) printf("Closing....\n");
    }
