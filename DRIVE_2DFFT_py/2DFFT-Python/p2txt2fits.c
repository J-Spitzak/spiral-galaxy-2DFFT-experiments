//
// P2TXT2FITS.C - This program will convert a FITS file in text format
//                (generated by IRAF/wtextimage into a floating point FITS file
//                with the same base name.   
//
//
// Version: 1.3  28-Aug-2017
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
// Usage: p2txt2fits [-r|--readsize] [-v|--verbose] <file list>
//
//      <file list> - List of FITS text files to bwe converted.
//
//      -r|--readsize - Force program to read the size of the FITS image from
//                      the first two bytes of the text file.  The default 
//                      behavior (recommended) is to have the program calculate
//                      size based on the image data in the text file.
//
//      -v|--verbose - Provides status messages as the program is working.
//
//
// Revision History:
//    1.3  28-Aug-2017 - Minor updates to eliminate compiler warnings on some
//                       Linux distributions
//    1.2  19-Feb-2017 - Fixed bug where files greater than 835 x 835 pixels
//                       would eb written incorrectly.
//    1.1  06-Dec-2016 - Add comments to arrayAlloc()
//    1.0  29-Nov-2016 - Initial Version:
//

#include    <math.h>
#include    <stdio.h>
#include    <stdlib.h>
#include    <string.h>
#include    <unistd.h>
#include    <getopt.h>
#include    "fitsio.h"

//  DIM_X and DIM_Y set the size of memory allocated for the data matrices, so
//   in effect are the maximum dimensions of the image file that can be
//   processed.  IMPORTANT NOTE:  These need to match the ones used in the
//   2DFFT program that produced the data files.

#define DIM_X	2048
#define DIM_Y	2048

//  DIM_RAD and DIM_THT 

#define DIM_RAD	2048
#define DIM_THT 1024

int     i,j;
int     opt;
int     ind;
int     x_dim;
int     y_dim;
int     status;
int     opterr = 0;
int     verbose = FALSE;
int     read_size = FALSE;
int     option_index = 0;

FILE    *fp_inp;
      
long    naxis=2;
long    naxes[2];

char    cmdstr[64];

float   **mat;
float   *data;

static  struct option long_options[] =
    {
    {"verbose", no_argument,     0, 'v'},
    {"readsize", no_argument,    0, 'r'},
    {0, 0, 0, 0}
    };

fitsfile    *fptr;

// Subroutines

//
// arrayAlloc() - This function will dynamically allocate a 2D array.  This is
//                needed because we need to dynamically allocate an aray to
//                generate the FITS file.  It has to be the exact size we need
//                and it needs to be a monolithic block of data (which is why
//                you can't use the tradiational series of mallocs per row to
//                create it).
//
// Arguments:
//      rows    - Number of rows (X dimension, slowest changing index)
//      columns - Number of columns (Y dimension, fastest changing index)
//
// Return Value:
//      rowptr  - pointer to base of 2D array
//

float** arrayAlloc(int rows, int columns) 
    {
    int a;

// Double index arrays have a header containing pointer to each row

    int header = rows * sizeof(float*);
    int body = 0;

    for(a=0; a<rows; a++) body+=columns;
    body*=sizeof(float);

// Allocate the total space needed

    float** rowptr = (float**)malloc(header + body);

    float* buf  = (float*)(rowptr + rows);

// Now populate all the row pointers with the correct values

    rowptr[0] = buf;
    for(a = 1; a < rows; ++a) 
        {
        rowptr[a] = rowptr[a-1] + columns;
        }
    return rowptr;
    }

//  Main Program Loop

int main(int argc, char **argv)
    {

// Get and process the command line arguments.

    while ((opt = getopt_long (argc, argv, "vr:", long_options, &option_index)) 
!= -1)
        {
        switch (opt)
            {
            case 'v':
                {
                verbose = TRUE;
                break;
                }
            case 'r':
                {
                read_size = TRUE;
                break;
                }
            default:
                {
                fprintf(stderr, "Usage: %s [-v] [-r] filenames\n", argv[0]);
                exit(EXIT_FAILURE);
                break;
                }
            }
        }

// Allocate the array for the text input file data

    data = (float *) malloc((DIM_RAD*DIM_THT*2+2) * sizeof(float *));
    if(NULL == data)
        {
        fprintf(stderr, "Memory allocation failed while allocating for data[]/n");
        exit(EXIT_FAILURE);
        }

// Process all the arguments provided. If there are no arguments, this will
//   simply fall through.
 
    if (verbose == TRUE) printf("Process Input Files\n");
    while (optind < argc)
        {
        if (verbose == TRUE) printf("File: %s\n", argv[optind]);

// Allocate the 2D matrix.  In order for the CFITSIO library functions to work,
//   the matrix needs to be same size as the FITS dimensions and contiguous.
//   This function will allocate and set up the matrix in mat.

        fp_inp=fopen(argv[optind],"r");

        printf("--- Reading Image: %s...",argv[optind]);

        ind=1;
        status=fscanf(fp_inp,"%f",&data[ind]);
        do
            {
            ind++;
            } while((fscanf(fp_inp,"%f",&data[ind]))!=EOF);

        ind--;

        printf("Done\n");

// Determine the size of the file so we can allocate the array

        if (read_size == TRUE)
            {
            if((data[1]==data[2]) && (data[1]>0.0) && (data[2]>0.0))
                {
                x_dim=(int) data[1];
                y_dim=(int) data[2];
                }
            else
                {
                printf("ERROR: File %s has dimensions %d,%d...Skipping\n",argv[optind],x_dim,y_dim);
                continue;
                }
            }
        else
            {
            x_dim=y_dim=sqrt(ind);
            }

        mat=arrayAlloc(x_dim, y_dim);

        if (verbose == TRUE) printf("%s --- dimensions : xdim=%d : ydim=%d\n",argv[optind],x_dim,y_dim);

        for(i=0;i<x_dim;i++)
            {
            for(j=0;j<y_dim;j++) 
                {
                mat[i][j]=data[(i*y_dim)+j];
                }
            }

// Now create the FITS file

// Set the size of the file

        naxes[0] = x_dim;
        naxes[1] = y_dim;

// The CFITSIO routines will fail if the file already exists

        sprintf(cmdstr,"rm %s.fits",argv[optind]);
        status=system(cmdstr);

        fits_create_file(&fptr,&cmdstr[3], &status);
    
        fits_create_img(fptr,FLOAT_IMG,naxis,naxes, &status);

        fits_write_img(fptr,TFLOAT,1,x_dim*y_dim,mat[0],&status);

        fits_close_file(fptr, &status);

        fits_report_error(stderr,status);

        optind++;
        free(mat);
        }
    return(0);
    }