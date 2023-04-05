//
// ASTRO.CPP - This class provides functions for astronomical analysis,
//             specifically the P2DFFT pitch angle analysis and hurricane
//             packages.
//
//
// Version 3.0: 12-Jun-2018
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
// Revision History:
//      3.0  12-Jun-2018: - Update FITS data read/write routines to use 2D
//                          functions and to compensate for row/col ordering
//                        - Fix fits_read() to allocate a buffer based on the 
//                          size specified and not a fixed nnumber
//      2.2  26-May-2018: - Add fits_write() function code
//                        - Add fits_header_write() function code
//                        - Add ASTRO_******* return codes
//                        - Update comments
//                        - Align error message output
//                        - Remove placeholders for fits_text_read/write
//                        - Fix bug  where multiple calls to a function
//                          might fail
//      2.1  16-Mar-2018: - Fix bug in debug statement (very meta)
//                        - Add warn flag and set_warn() function
//                        - Add warn flag check to print statements
//                        - Add prefixes to print statements (ERROR, WARNING)
//                        - Add astro_class.h version output to version()
//      2.0  05-Feb-2018: - Add version print function
//                        - Update comments
//                        - Use globals.h instead of local defs (e.g. DEBUG)
//                        - Add get_err function and remove extern var for
//                          error numbers
//                        - Removed deprecated cmd_path() function
//                        - Remove the astro_print variable and have all
//                          routines automatically print error messages
//                        - Remove requirement that images are square
//      1.3  28-Aug-2017: - Minor changes to eliminate compiler warnings on
//                          some Linux distributions
//      1.2  07-Jun-2017: - CRITICAL Fix:  Remove all C++11 code and replace
//                          it with traditional C++ code because some
//                          compiler/library combinations caused the C++11
//                          code to return erroneous results
//                        - Remove redefinition of BIN_DIR since it is 
//                          passed in on the compiler line from the makefile
//      1.1  21-Feb-2017: - Change fitsio.h include to system (<>) from
//                          from local ("")
//                        - All error message did not return the FITS text
//                          error message along with the code.  This was
//                          added where missing.
//                        - Remove unused variable
//      1.0  19-Feb-2017: - Initial version
//

#define ASTRO_VER   "3.0/20180612"

#include    <stdio.h>
#include    <string.h>
#include    <unistd.h>
#include    <fstream>
#include    <sstream>
#include    <magic.h>
#include    <sys/stat.h>
#include    <sys/types.h>

#include    "astro_class.h"
#include    <fitsio.h>

#include    "globals.h"

int         astro_warn=0;

//
// Define macro and variable for error handling
//

int     astro_errno=0;

#define set_astro_errno(err) (astro_errno = (err))

//
// FUNCTION BLOCK
//

//
// SET_WARN() - Set the value of the wrning flag to indicate if warnings 
//              should be printed to standard out
//
// Arguments: NONE
//      0 for no warnings.   Error value if there is a an error/warning.
//
// Return Value: NONE
//

void   astro::set_warn(int  value)
    {
    astro_warn=value;
    }


//
// GET_ERR() - This function will return the latest error number
//
// Arguments: NONE
//
// Return Value: Most recent error code defined in astro_class.h
//

int    astro::get_err()
    {
    return(astro_errno);
    }


//
// VERSION() - This function will print the current version.
//
// Arguments: NONE
//
// Return Value: NONE
//

void    astro::version()
    {
    printf("  -- Astro Class Include Version:  %s\n",ASTRO_H_VER);
    printf("  -- Astro Class Function Version:  %s\n",ASTRO_VER);
    }


//
// FILE_TYPE() - This function will return a value based on the file type
//               determined by the magic number.
//
// Arguments:
//      fname   - Text string of filename
//
// Return Value:
//      ASTRO_TXT_FILE    - Likely TXT FITS file
//      ASTRO_BIN_FILE    - Binary FITS file
//      ASTRO_UNK_FILE    - None of the Above
//

int    astro::file_type(std::string fname)
    {
    int         ret=-1;
    magic_t     handle;
    const char  *type;

    handle=magic_open(MAGIC_NONE|MAGIC_COMPRESS);
    magic_load(handle,NULL);
    type = magic_file(handle,fname.c_str());
    
    if (strstr(type,"FITS image data"))
        {
        ret=ASTRO_BIN_FILE;
        }
    else
        {
        if (strstr(type,"ASCII"))
            {
            ret=ASTRO_TXT_FILE;
            }
        else
            {
            ret=ASTRO_UNK_FILE;
            }
        }
    magic_close(handle);    
    return(ret);
    }


//
// FILE_EXISTS() - This function will return true if a give file exists,
//                 otherwise, it return false.
//
// Arguments:
//      fname   - Text string of filename
//
// Return Value:
//      TRUE    - File exists
//      FALSE   - File exists
//

bool    astro::file_exists(const char *fname)
    {
    std::ifstream infile(fname);
    return infile.good();
    }


//
//   READ_LINES() - Reads the contents of file and populates the file_rec
//                  structure with the results.
//
// Arguments:
//      fname   - String with file name to be read
//      rec     - Pointer to vector (array) of file_rec structs (astro_class.h)
//
// Return Value:
//      ASTRO_SUCCESS   - Success
//      ASTRO_FAILURE   - Failure (astro_errno will be set with detailed code)
//
// Errors:  Function will set astro_errno with return code (see astro_class.h)
//

int     astro::read_lines(std::string fname, std::vector<file_rec> *rec)
    {
    int         x, y;
    int         calc_rad;
    std::string token;

//
// Try to open file
//

    std::ifstream   fs;
    const char  * fname2 = fname.c_str();
    fs.open(fname2);

//
// Check if any error flags are set
//

    if (DEBUG) printf("DEBUG: astro::read_lines() Calling Open\n");

    if (!fs.good())
        {
        if (astro_warn) printf("WARNING: astro::read_lines: Filename Error\n");
        set_astro_errno(ASTRO_ERR_OPEN);
        return(ASTRO_FAILURE);
        }

    std::string     line;

    if (DEBUG) printf("DEBUG: astro::read_lines() Starting Loop\n");

    while (std::getline(fs, line))
        {
        calc_rad=0;
        if (DEBUG) std::cout << "DEBUG: Line: " << line << ":" << std::endl;

//
// Ignore any blank lines
//

        if (line.empty())
            {
            if (DEBUG) printf("DEBUG: Skipping Blank Line...\n");
            continue;
            }

//
// Set up variables to read in the line.  We stream the line read from
//   the file as a character stream, ss
//
 
        file_rec    f;
        std::istringstream  ss(line);

//
// Get the file name.  This should always be there.  If not, it's an error
//   and skip the line.
//

        std::getline(ss, f.name, ',');
        if (DEBUG) std::cout << "DEBUG: Name: " << f.name << std::endl;

//
// Try to get the keyword  entry, but it may not exist
//

        std::getline(ss, f.result, ',');

//
// The file may only contain filenames, so if we can't read a keyword, need
//   to determine the keyword and radius
//

        if (f.result.empty())
            {
//
// If there is EOF, we only have the filename, so we need to create a keyword
//

            calc_rad=1;

            std::string delim = ".";
            f.result = line.substr(0, line.find(delim));

            if (DEBUG) std::cout << "DEBUG: Calculated Result: " << f.result << std::endl;
            }
        else
            {
//       
// We have a valid keyword specified
//

            if (DEBUG) std::cout << "DEBUG: Read Result: " << f.result << std::endl;
            }

//
// Next we need to determine a radius.
//

        f.keyword ="outi";
        std::getline(ss, token, ',');

        if ((calc_rad==1) || (token.empty()))
            {
//
// There is no radius specified, so either need to read it or calculate it
//
// Determine if it's a text file by trying to read the header
//

            if (fits_dims(f.name, &x, &y))
                {
//
// It's not a binary FITS file, so set it to -1  and P2DFFT will calculate
//  the size
//

                if (DEBUG) std::cout << "DEBUG: Provisional Header Radius: -1" << std::endl;
                f.binary = 0;
                f.radius = -1;
                f.valid = 0;
                }

            else
                {
//
// Valid FITS header, now populate the structure.
//

                f.binary=1;

                if (DEBUG) std::cout << "DEBUG: Read Header Radius: " << x << " So " << (x-1)/2 << std::endl;
                f.radius=(x-1)/2;
                f.valid=1;
                }
            }
        else
            {
//
// Need to guess the format
//

            if ((f.name.substr(f.name.find_last_of(".") + 1)=="fits")||(f.name.substr(f.name.find_last_of(".") + 1)=="fts"))
                {
                f.binary=1;
                }
            else
                {
                f.binary=0;
                }
            const char  * token2 = token.c_str();
            f.radius=atoi(token2);
            f.valid=1;
            if (DEBUG) std::cout << "DEBUG: File Header Radius: " << f.radius << std::endl;
            }
        token.clear();
        rec->push_back(f);
        }

    return(ASTRO_SUCCESS);
    }


//
// FITS_DIMS() - Reads the FITS header of a file and returns the rows and
//               columns in the file.  This routine currently assumes only
//               two axes.
//
//               PLEASE NOTE: FITS and C/C++ have different senses of 
//               slowest and fastest varying dimensions (they are opposite).
//               This function preserves FITS ordering of the dimensions.
//
// Arguments:
//      fname   - Text string of filename to be read
//      rows    - Number of rows (X dimension, fastest changing index)
//      columns - Number of columns (Y dimension, slowest changing index)
//
// Return Value:
//      ASTRO_SUCCESS  - Success getting values
//      ASTRO_FAILURE  - Failure
//
// Errors:  Function will set astro_errno with return code (see astro_class.h)
//

int     astro::fits_dims(std::string fname, int *rows, int *cols)

    {
    int         status=0;
    char        err_text[81];
    char        *file;
    long        naxes[2];
    fitsfile    *dim_p=NULL;

    file=(char *) fname.c_str();
    if (fits_open_file(&dim_p, file, READONLY, &status))
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_dims:fits_open_file() Error %d: %s\n", status, err_text);
        set_astro_errno(ASTRO_ERR_OPEN);
        return(ASTRO_FAILURE);
        }

    if (fits_get_img_size(dim_p, 2, naxes, &status))
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_dims:fits_get_img_size() Error %d: %s\n", status, err_text);
        fits_close_file(dim_p, &status);
        set_astro_errno(ASTRO_ERR_GET_SIZE);
        return(ASTRO_FAILURE);
        }

    *rows=(int)naxes[0];
    *cols=(int)naxes[1];

    if (DEBUG) printf("DEBUG: rows =%d:cols=%d:\n",*rows,*cols);

    if (fits_close_file(dim_p, &status))
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_dims:fits_close_file() Error %d: %s\n", status, err_text);
        set_astro_errno(ASTRO_ERR_CLOSE);
        return(ASTRO_FAILURE);
        }

    return(ASTRO_SUCCESS);
    }


//
// FITS_HEADER_READ() - Reads all fields of the FITS header and returns a
//                      character array with the information (on record per
//                      string)
//
// Arguments:
//      fname   - Test string with FITS filename to be read
//      nkey    - Pointer to variable which will contain the number of keys
//                (fields) in the header
//
// Return Value:
//      char ** - Array of strings (one string per header record/field)
//
// Errors:  Function will return NULL andset astro_errno with return code
//          (see astro_class.h)
//

char    **astro::fits_header_read(char *fname, int *nkeys)
    {
    int         i, pos, status=0;
    char        card[2048];
    char        err_text[81];
    char        **hdr_blk;
    fitsfile    *hdr_p=NULL;

    if (DEBUG) printf("DEBUG: astro::fits_header_read:Call fits_open_file()\n");

    if (fits_open_file(&hdr_p, fname, READONLY, &status)) 
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_header_read:fits_open_file() Error %d: %s\n", status, err_text);
        set_astro_errno(ASTRO_ERR_OPEN);
        return(NULL);
        }
         
    if (DEBUG) printf("DEBUG: astro::fits_header_read:Call fits_get_hdrpos()\n");

    if (fits_get_hdrpos(hdr_p, nkeys, &pos, &status))
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_header_read:fits_get_hdrpos() Error %d: %s\n", status, err_text);
        fits_close_file(hdr_p, &status);
        set_astro_errno(ASTRO_ERR_HDR_POS);
        return(NULL);
        }

//
// Allocate memory for return value
//

    if (DEBUG)
        {
        printf("DEBUG: astro::fits_header_read:Call ArrayAlloc()\n");
        }

    if ((hdr_blk=CArrayAlloc(*nkeys, 2048))==NULL)
        {
        if (astro_warn) printf("WARNING: astro::fits_header_read:CArrayAlloc() Error\n");
        fits_close_file(hdr_p, &status);
        set_astro_errno(ASTRO_ERR_MALLOC);
        return(NULL);
        }

    if (DEBUG)
        {
        printf("DEBUG: astro::fits_header_read:Keys=%d, hdr_blk=%p,%p\n",*nkeys,hdr_blk,hdr_blk[0]);
        printf("DEBUG: astro::fits_header_read:Start read_record loop\n");
        }

    for (i=1; i <= *nkeys; i++)
        {
        if (fits_read_record(hdr_p, i, card, &status))
            {
            fits_get_errstatus(status,err_text);
            if (astro_warn) printf("WARNING: astro::fits_header_read:fits_read_record() Error %d: %s\n", status, err_text);
            free(hdr_blk);
            fits_close_file(hdr_p, &status);
            set_astro_errno(ASTRO_ERR_RD_REC);
            return(NULL);
            }

        if (DEBUG) printf("DEBUG: %d:**%s**\n",i,card);
        strncpy(hdr_blk[i-1],card,2048);
        }

    if (DEBUG) printf("DEBUG: astro::fits_header_read:Call fits_close_file()\n");

    if (fits_close_file(hdr_p, &status))
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_header_read:fits_close_file() Error %d: %s\n", status, err_text);
        free(hdr_blk);
        set_astro_errno(ASTRO_ERR_CLOSE);
        return(NULL);
        }

    return(hdr_blk);
    }


//
// FITS_HEADER_WRITE() - Routine to write one or more keys to a FITS file
//
// Arguments:
//      fname   - Text filename for FITS file to be written
//      keys    - Array of strings with the key names to be written
//      items   - Array of strings containing the value for each key
//      num     - number of keys to be written
//
// Return Value:
//      ASTRO_SUCCESS - All keys written successfully
//      ASTRO_FAILURE - Error writing one or more keys
//
// Errors:  Function will set astro_errno with return code (see astro_class.h)
//

int    astro::fits_header_write(char *fname, char keys[][32], char items[][80], int num)
    {
    int         i;
    int         status=0;
    char        err_text[81];
    fitsfile    *fptr=NULL;

//
// Open the file for read/write
//

    if (DEBUG) printf("DEBUG: astro::fits_header_write:fits_open_file\n");

    if (fits_open_file(&fptr, fname, READWRITE, &status))
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_header_write:open() Error %d: %s\n",status,err_text);
        set_astro_errno(ASTRO_ERR_OPEN);
        return(ASTRO_FAILURE);
        }

//
// Loop through the keys
//

    for (i=0; i<num; i++)
        {
        if (fits_write_key(fptr, TSTRING, keys[i], items[i], NULL, &status))
            {
            fits_get_errstatus(status,err_text);
            if (astro_warn) printf("WARNING: astro::fits_header_write:fits_write_key() Error %d: %s\n", status, err_text);
            set_astro_errno(ASTRO_ERR_KEY);
            return(ASTRO_FAILURE);
            }
        }

    if (fits_close_file(fptr, &status))
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_header_write:fits_close_file() Error %d: %s\n", status, err_text);
        set_astro_errno(ASTRO_ERR_CLOSE);
        return(ASTRO_FAILURE);
        }

    return(ASTRO_SUCCESS);
    }


//
// FITS_READ() - Reads the image data from a 2D binary FITS file and returns
//               it in one dimensional array.   All the data read will be
//               converted to floating point format regardless of it's 
//               encoding in the FITS file.   
//
//               PLEASE NOTE:  The sense of indices will be preserved.  The
//               fastest varying index will be rows and the slowest varying
//               index will be columns.  This is opposite of how C/C++ usually
//               manages these.
//
//               PLEASE ALSO NOTE:  FITS files will start data at index 1,
//               not zero.  Size will be set to the number of entries, but
//               this means the array in C will go from [0...size], instead
//               [0...size-1].
//
// Arguments:
//      fname   - Text filename for FITS file to be read
//      size    - Pointer to variable that will be set to the number of data
//                entries in the return array.
//
// Return Value:
//      float * - pointer to base of one dimensional array with image data
//                NOTE: <FIX>
//
// Errors:  Function will return NULL and set astro_errno with return code
//          (see astro_class.h)
//

float   *astro::fits_read(char *fname, int *size)
    {
    int         i, xnum, ynum, status=0;
    long        nelements, fpixel[2];
    char        err_text[81];
    float       *data;
    fitsfile    *p=NULL;

//
// Get Actual size of data in file.  Don't set astro_errno here because 
//   fits_dims() will set it appropriately if there is a problem.
//

    if (fits_dims(fname, &xnum, &ynum))
        {
        if (astro_warn) printf("WARNING: astro::fits_read:fits_dims() Error\n");
        return(NULL);
        }

        if (DEBUG) printf("DEBUG: astro::fits_read:dims xnum=%d, ynum=%d\n",xnum,ynum);

//
// Open the file for reading
//

    if (fits_open_file(&p, fname, READONLY, &status)) 
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_read:fits_open_file() Error %d: %s\n",status,err_text);
        set_astro_errno(ASTRO_ERR_OPEN);
        return(NULL);
        }
         
    if (DEBUG) printf("DEBUG: astro::fits_read:fits_file_open\n");

//
// Allocate memory for the data reading buffer based on the size dimensions
//

    if ((data=(float *)malloc((xnum*ynum)*sizeof(float *)))==NULL)
        {
        if (astro_warn) printf("WARNING: astro::fits_read:malloc() Error\n");
        fits_close_file(p, &status);
        set_astro_errno(ASTRO_ERR_MALLOC);
        return(NULL);
        }

    for (i=0; i< (xnum*ynum) ; i++) data[i]=0.0;

    if (DEBUG) printf("DEBUG: astro::fits_read:malloc %zu (%d x %d x float size) bytes\n",xnum*ynum*sizeof(float),xnum,ynum);

    fpixel[0]=fpixel[1]=(long) 1;
    nelements=(long)(xnum*ynum);

    if (fits_read_pix(p, TFLOAT, fpixel, nelements, NULL, data, NULL, &status))
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_read:fits_read_pix() Error %d: %s\n",status,err_text);
        free(data);
        fits_close_file(p, &status);
        set_astro_errno(ASTRO_ERR_READPIX);
        return(NULL);
        }

    if (fits_close_file(p, &status))
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_read:fits_close_file() Error %d: %s\n", status, err_text);
        free(data);
        set_astro_errno(ASTRO_ERR_CLOSE);
        return(NULL);
        }
    *size=xnum*ynum;
    return(data);
    }


//
// FITS_WRITE() - Write an image to a FITS file.  The image can be written to
//                an existing file or will create a new file.
//
// Arguments:
//      fname   - Text filename for FITS file to be written
//      data    - Image data to be written to the file
//      x_size  - Number of rows in the image
//      y_size  - Number of cloumns in the image
//      newfile - Flag to create a file.  If 1, create a new file (overwrite
//                any existing one.  If 0, will write to an existing file.
//      pname   - Text string for program name in PROGRAM key value
//      version - Test string with version information for PROGRAM key value
//
// Return Value:
//      ASTRO_SUCCESS - Write was successful
//      ASTRO_FAILURE - Write failed
//
// Errors:  Function will set astro_errno with return code (see astro_class.h)
//

int    astro::fits_write(char *fname, float *data, int x_size, int y_size, 
                         int newfile, const char *pname, const char *version)
    {
    int         status=0;
    int         checkx;
    int         checky;
    long        naxes[2];
    char        key[81];
    char        err_text[81];
    fitsfile    *fptr=NULL;
    long        fpixel[2]= { 1, 1};

    naxes[0]=(long) x_size;
    naxes[1]=(long) y_size;

    if (x_size < MIN_FITS || y_size < MIN_FITS || x_size > MAX_FITS || y_size > MAX_FITS)
        {
        if (astro_warn) printf("WARNING: astro::fits_write: Image Size Invalid\n");
        set_astro_errno(ASTRO_ERR_WRITE);
        return(ASTRO_FAILURE);
        }

    if (newfile)
        {
//
// Create a new FITS file
//

        if (fits_create_file(&fptr,fname, &status))
            {
            fits_get_errstatus(status,err_text);
            if (astro_warn) printf("WARNING: astro::fits_write:fits_create_file() Error %d: %s\n",status,err_text);
            set_astro_errno(ASTRO_ERR_CREATE);
            return(ASTRO_FAILURE);
            }

//
// Create a new FITS image unit in the file
//
 
       if (fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status))
            {
            fits_get_errstatus(status,err_text);
            if (astro_warn) printf("WARNING: astro::fits_write:fits_create_img() Error %d: %s\n",status,err_text);
            set_astro_errno(ASTRO_ERR_IMAGE);
            return(ASTRO_FAILURE);
            }
        }
    else
        {
//
// Check the size of the primary image unit
//

        if (fits_dims(fname, &checkx, &checky)) return(ASTRO_FAILURE);

        if ((checkx != x_size) || (checky != y_size))
            {
            if (astro_warn) printf("WARNING: astro::fits_write: Exist File Size Error %d:%d-%d:%d\n",x_size,checkx,y_size,checky);
            set_astro_errno(ASTRO_ERR_SIZE);
            return(ASTRO_FAILURE);
            }

//
// Open the file for read/write
//

        if (DEBUG) printf("DEBUG: astro::fits_read:fits_file_open\n");

        if (fits_open_file(&fptr, fname, READWRITE, &status))
            {
            fits_get_errstatus(status,err_text);
            if (astro_warn) printf("WARNING: astro::fits_write:fits_open_file() Error %d: %s\n",status,err_text);
            set_astro_errno(ASTRO_ERR_OPEN);
            return(ASTRO_FAILURE);
            }
        }

//
// Write the image to the FITS files
//

    if (fits_write_pix(fptr, TFLOAT, fpixel, long(x_size*y_size), data, &status))
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_write:fits_write_pix() Error %d: %s\n", status, err_text);
        set_astro_errno(ASTRO_ERR_WRITE);
        return(ASTRO_FAILURE);
        }

    sprintf(key,"HDU Created by %s/%s - %s", MAJOR_VERSION, pname, version);
    if (fits_write_key(fptr, TSTRING, "PROGRAM", key, NULL, &status))
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_write:fits_write_key() Error %d: %s\n", status, err_text);
        set_astro_errno(ASTRO_ERR_KEY);
        return(ASTRO_FAILURE);
        }

    if (fits_close_file(fptr, &status))
        {
        fits_get_errstatus(status,err_text);
        if (astro_warn) printf("WARNING: astro::fits_write:fits_close_file() Error %d: %s\n", status, err_text);
        set_astro_errno(ASTRO_ERR_CLOSE);
        return(ASTRO_FAILURE);
        }

    return(ASTRO_SUCCESS);
    }


//
// CarrayAlloc() - This function will dynamically allocate a 2D character 
//                 array.  This is needed because we need to dynamically
//                 allocate an array for many FITS functions.  It has to be 
//                 the exact size (required by many CFITSIO routines) and 
//                 needs to be a monolithic block of data (which is why you
//                 can't use the traditional series of mallocs per row to
//                 create it).
//
//                 PLEASE NOTE:  This function uses the C/C++ sense of
//                 indices and not the FITS sense (which is opposite).
//
// Arguments:
//      crows    - Number of rows (X dimension, slowest changing index)
//      ccolumns - Number of columns (Y dimension, fastest changing index)
//
// Return Value:
//      char **  - Pointer to the base of the allocated 2D array
//
// Errors:  Function will set astro_errno with return code (see astro_class.h)
//

char     **astro::CArrayAlloc(int crows, int ccolumns) 

    {
    int a;
    char **cptr;

//
// Double index arrays have a header containing pointer to each row
//

    int header = crows * sizeof(char*);
    int body = 0;

    if (DEBUG) printf("DEBUG: astro::ArrayAlloc (C): crows=%x, ccolumns=%x, header=%x\n",crows,ccolumns,header);

    for(a=0; a<crows; a++) body+=ccolumns;
    body*=sizeof(char);

    if (DEBUG) printf("DEBUG: astro::ArrayAlloc (C): body=%x\n",body);

//
// Allocate the total space needed
//

    if ((cptr = (char**)malloc(header + body))==NULL)
        {
        if (astro_warn) printf("WARNING: astro::ArrayAlloc:malloc() error\n");
        set_astro_errno(ASTRO_ERR_MALLOC);
        return(NULL);
        }

    if (DEBUG) printf("DEBUG: astro::ArrayAlloc (C): cptr=%p\n",cptr);

    char* buf  = (char*)(cptr + crows);

    if (DEBUG) printf("DEBUG: astro::ArrayAlloc (C): buf=%p\n",buf);

//
// Now populate all the row pointers with the correct values
//

    cptr[0] = buf;
    for(a = 1; a < crows; ++a) 
        {
        cptr[a] = cptr[a-1] + ccolumns;
        if (DEBUG) printf("DEBUG: astro::ArrayAlloc (C): cptr[%d]=%p\n",a,cptr[a]);
        }
    return(cptr);
    }


//
// arrayAlloc() - This function will dynamically allocate a 2D float 
//                array.  This is needed because we need to dynamically
//                allocate an array for many FITS functions.  It has to be 
//                the exact size (required by many CFITSIO routines) and 
//                needs to be a monolithic block of data (which is why you
//                can't use the traditional series of mallocs per row to
//                create it).
//
//                PLEASE NOTE:  This function uses the C/C++ sense of
//                indices and not the FITS sense (which is opposite).
//
// Arguments:
//      frows    - Number of rows (X dimension, slowest changing index)
//      fcolumns - Number of columns (Y dimension, fastest changing index)
//
// Return Value:
//      float ** - Pointer to the base of the allocated 2D array
//
// Errors:  Function will set astro_errno with return code (see astro_class.h)
//

float   **astro::ArrayAlloc(int frows, int fcolumns) 

    {
    int     a;
    float   **fptr;

//
// Double index arrays have a header containing pointer to each row
//

    int header = frows * sizeof(float*);
    int body = 0;

    for(a=0; a<frows; a++) body+=fcolumns;
    body*=sizeof(float);

//
// Allocate the total space needed
//

    if ((fptr = (float**)malloc(header + body))==NULL)
        {
        if (astro_warn) printf("WARNING: astro::ArrayAlloc:malloc() error\n");
        set_astro_errno(ASTRO_ERR_MALLOC);
        return(NULL);
        }

    float* buf  = (float*)(fptr + frows);

//
// Now populate all the row pointers with the correct values
//

    fptr[0] = buf;
    for(a = 1; a < frows; ++a) 
        {
        fptr[a] = fptr[a-1] + fcolumns;
        }
    return(fptr);
    }
