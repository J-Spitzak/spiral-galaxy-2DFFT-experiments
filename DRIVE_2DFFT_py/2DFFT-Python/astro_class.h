//
// ASTRO.CPP - This class provides functions for astronomical analysis,
//             specifically the 2DFFT pitch angle analysis and hurricane
//             packages.
//
//
// Version 2.0: 26-May-2018
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
//      2.0  26-May-2018: - Add fits_write() function
//                        - Add new error codes
//                        - Add return constants
//                        - Remove placeholders for fits_text_read/write
//                        - Add error code for invalud buffer size
//      1.2  16-Mar-2018: - Add set_warn() function
//                        - Add version constant
//      1.1  03-Feb-2018: - Update comments
//                        - Add version print function
//                        - Add get_err function
//                        - Remove astro_print declarations
//                        - Remove external astro_errno declarations (now it
//                          is a private variable)
//      1.0  17-Feb-2017: - Initial version
//

#define     ASTRO_H_VER     "2.0/20180526"

#include    <cstddef>
#include    <iostream>
#include    <string>
#include    <vector>

#include    <fitsio.h>

//
// Data structure used for the file parameters
//

struct  file_rec
    {
    int             valid;      /* Flag to indicate if data is valid         */
    std::string     name;       /* File name                                 */
    std::string     keyword;    /* Name of output file per radius            */
    std::string     result;     /* Prefix for overall output files           */
    int             radius;     /* Outer radius value                        */
    int             binary;     /* Is binary (1) FITS or ASCII text FITS (0) */
    };        

//
// Class definition values
//

class   astro   {
                public:
                    void    set_warn(int value);
                    int     get_err();
                    void    version();
                    int     file_type(std::string fname);
                    bool    file_exists(const char *fname);
                    int     fits_dims(std::string fname, int *rows, int *cols);
                    char    **fits_header_read(char *fname, int *keys);
                    int     fits_header_write(char *fname, char keys[][32], char items[][80], int num);
                    float  *fits_read(char *fname, int *size);
                    int    fits_write(char *fname, float *data, int x_size, int y_size, int newfile, const char *pname, const char *version);
                    char   **CArrayAlloc(int crows, int ccols);
                    float  **ArrayAlloc(int frows, int fcols);
                    int    read_lines(std::string fname, std::vector<file_rec> *rec);
                };

//
// astro_class file type return values
//

#define     ASTRO_TXT_FILE      0
#define     ASTRO_BIN_FILE      1
#define     ASTRO_UNK_FILE      -1

//
// astro_class error number definitions
//

#define     ASTRO_ERR_KEY       1025
#define     ASTRO_ERR_OPEN      1026
#define     ASTRO_ERR_SIZE      1027
#define     ASTRO_ERR_CLOSE     1028
#define     ASTRO_ERR_WRITE     1029
#define     ASTRO_ERR_IMAGE     1030
#define     ASTRO_ERR_CREATE    1031
#define     ASTRO_ERR_MALLOC    1032
#define     ASTRO_ERR_RD_REC    1033
#define     ASTRO_ERR_HDR_POS   1034
#define     ASTRO_ERR_DIMSIZE   1035
#define     ASTRO_ERR_READPIX   1036
#define     ASTRO_ERR_HOMEDIR   1037
#define     ASTRO_ERR_GET_SIZE  1038

//
// astro_class return codes
//

#define     ASTRO_SUCCESS       0
#define     ASTRO_FAILURE       1
