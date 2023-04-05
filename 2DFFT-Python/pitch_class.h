//
// PITCH_CLASS.H - This class provides functions for interpreting the FFT
//                 results from P2DFFT.
//
//
// Version 1.3: 07-Apr-2018
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
//      1.3  07-Apr-2018: - Change to put snr and fwm in the result_pa struct
//      1.2  16-Mar-2018: - Add get_warn() function
//      1.1  17-Feb-2018: - Add snr() and fwhm() functions
//                        - Update error code
//                        - Add return code definitions
//                        - Add include file version number
//      1.0  04-Feb-2018: - Initial version
//

#define     PITCH_H_VER   "1.3/20180407"

#include    <cstddef>
#include    <iostream>
#include    <string>

//
// Data structure used for the parsed FFT output
//

struct  fft_out
    {
    double      real;      /* Real component of FFT output              */
    double      imag;      /* Imaginary component of FFT output         */
    double      abs;       /* Absolute value of the complex number      */
    double      freq;      /* Frequency related to this complex number  */
    };

//
// Data structure for FFT basic output analysis
//

struct  result_pa
    {
    int         index;     /* Index of highest amplitude                */
    double      freq;      /* Frequency of highest amplitude            */
    double      amp;       /* Highest amplitude                         */
    double      avg_amp;   /* Average amplitude (noise level)           */
    double      pa;        /* Calculated pitch angle                    */
    double      phase;     /* calculated phase angle                    */
    double      snr;       /* calculated signal-to-noise ratio          */
    double      fwhm;      /* calculated full width half maximum        */
    };

//
// Class definition values
//

class   pitch {
              public:
                 void    set_warn(int value);
                 void    version();
                 int     get_err();
                 int     pitch_phase(fft_out *fft, int mode, result_pa *res);
                 int     snr(fft_out *fft, result_pa *res);
                 int     fwhm(fft_out *fft, result_pa *res);
              };

//
// Return codes.  The tan(2) function can return NaN for items with low to no
//                signal, which is not a real error.
//

#define     PITCH_RET_OK        1
#define     PITCH_RET_NAN       0
#define     PITCH_RET_ERR      -1

//
// Error Values
//

#define     PITCH_ERR_INVALID   2049
#define     PITCH_ERR_MAX_AMP   2050
#define     PITCH_ERR_ALLNANS   2051
#define     PITCH_ERR_SIGMA     2052
#define     PITCH_ERR_SCANFWHM  2053

