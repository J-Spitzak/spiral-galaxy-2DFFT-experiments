//
// PITCH_CLASS.CPP - This class provides functions for analysis and ordering
//                   of the FFT output data from P2DFFT.
//
//
// Version 1.3  07-Apr-2018
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
//      1.3  07-Apr-2018: - Change snr() and fwhm() to set calculated values in
//                          the return structure and just return a staus code
//      1.2  16-Mar-2018: - Fix bug due to FP rounding error in SNR
//                        - Add warn flag and set_warn() function, and add
//                          checks for this flag on print statements
//                        - Add prefix to print statements (ERROR, WARNING)
//      1.1  18-Feb-2018: - Add snr() and fwhm() functions
//                        - Fix small bug in pitch_phase()
//                        - Remove DC component from calculations
//                        - Add return codes from NaN and errors
//                        - Add include file version number in print
//                        - Clean up error messaging
//      1.0  05-Feb-2018: - Initial version
//

#define     PITCH_VER   "1.3/20180407"

#include    <stdio.h>
#include    <string.h>
#include    <unistd.h>
#include    <fstream>
#include    <sstream>
#include    <magic.h>
#include    <math.h>
#include    <sys/stat.h>
#include    <sys/types.h>

#include    "pitch_class.h"
#include    "globals.h"

int         pitch_warn=0;

//
// Define the error variable macro to set the error
//

int     pitch_errno=0;

#define set_pitch_errno(err) (pitch_errno = (err))

//
// CONSTANTS -- These must match the values used in the P2DFFT algorithms
//

#define LO_INDEX    824
#define HI_INDEX    1226

//
// FUNCTION BLOCK
//


//
// SET_WARN() - Sets the value of the warning flag which controls the
//              printing of warning messages
//
// Arguments:
//      value   - 0 for no warnings, non-zero for warnings
//
// Return Value: NONE
//

void    pitch::set_warn(int value)
    {
    pitch_warn=value;
    }


//
// VERSION() - This function will print the current version.
//
// Arguments: NONE
//
// Return Value: NONE
//

void    pitch::version()
    {
    printf("  -- Pitch Class Include Version:  %s\n",PITCH_H_VER);
    printf("  -- Pitch Class Function Version:  %s\n",PITCH_VER);
    }


//
// GET_ERR() - This function will return the lastest error number
//
// Arguments: NONE
//
// Return Value: Most recent error code defined in pitch_class.h
//

int    pitch::get_err()
    {
    return(pitch_errno);
    }


//
// PITCH_PHASE() - This public function will parse the FFT data and determine
//                 the highest frequency and amplitude.  It will also 
//                 calculate the corresponding pitch and phase angles.
//
//                 This routine replace the give_maximum_pitch_phase.f program
//
// Arguments:
//      fft     - Pointer to array of FFT output data
//      res     - Structure For return information
//
// Return Value:
//      PITCH_RET_OK       - Processing ended normally, results in res valid`
//      PITCH_RET_NAN      - Processing returned NaN (low to no signal)
//      PITCH_RET_ERR      - Error encountered, no results returned
//

int    pitch::pitch_phase(fft_out *fft, int mode, result_pa *res)
    {
    int     i;
    int     nan=1;
    int     index;

    double  a_max;

//
// Find Highest Amplitude/Frequency
//

    a_max=-255.0;
    index=-1;

    for (i=LO_INDEX; i <= HI_INDEX; i++)
        {
//
// First check if any number is not a NaN.  If the floating point returned NaN
//   values, it will not work.   We check by simply testing the number.  This
//   uses the fact that in C/C++ any comparison on a NaN will always be false
//   (even f == f).
//

        if (fft[i].abs == fft[i].abs) nan=0;
        
        if ((fft[i].abs > a_max) && (i != 1025))
            {
            index=i;
            a_max=fft[i].abs;
            }
        }

    if (DEBUG) printf("DEBUG: Max Amp %f, Index=%d\n", a_max, index);

    if (nan) return(PITCH_RET_NAN);
    
    if (index < 0)
        {
        if (pitch_warn) printf("WARNING: Can't locate maximum amplitude\n");
        set_pitch_errno(PITCH_ERR_MAX_AMP);
        return(0);
        }

    res->amp=fft[index].abs;
    res->freq=fft[index].freq;
    res->index=index;

//
// Calculate the pitch and phase angles
//

    res->pa=atan2((double)mode,fft[index].freq)*(1.0/GR_RAD);

    if (fabs(res->pa) > 90.0) res->pa-=180.0;

    res->phase=atan2(fft[index].imag,fft[index].real)*(1.0/GR_RAD)/mode;
    return(1);
    }


//
// SNR() - Returns the signal to noise value for a FFT data block.
//
// Arguments:
//      fft     - Pointer to array of FFT output data
//      res     - Return information from pitch phase
//
// Return Value:
//      PITCH_RET_OK       - Processing ended normally, results in res valid`
//      PITCH_RET_NAN      - Processing returned NaN (low to no signal)
//      PITCH_RET_ERR      - Error encountered, no results returned
//

int    pitch::snr(fft_out *fft, result_pa *res)
    {
    int     i;

    double  L=0.0;
    double  sigma;
    double  avg=0.0;
    double  vals=0.0;

//
// Find the average value, L, and set it in the return structure
//

    for (i=LO_INDEX; i <= HI_INDEX; i++)
        {
        if ((i != 1025) && (fft[i].abs == fft[i].abs))
            {
            avg+=fft[i].abs;
            vals+=1.0;
            }
        }

    if (vals < 0.5)
        {
        set_pitch_errno(PITCH_ERR_ALLNANS);
        return(PITCH_RET_ERR);
        }

    L=avg/vals;
    res->avg_amp=L;

    avg=0.0;

//
// Determine the SNR/sigma
//

    for (i=LO_INDEX; i <= HI_INDEX; i++)
        {
        if ((i != 1025) && (fft[i].abs == fft[i].abs)) avg+=pow((fft[i].abs-L),2.0);
        }
 
    sigma=pow((avg/vals),0.5);

    if (sigma <= 1e-10)
        {
        set_pitch_errno(PITCH_ERR_SIGMA);
        return(PITCH_RET_ERR);
        }

    res->snr=(res->amp-L)/sigma;
 
    if (DEBUG) printf("DEBUG: SNR=%g, Sigma=%g, L=%g\n", res->snr, sigma, L);

    if (res->snr == res->snr) return(PITCH_RET_OK);
    else return(PITCH_RET_NAN);
    }


//
// FWHM() - Returns the FWHM for an FFT data block.
//
//         NOTE: pitch_phase() and snr() must have been called and populated the
//               res data structure (with valid data) before calling this
//               routine.
//
// Arguments:
//      fft     - Pointer to array of FFT output data
//      res     - Return information from pitch phase
//
// Return Value:
//      PITCH_RET_OK       - Processing ended normally, results in res valid`
//      PITCH_RET_NAN      - Processing returned NaN (low to no signal)
//      PITCH_RET_ERR      - Error encountered, no results returned
//


int    pitch::fwhm(fft_out *fft, result_pa *res)
    {
    int     i;
    int    lo=0;
    int    hi=0;

    double  limit;

//
// Check the data in res block.  Ideally this should never happen, but we check
//   to make sure.
//

    if ((res->index < LO_INDEX) || (res->index > HI_INDEX))
        {
        if (pitch_warn) printf("WARNING: Invalid data in res block\n");
        set_pitch_errno(PITCH_ERR_INVALID);
        return(PITCH_RET_ERR);
        }

//
// Track the slope of the right side (positive direction) of the peak.  Stop
//   when the slope is no longer negative or noise floor is encountered.
//

    limit=res->amp - ((res->amp - res->avg_amp)/2.0);
    i=res->index+1;

    while (i <= HI_INDEX)
        {
        if (DEBUG) printf("DEBUG: Process Index %d, ABS=%f, LIMIT=%f\n",i,fft[i].abs,limit);
        if (i != 1025)
            {
            if (fft[i].abs < limit)
                {
                hi=i-1;
                break;
                }
            }
        i++;
        }

    i=res->index-1;

    while (i >= LO_INDEX)
        {
        if (DEBUG) printf("DEBUG: Process Index %d, ABS=%f, LIMIT=%f\n",i,fft[i].abs,limit);
        if (i != 1025)
            {
            if (fft[i].abs < limit)
                {
                lo=i+1;
                break;
                }
            }
        i--;
        }

    if (DEBUG) printf("DEBUG: Hi=%d, Lo=%d\n",hi,lo);

    if ((hi==0) || (lo==0))
        {
        set_pitch_errno(PITCH_ERR_SCANFWHM);
        return(PITCH_RET_ERR);
        }

    res->fwhm=(double)(hi - lo + 1);

    if (res->fwhm == res->fwhm) return(PITCH_RET_OK);
    else return(PITCH_RET_NAN);
    }

