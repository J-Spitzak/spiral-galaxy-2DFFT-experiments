#
#  Makefile for P2DFFT Suite of Programs
#
#           The program has been extensivley modified at the NC Museum of
#           Natural Sciences from previous 2DFFT versions so that it will
#           support all the changes in the P2DFFT project.
#
#  Version 5.1  20-Jun-2019
#
#  2DFFT (original) Author: Dr. Ivanio Puerari
#                           Instituto Nacional de Astrofisica,
#                           Optica y Electronica, 
#                           Santa Maria Tonantzintla,
#                           Puebla, Mexico
#
#  2DFFT (revised) Lead Author: Dr. Marc Seigar
#                               University of Minnesota Duluth,
#                               Duluth, MN USA
#
#  2DFFT (progenitor of P2DFFT) Lead Author: Dr. Benjamin Davis
#                                            Swinburne University of Technology.
#                                            Centre for Astrophysics and
#                                            Supercomputing
#                                            Melbourne, Victoria, Australia
#                                        http://d.umn.edu/~msseigar/2DFFT.html
#
#  P2DFFT By: Ian Hewitt & Dr. Patrick Treuthardt, 
#             NC Museum of Natural Sciences, 
#             Astronomy & Astrophysics Lab,
#             Raleigh, NC USA.
#             http://github.com/treuthardt/P2DFFT
#
#
#  LICENSE
#
#  P2DFFT Spiral Galaxy Arm Pitch Angle Analysis Suite
#  Copyright (c) 2016-2019  Ian B. Hewitt & Dr. Patrick Treuthardt
#
#  The program is free software:  you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by the Free
#  Software Foundation, version 3.
#
#  This program is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
#  details.
#
#  You should have received a copy of the GNU General Public License along with
#  this program.  If not, see < https://www.gnu.org/licenses >.
#
#  The authors can be contacted at:
#
#      North Carolina Museum of Natural Sciences
#      Astronomy & Astrophysics Laboratory
#      11 West Jones Street
#      Raleigh, NC, 27601  USA
#      +1.919.707.9800
#
#      -- or --
#
#      patrick.treuthardt@naturalsciences.org
#
#
#  Revision History:
#
#       5.1 20-Jun-2019 - Update comments to identify correct original author
#                       - Update for filename changes
#                       - Clarify licensing/contact information
#       5.0 21-Feb-2019 - Remove textify.cl and add zname.py
#                       - Add/remove filenames for v5 dist
#       4.9 29-Dec-2018 - Fix dist rule to remove deprected files
#       4.8 30-Oct-2018 - Minor change to add concept of patch leve vi RM number
#                       - Remove warning flag
#                       - Update install list for merger of files into pa.py
#                       - Remove deprecated files
#       4.7 20-Oct-2018 - Add rules for zoo.py
#       4.6 14-Sep-2018 - Change rules for spiral.c to be C++ and use 
#                         astro_class
#                       - Remove deprecated rules
#                       - Add generic variables for items needing classes
#                       - Updated rules for ifft (now C++)
#                       - Added rules for map
#                       - Added curl libraries to support new files
#                       - Remove deprecated programs
#                       - Add build of optional files seperately via "opt" and
#                         "optinstall" rules
#       4.5 26-Mar-2018 - Add dependencies on globals.h
#                       - Change some parameters to improve multithreading
#                         performance of FFTW
#       4.4 03-Feb-2018 - Change rules to build p2dfft.cpp instead of 2dfft.cpp
#                       - Remove rules for deprecated files
#                       - Add pitch_class.cpp into build rules
#                       - Add include files to build prerequisite checks
#       4.3 21-Jan-2018 - Add prequisite on install rule to automatically cause
#                         a build if needed.
#       4.2 01-Dec-2017 - Add dist rule to package the files into a tgz file
#       4.1 28-Aug-2017 - Update compiler flags and libraries to use FFTW3 and
#                         remove fourn.c references
#                       - Change ordering of libraries to prevent some
#                         distributions having compiler problems (like Ubuntu)
#                       - Update to change ifft to C++
#       4.0 23-Jul-2017 - Add blank gif files to install step
#       3.9 09-Jun-2017 - Add pthread library to linking step to support 
#                         multithreaded/re-entrant mode in cfitsio
#                       - Remove functions related to test package (now a 
#                         seperate package)
#       3.8 07-Jun-2017 - CRITICAL Fix:  Remove -std=c++11 compiler flag
#                         because it can lead to erroneous results on some
#                         Linux installations
#       3.7 21-Feb-2017 - Change default compilation option from debug (-g) to
#                         optimize (-O)
#                       - Add comment about not commenting BIN_DIR
#                       - Add missing BIN_DIR variable pass for C++ files
#       3.6 08-Feb-2017 - Update to remove scripter and change 2dfft to C++
#       3.5 05-Jan-2017 - Fix bug in CFLAGS definition in makefile
#       3.4 14-Dec-2016 - add avg.py in install portion of makefile
#       3.3 07-Dec-2016 - Fix bug in build order which would cause the first
#                         time make was typed to fail, but work the second
#       3.2 06-Dec-2016 - Add spiral (spiral generation program) to makefile
#                       - Add txt2fits program to makefile
#                       - Add ifft program to makefile
#       3.1 27-Nov-2016 - Change CFLAGS so the value for BIN_DIR is passed to
#                         the program being compiled
#       3.0 20-Sep-2016 - Modified to add testing data and test scripts to
#                         validate any changes in the code from the original
#       2.2 06-Sep-2016 - Move from fort77 to gfortran compiler
#       2.1    May-2016 - Add comments for all compiler flags and move them to
#                         a single variable
#       2.0    May-2016 - Modified to add parallel processing version and unify
#                         all build files into a single makefile
#

#
# Change BIN_DIR to the location where you will be installing the binary files.
#   This line should a single directory (no comments).
#

BIN_DIR = ~/bin

#
# Change VERSION if you will be making code changes and want to packge them
#   for distribution.  Recommend doing something like 4.2.x to identify the
#   base version that was used before changes.
#

VERSION = 5.2.2

CFLAGS = -O -DBIN_DIR='"$(BIN_DIR)"'
CCFLAGS = -O -DBIN_DIR='"$(BIN_DIR)"' -fopenmp
LIBS = -lmagic -lcfitsio -lfftw3 -lcurl -lpthread -lm
ASTRO = astro_class.cpp astro_class.h
PITCH = pitch_class.cpp pitch_class.h

all: p2ifft p2dfft p2spiral 

opt: p2txt2fits p2map

install: all
	mkdir -p $(BIN_DIR)
	cp p2dfft p2spiral p2ifft p2pa p2zoo p2zname $(BIN_DIR)

optinstall: opt
	mkdir -p $(BIN_DIR)
	cp p2boost p2logsp p2txt2fits p2map p2filter p2chart_freq $(BIN_DIR)

clean:
	rm -f *.o *.a core p2dfft p2spiral p2txt2fits p2ifft p2map

dist:
	tar czvf ../p2dfft-$(VERSION).tgz README.* CHANGES makefile* *.cpp *.h *.c GNU*  PA_Notes.* input.txt sp_input.txt p2pa p2zoo p2zname p2boost p2logsp p2filter p2chart_freq

p2dfft: p2dfft.cpp $(ASTRO) $(PITCH) globals.h
	g++ $(CCFLAGS) -o p2dfft p2dfft.cpp astro_class.cpp pitch_class.cpp $(LIBS) -fopenmp
	rm -f *.o

p2spiral: p2spiral.cpp $(ASTRO) globals.h
	g++ $(CCFLAGS) -o p2spiral p2spiral.cpp astro_class.cpp $(LIBS)
	rm -f *.o

p2ifft: p2ifft.cpp astro_class.cpp astro_class.h globals.h
	g++ $(CCFLAGS) -o p2ifft p2ifft.cpp astro_class.cpp $(LIBS)

p2map: p2map.cpp astro_class.cpp astro_class.h globals.h
	g++ $(CCFLAGS) -o p2map p2map.cpp astro_class.cpp $(LIBS)

.c: globals.h
	cc -o $* $(CFLAGS) $*.c $(LIBS)
	rm -f *.o

.cpp: globals.h
	g++ $(CFLAGS) -o $* $*.cpp $(LIBS)
	rm -f *.o

