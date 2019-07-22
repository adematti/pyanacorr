///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of CUTE.                                        //
//                                                                   //
// CUTE is free software: you can redistribute it and/or modify it   //
// under the terms of the GNU General Public License as published by //
// the Free Software Foundation, either version 3 of the License, or //
// (at your option) any later version.                               //
//                                                                   //
// CUTE is distributed in the hope that it will be useful, but       //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU  //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CUTE.  If not, see <http://www.gnu.org/licenses/>.     //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#ifndef _WINDOW_DEFINE_
#define _WINDOW_DEFINE_

#define MAX(a, b) (((a) > (b)) ? (a) : (b)) //Maximum of two numbers
#define MIN(a, b) (((a) < (b)) ? (a) : (b)) //Minimum of two numbers
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x))) //min(max(a,low),high)
//#define ABS(a)   (((a) < 0) ? -(a) : (a)) //Absolute value
#define MAX_ELLS 9
#define MAX_RADIALS 4
#define MAX_POLES 2
#define MAX_LOS 2

#ifdef _FLOAT32
typedef float histo_t;
#else
typedef double histo_t;
#endif //_FLOAT32

typedef enum {QUIET, INFO, DEBUG} VERBOSITY;
typedef enum {LIN, POLY} INTERPOL;
typedef enum {LOS_MIDPOINT, LOS_ENDPOINT, LOS_FIRSTPOINT} LOS_TYPE;
typedef enum {MULTI_ALL, MULTI_EVEN, MULTI_ODD} MULTI_TYPE;

//Selection
typedef struct {
	histo_t *x;
	histo_t *y;
	size_t *shape;
	size_t n_dim;
	size_t size;
} Selection;

typedef struct {
	size_t n_ells;
	MULTI_TYPE type;  
} Pole;

typedef struct {
	LOS_TYPE type;
	size_t n;
} LOS;

VERBOSITY verbose;

#endif //_WINDOW_DEFINE_
