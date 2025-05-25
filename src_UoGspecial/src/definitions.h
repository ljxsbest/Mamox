#pragma once

#include "Constants.h"
// in this file where the definitions for compile time are included

//#define QUICK_TEST

#define SIZEMAX 20

#define EQUILIBRIUM_MECH
//#define OPTIMIZED_RATE_RULES

#define MAX_TIME 100.0
#define TOL_GLOB_MULTI	1.0e-8
#define TOL_GLOB_SIN	1.0e-8
#define TOL_LOC_MULTI	1.0e-15
#define TOL_LOC_SIN		1.0e-9

#define EQUIVALENCE_RATIO 1.0

#ifndef QUICK_TEST
	//#define TEMPERATURES    600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800//, 1000, 1200
	#define TEMPERATURES    600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050 , 1100, 1150, 1200
	#define PRESSURES		1.0, 2.0, 4.0, 8.0, 16.0, 32.0
#endif // !QUICK_TEST
#ifdef QUICK_TEST
	#define TEMPERATURES    1000, 1200
	#define PRESSURES		1.0
#endif // QUICK_TEST



#define MAX_NUM_CORES 6