/*****
 *
 *
 *  common.h common definitions header file
 *
 *  Author: Johan Lofgren
 *
 *
 **/

#ifndef _COMMON_H
#define _COMMON_H

#include "arb.h"
//#include "flint/profiler.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <algorithm>
//#include <cstdint>
#include <cstring>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include "lookup.h"

//---------------------Types---------------------------------------
#define float_t double //change to long double if necessary
//#define size_t unsigned int

const float_t kappa = 0.00633257397764611; // 1/(64 pi^2)
const float_t EulerGamma = 0.5772156649015328606;
//const float_t EulerNumber = exp(1);
const float_t VEV = 246.0;
const float_t VEV2 = VEV*VEV;

//typedef int bool;
#define true 1
#define false 0
//-----------------------------------------------------------------

struct minOptions_t{
	minOptions_t() : stepSize(1.0), maxTemp(500), minTemp(0.0), tol(1e-2), locTol(1e-0), tempTol(1e-2), fieldMax(500.0), maxIter(200) {}

	float_t stepSize;
	float_t maxTemp;
	float_t minTemp;
	float_t tol;
	float_t locTol;
	float_t tempTol;
	float_t fieldMax;

	size_t maxIter;
};

void printtime(clock_t s, clock_t e);

//Integrals:
float_t integralInterpol(const float_t x2, const table_t* lookUpTable, const int n);

size_t Minimizer(minOptions_t options, gsl_multimin_fminimizer *s,gsl_vector* min, size_t numvar);

#endif
