/*****
 *
 *
 *  funcs.h functions header file
 *
 *  Author: Johan Löfgren
 *
 *
 **/

#ifndef _FUNCS_H
#define _FUNCS_H

#include "common.h"
#include "AbstractModel.h"

//Numerical functions:
double my_totPotFunction (const gsl_vector *v, void *params);
double my_deltaPotFunction (double T, void *params);
double my_deltaPotResummedFunction (double T, void *params);


//PTA methods:
void whenPT_PRM(Model* model, float_t& Tc, float_t& vc, size_t& errorFlag);
void whenPT_PRM_R(Model* model, float_t& Tc, float_t& vc, size_t& errorFlag);
void whenPT_L(Model* model, float_t& Tc, float_t& vc, size_t& errorFlag, minOptions_t options);
void whenPT_L_FAST(Model* model, float_t& Tc, float_t& vc, size_t& errorFlag, minOptions_t options);
int whichPhase(gsl_multimin_fminimizer *s,gsl_vector* ss, gsl_vector* min, gsl_multimin_function my_func, Model* model, size_t& errorFlag, minOptions_t options);
int whichPhase2(gsl_multimin_fminimizer *s,gsl_vector* ss, gsl_vector* min, gsl_multimin_function my_func, Model* model, size_t& errorFlag, minOptions_t options);



#endif
