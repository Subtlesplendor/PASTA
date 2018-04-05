#include "common.h"

//Used to check performance (remove before making public):
void printtime(clock_t s, clock_t e)
{
  printf("runtime: %f seconds \n", (double)(e-s)/CLOCKS_PER_SEC);
}

float_t integralInterpol(const float_t x, const table_t* lookUpTable, const int n)
{
    //int n = length;

    float_t xmin = lookUpTable[0].x;
    float_t xmax = lookUpTable[n-1].x;
    float_t h = (xmax-xmin)/((double)n-1);

    //float_t x = sqrt(x2);
    //printf("x = %f\n", x);

    int i;

    if (x==xmin)
    {
        i = 0;
    }
    else
    {
        float_t index = (x-xmin)/h;
        i = floor(index);
    }

    //printf("i = %i\n",i);
    float_t xdiff= x - lookUpTable[i].x;
    //printf("xdiff = %f,\n",xdiff);



    float_t result = lookUpTable[i].y + lookUpTable[i].yp * xdiff + lookUpTable[i].a * xdiff * xdiff + lookUpTable[i].b * xdiff * xdiff * xdiff;

    //printf("integral = %f \n", result);

    return result;

}


size_t Minimizer(minOptions_t options, gsl_multimin_fminimizer *s,gsl_vector* min, size_t numvar)
{
	//Minimizer minimizes a function using a simplex algorithm from GSL.

	int status;
	double size;
	size_t iter = 0;
	size_t itermax = options.maxIter;
	float_t fieldMax = options.fieldMax;
	float_t tol = options.tol;
    float_t locTol = options.locTol;

	float_t field;

	do
	{
		iter++;

		status = gsl_multimin_fminimizer_iterate (s);
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, tol);
		field = gsl_blas_dnrm2 (s->x);

       // printf ("%5lu %10.7e f() = %7.3f size = %.3f, status = %i\n",
       // 	      iter,
       //           gsl_vector_get (s->x, 0),
       //           s->fval, size, status);
	}
	while ( (status == GSL_CONTINUE) && (iter<=itermax) && (field<=fieldMax));

     // printf ("%5lu %10.3e f() = %7.3f size = %.3f, status = %i\n",
     // 	      iter,
     //           gsl_vector_get (s->x, 0),
     //           s->fval, size, status);

	//The following things check if all is complete.
	if ( (status!=GSL_SUCCESS) && (!((iter>itermax) && (field<locTol) ) ) )
	{
        //printf("fail: iter = %i, field = %f\n",iter,field);
	   	return 1;
	}
	else
	{
		gsl_vector_memcpy (min, s->x);
		return 0;
	}

	return 2;

}
