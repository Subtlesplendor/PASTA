/*****
 *
 *
 *  functions implementation file
 *
 *  Author: Johan Löfgren
 *
 *
 **/

#include "PTAfunctions.h"

//--------------------Numerical functions---------------------

//These are used by GSL methods.

//The total potential
double my_totPotFunction (const gsl_vector *v, void *params)
{

	Model* model = (Model *)params;
	float_t V;
	const size_t numvar = model->getNumvar();
	float_t *fields = new float_t[numvar];

    for (size_t i = 0; i < numvar; ++i)
    {
        fields[i] = gsl_vector_get(v, i);
		//printf("field = %f\n",fields[i]);
    }
	//printf("T = %f\n",model->getTemperature() );
	//printf("V = %f\n",model->totalPot(fields));
	V=model->totalPot(fields);
    delete[] fields;
	//delete model;
    return V;
}

//The difference between the potential evalueated at the EW min and the sym. min.
double my_deltaPotFunction (double T, void *params)
{

	Model* model = (Model *)params;
	model->setTemperature(T);
	float_t V;

/*    float_t treeMin[numVar] = {};
    for (size_t i = 0; i < numVar; ++i)
    {
        treeMin[i] = p->getTreemin()[i];
    }*/
	V = model->deltaPot(model->getTreemin());
    return V;
}

double my_deltaPotResummedFunction (double T, void *params)
{

	Model* model = (Model *)params;
	model->setTemperature(T);
	float_t V;

/*    float_t treeMin[numVar] = {};
    for (size_t i = 0; i < numVar; ++i)
    {
        treeMin[i] = p->getTreemin()[i];
    }*/
	V = model->deltaPotResummed(model->getTreemin());
    return V;
}

//------------------------------------------------------------------

//-------------Methods which are used to study the PT---------------

//Resummed:
void whenPT_PRM_R(Model* model, float_t& Tc, float_t& vc, size_t& errorFlag)
{

	float_t T0 = model->getT0();

	if (T0==0.0)
	{
		T0 = 1000;
	}

	model->setTemperature(0.0);
	//float_t dV = model->deltaPot(model->getTreemin());
	float_t dV = model->deltaPotResummed(model->getTreemin());

	model->setTemperature(T0);
	//float_t dVT0 = model->deltaPot(model->getTreemin());
	float_t dVT0 = model->deltaPotResummed(model->getTreemin());

	if (dVT0 * dV >=0.0f)
	{
		errorFlag = 1;
		return;
	}

    //--------------------Setting up the solver:-------------
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double rootDV = 0.0;
    double minTemp = 0.0, maxTemp = T0;
    double absTol = 0.0, relTol = 0.001;

    gsl_function F;
    F.function = &my_deltaPotResummedFunction;
    F.params = model;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, minTemp, maxTemp);
    //--------------------------------------------------------

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        rootDV = gsl_root_fsolver_root (s);
        minTemp = gsl_root_fsolver_x_lower (s);
        maxTemp = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (minTemp, maxTemp, absTol, relTol);

    }
    while (status == GSL_CONTINUE && iter < max_iter);

    if (status != GSL_SUCCESS)
    {
        Tc = 0.0;
        errorFlag = 1;
    }
    else
    {
        Tc = rootDV;
        model->findSphscale(Tc, vc, errorFlag);
    }

    gsl_root_fsolver_free (s);
    return;
}

void whenPT_PRM(Model* model, float_t& Tc, float_t& vc, size_t& errorFlag)
{

	float_t T0 = model->getT0();

	if (T0==0.0)
	{
		T0 = 1000;
	}

	model->setTemperature(0.0);
	float_t dV = model->deltaPot(model->getTreemin());
	//float_t dV = model->deltaPotResummed(model->getTreemin());

	model->setTemperature(T0);
	float_t dVT0 = model->deltaPot(model->getTreemin());
	//float_t dVT0 = model->deltaPotResummed(model->getTreemin());

	if (dVT0 * dV >=0.0f)
	{
		errorFlag = 1;
		return;
	}

    //--------------------Setting up the solver:-------------
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double rootDV = 0.0;
    double minTemp = 0.0, maxTemp = T0;
    double absTol = 0.0, relTol = 0.001;

    gsl_function F;
    F.function = &my_deltaPotFunction;
    F.params = model;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, minTemp, maxTemp);
    //--------------------------------------------------------

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        rootDV = gsl_root_fsolver_root (s);
        minTemp = gsl_root_fsolver_x_lower (s);
        maxTemp = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (minTemp, maxTemp, absTol, relTol);

    }
    while (status == GSL_CONTINUE && iter < max_iter);

    if (status != GSL_SUCCESS)
    {
        Tc = 0.0;
        errorFlag = 2;
    }
    else
    {
        Tc = rootDV;
        model->findSphscale(Tc, vc, errorFlag);
    }

    gsl_root_fsolver_free (s);
    return;
}

//This next method uses the traditional, gauge-dependent, method in Landau gauge.


//Uses a linear temperature search
void whenPT_L(Model* model, float_t& Tc, float_t& vc, size_t& errorFlag, minOptions_t options)
{

    //------------------initializing minimizer--------------------
	float_t maxTemp = options.maxTemp;
	float_t minTemp = options.minTemp;
	float_t tempTol = options.tempTol;
	float_t stepSize = options.stepSize;

    const gsl_multimin_fminimizer_type *minType = gsl_multimin_fminimizer_nmsimplex2;

    gsl_multimin_fminimizer *s;
    gsl_vector *ss, *min;
    gsl_multimin_function my_func;

    const size_t numVar = model->getNumvar();
    min = gsl_vector_alloc (numVar);
    ss = gsl_vector_alloc (numVar);
    gsl_vector_set_all (ss, stepSize); //initial step size
    s = gsl_multimin_fminimizer_alloc(minType, numVar);

    my_func.n = numVar;
    my_func.f = my_totPotFunction;
    my_func.params = model;


	int status;

	//check that T=Tmax is ok
 	gsl_vector_set_all(min, 0.0);
	model->setTemperature(maxTemp);
    gsl_multimin_fminimizer_set(s, &my_func, min, ss);
	status = whichPhase(s, ss, min, my_func, model, errorFlag, options);
	if (status==1)
	{
		errorFlag = 3;
		Tc = 0.0;
		vc = 0.0;
	}
	else if (errorFlag!=0)
	{
		Tc = 0.0;
		vc = 0.0;
	}


	gsl_vector_set_all(min, 246.234133);
	vc = gsl_blas_dnrm2 (min);
	float_t tempsize = gsl_blas_dnrm2 (min)/2.0;

	float_t intervalLength = 5;
	float_t newTemp = 0.0;

	while( (intervalLength >= tempTol) && (errorFlag==0) )
	{
 		gsl_vector_set_all(min,vc);
		gsl_vector_set_all (ss, stepSize);

		model->setTemperature(newTemp);
		//printf("-------T = %f-------\n",model->getTemperature());

	    gsl_multimin_fminimizer_set(s, &my_func, min,ss);

		status = whichPhase(s, ss, min, my_func, model, errorFlag, options);
        //printf("vc = %f\n",gsl_blas_dnrm2 (min));
		if (status==0)
		{
			newTemp-=intervalLength;
			intervalLength/=2.0;
		}
		else if (status==1)
		{
			vc = gsl_blas_dnrm2 (min);
			//printf("vc = %f \n",vc );
			//gsl_vector_set_all (ss, vc);
		}

		newTemp += intervalLength;

	}

	//printf("vc = %f\n",vc);

    gsl_vector_free (min);
    gsl_vector_free (ss);
    gsl_multimin_fminimizer_free (s);
	if (errorFlag==0)
	{
		Tc = newTemp;
		return;
	}
	else
	{
		Tc = 0.0;
		vc = 0.0;
		return;
	}

}


//Uses a bijection temperature search:
void whenPT_L_FAST(Model* model, float_t& Tc, float_t& vc, size_t& errorFlag, minOptions_t options)
{

    //------------------initializing minimizer--------------------
	float_t maxTemp = options.maxTemp;
	float_t minTemp = options.minTemp;
	float_t tempTol = options.tempTol;
	float_t stepSize = options.stepSize;

    const gsl_multimin_fminimizer_type *minType = gsl_multimin_fminimizer_nmsimplex2;

    gsl_multimin_fminimizer *s;
    gsl_vector *ss, *min;
    gsl_multimin_function my_func;

    const size_t numVar = model->getNumvar();
    min = gsl_vector_alloc (numVar);
    ss = gsl_vector_alloc (numVar);
    gsl_vector_set_all (ss, stepSize); //initial step size
    s = gsl_multimin_fminimizer_alloc(minType, numVar);

    my_func.n = numVar;
    my_func.f = my_totPotFunction;
    my_func.params = model;

	int status;

	//check that T=Tmax is ok
 	gsl_vector_set_all(min, 0.0);
	model->setTemperature(maxTemp);
    gsl_multimin_fminimizer_set(s, &my_func, min, ss);
	status = whichPhase(s, ss, min, my_func, model, errorFlag, options);
	if (status==1)
	{
		errorFlag = 3;
		Tc = 0.0;
		vc = 0.0;
	}
	else if (errorFlag!=0)
	{
		Tc = 0.0;
		vc = 0.0;
	}

	float_t intervalLength = maxTemp - minTemp;
	float_t newTemp = minTemp + intervalLength/2.0;
	gsl_vector_set_all(min, 246.234133);
	vc = gsl_blas_dnrm2 (min);
	float_t tempsize = gsl_blas_dnrm2 (min)/2.0;
	//printf("tempsize = %f\n",tempsize );

	while( (intervalLength >= tempTol) && (errorFlag==0) )
	{
		tempsize = vc/2;
 		gsl_vector_set_all(min,vc);
		gsl_vector_set_all (ss, stepSize);
		//gsl_vector_set_all(min, 0.0);
		//gsl_vector_set_all (ss, stepSize);
		//printf("tempsize = %f\n",tempsize);
		model->setTemperature(newTemp);
		//printf("-------T = %f-------\n",model->getTemperature());
		//printf("vc = %f\n", vc);
    	//gsl_vector_set_all (ss, options.stepSize); //initial step size
	    gsl_multimin_fminimizer_set(s, &my_func, min,ss);
		//printf("ss = %f\n",gsl_vector_get (ss, 0));

		status = whichPhase(s, ss, min, my_func, model, errorFlag, options);
        //printf("vc = %f\n",gsl_blas_dnrm2 (min));
		if (status==0)
		{
			maxTemp = newTemp;
			//gsl_vector_set_all (ss, stepSize);
		}
		else if (status==1)
		{
			minTemp = newTemp;
			vc = gsl_blas_dnrm2 (min);
			//printf("vc = %f \n",vc );
			//gsl_vector_set_all (ss, vc);
		}

		intervalLength /= 2.0;
		newTemp = minTemp + intervalLength/2.0;
	}

	//printf("vc = %f\n",vc);

    gsl_vector_free (min);
    gsl_vector_free (ss);
    gsl_multimin_fminimizer_free (s);
	if (errorFlag==0)
	{
		Tc = newTemp;
		return;
	}
	else
	{
		Tc = 0.0;
		vc = 0.0;
		return;
	}

}

//Used by the traditional method:
int whichPhase(gsl_multimin_fminimizer *s, gsl_vector* ss, gsl_vector* min, gsl_multimin_function my_func, Model* model, size_t& errorFlag, minOptions_t options)
{
	size_t numvar = model->getNumvar();
	double locTol = options.locTol;

	float_t *temp = new float_t[numvar];
	std::fill_n(temp,numvar,0);
	float_t V0 = model->totalPot(temp);
	delete[] temp;

	errorFlag = Minimizer(options, s, min, numvar);
	if (errorFlag!=0)
	{
		return 2;
	}

	double minLength = gsl_blas_dnrm2 (min);
	double fval = s->fval;
        // printf ("minLength = %7.3f, T =%7.3f,  f() = %7.4f, f(0) =%7.4f \n",
        //           minLength,
        //           model->getTemperature(),
        //           fval,V0);


	if ( (fval<=V0) && (minLength > locTol) )
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

int whichPhase2(gsl_multimin_fminimizer *s, gsl_vector* ss, gsl_vector* min, gsl_multimin_function my_func, Model* model, size_t& errorFlag, minOptions_t options)
{
	size_t numvar = model->getNumvar();
	double locTol = options.locTol;

	float_t *temp = new float_t[numvar];
	std::fill_n(temp,numvar,0);
	float_t V0 = model->totalPot(temp);
	delete[] temp;

	errorFlag = Minimizer(options, s, min, numvar);
	if (errorFlag!=0)
	{
		return 2;
	}

	double minLength = gsl_blas_dnrm2 (min);
	double fval = s->fval;
        printf ("minLength = %7.3f, T =%7.3f,  f() = %7.4f, f(0) =%7.4f \n",
                  minLength,
                  model->getTemperature(),
                  fval,V0);


	if ( (std::abs(1- fval/V0) < 1e-3) && (minLength < 30) )
	{
		printf("fval/V0 = %f, minLength =%f \n",fval/V0,minLength );
		return 0;
	}

	if ((fval>V0)&&(errorFlag==0))
	{
		size_t errorFlag2=0;
		gsl_vector *min2;
		min2 = gsl_vector_alloc (numvar);
	    gsl_vector_set_all (min2, 0.0);

		errorFlag2 = Minimizer(options, s, min2, numvar);
		double min2Length = gsl_blas_dnrm2 (min2);
		gsl_vector_free (min2);
		if (min2Length>locTol)
		{
			errorFlag = errorFlag2;
			return 1;
		}
		else
		{
			return 0;
		}
	}
	else if ( (fval<=V0) && (minLength > locTol) )
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

//------------------------------------------------------------------
