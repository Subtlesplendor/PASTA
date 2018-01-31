#include "AbstractModel.h"
//int counter = 0;
//-----------------------------MISC----------------------------

//Used by the numerical dM(T) function:
double my_RootFunction (double x, void *params)
{
    Model* model = (Model *)params;

    float_t detM = model->detM(x);
    return detM;
}

//Used by findSphscale:
double my_htPotFunction (const gsl_vector *v, void *params)
{

    Model* model = (Model *)params;

    const size_t numvar = model->getNumvar();
    float_t *fields = new float_t[numvar];

    for (size_t i = 0; i < numvar; ++i)
    {
        fields[i] = gsl_vector_get(v, i);
    }

    float_t htV = model->htPot(fields);
    delete[] fields;
    return htV;
}
//------------------------------------------------------------

//-------------------------Model--------------------------

            /*----Constructor & Destructor----*/
Model::Model(float_t* par, const size_t numParameters, Particle* particles[], const size_t numParticles,const float_t renScale2, const size_t scheme, const size_t numVariables): m_scheme(scheme), m_numVar(numVariables),m_numParameters(numParameters), m_numParticles(numParticles)
{
    this->m_pars = new float_t[m_numParameters];
    this->m_ctPars = new float_t[m_numParameters];
    this->m_renPars = new float_t[m_numParameters];
    this->m_treeMin = new float_t[m_numVar];
    this->m_totMin = new float_t[m_numVar];
    this->m_treeDeriv = new float_t[m_numVar];
    this->m_loopDeriv = new float_t[m_numVar];
    this->m_totDeriv = new float_t[m_numVar];
    this->m_refScales2 = new float_t[m_numParameters];

    setParticles(particles);
    setParameters(par);
    setRenpars(par);
    setTemperature(0.0);
    setGaugeXi(0.0);
    setRenormscale(renScale2);

    float_t *initialRefScales2 = new float_t[m_numParameters];
    std::fill_n(initialRefScales2, m_numParameters, 100*100);
    setRefscales(initialRefScales2);
    delete[] initialRefScales2;

    float_t* tempDeriv = new float_t[m_numVar];
    std::fill_n(tempDeriv, m_numVar,0.0);
    setTreeDeriv(tempDeriv);
    setLoopDeriv(tempDeriv);
    setTotDeriv(tempDeriv);
    delete[] tempDeriv;

}

Model::~Model()
{
    delete m_pars;
    delete m_ctPars;
    delete m_renPars;
    delete m_treeMin;
    delete m_totMin;
    delete m_treeDeriv;
    delete m_loopDeriv;
    delete m_totDeriv;
    delete m_refScales2;

}

            /*-------------------------------*/

            /*-------------Setters-----------*/

void Model::setParticles(Particle* particles[])
{

    m_particles = particles;

}
void Model::setParameters(float_t* par)
{
    size_t numparams = getNumparams();
    for (size_t i = 0; i < numparams; ++i)
    {
        m_pars[i] = par[i];
    }

}

void Model::setRenpars(float_t* renpars)
{
    size_t numparams = getNumparams();
    for (size_t i = 0; i < numparams; ++i)
    {
        m_renPars[i] = renpars[i];
    }

}

void Model::setCtparams(float_t* ctpars)
{
    size_t numparams = getNumparams();
    for (size_t i = 0; i < numparams; ++i)
    {
        m_ctPars[i] = ctpars[i];
    }

}

void Model::setTreemin(float_t* treemin)
{
    size_t numvar = getNumvar();
    for (size_t i = 0; i < numvar; ++i)
    {
        m_treeMin[i] = treemin[i];
    }

}

void Model::setTotmin(float_t* totmin)
{
    size_t numvar = getNumvar();
    for (size_t i = 0; i < numvar; ++i)
    {
        m_totMin[i] = totmin[i];
    }

}

void Model::setTreeDeriv(float_t* treeDeriv)
{
    size_t numvar = getNumvar();
    for (size_t i = 0; i < numvar; ++i)
    {
        m_treeDeriv[i] = treeDeriv[i];
    }
}

void Model::setLoopDeriv(float_t* loopDeriv)
{
    size_t numvar = getNumvar();
    for (size_t i = 0; i < numvar; ++i)
    {
        m_loopDeriv[i] = loopDeriv[i];
    }

}

void Model::setTotDeriv(float_t* totDeriv)
{
    size_t numvar = getNumvar();
    for (size_t i = 0; i < numvar; ++i)
    {
        m_totDeriv[i] = totDeriv[i];
    }

}

void Model::setRefscales(float_t* refScales2)
{
    size_t numpar = getNumparams();
    for (size_t i = 0; i < numpar; ++i)
    {
        m_refScales2[i] = refScales2[i];
    }

}

void Model::setCounterterms()
{
    size_t scheme = getScheme();
    size_t numParams = getNumparams();
    float_t *temp= new float_t[numParams];
    std::fill_n(temp, numParams, 0);


    if ((scheme == 1)||(scheme == 0))
    {
        setCtparams(temp);
        setRenpars(getParams());
    }
    else if (scheme == 2)
    {
        setCtparams(temp);

        float_t Q2 = getRenormscale();
        float_t* par = getParams();

        float_t *renPar = new float_t[numParams];
        //std::fill_n(renPar, numParams, 0);
        float_t* refScales2 = getRefscales();


        for (size_t i = 0; i < numParams; ++i)
        {
            renPar[i] = par[i] + 0.5*betaFunctions(i)*log(Q2/refScales2[i]);
        }

        setRenpars(renPar);


        delete[] renPar;

    }

    else
    {
        printf("Scheme nr %lu is not defined in this model\n", scheme);
        exit(EXIT_FAILURE);
    }

    delete[] temp;

}

            /*-------------------------------*/

            /*------------Modifiers----------*/

void Model::updateModel(float_t* par)
{
    setParameters(par);
    setTemperature(0.0);
    setCounterterms();

    this->findTreemin();
    size_t errorFlag=this->findT0();
    if (errorFlag != 0) {
        printf("Error finding T0 during update.\n");
    }
}

size_t Model::findT0()
{

    size_t T0errorFlag;

    //This is a numerical solver which is used in the standard case.

    //--------------------Setting up the solver:-------------
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double rootT0 = 0.0;
    double minTemp = 0.0, maxTemp = 1000;
    double absTol = 0.0, relTol = 0.001;

    double detMMin=detM(minTemp);
    double detMMax=detM(maxTemp);

    if (detMMin*detMMax >= 0.0f)
    {
        T0errorFlag = 1;
        return T0errorFlag;
    }

    gsl_function F;
    F.function = &my_RootFunction;
    F.params = this;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, minTemp, maxTemp);
    //--------------------------------------------------------

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        rootT0 = gsl_root_fsolver_root (s);
        minTemp = gsl_root_fsolver_x_lower (s);
        maxTemp = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (minTemp, maxTemp, absTol, relTol);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    if (status != GSL_SUCCESS)
    {
        setT0(0.0);
        T0errorFlag = 1;
    }
    else
    {
        setT0(rootT0);
        T0errorFlag = 0;
    }


    gsl_root_fsolver_free (s);

    return T0errorFlag;

}

void Model::findSphscale(float_t Tc, float_t& vc, size_t& errorFlag)
{

    this->setTemperature(Tc);

    const size_t numVar = this->getNumvar();
    float_t *treeMin = new float_t[numVar];
    //std::fill_n(treeMin, numVar, 0);
    for (size_t i = 0; i < numVar; ++i)
    {
        treeMin[i] = this->getTreemin()[i];
    }

    //------------------initializing minimizer--------------------
    minOptions_t options;

    const gsl_multimin_fminimizer_type *minType = gsl_multimin_fminimizer_nmsimplex2rand;

    gsl_multimin_fminimizer *s;
    gsl_vector *ss, *min;
    gsl_multimin_function my_func;

    min = gsl_vector_alloc (numVar);
    ss = gsl_vector_alloc (numVar);
    gsl_vector_set_all (ss, options.stepSize); //initial step size
    s = gsl_multimin_fminimizer_alloc(minType, numVar);

    for (size_t i = 0; i < numVar; ++i)
    {
        gsl_vector_set (min, i, treeMin[i]);
    }

    my_func.n = numVar;
    my_func.f = my_htPotFunction;
    my_func.params = this;

    gsl_multimin_fminimizer_set(s, &my_func, min,ss);
    //-------------------------------------------------------

    errorFlag = Minimizer(options, s, min, numVar);

    //gsl_vector_set_all (ss, vc);


    if (errorFlag==0)
    {
        vc = gsl_blas_dnrm2 (min);
    }
    else
    {
        vc = 0.0;
    }

    gsl_vector_free (min);
    gsl_vector_free (ss);
    gsl_multimin_fminimizer_free (s);

    delete[] treeMin;

    return;

}

            /*-------------------------------*/

            /*------Numerical functions------*/

float_t Model::betaFunctions(size_t parIndex) const
{
    return 0;
}

double Model::detM(float_t T) const
{
    return 0; //This is just a placeholder. Implement this in your model if you wish to use it.
}

void Model::findTreePotDeriv(float_t* fields)
{
    float_t* tempDeriv;
    std::fill_n(tempDeriv, m_numVar,0.0);
    setTreeDeriv(tempDeriv);
    delete[] tempDeriv;
}

float_t Model::totalPot(float_t* fields) const
{
    return  treeLevel(fields) + oneLoopLevel(fields);
}

void Model::findTotalPotDeriv(float_t* fields)
{
    float_t* temp = new float_t[m_numVar];
    findTreePotDeriv(fields);
    findOneLoopLevelDeriv(fields);
    for (size_t i = 0; i < m_numVar; i++) {
        temp[i]=m_treeDeriv[i]+m_loopDeriv[i];
    }
    setTotDeriv(temp);
    delete[] temp;
}

float_t Model::deltaPot(float_t* fields) const
{
    size_t numvar = getNumvar();
    float_t *temp = new float_t[numvar];
    std::fill_n(temp,numvar,0);
    float_t dv = treeLevel(fields) + oneLoopLevel(fields) - oneLoopLevel(temp);
    //printf("field = %f\n",fields[0]);
    //printf("phi0 = %f\n",temp[0]);
    //printf("V1phi = %f \n", oneLoopLevel(fields));
    //printf("V1_0 = %f \n", oneLoopLevel(temp));
    delete[] temp;
    return  dv;
}

float_t Model::deltaPotResummed(float_t* fields)
{
    size_t numvar = getNumvar();
    float_t *temp = new float_t[numvar];
    float_t *min = new float_t[numvar];
    float_t phi1 = getphi1();
    std::fill_n(temp,numvar,0);
    for (size_t i = 0; i < numvar; i++) {
        min[i] = fields[i] + phi1;
    }
    float_t dv = treeLevel(min) + oneLoopLevelResummed(min) - oneLoopLevelResummed(temp);
    //printf("min = %f\n",min[0]);
    //printf("temp = %f\n",temp[0]);
    //printf("field = %f\n",fields[0]);
    //printf("phi0 = %f\n",temp[0]);
    //printf("V1phi = %f \n", oneLoopLevel(fields));
    //printf("V1_0 = %f \n", oneLoopLevel(temp));

    printf("treeMin = %f, phi1 = %f, min = %f\n",fields[0], phi1, min[0]);
    //printf("deltaPotResummed\n");
    delete[] temp;
    delete[] min;
    return  dv;
}

float_t Model::tempPot(const float_t x2, bool fermionFlag) const
{
    //printf("x2 = %f\n",x2);
    float_t x = std::sqrt(std::abs(x2));
    //printf("x2 = %f\n",x2);
    //printf("x = %f\n",x);

    if (x2<-175.0) {
        slong prec;
        arb_t w1, w2, s1, s2, a;
        float_t arg;
        float_t result=0.0;
        float_t wavl = 2 * M_PI;

        arb_init(w1); arb_init(w2); arb_init(s1); arb_init(s2); arb_init(a);


        arg = 1.0 - fmod(x,wavl)/wavl;
        //printf("arg = %.5f\n", arg);

        arb_set_d(s1,-1.5);
        arb_set_d(s2,-2.5);
        arb_set_d(a,arg);

        prec=59; // ca 16 digit precision

        arb_hurwitz_zeta(w1, s1, a, prec);
        arb_hurwitz_zeta(w2, s2, a, prec);
        //printf("x2 = %f,arg = %.5f\n",x2, arg);
        //printf("w1 = %1.15f\n",arf_get_d(arb_midref(w1),ARF_RND_UP));
        //printf("w2 = %1.15f\n",arf_get_d(arb_midref(w2),ARF_RND_UP));
        //counter+=1;
        //printf("counter = %i\n",counter);
        //arb_printn(w1, digits, ARB_STR_CONDENSE * condense);
        //flint_printf("\n");
        double W1, W2;
        W1 = arf_get_d(arb_midref(w1),ARF_RND_UP);
        W2 = arf_get_d(arb_midref(w2),ARF_RND_UP);

        result = -46.64911554033297*pow(x,1.5)*W1 -219.8287780169573*pow(x,0.5)*W2 ;

        //flint_printf("result = %1.15f\n",result);

        arb_clear(s1); arb_clear(s2); arb_clear(w1); arb_clear(w2); arb_clear(a);
        flint_cleanup();

        return result;

    }
    else if (x2>400.0)
    {
        return exp(-x)*(0.646082328590889*pow(x,-3.5) -0.397589125286701*pow(x,-2.5) +0.385540969974983*pow(x,-1.5) -1.02810925326662*pow(x,-0.5) -2.34996400746656*pow(x,0.5) -1.25331413731550 *pow(x,1.5));
    }
    else if ( (x2>=-175.0) && (x2<=400.0) )
    {
        int length;
        float_t result;
        if (!fermionFlag)
        {
            length = 7500;
            result = integralInterpol(x2, JBosonTable, length);
        }
        else
        {
            length = 5000;
            result = integralInterpol(x2, JFermionTable, length);
        }
            return result;
    }
    else
    {
        printf("x = %f\n", x);
        //counter=0;
        exit(EXIT_FAILURE);// x not in range
    }


}

float_t Model::tempPotDeriv(const float_t x2, bool fermionFlag) const
{

    const float_t h = 0.00001;

    if ( (fermionFlag) && ((x2-2 * h)<0) ) {
        return M_PI * M_PI/24.0;
    }
    else
    {
        float_t result = (-tempPot(x2+2.0*h,fermionFlag) +8.0*tempPot(x2+h,fermionFlag) -8.0*tempPot(x2-h,fermionFlag) +tempPot(x2-2.0*h,fermionFlag))/(12.0 * h);
        return result;
    }


}

/*float_t Model::tempPotDeriv(const float_t x2, bool fermionFlag) const
{

    float_t x = std::sqrt(std::abs(x2));
    //printf("x = %f\n",x);

    //The low T expansion is common for fermions and bosons:
    if (x2>400.0)
    {
        return exp(-x)*(-0.423991528137771*pow(x,-5.5) +0.173945242312932*pow(x,-4.5) -0.0903611648378866*pow(x,-3.5) +0.0642568283291638*pow(x,-2.5)-0.0734363752333301*pow(x,-1.5) +0.234996400746656*pow(x,-0.5)+0.626657*pow(x,0.5));
    }
    else if ( (!fermionFlag) && (x2>0.0) && (x<=std::sqrt(0.5) ) )
    {
        return 0.822467033424113 - 0.785398163397448 * x + x2 * (0.306726072758470 -0.125 * log(x)) +0.000951514283074790 * x2 * x2 -0.0000103955878258706 * x2 * x2 * x2;
    }
    else if ( (fermionFlag) && (x2>0.0) && (x<=std::sqrt(0.2)) )
    {
                return 0.411233516712057 + x2 * (-0.133439277618483 + 0.0625 * log(x2) ) -0.00666059998152353 * x2 * x2 + 0.000322263222601989 * x2 * x2 * x2;
    }
    else if ( (x2>=0.0) && (x2<=400.0) )
    {

        int length;
        float_t result;
        if (!fermionFlag)
        {
            length = 5000;
            result = integralInterpol(x2, h3Table, length);
        }
        else
        {
            length = 5000;
            result = integralInterpol(x2, f3Table, length);
        }
            return result;
    }
    else if ( (!fermionFlag) && (x2<0.0) && x2>=-16.6){
                return 0.822467033424113 + x2 * (0.306726072758470 -0.125 * log(x)) +0.000951514283074790 * x2 * x2 -0.0000103955878258706 * x2 * x2 * x2 +1.60041209405272e-7 * x2 * x2 * x2 * x2 -2.81987914550857e-9 * x2 * x2 * x2 * x2 * x2 +5.34903246226569e-11* x2 * x2 * x2 * x2 * x2 * x2;
    }
    else if ( (!fermionFlag) && (x2<-16.6) ) {

        slong prec;
        arb_t w1, w2, w3, w4, w5;
        arb_t s1, s2, s3, s4, s5;
        arb_t a;
        float_t arg;
        float_t result;
        float_t wavl = 2 * M_PI;

        arb_init(w1); arb_init(w2); arb_init(w3); arb_init(w4); arb_init(w5);
        arb_init(s1); arb_init(s2); arb_init(s3); arb_init(s4); arb_init(s5);
        arb_init(a);

        arg = 1.0 - fmod(x,wavl)/wavl;;
        //printf("arg = %.5f\n", arg);

        arb_set_d(s1,-0.5);
        arb_set_d(s2,-1.5);
        arb_set_d(s3,-2.5);
        arb_set_d(s4,-3.5);
        arb_set_d(s5,-4.5);

        arb_set_d(a,arg);

        prec=59; // ca 16 digit precision

        arb_hurwitz_zeta(w1, s1, a, prec);
        arb_hurwitz_zeta(w2, s2, a, prec);
        arb_hurwitz_zeta(w3, s3, a, prec);
        arb_hurwitz_zeta(w4, s4, a, prec);
        arb_hurwitz_zeta(w5, s5, a, prec);
        //arb_printn(w, digits, ARB_STR_CONDENSE * condense);
        //flint_printf("\n");

        result =-5.56832799683171*pow(x,0.5)*arf_get_d(arb_midref(w1),ARF_RND_UP)
        -8.74670916381243*pow(x,-0.5)*arf_get_d(arb_midref(w2),ARF_RND_UP)
        +6.86964931302991*pow(x,-1.5)*arf_get_d(arb_midref(w3),ARF_RND_UP)
        -10.7908199072765*pow(x,-2.5)*arf_get_d(arb_midref(w4),ARF_RND_UP)
        +67.8007210938205*pow(x,-3.5)*arf_get_d(arb_midref(w5),ARF_RND_UP);

        flint_printf("result = %1.15f\n",result);

        arb_clear(s1); arb_clear(s2); arb_clear(s3); arb_clear(s4); arb_clear(s5);
        arb_clear(w1); arb_clear(w2); arb_clear(w3); arb_clear(w4);  arb_clear(w5);
        arb_clear(a);
        flint_cleanup();
        return result;

    }
    else
    {
        printf("x = %f\n", x);
        exit(EXIT_FAILURE);// x not in range
    }

} */

float_t Model::oneLoopLevel(float_t* fields) const
{
    //Holds info for each particle:
    float_t pMass2;
    int pDof;
    int pSign;
    float_t pSpin;

    //Holds values
    float_t c;
    float_t Vloop = 0;
    float_t x2 = 0.0;
    //printf("T = %f\n", m_temperature);
    float_t T2 = m_temperature * m_temperature;

    int ghostSign;
    bool fermionFlag;


    for (size_t i = 0; i < m_numParticles; ++i)
    {
        //printf("field = %f\n", fields[0]);
        pMass2 = m_particles[i]->getFieldmass(m_pars, m_gaugeXi, fields, 0.0);
        pDof = m_particles[i]->getDof();
        pSpin = m_particles[i]->getSpin();
        pSign = m_particles[i]->getSign();

        fermionFlag = (pSpin!=floor(pSpin));


        //ghosts have boson statistics but opposite sign.
        if ( (!fermionFlag) && (pSign<0) )
        {
            ghostSign = -1;
        }
        else
        {
            ghostSign = 1;
        }



        if ((pMass2!=0)&&(m_scheme!=0))
        {

            if (pSpin!=1)
            {
                c = 1.5;
            }
            else
            {
                c = 0.833333333333333;// 5/6 for spin-1 bosons
            }

            Vloop += 0.25 * kappa * pSign * pDof * pMass2 * pMass2 * (log(std::abs(pMass2)/m_renormScale2) - c);
        }
        /*
        This functions uses high temperature and low temperature approximations,
        with a smooth interpolation in between:
        */

        if (m_temperature==0)
        {
            //printf("T == 0\n");
            continue;
        }

        x2 = pMass2 / T2;
        //printf("particle # %lu has mass %f at phi = %f, T = %f, xi = %f\n",i,pMass2,fields[0],m_temperature,m_gaugeXi);

        Vloop += (0.5/(M_PI * M_PI)) * T2 * T2 * pDof*ghostSign*tempPot(x2, fermionFlag);

        //printf("particle # %lu has mass %f and Pot = %f at T = %f\n",i,pMass2,(0.5/(M_PI * M_PI)) * m_temperature*m_temperature*m_temperature*m_temperature * pDof*ghostSign*tempPot(x2, fermionFlag),m_temperature);

    }

    return Vloop;
}

float_t Model::getphi1()
{
    float_t* pMass2Deriv = new float_t[m_numVar];
    float_t* treeMin = new float_t[m_numVar];
    treeMin = getTreemin();
    for (size_t i = 0; i < m_numVar; i++) {
    	pMass2Deriv[i]=0;
    }


    pMass2Deriv = m_particles[1]->getFieldmassDeriv(m_pars, 0.0, treeMin, 0.0);

    float_t gbSEmin = gbSelfEnergy(treeMin,0);
    float_t phi1 = 0.0;
    float_t mgpmin = pMass2Deriv[0];

    phi1 = - gbSEmin/mgpmin;

    //printf("getphi1\n");
    //delete[] treeMin;
    //delete[] pMass2Deriv;
    return phi1;
}


float_t Model::oneLoopLevelResummed(float_t* fields)
{
    //Holds info for each particle:
    float_t pMass2;
    int pDof;
    int pSign;
    float_t pSpin;
    bool pGB;

    int ghostSign;
    bool fermionFlag;

    //Holds values
    float_t c;
    float_t Vloop = 0;
    float_t x2 = 0.0;
    //printf("T = %f\n", m_temperature);
    float_t T2 = m_temperature * m_temperature;

    float_t* treeMin = new float_t[m_numVar];
    float_t* origin = new float_t[m_numVar];
    treeMin = getTreemin();

    float_t* pMass2Deriv = new float_t[m_numVar];

    float_t fieldLength = 0.0;
    for (size_t i = 0; i < m_numVar; i++) {
        origin[i] = 0.0;
        fieldLength += fields[i];
    	pMass2Deriv[i]=0;
    }

    pMass2Deriv = m_particles[1]->getFieldmassDeriv(m_pars, 0.0, treeMin, 0.0);


    float_t gbSEmin = gbSelfEnergy(treeMin,0);
    float_t gbSE0 = gbSelfEnergy(origin,0);
    float_t phi1 = 0.0;
    float_t mgpmin = pMass2Deriv[0];
    //float_t gbSE = gbSE0;

    if (fieldLength!=0.0) {
        phi1 = getphi1();
        //gbSE = gbSEmin;
    }




    for (size_t i = 0; i < m_numParticles; ++i)
    {
        //printf("field = %f\n", fields[0]);
        pMass2 = m_particles[i]->getFieldmass(m_pars, m_gaugeXi, fields, 0.0);
        pDof = m_particles[i]->getDof();
        pSpin = m_particles[i]->getSpin();
        pSign = m_particles[i]->getSign();
        pGB = m_particles[i]->getGB();

        fermionFlag = (pSpin!=floor(pSpin));

        if (pGB)
        {
            if ((fieldLength!=0.0)) {
                pMass2 = m_particles[i]->getFieldmass(m_pars, m_gaugeXi, treeMin, 0.0) + phi1 * mgpmin + gbSEmin;
            }
            else
            {
            pMass2 = m_particles[i]->getFieldmass(m_pars, m_gaugeXi, fields, 0.0) + gbSE0;
            }
        }



        //ghosts have boson statistics but opposite sign.
        if ( (!fermionFlag) && (pSign<0) )
        {
            ghostSign = -1;
            if ((fieldLength!=0.0)) {
                pMass2 = m_particles[i]->getFieldmass(m_pars, m_gaugeXi, treeMin, 0.0);
            }
        }
        else
        {
            ghostSign = 1;
        }



        if ((pMass2!=0)&&(m_scheme!=0))
        {

            if (pSpin!=1)
            {
                c = 1.5;
            }
            else
            {
                c = 0.833333333333333;// 5/6 for spin-1 bosons
            }

            Vloop += 0.25 * kappa * pSign * pDof * pMass2 * pMass2 * (log(std::abs(pMass2)/m_renormScale2) - c);
        }
        /*
        This function uses high temperature and low temperature approximations,
        with a smooth interpolation in between:
        */

        if (m_temperature==0)
        {
            //printf("T == 0\n");
            continue;
        }

        x2 = pMass2 / T2;
        //printf("particle # %lu has mass %f at phi = %f, T = %f, xi = %f\n",i,pMass2,fields[0],m_temperature,m_gaugeXi);

        Vloop += (0.5/(M_PI * M_PI)) * T2 * T2 * pDof*ghostSign*tempPot(x2, fermionFlag);

        //printf("particle # %lu has mass %f and Pot = %f at T = %f\n",i,pMass2,(0.5/(M_PI * M_PI)) * m_temperature*m_temperature*m_temperature*m_temperature * pDof*ghostSign*tempPot(x2, fermionFlag),m_temperature);

    }
    //printf("oneLoopLevelResummed\n");
    //delete[] pMass2Deriv;
    //delete[] treeMin;
    //delete[] origin;
    return Vloop;
}

void Model::findOneLoopLevelDeriv(float_t* fields)
{

    //Holds info for each particle:
    float_t pMass2;
    int pDof;
    int pSign;
    float_t pSpin;


    float_t* pMass2Deriv = new float_t[m_numVar];
    float_t* loopDeriv = new float_t[m_numVar];

    for (size_t i = 0; i < m_numVar; ++i)
    {
    	pMass2Deriv[i]=0;
        loopDeriv[i]=0;
    }
    //Holds values
    float_t cp;
    float_t x2;
    float_t T2 = m_temperature * m_temperature;


    int ghostSign;
    bool fermionFlag;

    for (size_t i = 0; i < m_numParticles; ++i)
    {

        pMass2 = m_particles[i]->getFieldmass(m_pars, m_gaugeXi, fields, 0.0);
        pMass2Deriv = m_particles[i]->getFieldmassDeriv(m_pars, m_gaugeXi, fields, 0.0);
        pDof = m_particles[i]->getDof();
        pSpin = m_particles[i]->getSpin();
        pSign = m_particles[i]->getSign();

        fermionFlag = (pSpin!=floor(pSpin));

        //ghosts have boson statistics but opposite sign.
        if ( (!fermionFlag) && (pSign<0) )
        {
            ghostSign = -1;
        }
        else
        {
            ghostSign = 1;
        }
        for (size_t j = 0; j < m_numVar; ++j)
        {
        	if ((pMass2!=0)&&(m_scheme!=0))
        	{

            	if (pSpin!=1)
            	{
                	cp = 1;
            	}
            	else
            	{
                	cp = 0.333333333333333; // 1/3 for spin-1 bosons
            	}

        		loopDeriv[j]+= 0.5 * kappa * pSign * pDof  * pMass2Deriv[j] * pMass2 * (log(std::abs(pMass2)/m_renormScale2) - cp);

        	}

            /*
        	This functions uses high temperature and low temperature approximations,
        	with a smooth interpolation in between:
        	*/

            if (m_temperature==0)
            {
                //printf("T == 0\n");
                continue;
            }

            x2 = pMass2 / T2;
        	loopDeriv[j] += (0.5/(M_PI * M_PI)) * T2 * pMass2Deriv[j] * pDof*ghostSign*tempPotDeriv(x2, fermionFlag);

            printf("x2 = %f, ghostSign = %i, fermionFlag = %i, tDeriv = %f\n", x2, ghostSign, fermionFlag,(0.5/(M_PI * M_PI)) * T2 * pMass2Deriv[j] * pDof*ghostSign*tempPotDeriv(x2, fermionFlag));
            //printf("ghostSign = %i\n",ghostSign );


        }

    }
    setLoopDeriv(loopDeriv);
    //printf("findOneLoopLevelDeriv\n");
    delete[] pMass2Deriv;
    delete[] loopDeriv;
    //return m_loopDeriv;
}

float_t Model::gbSelfEnergy(float_t* fields,size_t dirIndex)
{

    //Holds info for each particle:
    float_t pMass2;
    int pDof;
    int pSign;
    float_t pSpin;
    bool pGB;

    float_t gbSE = 0;

    float_t* pMass2Deriv = new float_t[m_numVar];
    float_t* fieldTemp = new float_t[m_numVar];
    float_t fieldTempNorm2 = 0;
    float_t pM2DoverF = 0;

    for (size_t i = 0; i < m_numVar; ++i)
    {
    	pMass2Deriv[i] = 0;
        fieldTemp[i] = 1.0;
        fieldTempNorm2 += fieldTemp[i];
    }

    float_t fieldTempNorm = std::sqrt(fieldTempNorm2);

    //Holds values
    float_t cp;
    float_t x2;
    float_t T2 = m_temperature * m_temperature;

    int ghostSign;
    bool fermionFlag;

    for (size_t i = 0; i < m_numParticles; ++i)
    {
        pGB = m_particles[i]->getGB();
        pSpin = m_particles[i]->getSpin();
        pSign = m_particles[i]->getSign();

        fermionFlag = (pSpin!=floor(pSpin));

        if ( (!fermionFlag) && (pSign<0) )
        {
            ghostSign = -1;
        }
        else
        {
            ghostSign = 1;
        }


        if ( !(pGB) && (ghostSign == 1) )
        {
            pMass2Deriv = m_particles[i]->getFieldmassDeriv(m_pars, m_gaugeXi, fieldTemp, 0.0);
            pM2DoverF = pMass2Deriv[dirIndex]/fieldTempNorm;

            pMass2 = m_particles[i]->getFieldmass(m_pars, m_gaugeXi, fields, 0.0);
            pDof = m_particles[i]->getDof();


            //ghosts have boson statistics but opposite sign.



            if ((pMass2!=0)&&(m_scheme!=0))
            {

                if (pSpin!=1)
                {
                	cp = 1;
                }
                else
                {
                	cp = 0.333333333333333; // 1/3 for spin-1 bosons
                }

            	gbSE += 0.5 * kappa * pSign * pDof  * pM2DoverF * pMass2 * (log(std::abs(pMass2)/m_renormScale2) - cp);

            }

            if (m_temperature==0)
            {
                //printf("T == 0\n");
                continue;
            }

            x2 = pMass2 / T2;

            gbSE += (0.5/(M_PI * M_PI)) * T2 * pM2DoverF * pDof*ghostSign*tempPotDeriv(x2, fermionFlag);

        }

    }
    //printf("gbSelfEnergy\n");
    //delete[] pMass2Deriv;
    //delete[] fieldTemp;
    return gbSE;
}

            /*-------------------------------*/

//------------------------------------------------------------

//---------------------------PARTICLE-------------------------

            /*-----------Constructor---------*/

Particle::Particle(float_t mass, float_t spin, int dof, size_t numVar, bool isGB):  m_dof(dof), m_numVar(numVar), m_spin(spin), m_isGB(isGB)
{
    this->m_fieldMassDeriv = new float_t[m_numVar];
    setParticle(mass);
}
            /*-------------------------------*/

            /*------------Setters------------*/

void Particle::setParticle(float_t mass)
{
    m_mass = mass;

    if (floor(m_spin)==m_spin)
    {
        m_sign = 1;
    }
    else
    {
        m_sign = -1;
    }
    for (size_t i = 0; i < m_numVar; ++i)
    {
        m_fieldMassDeriv[i] = 0.0;
    }

}

float_t* Particle::getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T)
{
	 //Should be overwritten if used:
    for (size_t i = 0; i < m_numVar; ++i)
    {
        m_fieldMassDeriv[i] = 0.0;
    }

    return m_fieldMassDeriv;

}

            /*-------------------------------*/

//------------------------------------------------------------
