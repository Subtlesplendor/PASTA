#include "SM.h"

//-------------------------SM-----------------------------

            /*-----------Setters------------*/

/*void Model::setCounterterms()
{
    size_t scheme = getScheme();
    size_t numParams = getNumparams();
    float_t temp[numParams]= {};

    printf("hiPOT\n");

    if ((scheme == 1)||(scheme == 0))
    {
        setCtparams(temp);
        setRenpars(getParams());
    }
    else
    {
        printf("Scheme nr %i is not defined in this model\n", scheme);
        exit(EXIT_FAILURE);
    }
}*/

            /*-----------Modifiers-----------*/
void SM::findTreemin()
{
    float_t mu2 = getParams()[0];
    float_t lambda = getParams()[1];

    float_t *temp = new float_t[getNumvar()];

    temp[0] = sqrt( mu2 / lambda);

    setTreemin(temp);

    delete[] temp;
}

size_t SM::findT0()
{

    float_t lambda = getParams()[1];

    float_t D = (3* g2 + gprime2 + 8 * lambda + 4 * yt * yt);

    float_t temp = sqrt(16 * lambda * VEV2 / D);

    setT0(temp);

    return 0;
}

void SM::findSphscale(float_t Tc, float_t& vc, size_t& errorFlag)
{

    float_t T0 = getT0();

    vc = VEV * sqrt(1-Tc*Tc/(T0 * T0));

}

            /*-------------------------------*/

            /*------Numerical functions------*/

float_t SM::treeLevel(float_t* fields) const
{

    float_t mu2R = getRenparams()[0];
    float_t lambdaR = getRenparams()[1];

    return -0.5*mu2R * fields[0]*fields[0] + 0.25*lambdaR*fields[0]*fields[0]*fields[0]*fields[0];
}

void SM::findTreePotDeriv(float_t* fields)
{
    float_t mu2R = getRenparams()[0];
    float_t lambdaR = getRenparams()[1];

    float_t* treeDeriv = new float_t[1];

    treeDeriv[0]= -mu2R * fields[0] + lambdaR*fields[0]*fields[0]*fields[0];;
    setTreeDeriv(treeDeriv);
    delete[] treeDeriv;
}

float_t SM::htPot(float_t* fields) const
{

    float_t mu2 = getParams()[0];
    float_t lambda = getParams()[1];

    float_t temperature = getTemperature();

    return treeLevel(fields) - (temperature * temperature/96.0) *(16 * mu2 - 3 * fields[0] * fields[0] * (3 * g2 + gprime2 + 8 * lambda + 4 * yt * yt));
}

float_t SM::betaFunctions(size_t parIndex) const
{

    float_t* par = getParams();
    float_t mu2 = par[0];
    float_t lambda = par[1];

    float_t betaMu2 = 0.0063325739776461107 * mu2*( 12*lambda + 6*yt*yt- (3.0/2.0)*gprime2 - (9.0/2.0)*g2 );

    float_t betaLambda = 0.5*0.0063325739776461107 * (48 * lambda * lambda + (12 * yt*yt - 3*gprime2 - 9*g2 )*lambda*2 + (3.0/4.0)*gprime2*gprime2 +(3.0/2.0)*g2*gprime2 + (9.0/4.0)*g2*g2 - 12*yt*yt*yt*yt);

    if (parIndex == 0)
    {
        return betaMu2;
    }
    else if (parIndex ==1)
    {
        return betaLambda;
    }
    else
    {
        printf("parameter index out of bounds in betaFunctions, parIndex = %lu\n",parIndex);
        exit(EXIT_FAILURE);
    }
}

float_t SM::sphScale() const
{
    float_t temperature = getTemperature();
    float_t T0 = getT0();

    return VEV * sqrt(1-temperature*temperature/(T0 * T0));

}

/*double SM::detM(float_t temperature) const
{

    float_t mu2 = getParams()[0];
    float_t lambda = getParams()[1];

    return temperature*temperature*(0.1875*g2 + 0.0625 * gprime2 + 0.25*yt * yt + 0.5*lambda) - mu2;
}*/

            /*-------------------------------*/

//------------------------------------------------------------

//----------------------------Particles------------------------

float_t hboson :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return 3*par[1]*fields[0]*fields[0]  - par[0];
}

float_t* hboson :: getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T)
{
    m_fieldMassDeriv[0]=6*par[1]*fields[0];
    return m_fieldMassDeriv;
}

float_t chibosonN :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return par[1]*fields[0]*fields[0] - par[0] + xi*0.25*(g2+gprime2)*fields[0]*fields[0];

}

float_t* chibosonN :: getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T)
{
    m_fieldMassDeriv[0]=2*par[1]*fields[0] + xi*0.5*(g2+gprime2)*fields[0];
    return m_fieldMassDeriv;

}

float_t chibosonC :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return par[1]*fields[0]*fields[0] - par[0] + xi*0.25*g2*fields[0]*fields[0];

}

float_t* chibosonC :: getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T)
{
    m_fieldMassDeriv[0]=2*par[1]*fields[0] + xi*0.5*g2*fields[0];
    return m_fieldMassDeriv;

}

float_t wboson :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return 0.25*g2*fields[0]*fields[0];
}

float_t* wboson :: getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T)
{
    m_fieldMassDeriv[0]=0.5*g2*fields[0];
    return m_fieldMassDeriv;
}

float_t zboson :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return 0.25*(g2+gprime2)*fields[0]*fields[0];
}

float_t* zboson :: getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T)
{
    m_fieldMassDeriv[0]=0.5*(g2+gprime2)*fields[0];
    return m_fieldMassDeriv;
}

float_t topquark :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return 0.5 * yt * yt * fields[0] * fields[0];
}

float_t* topquark :: getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T)
{
    m_fieldMassDeriv[0]=yt * yt * fields[0];
    return m_fieldMassDeriv;
}

float_t wghost :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return xi*0.25*g2*fields[0]*fields[0];
}

float_t* wghost :: getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T)
{
    m_fieldMassDeriv[0]=xi*0.5*g2*fields[0];
    return m_fieldMassDeriv;
}

float_t zghost :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return xi*0.25*(g2+gprime2)*fields[0]*fields[0];
}

float_t* zghost :: getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T)
{
    m_fieldMassDeriv[0]=xi*0.5*(g2+gprime2)*fields[0];
    return m_fieldMassDeriv;
}

//------------------------------------------------------------
