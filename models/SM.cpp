#include "SM.h"

//-------------------------SM-----------------------------

            /*-----------Setters------------*/


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

    float_t mu2 = getParams()[0];
    float_t lambda = getParams()[1];

    return -0.5*mu2 * fields[0]*fields[0] + 0.25*lambda*fields[0]*fields[0]*fields[0]*fields[0];
}

float_t SM::htPot(float_t* fields) const
{

    float_t mu2 = getParams()[0];
    float_t lambda = getParams()[1];

    float_t temperature = getTemperature();

    return treeLevel(fields) - (temperature * temperature/96.0) *(16 * mu2 - 3 * fields[0] * fields[0] * (3 * g2 + gprime2 + 8 * lambda + 4 * yt * yt));
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



float_t chibosonN :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return par[1]*fields[0]*fields[0] - par[0] + xi*0.25*(g2+gprime2)*fields[0]*fields[0];

}

float_t chibosonC :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return par[1]*fields[0]*fields[0] - par[0] + xi*0.25*g2*fields[0]*fields[0];

}


float_t wboson :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return 0.25*g2*fields[0]*fields[0];
}


float_t zboson :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return 0.25*(g2+gprime2)*fields[0]*fields[0];
}

float_t topquark :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return 0.5 * yt * yt * fields[0] * fields[0];
}

float_t wghost :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return xi*0.25*g2*fields[0]*fields[0];
}

float_t zghost :: getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const
{

    return xi*0.25*(g2+gprime2)*fields[0]*fields[0];
}

//------------------------------------------------------------
