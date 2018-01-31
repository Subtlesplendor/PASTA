#ifndef SM_H
#define SM_H

#include "../src/AbstractModel.h"

/*This is the model file for the standard model (SM). It contains, in the following order,

    1. Definitions of physical constants
    2. Useful flags used in the model file
    3. The SM class
    4. The particle classes

*/

//----------------Physical constants--------------------------------
const float_t mtPhys2 = 175 * 175;
const float_t mw = 80.4;
const float_t mz  = 91.2;
const float_t mw2 = mw*mw;
const float_t mz2  = mz*mz;
const float_t mH = 125.09;
const float_t mH2 = mH*mH;
const float_t yt = sqrt(2 * mtPhys2/VEV2);
const float_t g2 = 0.653659*0.653659;
const float_t gprime2 = 0.349998*0.349998;
//----------------------------------------------------------------

//-----------------------Useful flags------------------------
const float_t SCALAR = 0;
const float_t FERMION = 0.5;
const float_t VECTOR = 1;

const size_t NUMVARSM = 1;
//-----------------------------------------------------


//----------------------The model---------------------
class SM : public Model
{
private:

public:

	SM(float_t* par, const size_t numParameters, Particle* particles[], const size_t numParticles, const float_t reScale2, const size_t scheme) : Model(par, numParameters, particles, numParticles, reScale2, scheme, NUMVARSM)
	{
		findTreemin();
        size_t errorFlag = findT0();
        setCounterterms();
	};

	virtual ~SM(){};

    //Setters:
    //void setCounterterms();

    //Modifiers:
    void findTreemin();
    size_t findT0();
    void findSphscale(float_t Tc, float_t& vc, size_t& errorFlag);
	void findTreePotDeriv(float_t* fields);

    //Numerical functions:
    float_t treeLevel(float_t* fields) const;
    float_t htPot(float_t* fields) const;
    //double detM(float_t T) const;
    float_t betaFunctions(size_t parIndex) const;
    float_t sphScale() const;

};

//---------------------------------------------------------

//----------------------The particles----------------------

class hboson : public Particle
{
private:
	//char* m_name;

public:
    hboson(float_t mass = 0.0): Particle(mass, SCALAR, 1, NUMVARSM, false)
    {
    }

    //virtual ~hboson(){};

 	virtual float_t getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const;
    virtual float_t* getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T);
};


class chibosonN : public Particle
{
private:
	//char* m_name;

public:
    chibosonN(float_t mass = 0.0): Particle(mass, SCALAR, 1, NUMVARSM, true)
    {
    }

    //virtual ~chibosonN(){};

 	virtual float_t getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const;
    virtual float_t* getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T);
};

class chibosonC : public Particle
{
private:
    //char* m_name;

public:
    chibosonC(float_t mass = 0.0): Particle(mass, SCALAR, 2, NUMVARSM, true)
    {
    }

    //virtual ~chibosonC(){};

    virtual float_t getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const;
    virtual float_t* getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T);
};


class wboson : public Particle
{
private:
	//char* m_name;

public:
    wboson(float_t mass = 0.0): Particle(mass, VECTOR, 6, NUMVARSM, false)
    {
    }

    //virtual ~wboson(){};

 	virtual float_t getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const;
    virtual float_t* getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T);
};


class zboson : public Particle
{
private:
	//char* m_name;

public:
    zboson(float_t mass = 0.0): Particle(mass, VECTOR, 3, NUMVARSM, false)
    {
    }

   	//virtual ~zboson(){};

 	virtual float_t getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const;
    virtual float_t* getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T);
};


class topquark : public Particle
{
private:
	//char* m_name;

public:
    topquark(float_t mass = 0.0): Particle(mass, FERMION, 12, NUMVARSM, false)
    {
    }

    //virtual ~topquark(){};

 	virtual float_t getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const;
    virtual float_t* getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T);
};

class wghost : public Particle
{
public:
    wghost(float_t mass = 0.0): Particle(mass, SCALAR, 2, NUMVARSM, false)
    {
        this->setSign(-1);
    }

    //virtual ~wghost(){};

    virtual float_t getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const;
    virtual float_t* getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T);
};

class zghost : public Particle
{
private:
    //char* m_name;

public:
    zghost(float_t mass = 0.0): Particle(mass, SCALAR, 1, NUMVARSM, false)
    {
        this->setSign(-1);
    }

    //virtual ~zghost(){};

    virtual float_t getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const;
    virtual float_t* getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T);
};

//-----------------------------------------------------

#endif
