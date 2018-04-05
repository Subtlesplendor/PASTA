#ifndef ABSTRACTMODEL_H
#define ABSTRACTMODEL_H


#include "common.h"
#include <cmath>

/*This file contains the declarations of the abstract classes from which
the actual model implementations inherit. In order, it contains:

    1. The particle class
    2. The potential class

*/


class Particle
{
private:
    const int m_dof;
    const size_t m_numVar;
    const float_t m_spin;
    const bool m_isGB;

    int m_sign;
    float_t m_mass;

public:
    Particle(float_t mass, float_t spin, int dof, size_t numVar, bool isGB);
    virtual ~Particle(){delete m_fieldMassDeriv;};

    //Setters:
    void setParticle(float_t mass);
    void setSign(int sign) {m_sign = sign;};

    //Getters:
    float_t getMass() const { return m_mass; }
    float_t getSpin() const { return m_spin; }
    int getSign() const { return m_sign; }
    int getDof() const { return m_dof; }
    bool getGB() const {return m_isGB; }

    //Modifiers:
    //virtual void findFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T, const size_t numVar); //CAN be overwritten in the inheriting class

    //Public variable (beware...)
    float_t* m_fieldMassDeriv;

    //Numerical functions:
    virtual float_t getFieldmass(const float_t* par, const float_t xi, const float_t* fields, const float_t T) const = 0; //MUST be overwritten in the inheriting class
    virtual float_t* getFieldmassDeriv(const float_t* par, const float_t xi, const float_t* fields, const float_t T);
};


class Model
{
private:
    //Initialized when a Model-object is created and cannot be changed while a program is running:
    const size_t m_scheme; // the chosen renormalization scheme, see manual.
    const size_t m_numVar; // the number of fields which can attain an EW vev.
    const size_t m_numParameters; // the number of parameters of the potential
    const size_t m_numParticles; // the number of particles included in the loops


    //Can be changed:
    float_t m_renormScale2;
    float_t m_gaugeXi;
    float_t* m_refScales2; //reference scales at which parameters are measured
    Particle** m_particles; // holds all the particles to be included in the loops
    float_t* m_pars; // parameters
    float_t* m_ctPars; // counterterms 'parameters'
    float_t* m_renPars; // parameters + counterterms
    float_t m_temperature; // temperature of the potential
    float_t m_T0;           //Temperature of second order PT in HT eff pot
    float_t* m_treeMin; // the minimum at tree-level
    float_t* m_totMin;  // the minimum at loop-level (T = 0)
    float_t* m_treeDeriv;
    float_t* m_loopDeriv;
    float_t* m_totDeriv;

    //float_t* m_potDeriv; // The derivative of the one loop potential -- for the self energy

public:
    Model(float_t* par, const size_t numParameters, Particle* particles[], const size_t numParticles,const float_t reScale2, const size_t scheme, const size_t numVariables);

    virtual ~Model();

    //Getters:
    const size_t getScheme() const {return m_scheme;};
    const size_t getNumvar() const {return m_numVar;};
    const size_t getNumparams() const {return m_numParameters;};
    const size_t getNumparticles() const {return m_numParticles;};
    const float_t getRenormscale() const {return m_renormScale2;};
    Particle** getParticles() const {return m_particles;};
    float_t* getParams() const {return m_pars;};
    float_t* getCtparams() const {return m_ctPars;};
    float_t* getRenparams() const {return m_renPars;};
    float_t getTemperature() const {return m_temperature;};
    float_t getGaugeXi() const {return m_gaugeXi;};
    float_t getT0() const {return m_T0;};
    float_t* getTreemin() const {return m_treeMin;};
    float_t* getTotmin() const {return m_treeMin;};
    float_t* getRefscales() const {return m_refScales2;};
    float_t* getTreeDeriv() const {return m_treeDeriv;};
    float_t* getLoopDeriv() const {return m_loopDeriv;};
    float_t* getTotDeriv() const {return m_totDeriv;};

    //Setters:
    void setParticles(Particle* particles[]);
    void setParameters(float_t* par);
    void setCtparams(float_t* ctpars);
    void setRenpars(float_t* renpars);
    void setTemperature(const float_t T){m_temperature = T;};
    void setGaugeXi(const float_t xi){m_gaugeXi = xi;};
    void setT0(const float_t T0){m_T0 = T0;};
    void setTreemin(float_t* treemin);
    void setTotmin(float_t* totmin);
    void setTreeDeriv(float_t* treeDeriv);
    void setLoopDeriv(float_t* loopDeriv);
    void setTotDeriv(float_t* totDeriv);
    void setRenormscale(const float_t reScale2){m_renormScale2 = reScale2;};
    void setRefscales(float_t* refScales2);
    virtual void setCounterterms();

    //Modifiers:
    void updateModel(float_t* par);
    virtual size_t findT0();
    virtual void findSphscale(float_t Tc, float_t& vc, size_t& errorFlag);
    //virtual void findOneLoopLevelDeriv(float_t* fields);
    virtual void findTreemin() = 0; //MUST be overwritten in the inheriting class
    virtual void findTreePotDeriv(float_t* fields);
    void findOneLoopLevelDeriv(float_t* fields);
    void findTotalPotDeriv(float_t* fields);

    //Numerical functions:
    float_t gbSelfEnergy(float_t* fields,size_t dirIndex);
    float_t oneLoopLevel(float_t* fields) const;
    float_t oneLoopLevelResummed(float_t* fields); //Needs work.
    float_t tempPot(const float_t x, bool fermionFlag) const;
    float_t tempPotDeriv(const float_t x, bool fermionFlag) const;
    float_t totalPot(float_t* fields) const;
    float_t deltaPot(float_t* fields) const;
    float_t deltaPotResummed(float_t* fields);   //Needs work.
    virtual float_t betaFunctions(size_t parIndex) const;
    virtual double detM(float_t T) const;
    virtual float_t treeLevel(float_t* fields) const = 0; //MUST be overwritten in the inheriting class
    virtual float_t htPot(float_t* fields) const = 0; //MUST be overwritten in the inheriting class
    float_t getphi1();

};


#endif
