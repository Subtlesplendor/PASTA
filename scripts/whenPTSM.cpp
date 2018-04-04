/*****
 *  Phase Transition Analyzer (PTA)
 *
 *
 *
 *  Author: Johan Lofgren
 *
 *
 **/

//Needed for analyzing the phase transition:
#include "../bin/PTAfunctions.h"
//Include your model here:
#include "../models/SM.h"

int main(int argc, char **argv)
{
    /*
        This is an example program that calculates when the electroweak phase transition occurs in the standard model. It can be modified to work for the model of your choice, or you can tweak the settings to see how the result is affected.

        Author: Johan Löfgren

    */
    //flint_cleanup();
    //flint_set_num_threads(1);
    //For Errorhandling of input:
    if(argc > 1)
    {
        printf("usage: no input\n");
        return 0;
    }

    //For measuring runtime:
    clock_t start, end;
    //const int N = 1e4;
    const int N = 1;

    //-------------------Options for the potential: -----------------
    const size_t scheme = 1;
    const float_t scale = 100*100;

    //The parameters of the potential of the SM, (mu2, lambda):
    float_t mh = 65.0855;
    float_t mh2 = mh * mh;
    float_t mu2 = 0.5 * mh2;
    float_t lambda = 0.5 * mh2 / VEV2;
    printf("mu2 = %f, lambda = %f\n",mu2,lambda);
    float_t par[2] = {mu2, lambda};
    const size_t numParameters = sizeof(par)/sizeof(par[0]);

    const size_t numVar = 1; // Only one field obtains a vev

    // The particles included in the effective potential
    const size_t numParticles = 8;
    Particle* particles[numParticles];

    particles[0] = new hboson(mh);
    particles[1] = new chibosonN(0.0);
    particles[2] = new chibosonC(0.0);
    particles[3] = new wboson(mw);
    particles[4] = new zboson(mz);
    particles[5] = new topquark(sqrt(mtPhys2));
    particles[6] = new wghost(0.0);
    particles[7] = new zghost(0.0);

    //Creating the potential:
    Model* mod;
    mod = new SM(par, numParameters, particles, numParticles, scale, scheme);
    float_t treeMin[numVar] = {};
    for (size_t i = 0; i < numVar; ++i)
    {
        treeMin[i] = mod->getTreemin()[i];
    }

    //-------------------------------------------------------------


    //------------------Options for minimizer--------------------
    minOptions_t options;
    //Uncomment a line and change it to modify an option.
    //options.stepSize = 1;
    //options.maxTemp = 1000;
    //options.minTemp = 0;
    //options.tol = 1e-3;
    //options.locTol = 4 * options.tol;
    //options.tempTol = options.tol;
    //options.maxIter = 1000;
    //----------------------------------------------------------


    //If you wish to study gauge dependence, you can try changing xi here.
    mod->setGaugeXi(0.0);
    float_t xi = mod->getGaugeXi();
    printf("xi = %f\n",xi);

    //-----------------Running the actual calculation--------
    size_t errorFlagL = 0;
    size_t errorFlagPRM = 0;
    float_t TcL;
    float_t vcL;
    float_t TcPRM;
    float_t vcPRM;
    start = clock();
    float_t T0=0;
    for (int i = 0; i <N; ++i)
    {
        T0 = mod->getT0();
        whenPT_L_FAST(mod, TcL, vcL, errorFlagL, options);
        whenPT_PRM(mod, TcPRM, vcPRM, errorFlagPRM);
    }
    end = clock();
    printtime(start, end);
    //-------------------------------------------------------

    printf("T0 = %f\n",T0);

    size_t T0errorFlag = mod->findT0();
    printf("T0errorflag = %lu\n",T0errorFlag );

    printf("------------------------------------------------\n");
    printf("LANDAU errorflag = %lu\n",errorFlagL );
    printf("LANDAU phase transition temperature = %f GeV\n",TcL );
    printf("LANDAU phase transition VEV = %f GeV\n",vcL );
    printf("LANDAU vc / Tc = %f\n",vcL/TcL );
    printf("------------------------------------------------\n");
    printf("PRM errorflag = %lu\n",errorFlagPRM );
    printf("PRM phase transition temperature = %f GeV\n",TcPRM );
    printf("PRM phase transition VEV = %f GeV\n",vcPRM );
    printf("PRM vc / Tc = %f\n",vcPRM/TcPRM );
    printf("------------------------------------------------\n");


    //--------------------Cleaning up---------------------------
    delete mod;
    for (size_t i = 0; i < numParticles; ++i)
    {
        delete particles[i];
    }
    //-------------------------------------------------------
    flint_cleanup();

    return 0;
}
