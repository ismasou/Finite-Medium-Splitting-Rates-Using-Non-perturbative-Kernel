
#include "Parameters.cpp"

// SETUP SIGMA P //
namespace Sigma{
    const static int Ns = 1024;
    std::array<double,Ns> pVals,SigmaVals;
    double pMin = 1e-4, pMax = 1e3;
    std::unique_ptr<gsl_spline,decltype(&gsl_spline_free)> SigmaSpline(gsl_spline_alloc(gsl_interp_cspline,Ns),&gsl_spline_free);
    gsl_interp_accel **Acc;


    double SigmaIntegrand(double x, void *p){
        double &z = *(double *) p;
        return x*(COne*Gamma(x) + Cz*Gamma(x/z)/(z*z) + CzB*Gamma(x/(1.0-z))/( (1.0-z)*(1.0-z)) )/(2.0*M_PI);
    }

    double GetSigma(double p){
        if( p>pMin && p < pMax ){
            return gsl_spline_eval(SigmaSpline.get(),p,Acc[omp_get_thread_num()]);
        }
        else if( p >= pMax ) {
            return 0.0;
        }
        
        return 0.0;
    }

    void Setup(){
        // std::cerr << "mDSqr " << mDSqr << std::endl;
        // std::cerr << "g " << g << std::endl;
        gsl_set_error_handler_off();
        gsl_function Fct; Fct.function = &SigmaIntegrand;
        Fct.params = &zVal;
        #pragma omp parallel for
        for(int ip=0;ip<Ns;ip++){
            double res,err;
            double p = pMin*std::exp( std::log(pMax/pMin) * double(ip)/double(Ns-1) );
            gsl_integration_qagiu(&Fct,p,1e-5,1e-8,Integration::NumberOfIntegrationPoints,Integration::wsp[omp_get_thread_num()],&res,&err);

            pVals[ip] = p;
            SigmaVals[ip] = res;
            // #pragma omp critical
            // std::cout << pVals[ip] << " " << SigmaVals[ip] << std::endl;
        }
        gsl_spline_init(SigmaSpline.get(),pVals.data(),SigmaVals.data(),Ns);
        Acc = new gsl_interp_accel*[NumberOfOpenMPThreads];
        #pragma omp parallel 
        {
            Acc[omp_get_thread_num()] = gsl_interp_accel_alloc();
        }
        gsl_set_error_handler(NULL);
    }


}


// Opacity Integration //
namespace Opacity{

    static int InitialIntegrand(const int *ndim, const double xx[], 
    const int *ncomp, double fval[], void *userdata){
        GSLVARIABLES &parms = *(GSLVARIABLES *) userdata;
        double z = parms.z;


        double xprime = (M_PI/2.0)*xx[0];
        double yprime = (M_PI/2.0)*xx[1];
        
        double p = mD*std::tan(xprime);
        double q = mD*std::tan(yprime);
        double p2 = p*p;
        double q2 = q*q;
        double jacobian = (mD*M_PI/2.0)*(( mDSqr + q2)/mDSqr)*(mD*M_PI/2.0)*(( mDSqr + p2)/mDSqr);

        double EpsilonP = Epsilon(p2,z);
        double SigmaP   = Sigma::GetSigma(p)*(2*P*z*(1.0-z));
        // double SigmaP = (2*P*z*(1.0-z))*(-1.2461732172989273*std::log(p) + 0.6230866086494636*std::log(4.145269999999999 + 1.*p*p));
        double Prefactor = -std::imag( std::complex<double>( 1.0 - std::cos(-EpsilonP*t)*std::exp( -SigmaP*t ) , -std::sin(-EpsilonP*t)*std::exp( -SigmaP*t ))/std::complex<double>( SigmaP, EpsilonP) );
        double pSqrtOverDeltaE = p2/EpsilonP;
        if(xx[0] >= (1.0 - 1e-5) || xx[1] >= (1.0 - 1e-5)){
            fval[0] = 0.0;

        }
        else{
            fval[0] = 0.5*Prefactor*p*jacobian*q*Gamma(q)*( COne*IntegrandFct(p2,q2,z)+ Cz*IntegrandFct(p2,z*z*q2,z)+ CzB*IntegrandFct(p2,(1.0-z)*(1.0-z)*q2,z) )/(4.0*M_PI*M_PI);

        }
        return 1;
    }

    double InitialForceTerm(double z){
        GSLVARIABLES Vars; Vars.z = z;

        double integral[1],error[1],prob[1];
        int nregions, neval,fail;
        double epsrel = 1e-5; double epsabs = 1e-5;

        Cuhre(2, 1, InitialIntegrand, &Vars, 1,
                epsrel, epsabs, 0,
                MINEVAL, MAXEVAL, 0,
                STATEFILE, SPIN,
                &nregions, &neval, &fail, integral, error, prob);

        if(fail == 1){
            std::cerr << "FAILED CUHRE INTEGRAL AT t = " << t << "\n";
            exit(0);
        }
        
        return (g*g/M_PI)*integral[0];
    }

    void OUTPUTRate(){
        double Result = InitialForceTerm(zVal);
        std::stringstream RateFolder;
        RateFolder << "OpacityImproved/Rate-P" << P << "-z" << zVal << ".txt";

        std::ofstream RateFile; 
        if( t == tmin ){
            RateFile << "# Parameter : " << Process << " gs = " << g << " " << "P = " << P << " " << "T = " << Temp << " " << "z = " << zVal << " mDSqr= " << mDSqr << "\n";
            RateFile << "# Columns: 1- Time: T^2 t/(P*z*(1-z))    2- 1/(T P(z)) dGamma/dz \n";
            RateFile.open(RateFolder.str().c_str(),std::ofstream::out);
        }
        else{
            RateFile.open(RateFolder.str().c_str(),std::ofstream::app);
        }
        RateFile << t << " " << Result << "\n";
        RateFile.close();
    }

    void Evolve(){
        int i = 0;
        dt *= 5;
        tmax *= 1.0; 
        std::cout << "gs = "<< g << "\n";
        std::cout << "mDSqr = "<< mDSqr << "\n";
        while( t <= tmax ){

            OUTPUTRate();
            if( i % 1 ==  0){
                std::cerr << "# Done t = " << t << "\n";
            }
            t+=dt;
            i++;
        }

    }
}

int main(int argc, char **argv){
  

    // COMMANDLINE ARGUMENTS //
    Konfig CommandlineArguments(argc,argv);

    CommandlineArguments.Getval("z",zVal);
    CommandlineArguments.Getval("P",P);

    int n =0,p=10000;
    cubacores(&n,&p);
    NumberOfOpenMPThreads = omp_get_max_threads();
    std::cerr<< "OMPThreads = " << NumberOfOpenMPThreads << "\n";

    Integration::Setup();
    #if(COLLISION_KERNEL == DEBYE_SCREENED_COULOMB_KERNEL )
        std::cerr << "# Using " << "DEBYE_SCREENED_COULOMB_KERNEL" << std::endl;
    #endif

    #if(COLLISION_KERNEL == LATTICE_EQCD_KERNEL )
        std::cerr << "# Using " << "LATTICE_EQCD_KERNEL" << std::endl;
        EQCD::Setup();
    #endif

    #if(COLLISION_KERNEL == LEADING_ORDER_KERNEL )
        std::cerr << "# Using " << "LEADING_ORDER_KERNEL" << std::endl;
        HardIntegralSpace::SetupHardIntegralSpline();
    #endif

    #if(COLLISION_KERNEL == NEXT_TO_LEADING_ORDER_KERNEL )
        std::cerr << "# Using " << "NEXT_TO_LEADING_ORDER_KERNEL" << std::endl;
        HardIntegralSpace::SetupHardIntegralSpline();
    #endif


    mDSqrOvergTSqr = (Nc/3.0 + Nf/6.0);
    mDOvergT = std::sqrt(mDSqrOvergTSqr);
    mDSqr = g*g*Temp*Temp*mDSqrOvergTSqr;
    mD    = std::sqrt(mDSqr);
    mQSqr = g*g*Temp*Temp*C_F/8.0;
    #if( PROCESS==GToGG )
        mz    = mDSqr/2;
        mzB   = mDSqr/2; 
        mOne  = mDSqr/2;
        COne  = C_A/2.0;  
        Cz    = C_A/2.0; 
        CzB   = C_A/2.0;
    #endif
    #if( PROCESS==QToQG )
        mz    = mDSqr/2;
        mzB   = mQSqr*2; 
        mOne  = mQSqr*2;
        COne  = C_A/2.0;  
        Cz    = (2.0*C_F - C_A)/2.0; 
        CzB   = C_A/2.0;
    #endif
    #if( PROCESS==GToQQ )
        mz    = mQSqr*2;
        mzB   = mQSqr*2; 
        mOne  = mDSqr/2;
        COne  = (2.0*C_F - C_A)/2.0;
        Cz    = C_A/2.0; 
        CzB   = C_A/2.0;
    #endif

    std::cerr << "# Parameter : ";
    std::cerr << "P = " << P << " ";
    std::cerr << "T = " << Temp << " ";
    std::cerr << "z = " << zVal << std::endl;
    std::cerr << "# Setup Sigma Spline " << std::endl;
    Sigma::Setup();
    std::cerr << "# Initialization Done " << std::endl;

    Opacity::Evolve();

    Integration::free();
    return 0;
}
