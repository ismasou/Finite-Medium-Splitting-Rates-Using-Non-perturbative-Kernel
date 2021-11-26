#include "Parameters.cpp"

// Opacity Integration //
namespace Opacity{
    static int IntegrandN1(const int *ndim, const double xx[], 
    const int *ncomp, double fval[], void *userdata){
        GSLVARIABLES &parms = *(GSLVARIABLES *) userdata;
        double z = parms.z;


        double xprime = (M_PI/2.0)*xx[0];
        double yprime = (M_PI/2.0)*xx[1];
        
        double p = mD*std::tan(xprime);
        double q = mD*std::tan(yprime);
        double p2 = p*p;
        double q2 = q*q;
        // double phi = (2.0*M_PI)*xx[2];
        // double cosphi = std::cos(phi);
        
        double jacobian = (mD*M_PI/2.0)*(( mDSqr + q2)/mDSqr)*(mD*M_PI/2.0)*(( mDSqr + p2)/mDSqr);

        double pSqrtOverDeltaE = p2/Epsilon(p2,z);
        if(xx[0] >= (1.0 - 1e-5) || xx[1] >= (1.0 - 1e-5)){
            fval[0] = 0.0;

        }
        else{
            // fval[0] = ((1.0 - std::cos(-Epsilon(p*p,z)*t))/Epsilon(p*p,z))*jacobian*p*q*Gamma(q)* ( COne*( pSqrtOverDeltaE - ( p*p - p*q*cosphi )/(Epsilon(( p*p - 2.0*p*q*cosphi + q*q ),z) ) ) 
            // + Cz*( pSqrtOverDeltaE - ( p*p + z*p*q*cosphi )/(Epsilon(( p*p + 2.0*p*z*q*cosphi + z*z*q*q ),z) ) )
            // + CzB*( pSqrtOverDeltaE -  ( p*p + (1.0-z)*p*q*cosphi )/(Epsilon(( p*p + 2.0*p*(1.0-z)*q*cosphi + (1.0-z)*(1.0-z)*q*q ),z) ) )  
            // )/(8.0*M_PI*M_PI*M_PI);

            //  - std::cos(Epsilon(p2,z)*t)
            fval[0] = 0.5*((1.0 - std::cos(Epsilon(p2,z)*t) )/Epsilon(p2,z))*p*jacobian*q*Gamma(q)*( COne*IntegrandFct(p2,q2,z)+ Cz*IntegrandFct(p2,z*z*q2,z)+ CzB*IntegrandFct(p2,(1.0-z)*(1.0-z)*q2,z) )/(4.0*M_PI*M_PI);
            // if(isnan(fval[0])){
            //     std::cerr << (xx[0] >= 1.0 - 1e-9 || xx[1] >= 1.0 - 1e-9) << "\n";
            //     std::cerr << xx[0] << "  " << xx[1] << "\n";
            //     exit(0);
            // }

        }
        return 1;
    }

    double OpacityN1(double z){
        GSLVARIABLES Vars; Vars.z = z;

        double integral[1],error[1],prob[1];
        int nregions, neval,fail;
        double epsrel = 1e-5; double epsabs = 1e-5;

        Cuhre(2, 1, IntegrandN1, &Vars, 1,
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

    static int IntegrandN2(const int *ndim, const double xx[], 
    const int *ncomp, double fval[], void *userdata){
        GSLVARIABLES &parms = *(GSLVARIABLES *) userdata;
        double z = parms.z;


        double xprime = (M_PI/2.0)*xx[0];
        double yprime = (M_PI/2.0)*xx[1];
        
        double p = mD*std::tan(xprime);
        double q = mD*std::tan(yprime);
        double p2 = p*p;
        double q2 = q*q;
        double phi = (2.0*M_PI)*xx[2];
        double cosphi = std::cos(phi);

        // double dt = t*xx[3];
        
        double jacobian = (2.0*M_PI)*(mD*M_PI/2.0)*(( mDSqr + q2)/mDSqr)*(mD*M_PI/2.0)*(( mDSqr + p2)/mDSqr);

        double pMq = std::sqrt(std::abs( p2 + q2 - 2.0*p*q*cosphi));

        auto I = std::complex<double>(0.0,1.0);
        auto Psi1p = std::complex<double>(0.0,0.0);
        auto Psi1pMq = std::complex<double>(0.0,0.0);
        
        int tID  = omp_get_thread_num();
        if( p >= xVals[0] && p <= xVals[Nx-1] ) { 
            Psi1p = I*gsl_spline_eval(PsiSplineImag.get(),p,Acc[tID]);
        }

        if( pMq >= xVals[0] && pMq <= xVals[Nx-1] ) { 
            Psi1pMq = I*gsl_spline_eval(PsiSplineImag.get(),pMq,Acc[tID]);
        }

        double Gam = COne*Gamma(q) + Cz*Gamma(q/z)/(z*z) + CzB*Gamma(q/z)/((1.-z)*(1.-z));

        auto EpsP   = Epsilon(p2,z);
        auto EpsPMQ = Epsilon(pMq*pMq,z);
        auto expP   = std::complex<double>( std::cos(-EpsP*t)  , std::sin(-EpsP*t));
        auto expPMQ = std::complex<double>( std::cos(-EpsPMQ*t), std::sin(-EpsPMQ*t));

        auto fact1  = ( -1.0 + expP*( 1.0 + I*EpsP*t))/( EpsP*EpsP );
        auto fact2  = ( EpsP*( expPMQ - 1.0 ) - EpsPMQ*( expP -1.0 ) )/( EpsP*EpsPMQ*( EpsP - EpsPMQ ) );

        // auto fact1  = expP*dt;
        // auto fact2  = ( expPMQ - expP)/( I*( EpsP - EpsPMQ) );

        if(xx[0] >= (1.0 - 1e-5) || xx[1] >= (1.0 - 1e-5) || pMq < 1e-10 || std::abs(EpsP - EpsPMQ) < 1e-10 || dt < 1e-15){
            fval[0] = 0.0;

        }
        else{
            fval[0] = -jacobian*p*q*Gam*std::real( Psi1p*fact1  - Psi1pMq*fact2*( p2-p*q*cosphi )/(pMq*pMq) )/(8.0*M_PI*M_PI*M_PI);
        }

        
        if(!std::isfinite(fval[0])){
            std::cerr << fval[0] << "\n";
            std::cerr << " q = " << q << " ";
            std::cerr << " p = " << p << "\n";
            exit(0);
        }

        return 1;
    }

    double OpacityN2(double z){
        GSLVARIABLES Vars; Vars.z = z;

        double integral[1],error[1],prob[1];
        int nregions, neval,fail;
        double epsrel = 1e-5; double epsabs = 1e-5;

        Cuhre(3, 1, IntegrandN2, &Vars, 1,
                epsrel, epsabs, 0,
                MINEVAL, MAXEVAL, 0,
                STATEFILE, SPIN,
                &nregions, &neval, &fail, integral, error, prob);

        if(fail == 1){
            std::cerr << "FAILED CUHRE INTEGRAL AT t = " << t << "\n";
            exit(0);
        }
        
        return (g*g/M_PI)*(2.0*P*z*(1.0-z))*integral[0];
    }


    void OUTPUTRate(){
        double Result  = OpacityN1(zVal);
        // double Result2 = OpacityN2(zVal);
        std::stringstream RateFolder;
        RateFolder << "Opacity/Rate-P" << P << "-z" << zVal << ".txt";

        std::ofstream RateFile; 
        if( t == tmin ){
            RateFile.open(RateFolder.str().c_str(),std::ofstream::out);
            RateFile << "# Parameter : "  << Process << " gs = " << g << " " << "P = " << P << " " << "T = " << Temp << " " << "z = " << zVal << "\n";
            RateFile << "# Columns: 1- Time: T^2 t/(P*z*(1-z))    2- 1/(T P(z)) dGamma/dz \n";
        }
        else{
            RateFile.open(RateFolder.str().c_str(),std::ofstream::app);
        }
        // RateFile << t << " " << Result << " " << Result2 << "\n";
        RateFile << t << " " << Result << "\n";
        RateFile.close();
    }

    void Evolve(){
        int i = 0;
        dt = 0.001; 
        tmax *= 10.0; 
        std::cout << "gs = "<< g << "\n";
        std::cout << "mDSqr = "<< mDSqr << "\n";
        while( t <= tmax ){

            OUTPUTRate();
            if( i % 1 ==  0){
                std::cerr << "# Done t = " << t << "\n";
            }
            t+= dt*std::exp( i/128.0*std::log(tmax/dt) ) ;
            i++;
        }

    }
}

int main(int argc, char **argv){
    Setup(argc,argv);

    for(int i=0;i<Nx;i++){
        int j = Nx-1-i;
        xVals[i] = xMin*std::exp(std::log(xMax/xMin)*i/(Nx-1));
        // xVals[i] = xMin + (xMax-xMin)*(std::cos(M_PI*j/(Nx-1)) + 1.0)/2.0;
        std::complex<double> Psi = PsiInit(xVals[i],zVal);
        PsiInitRealVals[i] = Psi.real();
        PsiInitImagVals[i] = Psi.imag();
        IntegrandSum[i]    = 0.0;
    }
    xVals[Nx-1] = xMax;
    
    gsl_spline_init(PsiSplineReal.get(),xVals.data(),PsiInitRealVals.data(),Nx);
    gsl_spline_init(PsiSplineImag.get(),xVals.data(),PsiInitImagVals.data(),Nx);
    gsl_spline_init(IntegrandSumSpline.get(),xVals.data(),IntegrandSum.data(),Nx);

    std::cerr << "# Initialization Done \n";
    Opacity::OUTPUTRate();
    std::cerr << "# OUTPUT Initial Condition \n";

    Opacity::Evolve();

    Clear();
    return 0;
}
