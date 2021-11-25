// include files
#include <memory>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <time.h>
#include <future>
#include <complex>
#include <vector>
#include <omp.h>

// gsl
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_bessel.h>

// CUBATATURE LIBRARY //
#include "cuba.h"


// DataFile Reader //
#ifndef CSVREADER
#define CSVREADER
#include "IO/CSV-Reader.cpp"
#endif 


#include "boost/math/quadrature/gauss_kronrod.hpp"
#include "boost/math/quadrature/sinh_sinh.hpp"
#include "boost/math/quadrature/tanh_sinh.hpp"
#include "boost/math/quadrature/exp_sinh.hpp"
#include "boost/math/special_functions/ellint_1.hpp"
#include "boost/math/special_functions/ellint_2.hpp"

// CUBA DEFINITIONS //
#define MINEVAL 3200
#define MAXEVAL 2560000
#define KEYY 0

// VEGAS PARAMETERS //
#define NSTART 10
#define NINCREASE 500
#define NBATCH 100
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

// SUAVE PARAMETERS //
#define NNEW 8000
#define NMIN 8
#define FLATNESS 15.

/// CUHRE PARAMETERS //
// #define KEY 0


#ifndef ZETA_2
#define ZETA_2 (M_PI*M_PI/6.0)
#endif

#ifndef ZETA_3
#define ZETA_3 double(1.202056903159594)
#endif

#ifndef EULERGAMMA
#define EULERGAMMA (0.5772156649015329)
#endif

// BESSEL FUNCTIONS //
#ifndef BesselJ0
#define BesselJ0(x) gsl_sf_bessel_Jn(0,x)
#endif

#ifndef BesselJ1
#define BesselJ1(x) gsl_sf_bessel_Jn(1,x)
#endif

#ifndef BesselK0
#define BesselK0(x) gsl_sf_bessel_Kn(0,x)
#endif

#ifndef BesselK1
#define BesselK1(x) gsl_sf_bessel_Kn(1,x)
#endif

#ifndef BesselK2
#define BesselK2(x) gsl_sf_bessel_Kn(2,x)
#endif

// DIVONNE PARAMETERS //
#define KEY1 11
#define KEY2 11
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0


#define DEBYE_SCREENED_COULOMB_KERNEL   666
#define LATTICE_EQCD_KERNEL 111
#define LEADING_ORDER_KERNEL 222
#define NEXT_TO_LEADING_ORDER_KERNEL 333
#define MULTIPLE_SOFT_SCATTERING_KERNEL 999

#ifndef COLLISION_KERNEL
#define COLLISION_KERNEL LATTICE_EQCD_KERNEL 
#endif

#define FULLRATE 666
#define OPACITY 999
#ifndef FULLRATE
#define RATETYPE FULLRATE
#endif

#define KERNEL 1

#define GToGG 11
#define QToQG 21
#define GToQQ 12

#ifndef PROCESS
#define PROCESS GToGG
#endif




#if( COLLISION_KERNEL != DEBYE_SCREENED_COULOMB_KERNEL )
#define SPLITSMALLQ 1
#endif

#if( COLLISION_KERNEL == DEBYE_SCREENED_COULOMB_KERNEL || COLLISION_KERNEL == LEADING_ORDER_KERNEL )
#undef SPLITSMALLQ
#define SPLITSMALLQ 1
#endif


#define EULER 666
#define RUNGE 999
#define IMEX 111
#if( SPLITSMALLQ == 0 )
    #define TIMESTEP EULER
#endif
#if( SPLITSMALLQ == 1 )
    #define TIMESTEP IMEX
#endif


#include "EllipticFct.cpp"
#include "IO/cfile.c"
// SIMULATION PARAMETERS //


double Nc = 3.0; double  C_F=(Nc*Nc-1.0)/(2.0*Nc);double C_A=Nc; double Nf=3.0; double d_A=Nc*Nc-1.0;
const double Tf =0.5;
const double MathcalN = (ZETA_3/ZETA_2)*(1.0+0.25*Nf);
double g = std::sqrt(2.763515753);//std::sqrt(0.3*4.0*M_PI);
double Temp = 1.0;
double COne,Cz,CzB;
double mDSqr,mD,mQSqr,mDSqrOvergTSqr,mDOvergT;
double mz,mzB,mOne;
double MediumLength = 2.0;
double P = 300.0; 
int NumberOfOpenMPThreads;
const int Nx = 128; double xMin = 1e-3;double xMax = 1e3;
double zVal = 0.1;
std::vector<double> xVals(Nx),PsiInitRealVals(Nx),PsiInitImagVals(Nx),IntegrandSum(Nx);
const gsl_interp_type* InterpType = gsl_interp_cspline;
std::unique_ptr<gsl_spline,decltype(&gsl_spline_free)> PsiSplineReal(gsl_spline_alloc(InterpType,Nx),&gsl_spline_free);
std::unique_ptr<gsl_spline,decltype(&gsl_spline_free)> PsiSplineImag(gsl_spline_alloc(InterpType,Nx),&gsl_spline_free);
std::unique_ptr<gsl_spline,decltype(&gsl_spline_free)> IntegrandSumSpline(gsl_spline_alloc(InterpType,Nx),&gsl_spline_free);
gsl_interp_accel **Acc,**Acc1;

double tmax = 1.0;
double dt = 1e-3;
double tmin = 0.0;
double t = tmin;

// Harmonic Oscillator Parameters
const double wBH = 1.5, GammaE = 0.577216, Xip = std::exp( GammaE + M_PI/4.0);
const double FityLO = std::exp(-2*GammaE + 2.0), FityNLO = 737.919, FityNP = 1./0.17025234962865;
std::string Process;

// Effective Mass //
inline double Meff(double z){
    return mz/(2.0*z*P) + mzB/(2.0*(1.0-z)*P) - mOne/(2.0*P);
}
// Normalized effective MAss //
inline double Meff2PzzB(double z){
    return (1.0-z)*mz + z*mzB - z*(1.0-z)*mOne;
}
// Kinetic Energy //
inline double Momentum(double x2, double z){
    return x2/(2.0*P*z*(1.0-z));
}

// Energy //
inline double deltaE(double x2, double z){
    return Momentum(x2,z) + Meff(z) ;
}

// Dimensionless energy //
inline double Epsilon(double x2, double z){
    return x2 + Meff2PzzB(z);
}

namespace Integration{
    int NumberOfIntegrationPoints = 2048*16;
    gsl_integration_workspace ** wsp;
    gsl_integration_qawo_table **QAWOTable;
    gsl_integration_qaws_table **QAWSTable;
    const int NClenshawCurtis = 256;
    double ClenshawCurtisWeights[NClenshawCurtis];

    inline double Points(int j){
        return std::cos(double(j+1)*M_PI/double(NClenshawCurtis+1));
    }

    // gsl_integration_cquad_workspace ** wspcquad;

    void Setup(){
        wsp = new gsl_integration_workspace*[NumberOfOpenMPThreads];
        QAWOTable = new gsl_integration_qawo_table*[NumberOfOpenMPThreads];
        QAWSTable = new gsl_integration_qaws_table*[NumberOfOpenMPThreads];

        // wspcquad = new gsl_integration_cquad_workspace*[NumberOfOpenMPThreads];
        #pragma omp parallel
        {
            int tID = omp_get_thread_num();
            wsp[tID] = gsl_integration_workspace_alloc(NumberOfIntegrationPoints);
            QAWOTable[tID] = gsl_integration_qawo_table_alloc(1.0,2.0*M_PI,GSL_INTEG_COSINE,NumberOfIntegrationPoints);
            QAWSTable[tID] = gsl_integration_qaws_table_alloc(-0.5,-0.5,0.0,0.0);

            // wspcquad[tID] = gsl_integration_cquad_workspace_alloc(NumberOfIntegrationPoints);
        }
        
        // double integral = 0.0;
        for(int i=0;i<NClenshawCurtis;i++){
            ClenshawCurtisWeights[i] = 0.0;
            double ti = (i+1)*M_PI/double(NClenshawCurtis+1); 
            for(int j=0;j<NClenshawCurtis;j+=2){
                ClenshawCurtisWeights[i] += (4.0/double(NClenshawCurtis+1))*std::sin(ti)*std::sin( (j+1)*ti)/double(j+1);
                
            }
            // integral += ClenshawCurtisWeights[i]*std::exp(Points(i));
        }
        // std::cerr << "ClenshawCurtis Integral of int_{-1}^{1} exp(x) = "<< integral << "";
        // std::cerr << " Error"<< 1.0-(integral/(std::exp(1.0)-std::exp(-1.0))) << "\n";
    }
    void free(){
        #pragma omp parallel
        {
            int tID = omp_get_thread_num();
            gsl_integration_workspace_free(wsp[tID]);
            gsl_integration_qawo_table_free(QAWOTable[tID]);
            gsl_integration_qaws_table_free(QAWSTable[tID]);

            // gsl_integration_cquad_workspace_free(wspcquad[tID]);
        }
    }
}

#if(COLLISION_KERNEL == LATTICE_EQCD_KERNEL )
    #include "EQCDInMomentumSpace.cpp"
#endif

#include "BroadeningKernel.cpp"

// Broadening Kernel //
inline double Gamma(double x){
    #if(COLLISION_KERNEL == DEBYE_SCREENED_COULOMB_KERNEL )
        //return g*g*Temp*mDSqr/(x*x*(x*x + mDSqr));
        return g*g*Temp*mDSqr/((x*x + mDSqr)*(x*x + mDSqr));
    #endif

    #if(COLLISION_KERNEL == LATTICE_EQCD_KERNEL )
        return EQCD::FullGamma(x);
    #endif

    #if(COLLISION_KERNEL == LEADING_ORDER_KERNEL )
        double qOvergT = x/(g*Temp);
        return BroadeningKernel::LO::CombinedTGammaq(qOvergT);
    #endif

    #if(COLLISION_KERNEL == NEXT_TO_LEADING_ORDER_KERNEL )
        double qOvergT = x/(g*Temp);
        return BroadeningKernel::NLO::CombinedTGammaq(qOvergT);
    #endif
}

struct GSLVARIABLES{
    double p,t,z,qCutoff,muSqr;
    gsl_spline *SplineReal,*SplineImag;
    int RealImag,Which;
};

// Initial condition 2D integrand //
static int InitialIntegrand(const int *ndim, const double xx[], 
    const int *ncomp, double fval[], void *userdata){
    GSLVARIABLES &parms = *(GSLVARIABLES *) userdata;
    double x = parms.p;
    double z = parms.z;

    double yprime = (mD*M_PI/2.0)*xx[0];
    
    double y = mD*std::tan(yprime/mD);
    double phi = (2.0*M_PI)*xx[1];
    double cosphi = std::cos(phi);
    
    double xMq2   = x*x - 2.0*x*y*cosphi + y*y;
    double xPzq2  = x*x + 2.0*x*z*y*cosphi + z*z*y*y;
    double xPzBq2 = x*x + 2.0*x*(1.0-z)*y*cosphi + (1.0-z)*(1.0-z)*y*y;

    double jacobian = (mD*M_PI/2.0)*(( mDSqr + y*y)/mDSqr)*(2.0*M_PI);

    double xSqrtOverDeltaE = x*x/Epsilon(x*x,z);

    fval[0] = jacobian*y*Gamma(y)* ( COne*( xSqrtOverDeltaE - ( x*x - x*y*cosphi )/(Epsilon(xMq2,z) ) ) 
        + Cz*( xSqrtOverDeltaE - ( x*x + z*x*y*cosphi )/(Epsilon(xPzq2,z) ) )
        + CzB*( xSqrtOverDeltaE -  ( x*x + (1.0-z)*x*y*cosphi )/(Epsilon(xPzBq2,z) ) )  
        )/(4.0*M_PI*M_PI);

    return 0;
}

// Helper function for the integrand //
inline double IntegrandFct(double p2, double q2,double z){
    double PzzB = P*z*(1.0-z);
    double MeffPzzB2 = Meff2PzzB(z);
    // return p2/deltaE(p2,z) - PzzB + PzzB*( q2 + MeffPzzB2 - p2 )/std::sqrt( p2*p2 + p2*(2.0*MeffPzzB2 - 2.0*q2) + std::pow( MeffPzzB2 + q2, 2 ) );
    return 2.0*p2/Epsilon(p2,z) - 1.0 + ( q2 + MeffPzzB2 - p2 )/std::sqrt( p2*p2 + p2*(2.0*MeffPzzB2 - 2.0*q2) + std::pow( MeffPzzB2 + q2, 2 ) );
}

// Initial condition 1D integrand //
double InitialIntegrand1D(double y, void *userdata){
    
    GSLVARIABLES &parms = *(GSLVARIABLES *) userdata;
    double p = parms.p;
    double z = parms.z;
    double p2 = p*p;
    
    // double q  = mD*std::tan(y/mD);
    double q  = y;
    double q2 = q*q;

    // double jacobian = (( mDSqr + q2)/mDSqr)/(2.0*M_PI);
    double jacobian = 1.0/(2.0*M_PI);
    
    // if( q > 1e2){
    //     return jacobian*q*Gamma(q)*(2.0*p2*P*z*(1.0-z)*(COne+Cz+CzB)/( p2 + 2.0*Meff(z)*P*z*(1.0-z)));
    // }
    // else if( q <  1e-2 ){
    //     return jacobian*q*Gamma(q)*(8.0*Meff(z)*p2*P*z*(1.0-z)*q2*( COne+Cz*z*z + CzB*(1.0-z)*(1.0-z) )/std::pow( p2+2.0*Meff(z)*P*z*(1.0-z) ,3.0)) ;
    // }

    // double Result =  jacobian*q*Gamma(q)*( COne*IntegrandFct(p2,q2,z)+ Cz*IntegrandFct(p2,z*z*q2,z)+ CzB*IntegrandFct(p2,(1.0-z)*(1.0-z)*q2,z) );
    double Result =  jacobian*q*( COne*Gamma(q) + Cz/(z*z)*Gamma(q/z) + CzB*Gamma(q/(1.0-z))/((1.0-z)*(1.0-z)) ) *IntegrandFct(p2,q2,z);

    return Result;
}

// Initial condition integrated //
double InitialForceTerm(double x, double z){
    GSLVARIABLES Vars; Vars.p = x; Vars.z = z;

    double integral[1],error[1],prob[1];
    int nregions, neval,fail;
    double epsrel = 1e-8; double epsabs = 1e-8;

    // Cuhre(2, 1, InitialIntegrand, &Vars, 1,
    //         epsrel, epsabs, 0,
    //         MINEVAL, MAXEVAL, 0,
    //         STATEFILE, SPIN,
    //         &nregions, &neval, &fail, integral, error, prob);
    
    using namespace Integration;
    gsl_function Fct;
    Fct.params = &Vars;
    Fct.function = &InitialIntegrand1D;
    // gsl_set_error_handler_off();

    auto g = [&](double s) {
        return InitialIntegrand1D(s,&Vars);
    };
    integral[0] = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(g, 0, std::numeric_limits<double>::infinity(), 5);

    // gsl_set_error_handler(NULL);
    // gsl_integration_qagiu(&Fct,0.0,epsabs,epsrel,NumberOfIntegrationPoints,wsp[omp_get_thread_num()],integral,error);
    return integral[0];
}

// Evaluating the Wave fct spline, with range conditions //
inline std::complex<double> SplineEval(double x, GSLVARIABLES & parms) { 
    int tID  = omp_get_thread_num();
    if( x >= xVals[0] && x <= xVals[Nx-1] ) { 
        return std::complex<double>( gsl_spline_eval(parms.SplineReal,x,Acc[tID]),gsl_spline_eval(parms.SplineImag,x,Acc1[tID]));
    }
    // else if( x < xVals[0] ){
    //     return std::complex<double>( 0.0, x*x)*(PsiInitImagVals[0]/(xVals[0]*xVals[0]));
    // }
    // else if( x > xVals[Nx-1] ){
    //     return std::complex<double>( 0.0, 1.0/(x*x) )*(PsiInitImagVals[Nx-1]*(xVals[Nx-1]*xVals[Nx-1]));
    // }
    return std::complex<double>(0.0,0.0);
}

// Derivative of the Spline wave fct interpolation ///
inline std::complex<double> SplineDerivEval(double x, GSLVARIABLES & parms) { 
    int tID  = omp_get_thread_num();
    if( x >= xMin && x <= xMax ) { 
        return std::complex<double>( gsl_spline_eval_deriv(parms.SplineReal,x,Acc[tID]),gsl_spline_eval_deriv(parms.SplineImag,x,Acc1[tID]));
        } 
    return std::complex<double>(0.0,0.0);
}

// Second Derivative of the Spline wave fct interpolation ///
inline std::complex<double> SplineDeriv2Eval(double x, GSLVARIABLES & parms) { 
    int tID  = omp_get_thread_num();
    if( x >= xMin && x <= xMax ) { 
        return std::complex<double>( gsl_spline_eval_deriv2(parms.SplineReal,x,Acc[tID]),gsl_spline_eval_deriv2(parms.SplineImag,x,Acc1[tID]));
        } 
    return std::complex<double>(0.0,0.0);
}

// Spline evaluation gain minus loss term // 
inline std::complex<double> SplinePMinusQ(double p, double q, double cosphi, GSLVARIABLES & parms){
    double Prefactor1 = 0.0;
    double Prefactor2 = 0.0;
    double pPq2       = std::fabs(p*p + 2.0*p*q*cosphi + q*q);
    double pPq        = std::sqrt(pPq2);
    double pMq2       = std::fabs(p*p - 2.0*p*q*cosphi + q*q);
    double pMq        = std::sqrt(pMq2);
    std::complex<double> res = std::complex<double>(0.0,0.0);
    if (p >= xVals[0] && p <= xVals[Nx - 1] )
    {
        int tID = omp_get_thread_num();
        if( pPq >= xVals[0] && pPq <= xVals[Nx - 1] )
        {
            Prefactor1 = 0.5 * (p * p + p * q * cosphi) / (pPq2);
            res = 0.5*std::complex<double>(gsl_spline_eval(parms.SplineReal, p, Acc[tID]), gsl_spline_eval(parms.SplineImag, p, Acc1[tID])) 
            - Prefactor1 * std::complex<double>(std::cos((p * p - pPq2) * parms.t), std::sin((p * p - pPq2) * parms.t)) * std::complex<double>(gsl_spline_eval(parms.SplineReal, pPq, Acc[tID]), gsl_spline_eval(parms.SplineImag, pPq, Acc[tID]));
        }
        if( pMq >= xVals[0] && pMq <= xVals[Nx - 1])
        {
            Prefactor2 = 0.5 * (p * p - p * q * cosphi) / (pMq2);
            res += 0.5*std::complex<double>(gsl_spline_eval(parms.SplineReal, p, Acc[tID]), gsl_spline_eval(parms.SplineImag, p, Acc1[tID]))
            - Prefactor2 * std::complex<double>(std::cos((p * p - pMq2) * parms.t), std::sin((p * p - pMq2) * parms.t)) * std::complex<double>(gsl_spline_eval(parms.SplineReal, pMq, Acc[tID]), gsl_spline_eval(parms.SplineImag, pMq, Acc[tID]));
        }
    }
    return res;
}

// Initial condition //
std::complex<double> PsiInit(double x, double z){
    // return std::complex<double>(0.0,1.0/(4.0*P*P*z*z*(1.0-z)*(1.0-z))*InitialForceTerm(x,z));
    return std::complex<double>(0.0,InitialForceTerm(x,z));
}


void Setup(int argc, char **argv){
  

    // COMMANDLINE ARGUMENTS //
    Konfig CommandlineArguments(argc,argv);

    CommandlineArguments.Getval("z",zVal);
    CommandlineArguments.Getval("P",P);
    CommandlineArguments.Getval("T",Temp);
    CommandlineArguments.Getval("g",g);

    int n =0,p=10000;
    cubacores(&n,&p);
    NumberOfOpenMPThreads = omp_get_max_threads();
    std::cerr<< "OMPThreads = " << NumberOfOpenMPThreads << "\n";

    Integration::Setup();
    #if(COLLISION_KERNEL == DEBYE_SCREENED_COULOMB_KERNEL )
    #endif

    #if(COLLISION_KERNEL == LATTICE_EQCD_KERNEL )
        EQCD::Setup();
    #endif

    #if(COLLISION_KERNEL == LEADING_ORDER_KERNEL || COLLISION_KERNEL == NEXT_TO_LEADING_ORDER_KERNEL )
        BroadeningKernel::HardIntegralSpace::SetupHardIntegralSpline();
    #endif
    mDSqrOvergTSqr = (Nc/3.0 + Nf/6.0);
    mDOvergT = std::sqrt(mDSqrOvergTSqr);
    mDSqr = g*g*Temp*Temp*mDSqrOvergTSqr;
    mD    = std::sqrt(mDSqr);
    mQSqr = g*g*Temp*Temp*C_F/8.0;
    #if( PROCESS==GToGG )
        Process.assign("G->GG");
        mz    = mDSqr/2;
        mzB   = mDSqr/2; 
        mOne  = mDSqr/2;
        COne  = C_A/2.0;  
        Cz    = C_A/2.0; 
        CzB   = C_A/2.0;
    #endif
    #if( PROCESS==QToQG )
        Process.assign("Q->GQ");
        mz    = mDSqr/2;
        mzB   = mQSqr*2; 
        mOne  = mQSqr*2;
        COne  = C_A/2.0;  
        Cz    = (2.0*C_F - C_A)/2.0; 
        CzB   = C_A/2.0;
    #endif
    #if( PROCESS==GToQQ )
        Process.assign("G->QQ");
        mz    = mQSqr*2;
        mzB   = mQSqr*2; 
        mOne  = mDSqr/2;
        COne  = (2.0*C_F - C_A)/2.0;
        Cz    = C_A/2.0; 
        CzB   = C_A/2.0;
    #endif

    Acc  = new gsl_interp_accel*[NumberOfOpenMPThreads];
    Acc1 = new gsl_interp_accel*[NumberOfOpenMPThreads];
    #pragma omp parallel 
    {
        int tID = omp_get_thread_num();
        Acc[tID]  = gsl_interp_accel_alloc();
        Acc1[tID] = gsl_interp_accel_alloc();
    }

    // Diffusion::Check();

    std::cerr << "# Parameter : ";
    std::cerr<< "gs = " << g << " ";
    std::cerr<< "P = " << P << " ";
    std::cerr<< "T = " << Temp << " ";
    std::cerr<< "z = " << zVal << "\n";

    tmax = 10.0;
    dt = 5e-3/Temp;//(2.0*P*zVal*(1.0-zVal));
    
}


void Clear(){
    Integration::free();
}
