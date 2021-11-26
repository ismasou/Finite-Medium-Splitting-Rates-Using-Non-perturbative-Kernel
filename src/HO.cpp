
inline double GammaPerturb( double b );
#include "Parameters.cpp"
auto I = std::complex<double>(0.0,1.0);

namespace LeadingOrder{

    double muSqr;
    double Qs,q0,qhat;
    std::complex<double> Omega;
    double qhatFct(double QSqr, double z){
        return q0*(COne + Cz*z*z + CzB*(1.0-z)*(1.0-z) )*std::log( QSqr / muSqr);
    }
    double QsFct(double qhat, double P, double z){
        return std::sqrt( P*z*(1.0-z)*qhat );
    }

    double Spectrum(){

        double RealOmega = std::real(Omega);
        double ImagOmega = std::imag(Omega);
        // Cos[a + I b] = Cos[a] Cosh[b] - I*Sin[a] Sinh[b]
        // Sin[a + I b] = Cosh[b] Sin[a] + I*Cos[a] Sinh[b]
        double Coshb = std::cosh(ImagOmega*t), Sinhb = std::sinh(ImagOmega*t);
        double Cosa  = std::cos(RealOmega*t),  Sina  = std::sin(RealOmega*t);
        auto CosOmegat  = Cosa*Coshb - I*Sina*Sinhb;
        return g*g/(4.0*M_PI*M_PI)*std::log( std::abs(CosOmegat));
    }

    double Rate(double t){

        double RealOmega = std::real(Omega);
        double ImagOmega = std::imag(Omega);
        // Cos[a + I b] = Cos[a] Cosh[b] - I*Sin[a] Sinh[b]
        // Sin[a + I b] = Cosh[b] Sin[a] + I*Cos[a] Sinh[b]
        double Coshb = std::cosh(ImagOmega*t), Sinhb = std::sinh(ImagOmega*t);
        double Cosa  = std::cos(RealOmega*t),  Sina  = std::sin(RealOmega*t);
        auto CosOmegat  = Cosa*Coshb - I*Sina*Sinhb;
        auto SinOmegat  = Sina*Coshb + I*Cosa*Sinhb;
        return -g*g/(4.0*M_PI*M_PI)*std::real( Omega*(SinOmegat/CosOmegat));
    }

    void Setup(){
        // muSqr = mDSqr/4.0*std::exp(-2.0+2.0*GammaE);
        #if(COLLISION_KERNEL == DEBYE_SCREENED_COULOMB_KERNEL )
            muSqr = mDSqr/4.0*std::exp(-2.0+2.0*GammaE);
            q0 = g*g*Temp*mDSqr/(4.0*M_PI);
        #endif
        #if(COLLISION_KERNEL == LATTICE_EQCD_KERNEL )
            muSqr = mDSqr/4.0*std::exp(-2.0+2.0*GammaE);
            q0 = g*g*g*g*Temp*Temp*Temp*MathcalN/(4.0*M_PI);
        #endif

        #if(COLLISION_KERNEL == LEADING_ORDER_KERNEL )
            muSqr = mDSqr/4.0*std::exp(-2.0+2.0*GammaE);
            q0 = g*g*g*g*Temp*Temp*Temp*MathcalN/(4.0*M_PI);
        #endif

        #if(COLLISION_KERNEL == NEXT_TO_LEADING_ORDER_KERNEL )
            muSqr = mDSqr/4.0*std::exp(-2.0+2.0*GammaE);
            q0 = g*g*g*g*Temp*Temp*Temp*MathcalN/(4.0*M_PI);
        #endif

        double OldVal = QsFct(qhatFct(10.0,zVal),P,zVal);
        double NewVal = QsFct(qhatFct(OldVal,zVal),P,zVal);
        double tol = 1e-8;
        // std::cerr << mDSqr <<std::endl;
        // std::cerr << muSqr <<std::endl;
        // std::cerr << OldVal <<std::endl;
        // std::cerr << COne << " " << Cz << " " << CzB <<std::endl;
        while( std::abs( OldVal/NewVal - 1.0 ) >= tol ){
            OldVal = NewVal;
            NewVal = QsFct(qhatFct(NewVal,zVal),P,zVal);
        }

        Qs = NewVal;
        qhat = qhatFct(Qs,zVal);
        std::cerr << "# Qs = " << Qs <<std::endl;

        Omega = std::complex<double>(0.5,-0.5)*std::sqrt( qhat/(P*zVal*(1.0-zVal)) );

    }

    double DeltaGamma(double b){
        // return -0.25*q0*b*b*std::log(Qs*b*b);

        #if(COLLISION_KERNEL == DEBYE_SCREENED_COULOMB_KERNEL )
            return  (g*g*Temp/(2.0*M_PI))*( GammaE + std::log( 0.5*mD * b ) + BesselK0(mD * b) ) - 0.25*q0*b*b*std::log(Qs/muSqr) ;
        #endif

        #if(COLLISION_KERNEL == LATTICE_EQCD_KERNEL )
            double gTb = g*Temp*b;
            return g*g*Temp*BroadeningKernel::NonPerturbative::GammabOvergSqrT(gTb)- 0.25*q0*b*b*std::log(Qs/muSqr) ;
        #endif

        #if(COLLISION_KERNEL == LEADING_ORDER_KERNEL )
            double gTb = g*Temp*b;
            return g*g*Temp*BroadeningKernel::LO::GammabOvergSqrT(gTb)- 0.25*q0*b*b*std::log(Qs/muSqr) ;
        #endif

        #if(COLLISION_KERNEL == NEXT_TO_LEADING_ORDER_KERNEL )
            double gTb = g*Temp*b;
            return g*g*Temp*BroadeningKernel::NLO::GammabOvergSqrT(gTb)- 0.25*q0*b*b*std::log(Qs/muSqr) ;
        #endif
    }
}

namespace NextToLeadingOrder{
    inline std::complex<double> F(std::complex<double> x){
        auto logX = std::complex<double>( std::log(std::abs(x)),std::arg(x));

        return (1.0/x)*(logX + GammaE);
    }
    inline std::complex<double> DerivF(std::complex<double> x){
        auto logX = std::complex<double>( std::log(std::abs(x)),std::arg(x));

        return (1.0/(x*x))*( 1.0 - logX - GammaE);
    }
    inline std::complex<double> kSqr(double s){

        double RealOmega = std::real(LeadingOrder::Omega);
        double ImagOmega = std::imag(LeadingOrder::Omega);

        // Cos[a + I b] = Cos[a] Cosh[b] - I*Sin[a] Sinh[b]
        // Sin[a + I b] = Cosh[b] Sin[a] + I*Cos[a] Sinh[b]
        double Coshb = std::cosh(ImagOmega*s), Sinhb = std::sinh(ImagOmega*s);
        double Cosa  = std::cos(RealOmega*s),  Sina  = std::sin(RealOmega*s);
        auto CosOmegas  = Cosa*Coshb - I*Sina*Sinhb;
        auto SinOmegas  = Sina*Coshb + I*Cosa*Sinhb;
        
        Coshb = std::cosh(ImagOmega*(t-s)); Sinhb = std::sinh(ImagOmega*(t-s));
        Cosa  = std::cos(RealOmega*(t-s));  Sina  = std::sin(RealOmega*(t-s));
        auto CosOmegatMs  = Cosa*Coshb - I*Sina*Sinhb;
        auto SinOmegatMs  = Sina*Coshb + I*Cosa*Sinhb;

        return I*P*zVal*(1.0-zVal)*0.5*LeadingOrder::Omega*( CosOmegas/SinOmegas - SinOmegatMs/CosOmegatMs );
    }


    inline std::complex<double> kSqrDeriv(double s){

        double RealOmega = std::real(LeadingOrder::Omega);
        double ImagOmega = std::imag(LeadingOrder::Omega);

        // Cos[a + I b] = Cos[a] Cosh[b] - I*Sin[a] Sinh[b]
        // Sin[a + I b] = Cosh[b] Sin[a] + I*Cos[a] Sinh[b]
        double Coshb = std::cosh(ImagOmega*(t-s)), Sinhb = std::sinh(ImagOmega*(t-s));
        double Cosa  = std::cos(RealOmega*(t-s)),  Sina  = std::sin(RealOmega*(t-s));
        auto CosOmegatMs  = Cosa*Coshb - I*Sina*Sinhb;

        return -I*P*zVal*(1.0-zVal)*0.5*LeadingOrder::Omega*LeadingOrder::Omega*( 1.0/(CosOmegatMs*CosOmegatMs) );
    }

    double Integrand(double s, void *userdata){

        auto ks = kSqr(s);

        auto k1  = ks/LeadingOrder::Qs;


        return 0.25*std::real( COne*F(-k1) + Cz*F(-k1/(zVal*zVal)) + CzB*F(-k1/((1.0-zVal)*(1.0-zVal))) );
    }

    double Spectrum(){
        
        double integral[1],error[1];
        double epsrel = 1e-3; double epsabs = 1e-2; 
        gsl_function Fct; Fct.function = Integrand;
        gsl_integration_qags(&Fct,0.0,t,epsabs,epsrel,Integration::NumberOfIntegrationPoints,Integration::wsp[omp_get_thread_num()],integral,error);
        return (g*g/(4.0*M_PI*M_PI))*LeadingOrder::q0/LeadingOrder::Qs*integral[0];
    }

    double IntegrandRate(double s, void *userdata){

        auto k1 = kSqr(s)/LeadingOrder::Qs;

        auto DerivKs = kSqrDeriv(s);

        return std::real((-DerivKs / LeadingOrder::Qs) * (DerivF(-k1) + (1.0 / (zVal * zVal)) * DerivF(-k1/(zVal*zVal)) + (1.0 / ((1.0 - zVal) * (1.0 - zVal))) * DerivF(-k1/((1.0-zVal)*(1.0-zVal)))));
    }

    double Rate(){
        
        double integral[1],error[1];
        double epsrel = 1e-3; double epsabs = 1e-2; 
        gsl_function Fct; Fct.function = IntegrandRate;
        gsl_integration_qags(&Fct,0.0,t,epsabs,epsrel,Integration::NumberOfIntegrationPoints,Integration::wsp[omp_get_thread_num()],integral,error);

        double Integrated = Integrand(t,NULL);
        return (g*g/(4.0*M_PI*M_PI))*LeadingOrder::q0/LeadingOrder::Qs*(Integrated + integral[0]);
    }

}

namespace NextToLeadingOrderIntegration{
    inline auto ComplexExp(std::complex<double> x){
        return std::exp( x.real() )*std::complex<double>( std::cos( x.imag() ) ,std::sin( x.imag() ) );
    }

    // Integrand //
    static int SpectrumIntegrand(const int *ndim, const double xx[], 
        const int *ncomp, double fval[], void *userdata){
        GSLVARIABLES &parms = *(GSLVARIABLES *) userdata;
        double z = zVal;
        double u = mD*std::tan((M_PI/2.0)*xx[0]);
        double s = t*xx[1];
        double u2 = u*u;

        double Jacobian = t*(mD*M_PI/2.0)*(( mDSqr + u2)/mDSqr);

        auto ks = NextToLeadingOrder::kSqr(s);

        auto Exp   = ComplexExp( ks*u2 );
        auto Expz  = ComplexExp( ks*u2/(z*z) );
        auto ExpzB = ComplexExp( ks*u2/((1.0-z)*(1.0-z)) );
        fval[0] = Jacobian*std::real( 2.0/u*LeadingOrder::DeltaGamma(u)*( COne*Exp + Cz*Expz + CzB*ExpzB ) );
        // fval[0] = std::real( -0.5*u*LeadingOrder::q0*std::log(LeadingOrder::Qs*u2)*( COne*Exp + Cz*Expz + CzB*ExpzB ) );
        return 1;
    }

    double Spectrum(){
        double integral[1],error[1],prob[1];
        int nregions, neval,fail;
        double epsrel = 1e-3; double epsabs = 1e-3;
        GSLVARIABLES Vars;

        Cuhre(2, 1, SpectrumIntegrand, &Vars, 1,
                    epsrel, epsabs, 0,
                    MINEVAL, MAXEVAL, 0,
                    STATEFILE, SPIN,
                    &nregions, &neval, &fail, integral, error, prob);
        return (g*g/(4.0*M_PI*M_PI))*integral[0];
    }

    // Integrand //
    static int RateIntegrand(const int *ndim, const double xx[], 
        const int *ncomp, double fval[], void *userdata){
        GSLVARIABLES &parms = *(GSLVARIABLES *) userdata;
        double z = zVal;
        double u = mD*std::tan((M_PI/2.0)*xx[0]);
        double s = t*xx[1];
        double u2 = u*u;

        double Jacobian = t*(mD*M_PI/2.0)*(( mDSqr + u2)/mDSqr);

        auto ks   = NextToLeadingOrder::kSqr(s);
        auto Deks = NextToLeadingOrder::kSqrDeriv(s);

        auto kt = NextToLeadingOrder::kSqr(t);


        auto Exp   = Deks*u2*ComplexExp( ks*u2 ) + ComplexExp( kt*u2 )/t;
        auto Expz  = Deks*u2/(z*z)*ComplexExp( ks*u2/(z*z) ) + ComplexExp( kt*u2/(z*z) )/t;
        auto ExpzB = Deks*u2/((1.0-z)*(1.0-z))*ComplexExp( ks*u2/((1.0-z)*(1.0-z)) ) + ComplexExp( kt*u2/((1.0-z)*(1.0-z)) )/t;

        fval[0] = Jacobian*std::real( 2.0/u*LeadingOrder::DeltaGamma(u)*( COne*Exp + Cz*Expz + CzB*ExpzB  ) );
        // fval[0] = std::real( -0.5*u*LeadingOrder::q0*std::log(LeadingOrder::Qs*u2)*( COne*Exp + Cz*Expz + CzB*ExpzB ) );
        return 1;
    }

    double OneDIntegrand(double x, void *p){
        double z = zVal;
        double u = mD*std::tan((M_PI/2.0)*x);
        double s = t;
        double u2 = u*u;

        double Jacobian = (mD*M_PI/2.0)*(( mDSqr + u2)/mDSqr);

        auto ks = NextToLeadingOrder::kSqr(s);

        auto Exp   = ComplexExp( ks*u2 );
        auto Expz  = ComplexExp( ks*u2/(z*z) );
        auto ExpzB = ComplexExp( ks*u2/((1.0-z)*(1.0-z)) );
        return Jacobian*std::real( 2.0/u*LeadingOrder::DeltaGamma(u)*( COne*Exp + Cz*Expz + CzB*ExpzB ) );
    }

    double Rate(){
        double integral[1],error[1],prob[1];
        int nregions, neval,fail;
        double epsrel = 1e-3; double epsabs = 1e-3;
        GSLVARIABLES Vars;

        Cuhre(2, 1, RateIntegrand, &Vars, 1,
                    epsrel, epsabs, 0,
                    MINEVAL, MAXEVAL, 0,
                    STATEFILE, SPIN,
                    &nregions, &neval, &fail, integral, error, prob);

        double Integrated = 0.0;
        gsl_function Fct; Fct.function = &OneDIntegrand;
        // gsl_integration_qags(&Fct,0.0,1.0,epsabs,epsrel,Integration::NumberOfIntegrationPoints,Integration::wsp[omp_get_thread_num()],&Integrated,error);

        return (g*g/(4.0*M_PI*M_PI))*(Integrated + integral[0]);
    }
}

void Evolve(){

    std::stringstream RateFolder;
    RateFolder << "HO/Rate-P" << P << "-z" << zVal << ".txt";
    std::ofstream RateFile; 
    RateFile.open(RateFolder.str().c_str(),std::ofstream::out);
    RateFile << "# Parameter : "  << Process << " gs = " << g << " " << "P = " << P << " " << "T = " << Temp << " " << "z = " << zVal  << " mDSqr= " << mDSqr << " Qs^2 = " << LeadingOrder::Qs << " q0= " << LeadingOrder::q0 << " q3= " << LeadingOrder::qhat << " muSqr= " << LeadingOrder::muSqr << "\n";

    int i = 0;
    dt = 1e-1; 
    tmax = 60.0; 


    int Nx = std::floor(tmax/dt)-1;
    double *tVals = new double[Nx],*SpectrumVals = new double[Nx],*SpectrumVals1 = new double[Nx],*SpectrumVals2 = new double[Nx];

    std::cerr << "gs = "<< g << "\n";
    std::cerr << "mDSqr = "<< mDSqr << "\n";
    t += dt;
    while( t <= tmax ){
        if( i % 100 ==  0){
            std::cerr << "# Done t = " << t << "\n";
        }

        tVals[i] = t;
        SpectrumVals[i]  = LeadingOrder::Spectrum();
        SpectrumVals1[i] = NextToLeadingOrder::Spectrum();
        SpectrumVals2[i] = NextToLeadingOrderIntegration::Spectrum();

        t+=dt;
        i++;
    }

    // gsl_spline *Spectrum  = gsl_spline_alloc(gsl_interp_linear,Nx);
    gsl_spline *Spectrum1 = gsl_spline_alloc(gsl_interp_linear,Nx);
    gsl_spline *Spectrum2 = gsl_spline_alloc(gsl_interp_linear,Nx);
    // gsl_spline_init(Spectrum,tVals,SpectrumVals,Nx);
    gsl_spline_init(Spectrum1,tVals,SpectrumVals1,Nx);
    gsl_spline_init(Spectrum2,tVals,SpectrumVals2,Nx);
    gsl_interp_accel *acc= gsl_interp_accel_alloc();
    // RateFile << " t Sum LODeriv NLODeriv NLODeriv2 LOSpect NLOSpect NLOSpect2 "<< std::endl;
    // RateFile << "# Parameter : " << Process << " gs = " << g << " " << "P = " << P << " " << "T = " << Temp << " " << "z = " << zVal << " mDSqr= " << mDSqr << "\n";
    RateFile << "# Columns: 1- Time: T^2 t/(P*z*(1-z))    2- 1/(g^4T P(z)) dGamma/dz \n";
    for(int i=0;i<Nx;i++){
        RateFile << tVals[i]/(2.*P*zVal*(1-zVal)) << " " 
                 << LeadingOrder::Rate(tVals[i]) + gsl_spline_eval_deriv(Spectrum2,tVals[i],acc) << " "
                //  << LeadingOrder::Rate(tVals[i]) << " "
                //  << gsl_spline_eval_deriv(Spectrum1,tVals[i],acc) << " "
                //  << gsl_spline_eval_deriv(Spectrum2,tVals[i],acc) << " "
                //  << SpectrumVals[i]  << " "
                //  << SpectrumVals1[i] << " "
                //  << SpectrumVals2[i] 
                 << std::endl;
    }
    RateFile.close();
}


int main(int argc, char **argv){
  
    Setup(argc,argv);
    #if(COLLISION_KERNEL == LATTICE_EQCD_KERNEL )
        BroadeningKernel::NonPerturbative::Setup(1);
    #endif
    #if(COLLISION_KERNEL == LEADING_ORDER_KERNEL )
        BroadeningKernel::HardIntegralSpace::SetupHardIntegralSpline();
        BroadeningKernel::LO::SetupBspaceSpline();
    #endif

    #if(COLLISION_KERNEL == NEXT_TO_LEADING_ORDER_KERNEL )
        BroadeningKernel::HardIntegralSpace::SetupHardIntegralSpline();
        BroadeningKernel::NLO::SetupBspaceSpline();
    #endif
    LeadingOrder::Setup();
    Evolve();

    Clear();
    return 0;
}
