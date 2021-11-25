#include "Parameters.cpp"

// 2D Force integrand //
static int ForceIntegrand(const int *ndim, const double xx[], 
    const int *ncomp, double fval[], void *userdata){
// int ForceIntegrand(unsigned ndim, const double *xx, void *userdata, unsigned fdim, double *fval){
    GSLVARIABLES &parms = *(GSLVARIABLES *) userdata;
    double p  = parms.p;
    double tau  = parms.t;
    double z  = parms.z;
    double xxMin  = std::atan(p/(mD*parms.qCutoff));
    double yprime = ( xxMin )*xx[0];
    double qprime = mD*std::tan(yprime);
    double q = p/qprime;
    double q2 = q*q;

    double phi = -M_PI + (2.0*M_PI)*xx[1];
    double cosphi = std::cos(phi);

    double jacobian = (mD*xxMin)*(( mDSqr + qprime*qprime)/mDSqr)*(2.0*M_PI)*q2/p;

    std::complex<double> res = std::complex<double>(0.0, 0.0);
    if (xx[0] > 1e-15 && xx[0] < (1.0 - 1e-5)){

        double Gammaq = ( COne*Gamma(q) + Cz/(z*z)*Gamma(q/z) + CzB/((1.0-z)*(1.0-z))*Gamma(q/(1.0-z)) );

        res = jacobian * q * Gammaq * ( SplinePMinusQ(p, q, cosphi, parms) ) / (4.0 * M_PI * M_PI);
    }

    fval[0] = res.real();
    fval[1] = res.imag();

    if(!isfinite(fval[0]) || !isfinite(fval[1])){
        std::cerr << fval[0] << " " << fval[1] << "\n";
        std::cerr << " q = " << q << " ";
        std::cerr << " p = " << p << "\n";
        exit(0);
    }

    return 1;
}



// 1D Force Integrand //
double ForceIntegrand1D( double xx, void *userdata){
    GSLVARIABLES &parms = *(GSLVARIABLES *) userdata;
    double p  = parms.p;
    double p2 = p*p;
    double tau  = parms.t;
    double z  = parms.z;

    double xxMin  = std::atan(p/(mD*parms.qCutoff));
    double yprime = ( xxMin )*xx;
    double qprime = mD*std::tan(yprime);
    double q = p/qprime;
    double q2 = q*q;

    double jacobian = (mD*xxMin)*(( mDSqr + qprime*qprime)/mDSqr)*(2.0*M_PI)*q2/p;

    // double yprime = (M_PI/2.0)*xx;
    // double qprime = mD*std::tan(yprime);
    // double q = p/qprime;
    // double q2 = q*q;
    // double jacobian = (M_PI)*(mD*M_PI/2.0)*( 1.0 + qprime*qprime/mDSqr)*(q2/p);

    // double q = ( 1.0 - xx)/xx;
    // double q2 = q*q;
    // double jacobian = (M_PI)/(xx*xx);

    std::complex<double> res = std::complex<double>(0.0,0.0);
    using namespace Integration;
    if( xx > 1e-15 && xx < (1.0-1e-5)  ){
        for(int k =0 ; k<NClenshawCurtis/2;k++){
            double phi = (M_PI)*Points(k)+M_PI;
            double cosphi = std::cos(phi);

            res += ClenshawCurtisWeights[k]*( SplinePMinusQ(p, q, cosphi, parms) );
        }

        double Gammaq = ( COne*Gamma(q) + Cz/(z*z)*Gamma(q/z) + CzB/((1.0-z)*(1.0-z))*Gamma(q/(1.0-z)) );
        res *= jacobian * q * Gammaq / (4.0 * M_PI * M_PI);
    }

    if(!isfinite(res.real()) || !isfinite(res.imag())){
        std::cerr << res.real() << " " << res.imag() << "\n";
        std::cerr << " q = " << q << " ";
        std::cerr << " p = " << p << "\n";
        exit(0);
    }
    if(parms.RealImag == 0){
        return res.real();
    }
    return res.imag();
}



// Small momentum exchange approximation //
namespace Diffusion{
    // Cutoff scale of momentum exchange // 
    inline double qCutFct(double p){
        return 0.05;
    }

    // Integrand Moments I_i helpers //
    double IntegrandMoment(double q, void * userdata){
        GSLVARIABLES &parms = *(GSLVARIABLES *) userdata;
        double p = parms.p;
        double tau = parms.t;
        double z = parms.z;
        double res;
        double q2 = q*q;

        double Gammaq = q*(COne*Gamma(q) + Cz/(z*z)*Gamma(q/z) + CzB/((1.0-z)*(1.0-z))*Gamma(q/(1.0-z)) );

        if(parms.Which == 0){
            if( parms.RealImag == 0){
                res = (1.0 - std::cos(-tau*q2) * gsl_sf_bessel_J0(tau * 2.0 * p*q));
            }
            else{
                res = - std::sin(-tau*q2) * gsl_sf_bessel_J0(tau * 2.0 * p*q);
            }

        }
        else if(parms.Which == 1){
            if( parms.RealImag == 0){
                res = -q2/(p*p)*std::cos(-tau*q2) * gsl_sf_bessel_Jn(2,tau * 2.0 * p*q);
            }
            else{
                res = -q2/(p*p)*std::sin(-tau*q2) * gsl_sf_bessel_Jn(2,tau * 2.0 * p*q);
            }

        }
        else if(parms.Which == 2){
            if( parms.RealImag == 0){
                res = -q/(p)*std::sin(-tau*q2) * gsl_sf_bessel_J1(tau * 2.0 * p*q);
            }
            else{
                res =  q/(p)*std::cos(-tau*q2) * gsl_sf_bessel_J1(tau * 2.0 * p*q);
            }
        }
        else if(parms.Which == 3){
            if( parms.RealImag == 0){
                if( tau * 2.0 * p*q > 1e-15 ){
                    res = q2/(p*p)*std::cos(-tau*q2) * ( 3.0*(1.0/(2.0*tau*p*q)) * gsl_sf_bessel_J1(tau * 2.0 * p*q) - 2.0*gsl_sf_bessel_J0(tau * 2.0 * p*q) );
                }
                else{
                    res = q2/(p*p)*std::cos(-tau*q2) * ( -0.5 );
                }
            }
            else{
                if( tau * 2.0 * p*q > 1e-15 ){
                    res = q2/(p*p)*std::sin(-tau*q2) * ( 3.0*(1.0/(2.0*tau*p*q)) * gsl_sf_bessel_J1(tau * 2.0 * p*q) - 2.0*gsl_sf_bessel_J0(tau * 2.0 * p*q) );
                }
                else{
                    res = q2/(p*p)*std::sin(-tau*q2) * ( -0.5 );
                }
            }

        }
        else if(parms.Which == 4){
            if( parms.RealImag == 0){
                if( tau * 2.0 * p*q > 1e-15 ){
                    res = q2/(p*p)*std::cos(-tau*q2) * ( (1.0/(2.0*tau*p*q)) * gsl_sf_bessel_J1(tau * 2.0 * p*q) - gsl_sf_bessel_Jn(2,tau * 2.0 * p*q) );
                }
                else{
                    res = q2/(p*p)*std::cos(-tau*q2) * ( 0.5 );
                }
            }
            else{
                if( tau * 2.0 * p*q > 1e-15 ){
                    res = q2/(p*p)*std::sin(-tau*q2) * ( (1.0/(2.0*tau*p*q)) * gsl_sf_bessel_J1(tau * 2.0 * p*q) - gsl_sf_bessel_Jn(2,tau * 2.0 * p*q) );
                }
                else{
                    res = q2/(p*p)*std::sin(-tau*q2) * ( 0.5 );
                }
            }

        }
        if(!isfinite(res)){
            std::cerr << "IntegranMomentI is not finite \n" ;
            exit(0);
        }

        return Gammaq*res/(2.0*M_PI);
    }

    // Integrating The moments //
    std::complex<double> IntegralI(double p, double t, double z, int Which){
        double qCutoff = qCutFct(p);
        using namespace Integration;
        double res1 =0.0 ,res2=0.0,err;
        double epsrel = 1e-3; double epsabs = 1e-3;
        GSLVARIABLES Vars; Vars.p = p; Vars.t = t; Vars.z = z;
        gsl_function Fct; Fct.function = &IntegrandMoment; Fct.params = &Vars;
        Vars.Which = Which;
        if(qCutoff > 1e-15){
            Vars.RealImag = 0;
            gsl_integration_qags(&Fct,0.0,qCutoff,epsabs,epsrel,NumberOfIntegrationPoints,wsp[omp_get_thread_num()],&res1,&err);

            Vars.RealImag = 1;
            gsl_integration_qags(&Fct,0.0,qCutoff,epsabs,epsrel,NumberOfIntegrationPoints,wsp[omp_get_thread_num()],&res2,&err);
        }
        return std::complex<double>(res1,res2);

    }
    

    // Forces //
    #if( SPLITSMALLQ == 0 )
            std::complex<double> DiffusionForce(gsl_spline* SplineReal, gsl_spline* SplineImag, double p, double z, double t, int k){
                return std::complex<double>(0.0,0.0);
            }
    #endif
    #if( SPLITSMALLQ == 1 )
        std::complex<double> DiffusionForce(gsl_spline* SplineReal, gsl_spline* SplineImag, double p, double z, double t, int k){
            std::complex<double> I10 = IntegralI(p,t,z,0);
            std::complex<double> I13 = IntegralI(p,t,z,1);
            std::complex<double> I2  = IntegralI(p,t,z,2);
            std::complex<double> I23 = IntegralI(p,t,z,3);
            std::complex<double> I3  = IntegralI(p,t,z,4);
            std::complex<double> P1  = I10 - I2 - I13;
            std::complex<double> P2  = p*( 2.0*I2 - I23 )*0.5;
            std::complex<double> P3  = -p*p*0.5*I3;

            std::complex<double> Psi       = std::complex<double>(0.0,0.0);
            std::complex<double> PsiDeriv  = std::complex<double>(0.0,0.0);
            std::complex<double> PsiDeriv2 = std::complex<double>(0.0,0.0);

            if( k >= 1 && k <= Nx-2 ) { 
                int tID   = omp_get_thread_num();
                Psi       = std::complex<double>( SplineReal->y[k],SplineImag->y[k]);

                int u = 3;
                if( k < u ){
                    // double dx1 = xVals[u] - xVals[u-1];
                    // double dx2 = xVals[u+1] - xVals[u];
                    // double dx = xVals[u+1] - xVals[u-1];
                    // // PsiDeriv  = std::complex<double>( (SplineReal->y[u]-SplineReal->y[u-1])/(dx1) ,(SplineImag->y[u]-SplineImag->y[u-1])/dx1);
                    // PsiDeriv  = (xVals[k]/xVals[u])*std::complex<double>( (SplineReal->y[u+1]-SplineReal->y[u-1])/(dx) ,(SplineImag->y[u+1]-SplineImag->y[u-1])/dx);
                    
                    // PsiDeriv2 = std::complex<double>( 
                    //     ( (SplineReal->y[u+1] - SplineReal->y[u])/dx2 - (SplineReal->y[u] - SplineReal->y[u-1])/dx1 ) /dx ,
                    //     ( (SplineImag->y[u+1] - SplineImag->y[u])/dx2 - (SplineImag->y[u] - SplineImag->y[u-1])/dx1 ) /dx );
                }
                else{

                    double dx1 = xVals[k] - xVals[k-1];
                    double dx2 = xVals[k+1] - xVals[k];
                    double dx = xVals[k+1] - xVals[k-1];
                    // PsiDeriv  = std::complex<double>( (SplineReal->y[k]-SplineReal->y[k-1])/(dx1) ,(SplineImag->y[k]-SplineImag->y[k-1])/dx1);
                    PsiDeriv  = std::complex<double>( (SplineReal->y[k+1]-SplineReal->y[k-1])/(dx) ,(SplineImag->y[k+1]-SplineImag->y[k-1])/dx);
                    
                    PsiDeriv2 = 2.0*std::complex<double>( 
                        ( (SplineReal->y[k+1] - SplineReal->y[k])/dx2 - (SplineReal->y[k] - SplineReal->y[k-1])/dx1 ) /dx ,
                        ( (SplineImag->y[k+1] - SplineImag->y[k])/dx2 - (SplineImag->y[k] - SplineImag->y[k-1])/dx1 ) /dx );
                }

                // Psi       = std::complex<double>( gsl_spline_eval(SplineReal,p,Acc[tID]),gsl_spline_eval(SplineImag,p,Acc1[tID]));
                // int u = 5;
                // if( p <= xVals[u] ){
                //     PsiDeriv = (p/xVals[u])*std::complex<double>( SplineReal->y[u],SplineImag->y[u]);
                // }
                // else{
                //     PsiDeriv = std::complex<double>( gsl_spline_eval_deriv(SplineReal,p,Acc[tID]),gsl_spline_eval_deriv(SplineImag,p,Acc1[tID]));
                // }

                // if( p <= xVals[u] ){
                //     PsiDeriv2 = std::complex<double>( SplineReal->y[u],SplineImag->y[u]);
                // }
                // else{
                //     PsiDeriv2 = std::complex<double>( gsl_spline_eval_deriv2(SplineReal,p,Acc[tID]),gsl_spline_eval_deriv2(SplineImag,p,Acc1[tID]));
                // }
                
                
            
            
            }

            return (2.0*P*z*(1.0-z))*(P1*Psi + P2*PsiDeriv + P3*PsiDeriv2 );
        }
    #endif
    

    #if( SPLITSMALLQ == 0 )
        void DiffusionCoeffs(double p, double z, double t, std::complex<double> &P1, std::complex<double> &P2, std::complex<double> &P3){
            P1  = std::complex<double>(0.0,0.0);
            P2  = std::complex<double>(0.0,0.0);
            P3  = std::complex<double>(0.0,0.0);
        }
    #endif
    #if( SPLITSMALLQ == 1 )
        void DiffusionCoeffs(double p, double z, double t, std::complex<double> &P1, std::complex<double> &P2, std::complex<double> &P3){
            std::complex<double> I10 = IntegralI(p,t,z,0);
            std::complex<double> I13 = IntegralI(p,t,z,1);
            std::complex<double> I2  = IntegralI(p,t,z,2);
            std::complex<double> I23 = IntegralI(p,t,z,3);
            std::complex<double> I3  = IntegralI(p,t,z,4);
            P1  = (2.0*P*z*(1.0-z))*(I10 - I2 - I13);
            P2  = (2.0*P*z*(1.0-z))*(p*0.5*(2.0*I2 - I23));
            P3  = (2.0*P*z*(1.0-z))*(-p*p*0.5*I3);
        }
    #endif
    
    // OUTPUT the Moments for checks //
    void Check(){

        for(int i=0; i<Nx;i++){
            std::cout << xVals[i] << " ";
            double p =  xVals[i];
            std::complex<double> I10 = IntegralI(p,0.01,zVal,0);
            std::complex<double> I13 = IntegralI(p,0.01,zVal,1);
            std::complex<double> I2  = IntegralI(p,0.01,zVal,2);
            std::complex<double> I23 = IntegralI(p,0.01,zVal,3);
            std::complex<double> I3  = IntegralI(p,0.01,zVal,4);
            std::complex<double> P1  = I10 - I2 - I13;
            std::complex<double> P2  = p*( 2.0*I2 - I23 )*0.5;
            std::complex<double> P3  = -p*p*0.5*I3;

            std::cout << P1.real() << " ";
            std::cout << P1.imag() << " ";
            std::cout << P2.real() << " ";
            std::cout << P2.imag() << " ";
            std::cout << P3.real() << " ";
            std::cout << P3.imag() << " ";
            std::cout << "\n";
        }
        exit(0);
    }
}

// Large momentum exchange force term //
std::complex<double> ForceTerm(gsl_spline* SplineReal, gsl_spline* SplineImag, double p, double z, double t, int k){
    GSLVARIABLES Vars; Vars.p = p; Vars.t = t; Vars.z = z;
    Vars.SplineReal = SplineReal;
    Vars.SplineImag = SplineImag;

    double integral[2],error[2],prob[2];
    int nregions, neval,fail;
    // double MidPoint = (xVals[Nx-1]+xVals[0])/2.0;
    // double DistanceToMidPoint = std::fabs(p - MidPoint)/(xVals[Nx-1]+xVals[0])*2.0;
    // double epsrel = ( 1e-3 + (1e-8 - 1e-3)*DistanceToMidPoint ); double epsabs = 1e-3 + (1e-8 - 1e-3)*DistanceToMidPoint;
    double epsrel = 1e-3; double epsabs = 1e-3;
    if( p > xVals[0] && p < xVals[Nx-1] ){

        #if( SPLITSMALLQ == 0 )
            double qCutoff = 0.0;
        #endif
        #if( SPLITSMALLQ == 1 )
            double qCutoff = Diffusion::qCutFct(p);
        #endif
        Vars.qCutoff = qCutoff;

        Cuhre(2, 2, ForceIntegrand, &Vars, 1,
                    epsrel, epsabs, 0,
                    MINEVAL, MAXEVAL, 0,
                    STATEFILE, SPIN,
                    &nregions, &neval, &fail, integral, error, prob);

        // Vars.p = xMin+1.0;
        // Vars.RealImag = 0;
        // integral[0]=ForceIntegrand1D((p-xMin)/(xMax-xMin),&Vars);
        // Vars.RealImag = 1;
        // integral[1]=ForceIntegrand1D((p-xMin)/(xMax-xMin),&Vars);
        
        // using namespace Integration;
        // gsl_set_error_handler_off();
        // double res,err;
        // gsl_function Fct; Fct.function = &ForceIntegrand1D; Fct.params = &Vars;
        // Vars.RealImag = 0;
        // gsl_integration_qags(&Fct,0.0,1.0,epsabs,epsrel,NumberOfIntegrationPoints,wsp[omp_get_thread_num()],&res,&err);
        // integral[0] = res;
        // //
        // Vars.RealImag = 1;
        // gsl_integration_qags(&Fct,0.0,1.0,epsabs,epsrel,NumberOfIntegrationPoints,wsp[omp_get_thread_num()],&res,&err);
        // integral[1] = res;

        // typedef typename std::complex<double>::value_type value_type;
        // auto f  = [&](double xx) { return ForceIntegrand1D1(xx,&Vars); };
        // boost::math::quadrature::sinh_sinh<value_type> integrator;
        // boost::math::quadrature::exp_sinh<value_type> integratorExp;
        // Vars.Which = 0;
        // std::complex<double> res1 = integratorExp.integrate(f,1.0,std::numeric_limits<double>::infinity());
        // integral[0] = res1.real();
        // integral[1] = res1.imag();
        // Vars.Which = 1;
        // res1 = integratorExp.integrate(f);
        // // integrator.integrate(f);
        // integral[0] += res1.real();
        // integral[1] += res1.imag();

    }
    else{
        integral[0] = 0.0; integral[1] = 0.0;
    }
    
    return  (2.0*P*z*(1.0-z))*(std::complex<double>(integral[0],integral[1]) );
}

// Updating the different Splines //
void UpdateSpline(){
    gsl_spline_init(PsiSplineReal.get(),xVals.data(),PsiInitRealVals.data(),Nx);
    gsl_spline_init(PsiSplineImag.get(),xVals.data(),PsiInitImagVals.data(),Nx);
    gsl_spline_init(IntegrandSumSpline.get(),xVals.data(),IntegrandSum.data(),Nx);
}

// Output the wave function //
void OUTPUTWave(){
    std::stringstream fname;
    fname <<"OUTPUT/WaveFct-" << t << ".txt";

    std::ofstream Wavefct; Wavefct.open(fname.str().c_str(),std::ofstream::out);
    Wavefct << "x PsiReal PsiImag ForceReal ForceImag DifReal DifImag IntegrandSum SplineReal SplineImag  SplineDervReal SplineDervImag  SplineDeriv2Real SplineDeriv2Imag \n";
    std::vector<std::complex<double>> Forc(Nx);
    std::vector<std::complex<double>> Forc1(Nx);

    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<Nx;i++){
        // std::cerr << "tID = "<< omp_get_thread_num() << " i= " << i <<"\n" ;
        Forc[i]  = ForceTerm(PsiSplineReal.get(),PsiSplineImag.get(),xVals[i],zVal,t,i);
        Forc1[i] = Diffusion::DiffusionForce(PsiSplineReal.get(),PsiSplineImag.get(),xVals[i],zVal,t,i);
    }

    for(int i=0;i<Nx;i++){
        Wavefct << xVals[i] << " ";
        Wavefct << PsiInitRealVals[i] << " ";
        Wavefct << PsiInitImagVals[i] << " ";
        Wavefct << std::real(Forc[i]) << " ";
        Wavefct << std::imag(Forc[i]) << " ";
        Wavefct << std::real(Forc1[i]) << " ";
        Wavefct << std::imag(Forc1[i]) << " ";
        Wavefct << IntegrandSum[i] << " ";
        Wavefct << gsl_spline_eval(PsiSplineReal.get(),xVals[i],Acc[omp_get_thread_num()]) << " ";
        Wavefct << gsl_spline_eval(PsiSplineImag.get(),xVals[i],Acc[omp_get_thread_num()]) << " ";
        Wavefct << gsl_spline_eval_deriv(PsiSplineReal.get(),xVals[i],Acc[omp_get_thread_num()]) << " ";
        Wavefct << gsl_spline_eval_deriv(PsiSplineImag.get(),xVals[i],Acc[omp_get_thread_num()]) << " ";
        Wavefct << gsl_spline_eval_deriv2(PsiSplineReal.get(),xVals[i],Acc[omp_get_thread_num()]) << " ";
        Wavefct << gsl_spline_eval_deriv2(PsiSplineImag.get(),xVals[i],Acc[omp_get_thread_num()]) << " ";
        Wavefct << "\n";
    }
    Wavefct.close();  
}

// Integrating to find the rate //
namespace RateIntegrand
{
    const int Nt = 21;
    gsl_interp_accel **XA, **TA;
    std::vector<double> tVals(Nt),PsiRealVals(Nt*Nx),PsiImagVals(Nt*Nx);
    std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> PsiReal(gsl_spline2d_alloc(gsl_interp2d_bicubic,Nx,Nt),&gsl_spline2d_free);
    std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> PsiImag(gsl_spline2d_alloc(gsl_interp2d_bicubic,Nx,Nt),&gsl_spline2d_free);
    // gsl_spline2d * PsiReal, * PsiImag;
    double tValsmin;
    double tValsmax;
    double CumulativeIntegral;

    inline void UpdateValue(int i, int tStep, double Real, double Imag){
        int j = tStep % (Nt - 1) ;
        // gsl_spline2d_set(PsiReal,PsiRealVals,i,j,Real);
        // gsl_spline2d_set(PsiImag,PsiImagVals,i,j,Imag);
        PsiRealVals[ j*Nx + i ] = Real;
        PsiImagVals[ j*Nx + i ] = Imag;
    }
    
    inline void UpdateTime(int tStep, double t){
        int i = tStep % (Nt-1) ;
        tVals[i] = t;
    }

    void Setup(double tmax_){
        tValsmin = tValsmax;
        tVals[0] = tValsmin;
        tValsmax = tmax_;
        tVals[Nt-1] = tValsmax;
        // gsl_spline2d_free(PsiReal.get());
        // gsl_spline2d_free(PsiImag.get());
        // PsiReal = gsl_spline2d_alloc(gsl_interp2d_bicubic,Nx,Nt);
        // PsiImag = gsl_spline2d_alloc(gsl_interp2d_bicubic,Nx,Nt);
        gsl_spline2d_init(PsiReal.get(),xVals.data(),tVals.data(),PsiRealVals.data(),Nx,Nt);
        gsl_spline2d_init(PsiImag.get(),xVals.data(),tVals.data(),PsiImagVals.data(),Nx,Nt);
        
    }

    inline std::complex<double> PsiEval(double x, double t){
        if(x >= xMin && x<= xMax && t >= tValsmin && t<= tValsmax){
            int tID = omp_get_thread_num();
            return std::complex<double>(gsl_spline2d_eval(PsiReal.get(),x,t,XA[tID],TA[tID]),gsl_spline2d_eval(PsiImag.get(),x,t,XA[tID],TA[tID]));
        }
        return std::complex<double>(0.0,0.0);
    }

    static int Integrand(const int *ndim, const double xx[], 
    const int *ncomp, double fval[], void *userdata){
        double t = tValsmin + (tValsmax - tValsmin)*xx[0];
        double x = xMin + (xMax - xMin)*xx[1];
        
        double jacobian = (xMax - xMin)*(tValsmax - tValsmin);

        double DeltaE = Epsilon(x*x,zVal);
        // std::cos(-deltaE(xVals[i]*xVals[i],zVal)*t)*PsiInitRealVals[i]-std::sin(-deltaE(xVals[i]*xVals[i],zVal)*t)*PsiInitImagVals[i]

        int tID = omp_get_thread_num();
        fval[0] = x*jacobian*( std::cos(-DeltaE*t)*gsl_spline2d_eval(PsiReal.get(),x,t,XA[tID],TA[tID]) - std::sin(-DeltaE*t)*gsl_spline2d_eval(PsiImag.get(),x,t,XA[tID],TA[tID]) )/(2.0*M_PI);
        return 0;
    }

    double Integral(){
        GSLVARIABLES Vars;
        double integral[1],error[1],prob[1];
        int nregions, neval,fail;
        double epsrel = 1e-5; double epsabs = 1e-5;

        Cuhre(2, 1, Integrand, &Vars, 1,
                epsrel, epsabs, 0,
                MINEVAL, MAXEVAL, 0,
                STATEFILE, SPIN,
                &nregions, &neval, &fail, integral, error, prob);

        CumulativeIntegral += (g*g/(2.0*M_PI))*integral[0];
        
        return CumulativeIntegral;
    }

    void Init(double tmin_, double tmax_){
        tValsmin = tmin_;
        tValsmax = tmax_;

        tVals[0] = tValsmin;

        CumulativeIntegral = 0.0;
        XA = new gsl_interp_accel*[NumberOfOpenMPThreads];
        TA = new gsl_interp_accel*[NumberOfOpenMPThreads];
        #pragma omp parallel 
        {
            int tID = omp_get_thread_num();
            XA[tID] = gsl_interp_accel_alloc();
            TA[tID] = gsl_interp_accel_alloc();
        }
    }
} // namespace name


// Output the rate //
void OUTPUTRate(){
    std::stringstream RateFolder;
    RateFolder << "OUTPUT/Rate-P" << P << "-z" << zVal << ".txt";

    double Result = 0.0;
    std::ofstream RateFile; 
    if( t == tmin ){
        RateFile.open(RateFolder.str().c_str(),std::ofstream::out);
        RateFile << "# Parameter : " << Process << " gs = " << g << " " << "P = " << P << " " << "T = " << Temp << " " << "z = " << zVal << " mDSqr= " << mDSqr << "\n";
    }
    else{
        Result = RateIntegrand::Integral();
        RateFile.open(RateFolder.str().c_str(),std::ofstream::app);
    }


    RateFile << t << " " << Result << "\n";
    RateFile.close();
}

// Do the Time evolution //
namespace TimeEvolution{
    std::vector<double> k1Real(Nx), k1Imag(Nx), k2Real(Nx), k2Imag(Nx), k3Real(Nx), k3Imag(Nx), k4Real(Nx), k4Imag(Nx);
    gsl_spline *k1SplineReal, *k1SplineImag, *k2SplineReal, *k2SplineImag, *k3SplineReal, *k3SplineImag, *k4SplineReal, *k4SplineImag;

    void Setup(){
        k1SplineReal = gsl_spline_alloc(InterpType,Nx);
        k2SplineReal = gsl_spline_alloc(InterpType,Nx);
        k3SplineReal = gsl_spline_alloc(InterpType,Nx);
        k4SplineReal = gsl_spline_alloc(InterpType,Nx);
        k1SplineImag = gsl_spline_alloc(InterpType,Nx);
        k2SplineImag = gsl_spline_alloc(InterpType,Nx);
        k3SplineImag = gsl_spline_alloc(InterpType,Nx);
        k4SplineImag = gsl_spline_alloc(InterpType,Nx);

    }

    void UpdateSpline( gsl_spline * Spline , std::vector<double> Values){
        gsl_spline_init(Spline,xVals.data(),Values.data(),Nx);
    }

    void Evolve(){
        int tStep = 1;

        dt = 5e-5; 
        tmax = 1.0;
        int Nt = 2048;
        // t = dt;
        while( t <= tmax ){
            // dt = 5e-5*std::exp( tStep/double(Nt)*std::log(tmax/dt) );

            #if( TIMESTEP == RUNGE )
                // Doing k1 integrand
                #pragma omp parallel for schedule(dynamic)
                for(int i=0;i<Nx;i++){
                    std::complex<double> Forc = ForceTerm(PsiSplineReal.get(),PsiSplineImag.get(),xVals[i],zVal,t,i);
                    double k1r = -Forc.real();
                    double k1i = -Forc.imag();
                    k1Real[i] = PsiInitRealVals[i] + dt*k1r/2.0;
                    k1Imag[i] = PsiInitImagVals[i] + dt*k1i/2.0;
                }
                // std::cerr << "# Done k1 ";
                
                UpdateSpline(k1SplineReal,k1Real);
                UpdateSpline(k1SplineImag,k1Imag);

                // Doing k2 integrand
                #pragma omp parallel for schedule(dynamic)
                for(int i=0;i<Nx;i++){
                    std::complex<double> Forc = ForceTerm(k1SplineReal,k1SplineImag,xVals[i],zVal,t+dt/2.0,i);
                    double k2r = -Forc.real();
                    double k2i = -Forc.imag();
                    k2Real[i] = PsiInitRealVals[i] + dt*k2r/2.0;
                    k2Imag[i] = PsiInitImagVals[i] + dt*k2i/2.0;
                }
                // std::cerr << " k2 ";

                UpdateSpline(k2SplineReal,k2Real);
                UpdateSpline(k2SplineImag,k2Imag);

                // Doing k3 integrand
                #pragma omp parallel for schedule(dynamic)
                for(int i=0;i<Nx;i++){
                    std::complex<double> Forc = ForceTerm(k2SplineReal,k2SplineImag,xVals[i],zVal,t+dt/2.0,i);
                    double k3r = -Forc.real();
                    double k3i = -Forc.imag();
                    k3Real[i] = PsiInitRealVals[i] + dt*k3r;
                    k3Imag[i] = PsiInitImagVals[i] + dt*k3i;
                }
                // std::cerr << " k3 ";

                UpdateSpline(k3SplineReal,k3Real);
                UpdateSpline(k3SplineImag,k3Imag);

                // Doing k4 integrand
                #pragma omp parallel for schedule(dynamic)
                for(int i=0;i<Nx;i++){
                    std::complex<double> Forc = ForceTerm(k3SplineReal,k3SplineImag,xVals[i],zVal,t+dt,i);
                    double k1r = (2.0)*(k1Real[i] - PsiInitRealVals[i] );
                    double k1i = (2.0)*(k1Imag[i] - PsiInitImagVals[i] );
                    double k2r = (2.0)*(k2Real[i] - PsiInitRealVals[i] );
                    double k2i = (2.0)*(k2Imag[i] - PsiInitImagVals[i] );
                    double k3r =       (k3Real[i] - PsiInitRealVals[i] );
                    double k3i =       (k3Imag[i] - PsiInitImagVals[i] );
                    double k4r = -dt*Forc.real();
                    double k4i = -dt*Forc.imag();
                    PsiInitRealVals[i] += (1.0/6.0)*( k1r + 2.0*k2r + 2.0*k3r + k4r );
                    PsiInitImagVals[i] += (1.0/6.0)*( k1i + 2.0*k2i + 2.0*k3i + k4i );
                    RateIntegrand::UpdateValue(i,tStep,PsiInitRealVals[i],PsiInitImagVals[i]);

                    // IntegrandSum[i]    += std::cos(-deltaE(xVals[i]*xVals[i],zVal)*t)*dt*PsiInitRealVals[i]-std::sin(-deltaE(xVals[i]*xVals[i],zVal)*t)*dt*PsiInitImagVals[i];
                }
                // std::cerr << " k4 ";
            #endif

            #if( TIMESTEP == EULER )
                #pragma omp parallel for schedule(dynamic)
                for(int i=0;i<Nx;i++){
                    std::complex<double> Forc = ForceTerm(PsiSplineReal.get(),PsiSplineImag.get(),xVals[i],zVal,t,i);
                    Forc += Diffusion::DiffusionForce(PsiSplineReal.get(),PsiSplineImag.get(),xVals[i],zVal,t,i);
                    PsiInitRealVals[i] += -dt*Forc.real();
                    PsiInitImagVals[i] += -dt*Forc.imag();
                    RateIntegrand::UpdateValue(i,tStep,PsiInitRealVals[i],PsiInitImagVals[i]);
                }
            #endif

            #if( TIMESTEP == IMEX )

                std::vector<std::complex<double>> ai(Nx),bi(Nx),ci(Nx),di(Nx),cip(Nx),dip(Nx);
                cip[0] =std::complex<double> (0.0,0.0);
                dip[0] =std::complex<double> (0.0,0.0);
                cip[Nx-1] =std::complex<double> (0.0,0.0);
                dip[Nx-1] =std::complex<double> (0.0,0.0);
                #pragma omp parallel for schedule(dynamic)
                for(int i=1;i <Nx-1;i++){
                    double dx  = xVals[i+1] - xVals[i-1];
                    double dx1 = xVals[i+1] - xVals[i];
                    double dx2 = xVals[i]   - xVals[i-1];
                    std::complex<double> P1,P2,P3;
                    Diffusion::DiffusionCoeffs(xVals[i],zVal,t+dt,P1,P2,P3);

                    ai[i] = dt*( -P2/dx + P3/(dx*dx2)) ;
                    bi[i] = 1.0 + dt*( P1 - P3*(1.0/dx1 + 1.0/dx2)/dx );
                    ci[i] = dt*( P2/dx + P3/(dx*dx1) );
                    di[i] = std::complex<double>(PsiInitRealVals[i],PsiInitImagVals[i]) - dt*ForceTerm(PsiSplineReal.get(),PsiSplineImag.get(),xVals[i],zVal,t,i);
                }

                for(int i=1;i <Nx-1;i++){
                    cip[i] = ci[i] / (bi[i] - ai[i] * cip[i - 1]);
                    dip[i] = (di[i] - ai[i] * dip[i - 1]) / (bi[i] - ai[i] * cip[i - 1]);
                }

                for(int i=Nx-2;i > 0;i--){

                    std::complex<double> Psi = dip[i] - cip[i]*std::complex<double>(PsiInitRealVals[i+1],PsiInitImagVals[i+1]);

                    PsiInitRealVals[i] = Psi.real();
                    PsiInitImagVals[i] = Psi.imag();
                    RateIntegrand::UpdateValue(i,tStep,PsiInitRealVals[i],PsiInitImagVals[i]);
                }
            #endif
            // Update Splines //
            gsl_spline_init(PsiSplineReal.get(),xVals.data(),PsiInitRealVals.data(),Nx);
            gsl_spline_init(PsiSplineImag.get(),xVals.data(),PsiInitImagVals.data(),Nx);
            // gsl_spline_init(IntegrandSumSpline.get(),xVals.data(),IntegrandSum.data(),Nx);
            
            t += dt;
            RateIntegrand::UpdateTime(tStep,t);
            if( tStep % (RateIntegrand::Nt-1) == 0){
                RateIntegrand::Setup(t);
                OUTPUTRate();
                std::cerr << "# Done t = " << t << "\n";
            }
            // if( tStep % 200 ==  0 ){
            //     OUTPUTWave();
            // }
            tStep++;
        }

    }
}

int main(int argc, char **argv){
    Setup(argc,argv);

    RateIntegrand::Init(0.0,0.0);
    // xMin = 0.0;
    for(int i=0;i<Nx;i++){
        int j = Nx-1-i;
        xVals[i] = xMin*std::exp(std::log(xMax/xMin)*i/(Nx-1));
        // xVals[i] = xMin + (xMax-xMin)*(std::cos(M_PI*j/(Nx-1)) + 1.0)/2.0;
        std::complex<double> Psi = PsiInit(xVals[i],zVal);
        PsiInitRealVals[i] = Psi.real();
        PsiInitImagVals[i] = Psi.imag();
        RateIntegrand::UpdateValue(i,0,PsiInitRealVals[i],PsiInitImagVals[i]);
        IntegrandSum[i]    = 0.0;
    }
    xVals[Nx-1] = xMax;
    
    gsl_spline_init(PsiSplineReal.get(),xVals.data(),PsiInitRealVals.data(),Nx);
    gsl_spline_init(PsiSplineImag.get(),xVals.data(),PsiInitImagVals.data(),Nx);
    gsl_spline_init(IntegrandSumSpline.get(),xVals.data(),IntegrandSum.data(),Nx);
    // Diffusion::Check();

    std::cerr << "# Initialization Done \n";
    OUTPUTRate();
    // OUTPUTWave();
    // std::cerr << "# OUTPUT Initial Condition \n";

    TimeEvolution::Setup();
    TimeEvolution::Evolve();
    

    Clear();
    return 0;
}
