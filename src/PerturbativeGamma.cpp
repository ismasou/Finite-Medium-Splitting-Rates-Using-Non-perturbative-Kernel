

///////////////////////////
// Leading-Order Kernel  //
///////////////////////////


class SplineSetup{
    public :

    double LowerCutOff,UpperCutOff;
    std::vector<double> xVals, FctVals;
    std::function<double(double)> FctToInterpolate, UVLimitFct, IRLimitFct;
    int SplineSize;

    gsl_spline * Spline;
    gsl_interp_accel ** Acc;

    void SetUpValuesLogGrid(){
        xVals.resize(SplineSize);
        FctVals.resize(SplineSize);
        // Computing the spline
        xVals[0] = LowerCutOff;
        FctVals[0] = FctToInterpolate(xVals[0]);
        LowerCutOff += 1e-10;
        #pragma omp parallel for schedule(dynamic)
        for(int i=1;i<SplineSize;i++){
            // xVals[i] = LowerCutOff + (UpperCutOff-LowerCutOff)*i/double(SplineSize-1);
            xVals[i] = LowerCutOff*std::exp(std::log(UpperCutOff/LowerCutOff)*i/double(SplineSize-1));
            FctVals[i] = FctToInterpolate(xVals[i]);
        }
        std::cerr << "\n";
    }

    void SetUpValues(){
        xVals.resize(SplineSize);
        FctVals.resize(SplineSize);
        // Computing the spline
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<SplineSize;i++){
            xVals[i] = LowerCutOff + (UpperCutOff-LowerCutOff)*i/double(SplineSize-1);
            FctVals[i] = FctToInterpolate(xVals[i]);
        }
        std::cerr << "\n";
    }

    inline double GetFct(double x){

        if(x<=LowerCutOff){
            return IRLimitFct(x);
        }
        else if(x>=UpperCutOff){
            return UVLimitFct(x);
        }
        else{
            return gsl_spline_eval(Spline,x,Acc[omp_get_thread_num()]);
        }
    }

    inline double EvalDerive(double x){
            return gsl_spline_eval_deriv(Spline,x,Acc[omp_get_thread_num()]);
    }

    SplineSetup(){

    }

    void Construct(std::function<double(double)> FctToInterpolate_, double LowerCutOff_, double UpperCutOff_, std::function<double(double)> IRLimitFct_, std::function<double(double)> UVLimitFct_, int SplineSize_ =512 ,const gsl_interp_type * InterpType = gsl_interp_linear, const int LogGrid = 1){
        this->FctToInterpolate = FctToInterpolate_;
        this->LowerCutOff = LowerCutOff_;
        this->UpperCutOff = UpperCutOff_;
        this->IRLimitFct = IRLimitFct_;
        this->UVLimitFct = UVLimitFct_;
        this->SplineSize = SplineSize_;
        
        if(LogGrid == 1){
            SetUpValuesLogGrid();
        }
        else{
            SetUpValues();
        }

        Spline = gsl_spline_alloc(InterpType,SplineSize);
        Acc = new gsl_interp_accel*[NumberOfOpenMPThreads];
        # pragma omp parallel
        {  
            Acc[omp_get_thread_num()] = gsl_interp_accel_alloc();
        }

        gsl_spline_init(Spline,xVals.data(),FctVals.data(),SplineSize);

    }

    ~SplineSetup(){
        delete Acc;
    }

};


namespace HardIntegralSpace{
    SplineSetup HardIntegralInterpolation;
    int SplineSize = 2048;
    gsl_spline * HardIntegralSpline;
    gsl_interp_accel ** Acc;
    double LowerCutoff = 1e-3;
    double UpperCutoff = 100.0;

    // Computing 'I' Helper function for LO kernel //
    inline double I(double qOverT, int n, int m){


        double Bes  = BesselK2(qOverT*std::sqrt(m*n));
        double Coef = (m*n)/(2.0*std::pow((m+n),3.0));

        return Coef*(qOverT*qOverT)*Bes; 
    }
    
    // HardIntegral Summed //
    double HardIntegralSummed( double qOvergT){

        double Res = 0.0; double OldRes = 0.0;
        double Tol = 1e-5;
        int SumLimit = 2048;

        

        for(int n=1;n<SumLimit;n++){

            for(int m =1; m<SumLimit;m++){
                double sign = 1.0;
                if( (m+n-1) % 2  ){
                    sign = -1.0;
                }

                Res += (2.0*C_A+sign*4.0*Nf*Tf)*I(g*qOvergT,n,m);
                
            }
            if( std::abs(1.0 - Res/OldRes) < Tol && n > 512 ){
                goto OutOfHere;
            }
            OldRes = Res;
        }
        // Error:
            std::cerr << "Finshed the sum without reaching tolerance on LO Gammaq \n";
            std::cerr << "qOvergT= " << qOvergT << "Gammaq= " << Res << std::endl;
            exit(0);
        OutOfHere:
            Res += (2.0*C_A+3.0*Nf*Tf)*ZETA_3;
            return Res/(M_PI*M_PI);
    }
    
    // Fermi Dirac distribution//
    inline double nF(double pOverT){
        return 1.0/(std::exp(pOverT)+1.0);
    }

    // Bose Einstein distribution//
    inline double nB(double pOverT){
        return 1.0/(std::exp(pOverT)-1.0);
    }

    // Hard Part Integrand //
    static int HardIntegrand(const int *ndim, const double xx[], 
        const int *ncomp, double fval[], void *userdata){

    // int HardIntegrand(unsigned int dim, const double xx[],
    //         void *userdata, unsigned int fdim, double fval[]){
        double qOvergT = *((double *) userdata);
        
        // Change Of Variable //
        double pperp  = std::tan(M_PI*0.5*xx[0]);
        double pz     = std::tan(M_PI*(xx[1]-0.5));
        double Phi    = xx[2]*(2*M_PI);

        // Jacobian //
        double Jacobian  = (M_PI/(1.0+std::cos(M_PI*xx[0]))); // pperp jacobian 
               Jacobian *= (2.0*M_PI/(1.0+std::cos(M_PI*(2.0*xx[1]-1.0)))); // pz jacobian
               Jacobian *= (2.0*M_PI); // Phi jacobian
               
        double pOverT = std::sqrt(pperp*pperp + pz*pz);

        double primeOverT = pOverT + g * (g* qOvergT*qOvergT + 2.0*qOvergT*pperp*std::cos(Phi))/(2.0*(pOverT-pz));

        double Bose = 2.0*C_A*nB(pOverT) * ( 1.0+nB(primeOverT) );
        double Fermi = 4.0*Nf*Tf*nF(pOverT) * ( 1.0-nF(primeOverT) );
        
        fval[0] = Jacobian*pperp*((pOverT-pz)/pOverT)*(Bose+Fermi)/(std::pow(2.0*M_PI,3.0));

        if(pOverT < 1e-8 ){
            // std::cerr << fval[0] << " " << qOvergT << " " << pperp << " " << pz << " " << Phi << "\n";
            fval[0] = 0.0 ;
        }
        return 1;
    }

    // Hard Part //
    double HardIntegralIntegration(double qOvergT){
        if(qOvergT == 0.0 ){
            return mDSqrOvergTSqr;
        }
        else{
            double epsrel = 1e-8;
            double epsabs = 0.0;

            double integral[1], error[1], prob[1];
            int nregions, neval, fail;

            Cuhre(3, 1, HardIntegrand, &qOvergT, 1,
                epsrel, epsabs, 0,
                MINEVAL, MAXEVAL, KEYY,
                STATEFILE, SPIN,
                &nregions, &neval, &fail, integral, error, prob);

            // double xmin[3] = {0.0,0.0,0.0};
            // double xmax[3] = {1.0,1.0,1.0};
            // hcubature(3,HardIntegrand,&qOvergT,3,xmin,xmax,MAXEVAL,epsabs,epsrel,ERROR_INDIVIDUAL,integral,error);
            // Divonne(3, 1, HardIntegrand, &qOvergT, 1,
            //     epsrel, epsabs, 0, 654,
            //     MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
            //     BORDER, MAXCHISQ, MINDEVIATION,
            //     NGIVEN, 3, NULL, NEXTRA, NULL,
            //     STATEFILE, SPIN,
            //     &nregions, &neval, &fail, integral, error, prob);
            
            return integral[0];
        }
    }


    // Hard Part Spline evaluation//
    double HardIntegralEval(double qOvergT){
        if(qOvergT < UpperCutoff && qOvergT>LowerCutoff){
            return gsl_spline_eval(HardIntegralSpline,qOvergT,Acc[omp_get_thread_num()]);
        }
        else{
            std::cerr << "HARD INTEGRAL OUT OF RANGE INERPOLATION AT Q=" << qOvergT << "\n";
            exit(0);
        }
    }

    // Full Hard Part //
    double HardIntegralCombined(double qOvergT){
        return HardIntegralInterpolation.GetFct(qOvergT);
    }

    void SetupHardIntegralSpline(){

        std::function<double(double)> IRLimit = [](double x){ return mDSqrOvergTSqr - g*C_A*x/16.0;};
        std::function<double(double)> UVLimit = [](double x){ return MathcalN;};

        HardIntegralInterpolation.Construct(HardIntegralIntegration,LowerCutoff,UpperCutoff,IRLimit, UVLimit,SplineSize,gsl_interp_linear,1);

        std::cerr << gsl_spline_eval_deriv(HardIntegralInterpolation.Spline,LowerCutoff,HardIntegralInterpolation.Acc[omp_get_thread_num()]) << "\n";
    }
}


namespace LO{
    SplineSetup LObSpaceSetup;
    int SplineSize = 2048;
    double LowerCutoff = 1e-1;
    double UpperCutoff = 1000.0;
    double SmallbSqrLogCoef = -MathcalN/(8.0*M_PI);
    double SmallbSqrCoef = MathcalN/(8.0*M_PI)*(1.0-EULERGAMMA-std::log(mDOvergT*0.5));
    double ConstantLargebLimit, OneOverbLargebLimit;

    double LargebConstant ;
    

    // Small q behavior //
    double SmallQ(double qOvergT){
        return (mDSqrOvergTSqr-qOvergT*g*C_A/16.0)/(qOvergT*qOvergT*(qOvergT*qOvergT+mDSqrOvergTSqr));
    }

    // Large q behavior //
    double LargeQ(double qOvergT){
        return MathcalN/(qOvergT*qOvergT*(qOvergT*qOvergT+mDSqrOvergTSqr));
    }


    // The leading order kernel in momentum space//
    double CombinedTGammaq(double qOvergT){
        return HardIntegralSpace::HardIntegralCombined(qOvergT)/(qOvergT*qOvergT*(qOvergT*qOvergT+mDSqrOvergTSqr));
    }

}

namespace NLO{
    SplineSetup NLObSpaceSetup;
    int SplineSize = 2048;
    double LowerCutoff = 1e-5;
    double UpperCutoff = 1e3;

    // Subtraction //
    inline double Subtraction(double qOvergT){
        return mDSqrOvergTSqr/std::pow(qOvergT,4.0) - g*C_A/(16.0*std::pow(qOvergT,3.0));
    }


    // Hard Part //
    double HardGammaWithSubtraction(double qOvergT){
        return (HardIntegralSpace::HardIntegralCombined(qOvergT)-mDSqrOvergTSqr+g*C_A/16.0*qOvergT)/std::pow(qOvergT,4.0) ;
    }

    // Soft LO Part //
    inline double SoftLOGamma( double qOvergT ){
        return mDSqrOvergTSqr/(qOvergT*qOvergT*(qOvergT*qOvergT+mDSqrOvergTSqr));
    }

    // Soft NLO Part //
    double SoftNLOGamma( double qOvergT ){
        if(qOvergT < 1e-3){
            return g*C_A*7.0/(32.0*std::pow(qOvergT,3.0));
        }
        else{
            double Res;
            Res =  7.0/(32.0*std::pow(qOvergT,3.0));
            Res += - ( mDOvergT + 2.0*(qOvergT-mDSqrOvergTSqr/qOvergT)*std::atan(qOvergT/mDOvergT) )/( 4.0*M_PI*std::pow((qOvergT*qOvergT+mDSqrOvergTSqr),2.0) );
            Res += ( mDOvergT - (  0.5*qOvergT + 2.0*mDSqrOvergTSqr/qOvergT )*std::atan(0.5*qOvergT/mDOvergT) )/( 8.0*M_PI*std::pow(qOvergT,4.0) );
            Res += -std::atan(qOvergT/mDOvergT)/( 2.0*M_PI*qOvergT*(qOvergT*qOvergT+mDSqrOvergTSqr) );
            Res += std::atan(0.5*qOvergT/mDOvergT)/( 2.0*M_PI*std::pow(qOvergT,3.0));
            Res += mDOvergT/(4.0*M_PI*(qOvergT*qOvergT+mDSqrOvergTSqr))*
            ( 3.0/(qOvergT*qOvergT+4.0*mDSqrOvergTSqr) -2.0/(qOvergT*qOvergT+mDSqrOvergTSqr) -1.0/(qOvergT*qOvergT));
            return g*C_A*Res;
        }
    }

    double CombinedTGammaq(double qOvergT){
        if (qOvergT <1e-3){
            return  SoftNLOGamma(qOvergT)+SoftLOGamma(qOvergT);
        }
        else{
            return SoftLOGamma(qOvergT)+SoftNLOGamma(qOvergT)+HardGammaWithSubtraction(qOvergT);
        }
    }
   
}
