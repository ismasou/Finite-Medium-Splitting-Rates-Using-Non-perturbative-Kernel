namespace EQCD
{

    double gi = std::sqrt(2.763515753);
    double blogbCoef = gi * gi * (0.027075) / C_F;

    // VALUES OF ASYMPTOTICS FOR b-> 0 //
    double bbLnbCoef = 3.0*ZETA_3/(4.0*pow(M_PI,3))*(1.0+Nf/4.0);
    double bbCoef = bbLnbCoef*(1.0-EULERGAMMA-std::log(mDOvergT*0.5)); //0.129*gSqr/4.0;
    double logCoef = C_F*g*g*( mDSqrOvergTSqr*(1.0/24.0 -1.0/(4*M_PI*M_PI))+Nc/(8*M_PI*M_PI ))/M_PI;


    double SigmaEQCD = (0.2867) * gi / C_F;
    int Ng;
    const int Npq = 512;
    double xyMin = 1e-3, xyMax = 1e3,qMin,qMax;
    double *GammaVals, *qVals;
    gsl_spline *GammaHatQInt;
    gsl_interp_accel **GammaHatQAcc, **xacc1, **xacc2;

    inline double LargeQ(double x)
    {
        return 8.0 * M_PI * blogbCoef / (x * x * x * x);
    }

    inline double SmallQ(double x)
    {
        return 2.0 * M_PI * SigmaEQCD / (x * x * x);
    }

    inline double GammaHatLargeBLimitTransform(double qOvergT){
        return  (2.0*M_PI)*EQCD::SigmaEQCD/(C_F)/(qOvergT*qOvergT*qOvergT)
               +(2.0*M_PI)*logCoef/C_F/(qOvergT*qOvergT);
    }


    inline double FullGamma(double q)
    {
        double qOvergT = q / (gi*Temp);
        if(qOvergT < qMin){
                return GammaHatLargeBLimitTransform(qOvergT);
            }
            else if (qOvergT > qMax){
                return MathcalN/pow(qOvergT,4.0);
            }
            else{
                int tID=omp_get_thread_num();
                if(!isfinite(gsl_spline_eval(GammaHatQInt,qOvergT,GammaHatQAcc[tID]))){
                    std::cerr << "NOT FINITE" << "\n";
                    exit(0);
                }
                return gsl_spline_eval(GammaHatQInt,qOvergT,GammaHatQAcc[tID]);
            }
    }
        
    void SetupFromFile(std::string fname){
        // Reading Spline //
        CSVReader ReadEQCDSpline(fname);
        std::vector<std::vector<double>> Data = ReadEQCDSpline.getData();

        int Nq = Data.size();

        // CREATE TABLES //
        std::vector<double> qVals(Nq);
        std::vector<double> GammaHatQVals(Nq);
        

        for(int i=0 ; i< Nq; i++){
            qVals[i] = Data[i][0];
            GammaHatQVals[i] = Data[i][2];
            // std::cerr << qVals[i] << " " << GammaHatQVals[i] << "\n";
        }

        qMin = qVals[0];
        qMax = qVals[Nq-1];

        // SETUP GSL INTERPOLATION //
        int NumberOfOpenMPThreads=omp_get_max_threads();
        
        GammaHatQAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        
        #pragma omp parallel for
        for(int i=0;i<NumberOfOpenMPThreads;i++){
            GammaHatQAcc[i]= gsl_interp_accel_alloc ();
        }
        
        GammaHatQInt = gsl_spline_alloc (gsl_interp_akima, Nq);
        gsl_spline_init(GammaHatQInt,qVals.data(),GammaHatQVals.data(),Nq);

        std::ofstream OutStream;
        OutStream.open((fname+"-Spline.txt").c_str());
        // OUTPUT SPLINE //
        int NqOut = 2048;
        for(int iq=0;iq<NqOut;iq++){
            
            double dq=(600-1e-4)/(NqOut-1); double q=1e-4+dq*iq;

            // CREATE OUTPUT //
            OutStream << q << " " << FullGamma(q) << " "; 
            OutStream << std::endl;
            
        }
        
        OutStream.close();
        
    }

    void Setup()
    {
        int DataPoint = 1;
        std::vector<std::vector<double>> Parameters;
        //                                             g**2       Temp(MeV)      Nf         A         Sigma    qhat    bbLogb
        Parameters.push_back(std::vector<double>{3.725027366, 250, 3, -0.67178761, 0.2836, 0.1465, -0.021106});
        Parameters.push_back(std::vector<double>{2.763515753, 500, 3, -0.48850742, 0.2867, 0.185, -0.027075});
        Parameters.push_back(std::vector<double>{2.210169201, 1000, 4, 0.11749608, 0.2901, 0.3136, -0.039115});
        Parameters.push_back(std::vector<double>{1.066560916, 100000, 5, 0.2323915, 0.2952, 0.5665, -0.077584});

        gi = std::sqrt(Parameters[DataPoint][0]);
        g  = gi;

        Nf = Parameters[DataPoint][2];
        double A  = Parameters[DataPoint][3];
        SigmaEQCD = Parameters[DataPoint][4] * g;
        bbLnbCoef = -(gi * gi / C_F) * Parameters[DataPoint][6];
        bbCoef    = (gi * gi / C_F) * Parameters[DataPoint][5] / 4.0 + (gi * gi / C_F) * Parameters[DataPoint][6] * std::log(gi);
        std::cerr << "# Lattice DataPoint = "<< DataPoint <<std::endl;
        std::stringstream EQCDSpline;
        EQCDSpline << "Data/CqTransform-T" << Parameters[DataPoint][1] << ".txt";
        
        SetupFromFile(EQCDSpline.str().c_str());
        xacc1 = new gsl_interp_accel *[NumberOfOpenMPThreads];
        xacc2 = new gsl_interp_accel *[NumberOfOpenMPThreads];
        #pragma omp parallel
        {
            int tID = omp_get_thread_num();
            xacc1[tID] = gsl_interp_accel_alloc();
            xacc2[tID] = gsl_interp_accel_alloc();
        }


        // SetupMoments();
        // std::ofstream RateFile;
        // RateFile.precision(15);
        // RateFile.open("OUTPUT/EQCDSpline.txt",std::ofstream::out);
        // for(int i=0; i<1024; i++){
        //     double x = 1e-4*std::exp(std::log(1e4/1e-4)*i/(1023.0));
        //     RateFile << x << " " << FullGamma(x) << "\n";
        // }
        // RateFile.close();
    }
} // namespace EQCD