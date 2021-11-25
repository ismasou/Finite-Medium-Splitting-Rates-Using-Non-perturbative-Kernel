#include <ctime>
#include <cstring>
#include <functional>
#include <cstdlib>
#include <complex>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <set>
#include <vector>
#include <algorithm>
#include <omp.h>

// GSL MATH ROUTINES //
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

// GSL INTERPOLATION //
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>

// GSL INTEGRATION ROUTINES //
#include <gsl/gsl_integration.h>

// Cubature integration //
#include "cubature/hcubature.c"

// DataFile Reader //
#ifndef CSVREADER
#define CSVREADER
#include "IO/CSV-Reader.cpp"
#endif 

#ifndef ZETA_2
#define ZETA_2 (M_PI * M_PI / 6.0)
#endif

#ifndef ZETA_3
#define ZETA_3 double(1.202056903159594)
#endif

#ifndef EULERGAMMA
#define EULERGAMMA (0.5772156649015329)
#endif

// BESSEL FUNCTIONS //
#ifndef BesselJ0
#define BesselJ0(x) gsl_sf_bessel_Jn(0, x)
#endif

#ifndef BesselJ1
#define BesselJ1(x) gsl_sf_bessel_Jn(1, x)
#endif

#ifndef BesselK0
#define BesselK0(x) gsl_sf_bessel_Kn(0, x)
#endif

#ifndef BesselK1
#define BesselK1(x) gsl_sf_bessel_Kn(1, x)
#endif

#ifndef BesselK2
#define BesselK2(x) gsl_sf_bessel_Kn(2, x)
#endif

namespace BroadeningKernel
{

    // COLOR, FLAVOR AND DEGENERACY FACTORS //
    static const double Nc = 3.0;
    double Nf = 3.0;
    static const double C_A = Nc;
    static const double C_F = (Nc * Nc - 1.0) / (2.0 * Nc);
    static const double Tf = 0.5;
    static const double nuG = 2.0 * (Nc * Nc - 1.0);
    static const double nuQ = 2.0 * Nc;
    static const double mDSqrOvergTSqr = (Nc / 3.0 + Nf / 6.0);
    static const double mDOvergT = std::sqrt(mDSqrOvergTSqr);
    static const double MathcalN = (ZETA_3 / ZETA_2) * (1.0 + 0.25 * Nf);

    double g = std::sqrt(2.763515753);
    double T = 1.0;
    double mDSqr = g * g * T * T * mDSqrOvergTSqr;
    double MeffQSqr = C_F * g * g * T * T / 4.0;

    int NumberOfOpenMPThreads;
    gsl_integration_workspace **wsp;
    int NumIntegrationPoints = 2048;

    void SetupGslIntegration()
    {
        wsp = new gsl_integration_workspace *[NumberOfOpenMPThreads];
#pragma omp parallel
        {
            wsp[omp_get_thread_num()] = gsl_integration_workspace_alloc(NumIntegrationPoints);
        }

        // Set gsl Error Off //
        gsl_set_error_handler_off();
    }

    void FreeGslIntegration()
    {
#pragma omp parallel
        {
            gsl_integration_workspace_free(wsp[omp_get_thread_num()]);
        }
    }

    namespace NonPerturbative
    {
        std::vector<std::vector<double>> Parameters;

        // COUPLING CONSTANT //
        double gSqr = 2.763515753;
        double g = sqrt(gSqr);

        // VALUES OF ASYMPTOTICS FOR b-> 0 //
        double bbLnbCoef = 3.0 * ZETA_3 / (4.0 * pow(M_PI, 3)) * (1.0 + Nf / 4.0);
        double bbCoef = bbLnbCoef * (1.0 - EULERGAMMA - std::log(mDOvergT * 0.5));
        double logCoef = C_F * g * g * (mDSqrOvergTSqr * (1.0 / 24.0 - 1.0 / (4 * M_PI * M_PI)) + Nc / (8 * M_PI * M_PI)) / M_PI;

        // VALUES OF ASYMPTOTICS FOR b->INFINITY //
        double A;
        double SigmaEQCD;

        // GSL INTERPOLATION OBJECTS //
        gsl_spline *GammaSpline;
        gsl_spline *GammaDerivativeSpline;
        int LinearInterpPoints;
        double *BperpVals_BSpline, *GammaVals_BSpline, *GammaDerivativeVals_BSpline;

        gsl_interp_accel **xacc, **yacc;

        inline double GammaHatSmallBLimit(double gTb)
        {
            return (-bbLnbCoef * (gTb) * (gTb)*std::log(gTb) + bbCoef * (gTb) * (gTb));
        }

        inline double GammaHatSmallBLimitDerivative(double gTb)
        {
            return (-bbLnbCoef * gTb * (1.0 + 2.0 * std::log(gTb)) + 2.0 * bbCoef * (gTb));
        }

        inline double GammaHatLargeBLimit(double gTb)
        {
            return (A + SigmaEQCD * (gTb) + logCoef * std::log(gTb)) / (C_F);
        }
        inline double GammaHatLargeBLimitDerivative(double gTb)
        {
            return (SigmaEQCD + logCoef / gTb) / (C_F);
        }

        inline double GammaHatInterpolated(double gTb)
        {
            if (gTb >= BperpVals_BSpline[0] && gTb <= BperpVals_BSpline[LinearInterpPoints - 1])
            {
                return gsl_spline_eval(GammaSpline, gTb, xacc[omp_get_thread_num()]);
            }
            else
            {
                std::cout << "INTERPOLATION OUT OF RANGE b =" << gTb << std::endl;
                exit(0);
                // return 0.0;
            }
        }

        inline double GammaHatInterpolatedDerivative(double gTb)
        {
            if (gTb >= BperpVals_BSpline[0] && gTb <= BperpVals_BSpline[LinearInterpPoints - 1])
            {
                return gsl_spline_eval(GammaDerivativeSpline, gTb, yacc[omp_get_thread_num()]);
            }
            else
            {
                std::cout << "INTERPOLATION OUT OF RANGE b =" << gTb << std::endl;
                exit(0);
                // return 0.0;
            }
        }

        inline double GammabOvergSqrT(double gTb)
        {

            // SMALL B-LIMIT //
            if (gTb <= BperpVals_BSpline[0])
            {
                return GammaHatSmallBLimit(gTb);
            }

            // SPLINE //
            else if (gTb < BperpVals_BSpline[LinearInterpPoints - 1])
            {
                return GammaHatInterpolated(gTb);
            }

            // LARGE B-LIMIT //
            else
            {
                return GammaHatLargeBLimit(gTb);
            }
        }

        inline double GammabOvergSqrTDerivative(double gTb)
        {

            // SMALL B-LIMIT //
            if (gTb <= BperpVals_BSpline[0])
            {
                return GammaHatSmallBLimitDerivative(gTb);
            }

            // SPLINE //
            else if (gTb < BperpVals_BSpline[LinearInterpPoints - 1])
            {
                return GammaHatInterpolatedDerivative(gTb);
            }

            // LARGE B-LIMIT //
            else
            {
                return GammaHatLargeBLimitDerivative(gTb);
            }
        }

        void Setup(int DataPoint)
        {

            // SETUP OMP //
            NumberOfOpenMPThreads = omp_get_max_threads();

            //                                             g**2       Temp(MeV)      Nf         A         Sigma    qhat    bbLogb
            Parameters.push_back(std::vector<double>{3.725027366, 250, 3, -0.67178761, 0.2836, 0.1465, -0.021106});
            Parameters.push_back(std::vector<double>{2.763515753, 500, 3, -0.48850742, 0.2867, 0.185, -0.027075});
            Parameters.push_back(std::vector<double>{2.210169201, 1000, 4, 0.11749608, 0.2901, 0.3136, -0.039115});
            Parameters.push_back(std::vector<double>{1.066560916, 100000, 5, 0.2323915, 0.2952, 0.5665, -0.077584});

            g = std::sqrt(Parameters[DataPoint][0]);
            Nf = Parameters[DataPoint][2];
            A = Parameters[DataPoint][3];
            SigmaEQCD = Parameters[DataPoint][4] * g;
            bbLnbCoef = -(g * g / C_F) * Parameters[DataPoint][6];
            bbCoef = (g * g / C_F) * Parameters[DataPoint][5] / 4.0 + (g * g / C_F) * Parameters[DataPoint][6] * std::log(g);

            std::stringstream EQCDSpline;
            EQCDSpline << "Data/EQCDSpline-T" << Parameters[DataPoint][1] << ".txt";

            // Reading Spline //
            CSVReader ReadEQCDSpline(EQCDSpline.str().c_str());
            std::vector<std::vector<double>> Data = ReadEQCDSpline.getData();
            LinearInterpPoints = Data.size();
            BperpVals_BSpline = new double[LinearInterpPoints];
            GammaVals_BSpline = new double[LinearInterpPoints];
            GammaDerivativeVals_BSpline = new double[LinearInterpPoints];

            std::cerr << "# Reading File with n = " << LinearInterpPoints << " Points \n";
            for (int i = 0; i < LinearInterpPoints; i++)
            {
                BperpVals_BSpline[i] = Data[i][0];
                GammaVals_BSpline[i] = Data[i][1];
                GammaDerivativeVals_BSpline[i] = Data[i][2];
            }

            GammaSpline = gsl_spline_alloc(gsl_interp_cspline, LinearInterpPoints);
            GammaDerivativeSpline = gsl_spline_alloc(gsl_interp_linear, LinearInterpPoints);
            gsl_spline_init(GammaSpline, BperpVals_BSpline, GammaVals_BSpline, LinearInterpPoints);
            gsl_spline_init(GammaDerivativeSpline, BperpVals_BSpline, GammaDerivativeVals_BSpline, LinearInterpPoints);

            xacc = new gsl_interp_accel *[NumberOfOpenMPThreads];
            yacc = new gsl_interp_accel *[NumberOfOpenMPThreads];
            for (int tID = 0; tID < NumberOfOpenMPThreads; tID++)
            {
                xacc[tID] = gsl_interp_accel_alloc();
                yacc[tID] = gsl_interp_accel_alloc();
            }
            // for (int i = 0; i < LinearInterpPoints; i++)
            // {
            //     double b = Data[i][0];
            //     std::cout << Data[i][0] << " ";
            //     std::cout << Data[i][1] << " ";
            //     std::cout << GammaPerturb(Data[i][0]) << " ";
            //     std::cout << MathcalN * b * b / (16.0 * M_PI) * std::log(Xip * mDSqrOvergTSqr * b * b / 4.0 / FityNP) << std::endl;
            // }
        }

    }

    ///////////////////////////
    // Leading-Order Kernel  //
    ///////////////////////////

    class SplineSetup
    {
    public:
        double LowerCutOff, UpperCutOff, *xVals, *FctVals;
        std::function<double(double)> FctToInterpolate, UVLimitFct, IRLimitFct;
        int SplineSize;

        gsl_spline *Spline;
        gsl_interp_accel **Acc;

        void SetUpValuesLogGrid()
        {
            xVals = new double[SplineSize];
            FctVals = new double[SplineSize];
            // Computing the spline
            xVals[0] = LowerCutOff;
            FctVals[0] = FctToInterpolate(xVals[0]);
            LowerCutOff += 1e-10;
#pragma omp parallel for schedule(dynamic)
            for (int i = 1; i < SplineSize; i++)
            {
                // xVals[i] = LowerCutOff + (UpperCutOff-LowerCutOff)*i/double(SplineSize-1);
                xVals[i] = LowerCutOff * std::exp(std::log(UpperCutOff / LowerCutOff) * i / double(SplineSize - 1));
                FctVals[i] = FctToInterpolate(xVals[i]);
            }
        }

        void SetUpValues()
        {
            xVals = new double[SplineSize];
            FctVals = new double[SplineSize];
// Computing the spline
#pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < SplineSize; i++)
            {
                xVals[i] = LowerCutOff + (UpperCutOff - LowerCutOff) * i / double(SplineSize - 1);
                FctVals[i] = FctToInterpolate(xVals[i]);
            }
        }

        inline double GetFct(double x)
        {

            if (x <= LowerCutOff)
            {
                return IRLimitFct(x);
            }
            else if (x >= UpperCutOff)
            {
                return UVLimitFct(x);
            }
            else
            {
                return gsl_spline_eval(Spline, x, Acc[omp_get_thread_num()]);
            }
        }

        inline double EvalDerive(double x)
        {
            return gsl_spline_eval_deriv(Spline, x, Acc[omp_get_thread_num()]);
        }

        SplineSetup()
        {
        }

        void Construct(std::function<double(double)> FctToInterpolate_, double LowerCutOff_, double UpperCutOff_, std::function<double(double)> IRLimitFct_, std::function<double(double)> UVLimitFct_, int SplineSize_ = 512, const gsl_interp_type *InterpType = gsl_interp_linear, const int LogGrid = 1)
        {
            this->FctToInterpolate = FctToInterpolate_;
            this->LowerCutOff = LowerCutOff_;
            this->UpperCutOff = UpperCutOff_;
            this->IRLimitFct = IRLimitFct_;
            this->UVLimitFct = UVLimitFct_;
            this->SplineSize = SplineSize_;

            if (LogGrid == 1)
            {
                SetUpValuesLogGrid();
            }
            else
            {
                SetUpValues();
            }

            Spline = gsl_spline_alloc(InterpType, SplineSize);
            Acc = new gsl_interp_accel *[NumberOfOpenMPThreads];

#pragma omp parallel
            {
                Acc[omp_get_thread_num()] = gsl_interp_accel_alloc();
            }

            gsl_spline_init(Spline, xVals, FctVals, SplineSize);
        }

        ~SplineSetup()
        {
            delete Acc;
            delete xVals;
            delete FctVals;
            delete Spline;
        }
    };

    namespace HardIntegralSpace
    {
        SplineSetup HardIntegralInterpolation;
        int SplineSize = 512;
        gsl_spline *HardIntegralSpline;
        gsl_interp_accel **Acc;
        double LowerCutoff = 1e-3;
        double UpperCutoff = 100.0;

        // Computing 'I' Helper function for LO kernel //
        inline double I(double qOverT, int n, int m)
        {

            double Bes = BesselK2(qOverT * std::sqrt(m * n));
            double Coef = (m * n) / (2.0 * std::pow((m + n), 3.0));

            return Coef * (qOverT * qOverT) * Bes;
        }

        // HardIntegral Summed //
        double HardIntegralSummed(double qOvergT)
        {

            double Res = 0.0;
            double OldRes = 0.0;
            double Tol = 1e-5;
            int SumLimit = 2048;

            for (int n = 1; n < SumLimit; n++)
            {

                for (int m = 1; m < SumLimit; m++)
                {
                    double sign = 1.0;
                    if ((m + n - 1) % 2)
                    {
                        sign = -1.0;
                    }

                    Res += (2.0 * C_A + sign * 4.0 * Nf * Tf) * I(g * qOvergT, n, m);
                }
                if (std::abs(1.0 - Res / OldRes) < Tol && n > 512)
                {
                    goto OutOfHere;
                }
                OldRes = Res;
            }
            // Error:
            std::cerr << "Finshed the sum without reaching tolerance on LO Gammaq \n";
            std::cerr << "qOvergT= " << qOvergT << "Gammaq= " << Res << std::endl;
            exit(0);
        OutOfHere:
            Res += (2.0 * C_A + 3.0 * Nf * Tf) * ZETA_3;
            return Res / (M_PI * M_PI);
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
        int HardIntegrand(unsigned int dim, const double xx[],
                          void *userdata, unsigned int fdim, double fval[])
        {
            double qOvergT = *((double *)userdata);

            // Change Of Variable //
            double pperp = std::tan(M_PI * 0.5 * xx[0]);
            double pz = std::tan(M_PI * (xx[1] - 0.5));
            double Phi = xx[2] * (2 * M_PI);

            // Jacobian //
            double Jacobian = (M_PI / (1.0 + std::cos(M_PI * xx[0])));               // pperp jacobian
            Jacobian *= (2.0 * M_PI / (1.0 + std::cos(M_PI * (2.0 * xx[1] - 1.0)))); // pz jacobian
            Jacobian *= (2.0 * M_PI);                                                // Phi jacobian

            double pOverT = std::sqrt(pperp * pperp + pz * pz);

            double primeOverT = pOverT + g * (g * qOvergT * qOvergT + 2.0 * qOvergT * pperp * std::cos(Phi)) / (2.0 * (pOverT - pz));

            double Bose = 2.0 * C_A * nB(pOverT) * (1.0 + nB(primeOverT));
            double Fermi = 4.0 * Nf * Tf * nF(pOverT) * (1.0 - nF(primeOverT));

            fval[0] = Jacobian * pperp * ((pOverT - pz) / pOverT) * (Bose + Fermi) / (std::pow(2.0 * M_PI, 3.0));

            if (pOverT < 1e-8)
            {
                // std::cerr << fval[0] << " " << qOvergT << " " << pperp << " " << pz << " " << Phi << "\n";
                fval[0] = 0.0;
            }
            return 1;
        }

        // Hard Part //
        double HardIntegralIntegration(double qOvergT)
        {
            if (qOvergT == 0.0)
            {
                return mDSqrOvergTSqr;
            }
            else
            {
                double epsrel = 1e-8;
                double epsabs = 0.0;

                double integral[1], error[1];
                double xmin[3] = {0.0,0.0,0.0};
                double xmax[3] = {1.0,1.0,1.0};


                hcubature(1,HardIntegrand,&qOvergT,3,xmin,xmax,10000,epsabs,epsrel,ERROR_INDIVIDUAL,integral,error);

                return integral[0];
            }
        }

        // Hard Part Spline evaluation//
        double HardIntegralEval(double qOvergT)
        {
            if (qOvergT < UpperCutoff && qOvergT > LowerCutoff)
            {
                return gsl_spline_eval(HardIntegralSpline, qOvergT, Acc[omp_get_thread_num()]);
            }
            else
            {
                std::cerr << "HARD INTEGRAL OUT OF RANGE INERPOLATION AT Q=" << qOvergT << "\n";
                exit(0);
            }
        }

        // Full Hard Part //
        double HardIntegralCombined(double qOvergT)
        {
            return HardIntegralInterpolation.GetFct(qOvergT);
        }

        void SetupHardIntegralSpline()
        {

            // SETUP OMP //
            NumberOfOpenMPThreads = omp_get_max_threads();

            std::function<double(double)> IRLimit = [](double x) { return mDSqrOvergTSqr - g * C_A * x / 16.0; };
            std::function<double(double)> UVLimit = [](double x) { return MathcalN; };

            HardIntegralInterpolation.Construct(HardIntegralSummed, LowerCutoff, UpperCutoff, IRLimit, UVLimit, SplineSize, gsl_interp_linear, 1);
        }
    }

    namespace FourierIntegration
    {
        struct IntegrationVariables
        {
            double gTb;
            double qOvergT;
            std::function<double(double)> FctToIntegrate;
        };

        double IntegrandBesselFromQToB(double x, void *userdata)
        {
            IntegrationVariables Vars = *((IntegrationVariables *)userdata);
            double gTb = Vars.gTb;
            double qOvergT = mDOvergT * std::tan(mDOvergT * x);
            double Gammaq = Vars.FctToIntegrate(qOvergT);
            // Jacobian //
            double Jacobian = (qOvergT * qOvergT + mDSqrOvergTSqr);

            if (qOvergT * gTb < 1E-5)
            {
                return Jacobian * qOvergT * (0.5 / M_PI) * (0.25 * qOvergT * gTb * qOvergT * gTb) * Gammaq;
            }
            else
            {
                return Jacobian * qOvergT * (0.5 / M_PI) * (1.0 - BesselJ0(qOvergT * gTb)) * Gammaq;
            }
        }

        double FourierIntegrationFromQtoB(double gTb, std::function<double(double)> Gammaq)
        {
            double BesselIntegral, error;
            IntegrationVariables Vars;
            Vars.gTb = gTb;
            Vars.FctToIntegrate = Gammaq;
            gsl_function Fct;
            Fct.params = &Vars;
            Fct.function = &IntegrandBesselFromQToB;

            gsl_integration_qags(&Fct, 0, (0.5 * M_PI / mDOvergT), 0, 1e-9, NumIntegrationPoints, wsp[omp_get_thread_num()], &BesselIntegral, &error);
            return BesselIntegral;
        }

        double TransformationOfIRSubtraction(double qOvergT)
        {
            return 0.0;
        }

        double TransformationOfUVSubtraction(double qOvergT)
        {
            return 0.0;
        }
    }

    namespace LO
    {
        SplineSetup LObSpaceSetup;
        int SplineSize = 2048;
        double LowerCutoff = 1e-1;
        double UpperCutoff = 1000.0;
        double SmallbSqrLogCoef = -MathcalN / (8.0 * M_PI);
        double SmallbSqrCoef = MathcalN / (8.0 * M_PI) * (1.0 - EULERGAMMA - std::log(mDOvergT * 0.5));
        double ConstantLargebLimit, OneOverbLargebLimit;

        double LargebConstant;

        // Small q behavior //
        double SmallQ(double qOvergT)
        {
            return (mDSqrOvergTSqr - qOvergT * g * C_A / 16.0) / (qOvergT * qOvergT * (qOvergT * qOvergT + mDSqrOvergTSqr));
        }

        // Large q behavior //
        double LargeQ(double qOvergT)
        {
            return MathcalN / (qOvergT * qOvergT * (qOvergT * qOvergT + mDSqrOvergTSqr));
        }

        // The leading order kernel in momentum space//
        double CombinedTGammaq(double qOvergT)
        {
            return HardIntegralSpace::HardIntegralCombined(qOvergT) / (qOvergT * qOvergT * (qOvergT * qOvergT + mDSqrOvergTSqr));
        }

        // The leading order kernel for large b in impact parameter space //
        inline double LargebGammabOvergSqrT(double gTb)
        {
            return (0.5 / M_PI) * (EULERGAMMA + std::log(0.5 * mDOvergT * gTb) + BesselK0(gTb * mDOvergT));
        }

        inline double LinearPartLargeGammabOvergSqrT(double gTb)
        {
            // return -g*C_A/16.0*( 0.25/mDOvergT - 0.5/(M_PI*gTb*mDSqrOvergTSqr) );
            return ConstantLargebLimit; //+OneOverbLargebLimit/gTb;
        }

        // The leading order kernel for small b in impact parameter space //
        inline double SmallbGammabOvergSqrT(double gTb)
        {
            return (MathcalN / mDSqrOvergTSqr) * (0.5 / M_PI) * (EULERGAMMA + std::log(0.5 * mDOvergT * gTb) + BesselK0(gTb * mDOvergT));
        }

        // The leading order kernel in impact parameter space //
        double GammabOvergSqrTIntegration(double gTb)
        {
            std::function<double(double)> Subtracted = [](double x) { return CombinedTGammaq(x) - LargeQ(x); };
            // std::function < double(double)> Subtracted = [] (double x){ return 1.0/(x*(x*x+mDSqrOvergTSqr));};
            // std::function < double(double)> Subtracted = [] (double x){ return CombinedTGammaq(x)-SmallQ(x);};
            return FourierIntegration::FourierIntegrationFromQtoB(gTb, Subtracted) + SmallbGammabOvergSqrT(gTb);
        }

        double GammabOvergSqrT(double gTb)
        {
            return LObSpaceSetup.GetFct(gTb);
        }

        void SetupBspaceSpline()
        {


            // SETUP OMP //
            NumberOfOpenMPThreads = omp_get_max_threads();

            // Setup GSL Integration space //
            SetupGslIntegration();

            ConstantLargebLimit = GammabOvergSqrTIntegration(1e4) - LargebGammabOvergSqrT(1e4);
            OneOverbLargebLimit = (GammabOvergSqrTIntegration(70) - LargebGammabOvergSqrT(70) - ConstantLargebLimit) * 70;
            // std::cerr << ConstantLargebLimit << "+" << OneOverbLargebLimit << "/b \n";

            std::function<double(double)> UVLimit = [](double gTb) { return LargebGammabOvergSqrT(gTb) + LinearPartLargeGammabOvergSqrT(gTb); };

            LObSpaceSetup.Construct(GammabOvergSqrTIntegration, LowerCutoff, UpperCutoff, SmallbGammabOvergSqrT, UVLimit, SplineSize, gsl_interp_cspline, 1);

            // Free Gsl Integration //
            FreeGslIntegration();
        }

    }

    namespace NLO
    {

        SplineSetup NLObSpaceSetup;
        int SplineSize = 2048;
        double LowerCutoff = 1e-5;
        double UpperCutoff = 1e3;

        // Subtraction //
        inline double Subtraction(double qOvergT)
        {
            return mDSqrOvergTSqr / std::pow(qOvergT, 4.0) - g * C_A / (16.0 * std::pow(qOvergT, 3.0));
        }

        // Hard Part //
        double HardGammaWithSubtraction(double qOvergT)
        {
            return (HardIntegralSpace::HardIntegralCombined(qOvergT) - mDSqrOvergTSqr + g * C_A / 16.0 * qOvergT) / std::pow(qOvergT, 4.0);
        }

        // Soft LO Part //
        inline double SoftLOGamma(double qOvergT)
        {
            return mDSqrOvergTSqr / (qOvergT * qOvergT * (qOvergT * qOvergT + mDSqrOvergTSqr));
        }

        // Soft NLO Part //
        double SoftNLOGamma(double qOvergT)
        {
            if (qOvergT < 1e-3)
            {
                return g * C_A * 7.0 / (32.0 * std::pow(qOvergT, 3.0));
            }
            else
            {
                double Res;
                Res = 7.0 / (32.0 * std::pow(qOvergT, 3.0));
                Res += -(mDOvergT + 2.0 * (qOvergT - mDSqrOvergTSqr / qOvergT) * std::atan(qOvergT / mDOvergT)) / (4.0 * M_PI * std::pow((qOvergT * qOvergT + mDSqrOvergTSqr), 2.0));
                Res += (mDOvergT - (0.5 * qOvergT + 2.0 * mDSqrOvergTSqr / qOvergT) * std::atan(0.5 * qOvergT / mDOvergT)) / (8.0 * M_PI * std::pow(qOvergT, 4.0));
                Res += -std::atan(qOvergT / mDOvergT) / (2.0 * M_PI * qOvergT * (qOvergT * qOvergT + mDSqrOvergTSqr));
                Res += std::atan(0.5 * qOvergT / mDOvergT) / (2.0 * M_PI * std::pow(qOvergT, 3.0));
                Res += mDOvergT / (4.0 * M_PI * (qOvergT * qOvergT + mDSqrOvergTSqr)) *
                       (3.0 / (qOvergT * qOvergT + 4.0 * mDSqrOvergTSqr) - 2.0 / (qOvergT * qOvergT + mDSqrOvergTSqr) - 1.0 / (qOvergT * qOvergT));
                return g * C_A * Res;
            }
        }

        double CombinedTGammaq(double qOvergT)
        {
            if (qOvergT < 1e-3)
            {
                return SoftNLOGamma(qOvergT) + SoftLOGamma(qOvergT);
            }
            else
            {
                return SoftLOGamma(qOvergT) + SoftNLOGamma(qOvergT) + HardGammaWithSubtraction(qOvergT);
            }
        }

        // Hard Integral Fourier transform //
        double HardGammaSubtractedbOvergSqrT(double gTb)
        {
            return FourierIntegration::FourierIntegrationFromQtoB(gTb, HardGammaWithSubtraction);
        }

        // Soft LO Fourier transform //
        double SoftLOGammabOvergSqrT(double gTb)
        {
            return (0.5 / M_PI) * (EULERGAMMA + std::log(0.5 * mDOvergT * gTb) + BesselK0(gTb * mDOvergT));
        }

        // Soft NLO Fourier transform //
        double SoftNLOGammabOvergSqrT(double gTb)
        {
            std::function<double(double)> SoftNLOGammaSubtracted = [](double x) { return SoftNLOGamma(x) + g * C_A / 16.0 / (x * x * x); };
            return FourierIntegration::FourierIntegrationFromQtoB(gTb, SoftNLOGammaSubtracted) - (0.5 / M_PI) * gTb * g * C_A / 16.0;
        }

        // The next to leading order kernel in impact parameter space //
        double GammabOvergSqrTIntegration(double gTb)
        {
            // std::function<double(double)> Subtracted = [] (double qOvergT){ return CombinedTGammaq(qOvergT) - g*C_A*7.0/32.0/(qOvergT*qOvergT*qOvergT);};
            // return FourierIntegration::FourierIntegrationFromQtoB(gTb,Subtracted) + (0.5/M_PI)*gTb*g*C_A*7/32.0;
            std::function<double(double)> Subtracted = [](double qOvergT) { return CombinedTGammaq(qOvergT) - MathcalN / (qOvergT * qOvergT * (qOvergT * qOvergT + mDSqrOvergTSqr)); };
            return FourierIntegration::FourierIntegrationFromQtoB(gTb, Subtracted) + (MathcalN / mDSqrOvergTSqr) * SoftLOGammabOvergSqrT(gTb);
        }

        // The next to leading order kernel in impact parameter space //
        double GammabOvergSqrT(double gTb)
        {
            return NLObSpaceSetup.GetFct(gTb);
        }

        double UVLimitInBSpace(double gTb)
        {
            return (0.5 / M_PI) * gTb * g * C_A * 7.0 / 32.0 + SoftLOGammabOvergSqrT(gTb);
        }

        double IRLimitInBSpace(double gTb)
        {
            return (MathcalN / mDSqrOvergTSqr) * SoftLOGammabOvergSqrT(gTb);
        }

        void SetupBspaceSpline()
        {

            // SETUP OMP //
            NumberOfOpenMPThreads = omp_get_max_threads();

            
            // Setup GSL Integration space //
            SetupGslIntegration();
            NLObSpaceSetup.Construct(GammabOvergSqrTIntegration, LowerCutoff, UpperCutoff, IRLimitInBSpace, UVLimitInBSpace, SplineSize, gsl_interp_cspline, 1);
            // Free Gsl Integration //
            FreeGslIntegration();
        }
    }

}