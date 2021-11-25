double EllipticE(double k){
    if( k < 1.0 ){
        const double y = 1.0 - k;
        gsl_sf_result rf,rd;
        const int rfstatus = gsl_sf_ellint_RF_e(0.0, y, 1.0, GSL_MODE_DEFAULT, &rf);
        const int rdstatus = gsl_sf_ellint_RD_e(0.0, y, 1.0, GSL_MODE_DEFAULT, &rd);
        return rf.val - k/3.0 * rd.val;
    }
    else{
        std::cerr << " Elliptic E out of range k = "<< k << "\n";
        exit(0);
    }
    return 0.0;

}

double EllipticK(double k){
    if( k < 1.0 ){
        const double y = 1.0 - k;
        gsl_sf_result rf;
        const int rfstatus = gsl_sf_ellint_RF_e(0.0, y, 1.0, GSL_MODE_DEFAULT, &rf);
        return rf.val;
    }
    else{
        std::cerr << " Elliptic K out of range k = "<< k << "\n";
        exit(0);
    }
    return 0.0;

}