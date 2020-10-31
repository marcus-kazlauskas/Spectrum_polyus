//
//  functions.h
//  Program_4
//
//  Created by Kozlov Mikhail on 10.11.16.
//  Copyright © 2016 Kozlov Mikhail. All rights reserved.
//

#ifndef functions_h
#define functions_h

#include <math.h>
#include <complex>

const double pi = 3.14159;

double beta(double teta){
    return pi*cos(teta)/2;
}

double p(double n, double teta){
    return n*cos(teta);
}
double q(double n, double teta){
    return cos(teta)/n;
}

double teta(double teta1, double n){
    return acos(sqrt(1-pow(sin(teta1/n), 2)));
}

class matrix2{
public:
    std::complex<double> mE[2][2], mM[2][2];
    std::complex<double> rTE, rTM;
    double RTE, RTM;
    
    void mulE(std::complex<double> m[2][2]){
        std::complex<double> mBuf[2][2];
        
        mBuf[0][0] = mE[0][0];
        mBuf[0][1] = mE[0][1];
        mBuf[1][0] = mE[1][0];
        mBuf[1][1] = mE[1][1];
        
        mE[0][0] = mBuf[0][0]*m[0][0]+mBuf[0][1]*m[1][0];
        mE[0][1] = mBuf[0][0]*m[0][1]+mBuf[0][1]*m[1][1];
        mE[1][0] = mBuf[1][0]*m[0][0]+mBuf[1][1]*m[1][0];
        mE[1][1] = mBuf[1][0]*m[0][1]+mBuf[1][1]*m[1][1];
    }
    
    void mulM(std::complex<double> m[2][2]){
        std::complex<double> mBuf[2][2];
        
        mBuf[0][0] = mM[0][0];
        mBuf[0][1] = mM[0][1];
        mBuf[1][0] = mM[1][0];
        mBuf[1][1] = mM[1][1];
        
        mM[0][0] = mBuf[0][0]*m[0][0]+mBuf[0][1]*m[1][0];
        mM[0][1] = mBuf[0][0]*m[0][1]+mBuf[0][1]*m[1][1];
        mM[1][0] = mBuf[1][0]*m[0][0]+mBuf[1][1]*m[1][0];
        mM[1][1] = mBuf[1][0]*m[0][1]+mBuf[1][1]*m[1][1];
    }
    
    void setTE(double n, double teta1){
        mE[0][0].real(cos(beta(teta(teta1, n))));
        mE[0][1].imag(-sin(beta(teta(teta1, n)))/p(n, teta(teta1, n)));
        mE[1][0].imag(-sin(beta(teta(teta1, n)))*p(n, teta(teta1, n)));
        mE[1][1].real(cos(beta(teta(teta1, n))));

    }
    
    void setTM(double n, double teta1){
        mM[0][0].real(cos(beta(teta(teta1, n))));
        mM[0][1].imag(-sin(beta(teta(teta1, n)))/q(n, teta(teta1, n)));
        mM[1][0].imag(-sin(beta(teta(teta1, n)))*q(n, teta(teta1, n)));
        mM[1][1].real(cos(beta(teta(teta1, n))));
        
    }
    
    matrix2(void){
        mE[0][0].real(1);
        mE[0][1].real(0);
        mE[1][0].real(0);
        mE[1][1].real(1);
        
        mE[0][0].imag(0);
        mE[0][1].imag(0);
        mE[1][0].imag(0);
        mE[1][1].imag(0);
        
        mM[0][0].real(1);
        mM[0][1].real(0);
        mM[1][0].real(0);
        mM[1][1].real(1);

        mM[0][0].imag(0);
        mM[0][1].imag(0);
        mM[1][0].imag(0);
        mM[1][1].imag(0);
    }
    
    void refTE(matrix2 mat, double n, double teta1);
    void refTM(matrix2 mat, double n, double teta1);
};

void matrix2::refTE(matrix2 mat, double n, double teta1){
    rTE = (mat.mE[0][0]*p(1, teta1)-mat.mE[1][1]*p(n, teta(teta1, n))+(mat.mE[0][1]*p(n, teta(teta1, n))*p(1, teta1)-mat.mE[1][0]))/
          (mat.mE[0][0]*p(1, teta1)+mat.mE[1][1]*p(n, teta(teta1, n))+(mat.mE[0][1]*p(n, teta(teta1, n))*p(1, teta1)+mat.mE[1][0]));
}

void matrix2::refTM(matrix2 mat, double n, double teta1){
    rTM = (mat.mM[0][0]*q(1, teta1)-mat.mM[1][1]*q(n, teta(teta1, n))+(mat.mM[0][1]*q(n, teta(teta1, n))*q(1, teta1)-mat.mM[1][0]))/
          (mat.mM[0][0]*q(1, teta1)+mat.mM[1][1]*q(n, teta(teta1, n))+(mat.mM[0][1]*q(n, teta(teta1, n))*q(1, teta1)+mat.mM[1][0]));
}

double mod_2(std::complex<double> num){
    return pow(num.real(), 2)+pow(num.imag(), 2);
}

#endif /* functions_h */
