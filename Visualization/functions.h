#ifndef functions_h
#define functions_h

#include <math.h>
#include <complex>

const double pi = 3.14159;

double beta(double teta, double lambda0, double lambda){    // элемент матрицы Абелеса для слоя с оптической длиной хода луча,
    return pi*cos(teta)*lambda0/(2*lambda);                 // равной четверти длины волны lambda0, на которую рассчитано зеркало
}

double p(double n, double teta){                    // для луча с TE-поляризацией (или s)
    return n*cos(teta);
}
double q(double n, double teta){                    // для луча с TM-поляризацией (или p)
    return cos(teta)/n;
}

double teta(double teta1, double n){                // угол падения луча в слое
    return acos(sqrt(1-pow(sin(teta1/n), 2)));
}

class matrix2{                                      // класс для характеристической матрицы Абелеса
public:
    std::complex<double> mE[2][2], mM[2][2];
    std::complex<double> rTE, rTM;                                          // коэф. отражения по амплитуде
    double RTE, RTM;                                                        // коэф. отражения по интенсивности
    
    void mulE(std::complex<double> m[2][2]){                                // умножение матрицы многослойника на матрицу следующего верхнего слоя для TE-поляризации
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
    
    void mulM(std::complex<double> m[2][2]){                                // умножение матрицы многослойника на матрицу следующего верхнего слоя для TM-поляризации
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
    
    void setTE(double n, double teta1, double lambda0, double lambda){      // задать матрицу для луча с TE-поляризацией и длиной волны lambda, в среде с преломлением n,
        mE[0][0].real(cos(beta(teta(teta1, n), lambda0, lambda)));          // под углом teta1 к нормали, падающего на слой, рассчитанный на длину волны labda0
        mE[0][0].imag(0);
        mE[0][1].real(0);
        mE[0][1].imag(-sin(beta(teta(teta1, n), lambda0, lambda))/p(n, teta(teta1, n)));
        mE[1][0].real(0);
        mE[1][0].imag(-sin(beta(teta(teta1, n), lambda0, lambda))*p(n, teta(teta1, n)));
        mE[1][1].real(cos(beta(teta(teta1, n), lambda0, lambda)));
        mE[1][1].imag(0);

    }
    
    void setTM(double n, double teta1, double lambda0, double lambda){      // задать матрицу для луча с TM-поляризацией и длиной волны lambda, в среде с преломлением n,
        mM[0][0].real(cos(beta(teta(teta1, n), lambda0, lambda)));          // под углом teta1 к нормали, падающего на слой, рассчитанный на длину волны lambda0
        mM[0][0].imag(0);
        mM[0][1].real(0);
        mM[0][1].imag(-sin(beta(teta(teta1, n), lambda0, lambda))/q(n, teta(teta1, n)));
        mM[1][0].real(0);
        mM[1][0].imag(-sin(beta(teta(teta1, n), lambda0, lambda))*q(n, teta(teta1, n)));
        mM[1][1].real(cos(beta(teta(teta1, n), lambda0, lambda)));
        mM[1][1].imag(0);
        
    }
    
    matrix2(void){                                                          // конструктор класса, задаёт единичную матрицу
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
    
    void refTE(matrix2 mat, double n, double teta1);        // коэф. отражения по амплитуде для TE-поляризации
    void refTM(matrix2 mat, double n, double teta1);        // коэф. отражения по амплитуде для TM-поляризации
};

void matrix2::refTE(matrix2 mat, double n, double teta1){
    rTE = (mat.mE[0][0]*p(1, teta1)-mat.mE[1][1]*p(n, teta(teta1, n))+(mat.mE[0][1]*p(n, teta(teta1, n))*p(1, teta1)-mat.mE[1][0]))/
          (mat.mE[0][0]*p(1, teta1)+mat.mE[1][1]*p(n, teta(teta1, n))+(mat.mE[0][1]*p(n, teta(teta1, n))*p(1, teta1)+mat.mE[1][0]));
}

void matrix2::refTM(matrix2 mat, double n, double teta1){
    rTM = (mat.mM[0][0]*q(1, teta1)-mat.mM[1][1]*q(n, teta(teta1, n))+(mat.mM[0][1]*q(n, teta(teta1, n))*q(1, teta1)-mat.mM[1][0]))/
          (mat.mM[0][0]*q(1, teta1)+mat.mM[1][1]*q(n, teta(teta1, n))+(mat.mM[0][1]*q(n, teta(teta1, n))*q(1, teta1)+mat.mM[1][0]));
}

double mod_2(std::complex<double> num){                     // квадрат модуля комплексного числа для подсчёта интенсивности
    return pow(num.real(), 2)+pow(num.imag(), 2);
}

#endif /* functions_h */
