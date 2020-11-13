//
//  main.cpp
//  Program_4
//
//  Created by Kozlov Mikhail on 10.11.16.
//  Copyright © 2016 Kozlov Mikhail. All rights reserved.
//

#include <iostream>

#include "functions.h"

using namespace std;

int main(int argc, const char * argv[]) {
//    A*(H*L)^N*H*G
    
    matrix2 High, Low, Matrix;      // матрица слоя с высоким показателем преломления (n2), с низким показателем преломления (n3), многослойника
    double teta1 = pi/6;            // угод падения луча
 //    double teta1 = pi/4;
    double ng = 1.457;              // показатель преломления подложки из кварца
//    double ng = 1.457;
    double n2 = 2.4;                // показатель преломления TiO2
//    double n2 = 2.3;
    double n3 = 1.457;              // показатель преломления кварца
//    double n3 = 1.39;
    int i = 0;
    int j = 0;
    
    High.setTE(n2, teta1);
    High.setTM(n2, teta1);
    Low.setTE(n3, teta1);
    Low.setTM(n3, teta1);
    
    Matrix.mulE(High.mE);
    Matrix.mulM(High.mM);
    
    Matrix.refTE(Matrix, ng, teta1);
    Matrix.RTE = mod_2(Matrix.rTE);
    Matrix.refTM(Matrix, ng, teta1);
    Matrix.RTM = mod_2(Matrix.rTM);

    // вывод коэф. пропускания для TE- и TM-поляризаций при N=0
    cout << 1-Matrix.RTE << " " << 1-Matrix.RTM << endl;
    
    while (1-Matrix.RTM > 0.002) {
        Matrix.mulM(Low.mM);
        Matrix.mulM(High.mM);
        Matrix.refTM(Matrix, ng, teta1);
        Matrix.RTM = mod_2(Matrix.rTM);
        j += 1;
    }
    
    for(i = 0; i < j; i++){
        Matrix.mulE(Low.mE);
        Matrix.mulE(High.mE);
        Matrix.refTE(Matrix, ng, teta1);
        Matrix.RTE = mod_2(Matrix.rTE);
    }

    // вывод коэф. пропускания и N для TE- и TM-поляризаций
    cout << 1-Matrix.RTE << " " << i << " " << 1-Matrix.RTM << " " << j << endl;
    
    return 0;
}
