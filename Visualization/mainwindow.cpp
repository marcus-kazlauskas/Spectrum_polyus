#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "functions.h"

#include <QGraphicsView>
#include <QDoubleSpinBox>
#include <QPen>
#include <QGraphicsTextItem>
#include <iostream>

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    QGraphicsScene *graph = new QGraphicsScene;
    QPen pen(Qt::red, 2);
    QPen pen1(Qt::black, 2);
    QPen pen2(Qt::green, 2);
    QPen pen3(Qt::blue, 2);
    QString str;
    matrix2 High, Low, Matrix;

    double teta1 = pi/6;            // угол падения луча
//    double teta1 = pi/4;
    double ng = 1.457;              // показатель преломления подложки из кварца
    double n2 = 2.4;                // показатель преломления TiO2
//    double n2 = 2.3;
    double n3 = 1.457;              // показатель преломления кварца
//    double n3 = 1.39;
    double lambda = 632;            // длина волны луча
    double lambdaI = 400;           // длина волны начала спектра
    double lambdaF = 900;           // длина волны конца спектра
    double l = lambdaI;
    double TTM, TTE;
    int N = 8;                      // число пар слоёв в периодической структуре зеркала
//    int N = 15;

    double scaleX = 3;              // параметры осей координат
    double scaleY = 800;
    double gap = 10;

    ui->graphicsView->setRenderHint(QPainter::Antialiasing, true);
    QGraphicsTextItem *textItem = new QGraphicsTextItem;

    graph->addLine(lambda*scaleX, 0, lambda*scaleX, -scaleY, pen2);                         // отрисовка осей координат
    str = QString("%1").arg(lambda);
    textItem = graph->addText(str);
    textItem->setPos(lambda*scaleX-gap, 0.5*gap);

    for(int i = 0; i <= 10; i++){
        graph->addLine(lambdaI*scaleX, -i*scaleY*0.1, lambdaF*scaleX, -i*scaleY*0.1, pen1);
        str = QString("%1").arg(i*0.1);
        textItem = graph->addText(str);
        textItem->setPos(lambdaI*scaleX-3*gap, -i*scaleY*0.1-gap);

        graph->addLine((lambdaI+i*50)*scaleX, 0, (lambdaI+i*50)*scaleX, -scaleY, pen1);
        str = QString("%1").arg(lambdaI+i*50);
        textItem = graph->addText(str);
        textItem->setPos((lambdaI+i*50)*scaleX-gap, 0.5*gap);
    }
    //------------------------------        // получение коэф. пропускания для начала спектра
    High.setTM(n2, teta1, lambda, l);
    Low.setTM(n3, teta1, lambda, l);
    Matrix.mulM(High.mM);

    for(int i = 0; i < N; i++){
        Matrix.mulM(Low.mM);
        Matrix.mulM(High.mM);
    }

    Matrix.refTM(Matrix, ng, teta1);
    Matrix.RTM = mod_2(Matrix.rTM);

    TTM = 1-Matrix.RTM;
    //------------------------------
    High.setTE(n2, teta1, lambda, l);
    Low.setTE(n3, teta1, lambda, l);
    Matrix.mulM(High.mE);

    for(int i = 0; i < N; i++){
        Matrix.mulE(Low.mE);
        Matrix.mulE(High.mE);
    }

    Matrix.refTE(Matrix, ng, teta1);
    Matrix.RTE = mod_2(Matrix.rTE);

    TTE = 1-Matrix.RTE;
    //----------------
    l += 0.1;

    while(l <= lambdaF){                    // вычисление коэф. пропускания до конца спектра для обеих поляризаций с отрисовкой графика
        Matrix.mM[0][0].real(1);
        Matrix.mM[0][1].real(0);
        Matrix.mM[1][0].real(0);
        Matrix.mM[1][1].real(1);

        Matrix.mM[0][0].imag(0);
        Matrix.mM[0][1].imag(0);
        Matrix.mM[1][0].imag(0);
        Matrix.mM[1][1].imag(0);

        High.setTM(n2, teta1, lambda, l);
        Low.setTM(n3, teta1, lambda, l);
        Matrix.mulM(High.mM);

        for(int i = 0; i < N; i++){
            Matrix.mulM(Low.mM);
            Matrix.mulM(High.mM);
        }

        Matrix.refTM(Matrix, ng, teta1);
        Matrix.RTM = mod_2(Matrix.rTM);

        graph->addLine((l-0.1)*scaleX, -TTM*scaleY, l*scaleX, -(1-Matrix.RTM)*scaleY, pen);
        TTM = 1-Matrix.RTM;
        //---------------------
        Matrix.mE[0][0].real(1);
        Matrix.mE[0][1].real(0);
        Matrix.mE[1][0].real(0);
        Matrix.mE[1][1].real(1);

        Matrix.mE[0][0].imag(0);
        Matrix.mE[0][1].imag(0);
        Matrix.mE[1][0].imag(0);
        Matrix.mE[1][1].imag(0);

        High.setTE(n2, teta1, lambda, l);
        Low.setTE(n3, teta1, lambda, l);
        Matrix.mulE(High.mE);

        for(int i = 0; i < N; i++){
            Matrix.mulE(Low.mE);
            Matrix.mulE(High.mE);
        }

        Matrix.refTE(Matrix, ng, teta1);
        Matrix.RTE = mod_2(Matrix.rTE);

        graph->addLine((l-0.1)*scaleX, -TTE*scaleY, l*scaleX, -(1-Matrix.RTE)*scaleY, pen3);
        TTE = 1-Matrix.RTE;
        //----------------
        l += 0.1;
    }

//    cout << Matrix.RTM << endl;

    ui->graphicsView->setScene(graph);
}

MainWindow::~MainWindow()
{
    delete ui;
}
