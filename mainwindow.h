#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QObject>
#include <QMainWindow>
#include <QWidget>
#include <QApplication>
#include <QDebug>
#include <iomanip>
#include <stdlib.h>
#include <utility>
#include <iostream>     // std::cout
#include <iterator>     // std::next
#include <list>         // std::list
#include <algorithm>    // std::for_each
#include "compressor.h"
#include "structs.h"
#include <QPainter>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QSizePolicy>
#include <QGroupBox>
#include <QRect>
#include <QScreen>
#include <QGuiApplication>
#include <QPushButton>
#include <QMouseEvent>
#include <QPointF>
#include <QLineEdit>
#include "RegionReduce.h"
#include "QTimer"
#include "utility/math/polygon.h"

#define BigNum 65535

class mainwindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit mainwindow(QWidget *parent = nullptr);
    ~mainwindow();

    void getArguments(int argc, char *argv[]);
    virtual void paintEvent(QPaintEvent *);
    virtual void mousePressEvent(QMouseEvent *event);
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void mouseReleaseEvent(QMouseEvent *event);

private:
    QPainter *painter;
    QWidget *canvas;
    QWidget *window;
    QGroupBox *paintArea;
    QGroupBox *controls;
    QPushButton *redraw;
    QPushButton *clear;
    QPushButton *setTh;
    QLineEdit *thEdit;
    QTimer t;
    compressor *c;

    int argc;
    char **argv;

    int pointX = 0;
    int pointY = 0;
    uint curPointIndex = BigNum;

    GroupFlight::Point  tempPoint;
    util::math::Point   tempTestPoint;

    std::vector<GroupFlight::Point> m_inRoute;
    std::vector<GroupFlight::Point> m_outRoute;
    util::math::Polygon inPoly;
    util::math::Polygon outPoly;
    double m_deltaTh = 0;
    bool toggleVal = false;

signals:
    void repaintSignal();

public slots:
    void redrawPoly();
    void clearCanvas();
    void setOffset();

};

#endif // MAINWINDOW_H
