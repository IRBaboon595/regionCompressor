#ifndef COMPRESSOR_H
#define COMPRESSOR_H

#include <stdlib.h>
#include <vector>
#include "structs.h"
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <QDebug>

using namespace std;
using namespace GroupFlight;

enum class TriangleParams{
    SideA = 0,
    SideB,
    SideC,
    AngleA,
    AngleB,
    AngleC
};

class compressor
{
public:
    compressor();
    ~compressor();

    void setInVals(std::vector<GroupFlight::Point> vals);
    void setOutVals(std::vector<GroupFlight::Point> vals);
    void setFileVals(std::vector<char> vals);
    void setRouteLen(int len);
    void setTriangleParameters(double sa, double sb, double sc, double aa, double ab, double ac);
    void setDeltaTh(double val);

    std::vector<GroupFlight::Point> inValues();
    std::vector<GroupFlight::Point> outValues();
    std::vector<char> fileValues();
    int routeLen();
    double deltaTh();
    double TriangleParameters(TriangleParams num);

    void addValue(GroupFlight::Point val);
    GroupFlight::Point getLastValue();

    std::vector<GroupFlight::Point> getRouteFromRAW(std::vector<char> vals);

    void parseAngle(int pointNum);
    void parseRoute();
    double calcAngle(GroupFlight::Point a, GroupFlight::Point b, GroupFlight::Point c);
    void clearDuplets();
    void getRouteFromFile();
    void filterOutVect();
    bool filterInVect();

    GroupFlight::Point findIntersectPoint(GroupFlight::Line line1, GroupFlight::Line line2);
    LineEquation calcCoefs(GroupFlight::Line line);
    void findOutIntercept();

    void clearAll();

private:
    std::vector<GroupFlight::Point> m_inRoute;
    std::vector<GroupFlight::Point> m_outRoute;
    double m_deltaTh = 0;

    ofstream m_fileLog;
    ifstream m_fileData;

    std::vector<char> m_fileVector;

    int m_routeLen = 0;
    double sideA = 0;
    double sideB = 0;
    double sideC = 0;
    double angleA = 0;
    double angleB = 0;
    double angleC = 0;

};

#endif // COMPRESSOR_H
