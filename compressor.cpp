#include "compressor.h"


compressor::compressor()
{
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    m_fileLog.open("Log001.txt");
    m_fileLog << "Log File Created" << " " << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
    m_fileLog.close();
}

void compressor::setInVals(std::vector<GroupFlight::Point> vals)
{
    this->m_inRoute = vals;
}

void compressor::setOutVals(std::vector<GroupFlight::Point> vals)
{
    this->m_outRoute = vals;
}

void compressor::setFileVals(std::vector<char> vals)
{
    this->m_fileVector = vals;
}

void compressor::setRouteLen(int len)
{
    this->m_routeLen = len;
}

void compressor::setTriangleParameters(double sa, double sb, double sc, double aa, double ab, double ac)
{
    this->sideA = sa;
    this->sideB = sb;
    this->sideC = sc;
    this->angleA = aa;
    this->angleB = ab;
    this->angleC = ac;
}

void compressor::setDeltaTh(double val)
{
    this->m_deltaTh = val;
}

std::vector<GroupFlight::Point> compressor::inValues()
{
    return this->m_inRoute;
}

std::vector<GroupFlight::Point> compressor::outValues()
{
    return this->m_outRoute;
}

std::vector<char> compressor::fileValues()
{
    return this->m_fileVector;
}

int compressor::routeLen()
{
    return this->m_routeLen;
}

double compressor::deltaTh()
{
    return this->m_deltaTh;
}

double compressor::TriangleParameters(TriangleParams num)
{
    switch(num)
    {
        case TriangleParams::SideA:
            return this->sideA;
        break;
        case TriangleParams::SideB:
            return this->sideB;
        break;
        case TriangleParams::SideC:
            return this->sideC;
        break;
        case TriangleParams::AngleA:
            return this->angleA;
        break;
        case TriangleParams::AngleB:
            return this->angleB;
        break;
        case TriangleParams::AngleC:
            return this->angleC;
        break;
        default:
        return 0;
        break;
    }
}

void compressor::addValue(GroupFlight::Point val)
{
    this->m_inRoute.push_back(val);
}

GroupFlight::Point compressor::getLastValue()
{
    return this->m_inRoute.back();
}

std::vector<GroupFlight::Point> compressor::getRouteFromRAW(std::vector<char> vals)
{
    int i = 0;
    double x = 0;
    double y = 0;
    double *coord = &x;
    GroupFlight::Point tempPoint;
    string tempString;
    tempString.clear();
    std::vector<GroupFlight::Point> route;

    for(int t = 0; t < m_routeLen; t++)
    {
        while(vals.at(i) != '/')
        {
            if(coord == &x)
            {
                if(m_fileVector.at(i) != ';')
                {
                    tempString.push_back(vals.at(i));
                    i++;
                }
                else
                {
                    x = stod(tempString);
                    tempString.clear();
                    coord = &y;
                    i++;
                }
            }
            else if(coord == &y)
            {
                tempString.push_back(vals.at(i));
                i++;

            }
        }
        y = stod(tempString);
        tempString.clear();
        coord = &x;
        i++;
        tempPoint.x = x;
        tempPoint.y = y;
        route.push_back(tempPoint);
    }
    return route;
}

void compressor::parseAngle(int pointNum)
{
    GroupFlight::Point pointA;
    GroupFlight::Point pointB(m_inRoute.at(pointNum));
    GroupFlight::Point pointC;
    GroupFlight::Point resultPoint;
    GroupFlight::Point candidante_1;
    GroupFlight::Point candidante_2;
    GroupFlight::Line a_side;
    GroupFlight::Line b_side;
    GroupFlight::Line c_side;
    double sideA = 0;
    double sideB = 0;
    double sideC = 0;
    double angleA = 0;
    double angleB = 0;
    double angleC = 0;
    double bisLen = 0;
    double bisAngle = 0;
    bool pointState = false;
    double testLine = 0;

    a_side.p2 = pointB;
    b_side.p1 = pointB;

    if(pointNum == 0)
    {
        pointA = m_inRoute.at(m_inRoute.size() - 1);
    }
    else
    {
        pointA = m_inRoute.at(pointNum - 1);
    }
    a_side.p1 = pointA;
    c_side.p2 = pointA;


    if(uint32_t(pointNum) == (m_inRoute.size() - 1))
    {

        pointC = m_inRoute.at(0);
    }
    else
    {
        pointC = m_inRoute.at(pointNum + 1);
    }
    b_side.p2 = pointC;
    c_side.p1 = pointC;

    sideA = GroupFlight::lineLength(a_side);
    sideB = GroupFlight::lineLength(b_side);
    sideC = GroupFlight::lineLength(c_side);
    angleA = acos((pow(sideB, 2) + pow(sideC, 2) - pow(sideA, 2))/(2 * sideB * sideC));
    angleB = acos((pow(sideA, 2) + pow(sideC, 2) - pow(sideB, 2))/(2 * sideA * sideC));
    angleC = GroupFlight::kPi - angleB - angleA;

    bisLen = m_deltaTh / (cos(GroupFlight::kHalfPi - GroupFlight::kHalf * angleC));
    testLine = 0.1 * m_deltaTh;

    bisAngle = calcAngle(pointA, pointB, pointC);

    candidante_1.x = testLine * cos(bisAngle) + pointB.x;
    candidante_1.y = testLine * sin(bisAngle) + pointB.y;

    candidante_2.x = testLine * cos(bisAngle - GroupFlight::kPi) + pointB.x;
    candidante_2.y = testLine * sin(bisAngle - GroupFlight::kPi) + pointB.y;

    pointState = GroupFlight::isRegionContainsPoint(m_inRoute, candidante_1);

    if(cos(GroupFlight::kHalfPi - GroupFlight::kHalf * angleC) < 0.2)
    {
        double tempLen = 1;

        if(pointState)
        {
            while(cos(GroupFlight::kHalfPi - GroupFlight::kHalf * angleC) < 0.2)
            {
                pointB.x += tempLen * cos(bisAngle);
                pointB.y += tempLen * sin(bisAngle);

                a_side.p2 = pointB;
                b_side.p1 = pointB;

                sideA = GroupFlight::lineLength(a_side);
                sideB = GroupFlight::lineLength(b_side);

                angleA = acos((pow(sideB, 2) + pow(sideC, 2) - pow(sideA, 2))/(2 * sideB * sideC));
                angleB = acos((pow(sideA, 2) + pow(sideC, 2) - pow(sideB, 2))/(2 * sideA * sideC));
                angleC = GroupFlight::kPi - angleB - angleA;
            }
            m_inRoute.erase(m_inRoute.begin() + pointNum);
            if(pointNum != 0)
                m_outRoute.erase(m_outRoute.end() - 1);
            if(angleC > 0.5)
            {
                m_inRoute.emplace(m_inRoute.begin() + pointNum, pointB);
                parseAngle(pointNum - 1);
                bisLen = m_deltaTh / (cos(GroupFlight::kHalfPi - GroupFlight::kHalf * angleC));
                bisAngle = calcAngle(pointA, pointB, pointC);

                candidante_1.x = bisLen * cos(bisAngle) + pointB.x;
                candidante_1.y = bisLen * sin(bisAngle) + pointB.y;

                resultPoint.x = candidante_1.x;
                resultPoint.y = candidante_1.y;

                m_outRoute.push_back(resultPoint);
            }
            else
            {
                if(pointNum != 0)
                    parseAngle(pointNum - 1);
                parseAngle(pointNum);
            }
        }
        else
        {
            GroupFlight::Point pointInter_1;
            GroupFlight::Point pointInter_2;
            GroupFlight::Point pointInter_3;
            GroupFlight::Coef testCoef_1;
            GroupFlight::Coef testCoef_2;
            GroupFlight::Coef testCoef_3;
            GroupFlight::Coef orthCoef;
            GroupFlight::Line testLine_1;
            GroupFlight::Line testLine_2;
            GroupFlight::Line testLine_3;
            GroupFlight::Line orthLine;
            candidante_2.x = bisLen * cos(bisAngle - GroupFlight::kPi) + pointB.x;
            candidante_2.y = bisLen * sin(bisAngle - GroupFlight::kPi) + pointB.y;

            pointInter_1.x = m_deltaTh * cos(bisAngle - GroupFlight::kPi) + pointB.x;
            pointInter_1.y = m_deltaTh * sin(bisAngle - GroupFlight::kPi) + pointB.y;

            testLine_1.p1 = pointB;
            testLine_1.p2 = pointInter_1;

            testCoef_1 = calcCoefs(testLine_1);
            orthCoef.k = - (1 / testCoef_1.k);
            orthCoef.b = pointInter_1.y - pointInter_1.x * orthCoef.k;

            testLine_2.p1 = m_outRoute.at(m_outRoute.end() - m_outRoute.begin() - 1);
            testLine_2.p2 = candidante_2;
            testCoef_2 = calcCoefs(testLine_2);

            pointInter_2.x = ((orthCoef.b - testCoef_2.b) / (testCoef_2.k - orthCoef.k));
            pointInter_2.y = orthCoef.k * pointInter_2.x + orthCoef.b;
            testLine_3.p1 = pointInter_1;
            testLine_3.p2 = pointInter_2;
            testCoef_3 = calcCoefs(testLine_3);

            pointInter_3.x = pointInter_1.x + (pointInter_1.x - pointInter_2.x);
            pointInter_3.y = testCoef_3.k * pointInter_3.x + testCoef_3.b;

            m_outRoute.push_back(pointInter_2);
            m_outRoute.push_back(pointInter_3);
        }
    }
    else
    {
        candidante_1.x = bisLen * cos(bisAngle) + pointB.x;
        candidante_1.y = bisLen * sin(bisAngle) + pointB.y;

        candidante_2.x = bisLen * cos(bisAngle - GroupFlight::kPi) + pointB.x;
        candidante_2.y = bisLen * sin(bisAngle - GroupFlight::kPi) + pointB.y;

        if(pointState)
        {
            resultPoint.x = candidante_1.x;
            resultPoint.y = candidante_1.y;
        }
        else
        {
            resultPoint.x = candidante_2.x;
            resultPoint.y = candidante_2.y;
        }

        m_outRoute.push_back(resultPoint);
    }
}

void compressor::parseRoute()
{
    if(filterInVect())
    {
        for(uint32_t i = 0; i < m_inRoute.size(); i++)
        {
            parseAngle(i);
        }
    }
    if(m_outRoute.size() > 3)
        findOutIntercept();
}

double compressor::calcAngle(GroupFlight::Point a, GroupFlight::Point b, GroupFlight::Point c)
{
    double result = 0;
    double angle_1 = 0;
    double angle_2 = 0;

    if((a.x - b.x) == 0)
    {
        if((a.y - b.y) > 0)
        {
            angle_1 = GroupFlight::kHalfPi;
        }
        else            //Точки у которых совпадают х и у фильтрованы ранее
        {
            angle_1 = -GroupFlight::kHalfPi;
        }
    }
    else if((a.y - b.y) == 0)
    {
        if((a.x - b.x) > 0)
        {
            angle_1 = 0;
        }
        else
        {
            angle_1 = -GroupFlight::kPi;
        }
    }
    else if((a.x - b.x) > 0)
    {
        if((a.y - b.y) > 0)
        {
            angle_1 = atan((a.y - b.y) / (a.x - b.x));
        }
        else
        {
            angle_1 = (GroupFlight::kHalfPi - abs(atan((a.y - b.y) / (a.x - b.x)))) + GroupFlight::kThreeDivTwoPi;
        }
    }
    else
    {
        if((a.y - b.y) > 0)
        {
            angle_1 = (GroupFlight::kHalfPi - abs(atan((a.y - b.y) / (a.x - b.x)))) + GroupFlight::kHalfPi;
        }
        else
        {
            angle_1 = atan((a.y - b.y) / (a.x - b.x)) + GroupFlight::kPi;
        }
    }

    if((c.x - b.x) == 0)
    {
        if((c.y - b.y) > 0)
        {
            angle_2 = GroupFlight::kHalfPi;
        }
        else            //Точки у которых совпадают х и у фильтрованы ранее
        {
            angle_2 = -GroupFlight::kHalfPi;
        }
    }
    else if((c.y - b.y) == 0)
    {
        if((c.x - b.x) > 0)
        {
            angle_2 = 0;
        }
        else
        {
            angle_2 = -GroupFlight::kPi;
        }
    }
    else if((c.x - b.x) > 0)
    {
        if((c.y - b.y) > 0)
        {
            angle_2 = atan((c.y - b.y) / (c.x - b.x));
        }
        else
        {
            angle_2 = (GroupFlight::kHalfPi - abs(atan((c.y - b.y) / (c.x - b.x)))) + GroupFlight::kThreeDivTwoPi;
        }
    }
    else
    {
        if((c.y - b.y) > 0)
        {
            angle_2 = (GroupFlight::kHalfPi - abs(atan((c.y - b.y) / (c.x - b.x)))) + GroupFlight::kHalfPi;
        }
        else
        {
            angle_2 = atan((c.y - b.y) / (c.x - b.x)) + GroupFlight::kPi;
        }
    }


    result = angle_1 + ((angle_2 - angle_1) * GroupFlight::kHalf);
    return result;
}

void compressor::clearDuplets()
{
    uint32_t size = m_inRoute.size();
    for(uint32_t i = 0; i < (size - 1); i++)
    {
        if(m_inRoute.at(i) == m_inRoute.at(i + 1))
        {
            m_inRoute.erase(m_inRoute.begin() + i);
            size--;
            i--;
        }
    }
}

void compressor::getRouteFromFile()
{
    char buff[10000];
    int i = 0;
    string tempString;
    tempString.clear();

    m_fileData.open("Data001.txt");
    while(!m_fileData.eof())
    {
        m_fileData >> buff[i];
        m_fileVector.push_back(buff[i]);
        if(m_fileVector.at(i) == '/')
                m_routeLen++;
        i++;
    }

    while(m_fileVector.at(0) != '/')
    {
        tempString.push_back(m_fileVector.at(0));
        m_fileVector.erase(m_fileVector.begin());
    }
    m_fileVector.erase(m_fileVector.begin());
    m_deltaTh = stod(tempString);
    m_routeLen--;
    m_inRoute = getRouteFromRAW(m_fileVector);
    clearDuplets();

    std::cout << std::distance(m_inRoute.begin(), m_inRoute.end()) << std::endl;
    std::cout << m_inRoute.size() << std::endl;

    /*std::cout << "Enter delta value: " << std::endl;
    std::cin >> m_deltaTh;
    std::cout << m_deltaTh << std::endl;*/
}

GroupFlight::Point compressor::findIntersectPoint(GroupFlight::Line line1, GroupFlight::Line line2)
{
    GroupFlight::Coef coef_1(calcCoefs(line1).k, calcCoefs(line1).b);
    GroupFlight::Coef coef_2(calcCoefs(line2).k, calcCoefs(line2).b);
    GroupFlight::Point interceptPoint;

    interceptPoint.x = ((coef_2.b - coef_1.b) / (coef_1.k - coef_2.k));
    interceptPoint.y = coef_2.k * interceptPoint.x + coef_2.b;

    return interceptPoint;
}

GroupFlight::Coef compressor::calcCoefs(GroupFlight::Line line)
{
    GroupFlight::Coef coef;

    coef.k = ((line.p2.y - line.p1.y) / (line.p2.x - line.p1.x));
    coef.b = line.p2.y - coef.k * line.p2.x;

    return coef;
}

void compressor::findOutIntercept()
{
    int num = m_outRoute.end() - m_outRoute.begin();
    GroupFlight::Point point1;
    GroupFlight::Point point2;
    GroupFlight::Point resultPoint;
    GroupFlight::Line line1;
    GroupFlight::Line line2;

    for(int i = 0; i < (num - 2); i++)
    {
        point1 = m_outRoute.at(i);
        line1.p1 = m_outRoute.at(i);
        line1.p2 = m_outRoute.at(i + 1);
        line2.p1 = m_outRoute.at(i + 2);
        if(uint32_t(i + 3) == m_outRoute.size())
        {
            point2 = m_outRoute.at(0);
            line2.p2 = m_outRoute.at(0);
        }
        else
        {
            point2 = m_outRoute.at(i + 3);
            line2.p2 = m_outRoute.at(i + 3);
        }

        if(GroupFlight::isIntersects(line1, line2))
        {
            resultPoint = findIntersectPoint(line1, line2);

            m_outRoute.erase(m_outRoute.begin() + i);
            m_outRoute.erase(m_outRoute.begin() + i);
            m_outRoute.erase(m_outRoute.begin() + i);
            m_outRoute.erase(m_outRoute.begin() + i);

            m_outRoute.emplace(m_outRoute.begin() + i, point2);
            m_outRoute.emplace(m_outRoute.begin() + i, resultPoint);
            m_outRoute.emplace(m_outRoute.begin() + i, point1);
            num--;
        }
    }
}

bool compressor::filterInVect()
{
    if(m_inRoute.size() < 3)                                   // Не менее трёх точек - иначе возврат
        return false;
    std::vector<GroupFlight::Line> testVector;
    uint32_t vectSize = m_inRoute.end() - m_inRoute.begin();
    int distance = m_inRoute.end() - m_inRoute.begin();

    for(int i = 0; i < distance; i++)                           // Выполняется заполнение тест-вектора линиями на базе точек исходного вектора
    {
        if(i == (distance - 1))
        {
            testVector.push_back(GroupFlight::Line(m_inRoute.at(i), m_inRoute.at(0)));
        }
        else
        {
            testVector.push_back(GroupFlight::Line(m_inRoute.at(i), m_inRoute.at(i + 1)));
        }
    }
    while(testVector.size() > 2)                                // Выполняется проверка на пересечение линии n с линиями n + 2 + i (где i - порядковый номер в векторе)
    {
        if(testVector.size() == vectSize)
            distance--;
        else
            distance = testVector.size();
        for(int i = 2; i < distance; i++)
        {
            if(GroupFlight::isIntersects(testVector.at(0), testVector.at(i)))               // Если обнаружено пересечение, то выходной поток сигнализирует ошибку - возвращается False
            {
                std::cout << "Route invalid - set another. Got border intercetion" << std::endl;
                m_inRoute.clear();
                return false;
            }
        }
        testVector.erase(testVector.begin());                   // Если в ходе проверки, не обнаржено пересечений ни с одной другой линией в массиве, то элемент удаляется из контейнера
    }
    return true;
}

void compressor::filterOutVect()
{
    m_inRoute.clear();
    m_inRoute = m_outRoute;
    m_outRoute.clear();
    parseRoute();
}

void compressor::clearAll()
{
    m_inRoute.clear();
    m_outRoute.clear();
    m_fileVector.clear();
    m_deltaTh = 0;
    m_routeLen = 0;
    setTriangleParameters(0, 0, 0, 0, 0, 0);
}

compressor::~compressor()
{

}
