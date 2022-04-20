#ifndef REGIONREDUCE_H
#define REGIONREDUCE_H

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

/*!
 *  \brief Данный файл содержит набор функций дял расчёта.
 *  Объявленные ниже функции работают с std:vector<GroupFlight::Point>.
 *  Выполняется расчёт полигона встроенного в искомый с учётом отстройки на заданное пользователем расстояние.
 *  Имеются дополнительные функции проверки и корректировки.
 *  Данный заголовочный файл будет дополняться новыми функциями для реализации корректной работы основного алгоритма.
 *  Для работы реализованного алгоритма необходим файл structs.h из библиотеки util.
*/
namespace RegionReduce
{

enum class Status : uint8_t
{
    Ok,
    ZeroInsideSharp,
    ZeroOutsideSharp
};

/*!
 *  \brief Функция работает со вектором исходных точек.
 *  Выполняется проверка на наличие пересечений в исходном полигоне.
 *  В ходе проверки не проверяется пересечение линии с ближайшими соседями по полигону.
 *  \param m_inRoute - Вектор искомых точек. Объект класса std::vector<GroupFlight::Point>.
 *  \return Булевая величина. True - пересечений в полигоне нет. False - персечения в полигоне есть.
*/
inline bool filterVect(std::vector<GroupFlight::Point> m_Route)
{
    if(m_Route.size() < 3)                                   // Не менее трёх точек - иначе возврат
        return false;
    std::vector<GroupFlight::Line> testVector;
    uint vectSize = m_Route.end() - m_Route.begin();
    uint distance = m_Route.end() - m_Route.begin();

    for(uint i = 0; i < distance; i++)                           // Выполняется заполнение тест-вектора линиями на базе точек исходного вектора
    {
        if(i == (distance - 1))
        {
            testVector.push_back(GroupFlight::Line(m_Route.at(i), m_Route.at(0)));
        }
        else
        {
            testVector.push_back(GroupFlight::Line(m_Route.at(i), m_Route.at(i + 1)));
        }
    }
    while(testVector.size() > 2)                                // Выполняется проверка на пересечение линии n с линиями n + 2 + i (где i - порядковый номер в векторе)
    {
        if(testVector.size() == vectSize)
            distance--;
        else
            distance = testVector.size();
        for(uint i = 2; i < distance; i++)
        {
            if(GroupFlight::isIntersects(testVector.at(0), testVector.at(i)))               // Если обнаружено пересечение, то выходной поток сигнализирует ошибку - возвращается False
            {
                std::cout << "Route invalid - set another. Got border intercetion" << std::endl;
                //m_inRoute.clear();
                return false;
            }
        }
        testVector.erase(testVector.begin());                   // Если в ходе проверки, не обнаржено пересечений ни с одной другой линией в массиве, то элемент удаляется из контейнера
    }
    return true;
}

/*!
 *  \brief Данная функция выполняет расчёт угла между биссектрисой угла текущей вершины и осью x.
 *  Выполняется расчёт угла на одну сторону треугольника, затем на вторую.
 *  Результатом является среднее арифметическое этих углов, плюс меньший угол.
 *  Важно! Данная функция не учитывает ситуацию в которой в массиве могут быть дублирующие точки.
 *  Пользователь должен исключить применение данной функции с массивом в котором есть дублирующие точки.
 *  \param a - вершина треугольника с индексом "-1" по отношению к искомой вершине. Обьект GroupFlight::Point.
 *  \param b - искомая вершина треуголника. Обьект GroupFlight::Point.
 *  \param c - вершина треугольника с индексом "+1" по отношению к искомой вершине. Обьект GroupFlight::Point.
 *  \return Угол между биссектрисой и осью Х. Число с плавающей точкой с двойной точностью (double).
*/
inline double calcAngle(const GroupFlight::Point a, const GroupFlight::Point b, const GroupFlight::Point c)
{
    double result = 0;
    double angle_1 = 0;
    double angle_2 = 0;

    if((a.x - b.x) == 0)                                    // Определяем лежит ли сторона треугольника на оси у
    {
        if((a.y - b.y) > 0)                                 // Где координата у при этом?
        {
            angle_1 = GroupFlight::kHalfPi;
        }
        else
        {
            angle_1 = -GroupFlight::kHalfPi;
        }
    }
    else if((a.y - b.y) == 0)                               // ...или на оси у
    {
        if((a.x - b.x) > 0)                                 // Где координата х при этом?
        {
            angle_1 = 0;
        }
        else
        {
            angle_1 = -GroupFlight::kPi;
        }
    }
    else if((a.x - b.x) > 0)      // Вычисляем угол стороны а относительно оси х с учётом квадранта (1 - 0 градусов, 2 - 90 градусов, 3 - 180 градусов, 4 - 270 градусов)
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

    if((c.x - b.x) == 0)            // То же самое для стороны b
    {
        if((c.y - b.y) > 0)
        {
            angle_2 = GroupFlight::kHalfPi;
        }
        else
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
    else if((c.x - b.x) > 0)     // Вычисляем угол стороны а относительно оси х с учётом квадранта (1 - 0 градусов, 2 - 90 градусов, 3 - 180 градусов, 4 - 270 градусов)
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

    result = angle_1 + ((angle_2 - angle_1) * GroupFlight::kHalf);  // Расчёт - среднее между двумя углами плюс сдвиг на меньший из углов сторон

    if(abs(angle_2 - angle_1) > GroupFlight::kPi)
    {
        result += GroupFlight::kPi;
    }

    return result;
}

/*!
 *  \brief Данная функция выполняет расчёт коэффициентов k и b для простого линейного уравнения.
 *  \param line - линия. Обьект GroupFlight::Line.
 *  \return Коэффициенты линейного уравнения. Обьект RegionReduce::Coef.
*/
inline GroupFlight::Coef calcCoefs(const GroupFlight::Line line)
{
    GroupFlight::Coef coef;

    coef.k = ((line.p2.y - line.p1.y) / (line.p2.x - line.p1.x));
    coef.b = line.p2.y - coef.k * line.p2.x;

    return coef;
}

/*! \brief Расчёт медианы. Данная функция выполняет расчёт угла между биссектрисой угла текущей вершины и осью x.
 *  Выполняется расчёт угла на одну сторону треугольника, затем на вторую.
 *  Результатом является среднее арифметическое этих углов, плюс меньший угол.
 *  Важно! Данная функция не учитывает ситуацию в которой в массиве могут быть дублирующие точки.
 *  Пользователь должен исключить применение данной функции с массивом в котором есть дублирующие точки.
 *  \param a - вершина треугольника с индексом "-1" по отношению к искомой вершине. Обьект GroupFlight::Point.
 *  \param b - искомая вершина треуголника. Обьект GroupFlight::Point.
 *  \param c - вершина треугольника с индексом "+1" по отношению к искомой вершине. Обьект GroupFlight::Point.
 *  \return Угол между медианой и осью Х. Число с плавающей точкой с двойной точностью (double).
*/
inline double calculAngle(const GroupFlight::Point a, const GroupFlight::Point b, const GroupFlight::Point c)
{
    double result = 0;
    GroupFlight::Line sideB(c, a);
    GroupFlight::Coef sideB_coef(calcCoefs(sideB));
    double x = abs((c.x - a.x) / 2);
    x += (c.x > a.x ? a.x : c.x);
    double y = x * sideB_coef.k + sideB_coef.b;
    GroupFlight::Point bisPoint(x, y);

    result = atan((y - b.y)/(x - b.x));

    return result;
}

/*!
 *  \brief Данная функция выполняет расчёт точки пересечения двух линий.
 *  Данную функцию можно использовать только с заведомо пересекающимися линиями.
 *  \param line1 - линия 1. Обьект GroupFlight::Line.
 *  \param line2 - линия 2. Обьект GroupFlight::Line.
 *  \return Точка пересечения. Обьект GroupFlight::Coef.
 */
inline GroupFlight::Point findIntersectPoint(GroupFlight::Line line1, GroupFlight::Line line2)
{
    GroupFlight::Coef coef_1(calcCoefs(line1).k, calcCoefs(line1).b);
    GroupFlight::Coef coef_2(calcCoefs(line2).k, calcCoefs(line2).b);
    GroupFlight::Point interceptPoint;

    interceptPoint.x = ((coef_2.b - coef_1.b) / (coef_1.k - coef_2.k));
    interceptPoint.y = coef_2.k * interceptPoint.x + coef_2.b;

    return interceptPoint;
}

/*!
 *  \brief Данная функция выполняет поиск пересекающихся линий в полигоне.
 *  В общем алгоритме она применяется для проверки выходного полигона на наличие пересекающихся сторон.
 *  На текущий момент выполняется проверка нет ли пересечений между линиями "через одну", чтобы убрать петли в полигоне.
 *  Петли в полигоне возникают, как артефакты расчёта особо острых внутренних углов.
 *  В случае обнаружения "петли", из вектора удаляются точки с индексом n, n + 1, n + 2, n + 3. И добавляются n, "обнаруженная точка персечения", n + 3.
 *  \param m_outRoute - вектор точек. Объект std::vector<GroupFlight::Point>.
 *  \return Проверенный вектор точек. Объект std::vector<GroupFlight::Point>.
 */
inline bool findOutIntercept(std::vector<GroupFlight::Point> *m_outRoute)
{
    int num = m_outRoute->size();
    GroupFlight::Point point1;
    GroupFlight::Point point2;
    GroupFlight::Point resultPoint;
    GroupFlight::Line line1;
    GroupFlight::Line line2;

    for(int i = 0; i < (num - 2); i++)
    {
        point1 = m_outRoute->at(i);                             // Проверяем пересекаются ли линии n и n + 2
        line1.p1 = m_outRoute->at(i);
        line1.p2 = m_outRoute->at(i + 1);
        line2.p1 = m_outRoute->at(i + 2);
        if(uint32_t(i + 3) == m_outRoute->size())
        {
            point2 = m_outRoute->at(0);
            line2.p2 = m_outRoute->at(0);
        }
        else
        {
            point2 = m_outRoute->at(i + 3);
            line2.p2 = m_outRoute->at(i + 3);
        }

        if(GroupFlight::isIntersects(line1, line2))                 // Если обнаружено пересечение...
        {
            resultPoint = findIntersectPoint(line1, line2);
            qDebug() << "message 1";
            m_outRoute->erase(m_outRoute->begin() + i);             // ...удаляем точки формирующие пересекающие кривые...
            m_outRoute->erase(m_outRoute->begin() + i);
            m_outRoute->erase(m_outRoute->begin() + i);
            m_outRoute->erase(m_outRoute->begin() + i);
            qDebug() << "message 2";

            m_outRoute->emplace(m_outRoute->begin() + i, point2);           // ...добавляем точки без персечения (см. бриф).
            m_outRoute->emplace(m_outRoute->begin() + i, resultPoint);
            m_outRoute->emplace(m_outRoute->begin() + i, point1);
            num--;
        }
    }
}

/*!
 *  \brief Основная функция.
 *  Выполняет расчёт точки для всписываемого полигона.
 *  Расчёт выполняется для вершины треугольника образованного точками входного полигона с индексами:
 *  n - 1, n, n + 1.
 *  \param pointNum - индекс текущей вершины во входном контейнере. Число int.
 *  \param deltaTh - желаемое сжатие расчётного полигона. Число с плавающей точкой с двойной точностью (double).
 *  \param *m_inRoute - вектор точек. Указатель на обьект std::vector<GroupFlight::Point>.
 *  \param *m_outRoute - вектор точек. Указатель на обьект std::vector<GroupFlight::Point>.
 *  \return Статус отработанной вершины. Обьект Enum-класса Status.
 */
inline Status parseAngle(uint pointNum, double deltaTh, std::vector<GroupFlight::Point> *m_inRoute, std::vector<GroupFlight::Point> *m_outRoute)
{
    GroupFlight::Point pointA;                                  // Вершины треугольника заданы точками, углами и рёбрами.
    GroupFlight::Point pointB(m_inRoute->at(pointNum));
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

    a_side.p1 = pointB;
    c_side.p2 = pointB;
    if(pointNum == 0)
    {
        pointA = m_inRoute->back();
    }
    else
    {
        pointA = m_inRoute->at(pointNum - 1);
    }
    b_side.p2 = pointA;
    c_side.p1 = pointA;

    if(pointNum == (m_inRoute->size() - 1))
    {
        pointC = m_inRoute->front();
    }
    else
    {
        pointC = m_inRoute->at(pointNum + 1);
    }
    a_side.p2 = pointC;
    b_side.p1 = pointC;

    sideA = GroupFlight::lineLength(a_side);
    sideB = GroupFlight::lineLength(b_side);
    sideC = GroupFlight::lineLength(c_side);
    angleA = acos((pow(sideB, 2) + pow(sideC, 2) - pow(sideA, 2))/(2 * sideB * sideC));
    angleB = acos((pow(sideA, 2) + pow(sideC, 2) - pow(sideB, 2))/(2 * sideA * sideC));
    angleC = GroupFlight::kPi - angleB - angleA;                   // Расчёт параметров треугольника закончен

    bisLen = deltaTh / (cos(GroupFlight::kHalfPi - GroupFlight::kHalf * angleB));       // Длина биссектрисы
    testLine = 0.1 * deltaTh;   // Длина тест-линии для определения типа угла внутренний/внешний

    bisAngle = calcAngle(pointA, pointB, pointC);   // Угол биссектрисы относительно оси Х

    candidante_1.x = testLine * cos(bisAngle) + pointB.x;   // Две тест-точки. Одна вне полигона, а другая внутри.
    candidante_1.y = testLine * sin(bisAngle) + pointB.y;

    candidante_2.x = testLine * cos(bisAngle - GroupFlight::kPi) + pointB.x;
    candidante_2.y = testLine * sin(bisAngle - GroupFlight::kPi) + pointB.y;

    pointState = GroupFlight::isRegionContainsPoint(*m_inRoute, candidante_1);  // Проверяем - находиться ли одна из точек внутри полигона.

    /******************************** Расчёт ***************************************/

    std::vector<GroupFlight::Point> testRoute(*m_inRoute);
    testRoute.erase(testRoute.begin() + pointNum);
    bool pointBState = GroupFlight::isRegionContainsPoint(testRoute, pointB);

    if(pointBState)
    {
        if(angleB < 0.5)
        {
            if(pointNum == 0)
            {
                return Status::ZeroInsideSharp;
            }
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
            candidante_2.x = bisLen * cos(bisAngle - GroupFlight::kPi) + pointB.x;
            candidante_2.y = bisLen * sin(bisAngle - GroupFlight::kPi) + pointB.y;

            pointInter_1.x = deltaTh * cos(bisAngle - GroupFlight::kPi) + pointB.x;
            pointInter_1.y = deltaTh * sin(bisAngle - GroupFlight::kPi) + pointB.y;

            testLine_1.p1 = pointB;
            testLine_1.p2 = pointInter_1;
            testCoef_1 = calcCoefs(testLine_1);

            orthCoef.k = - (1 / testCoef_1.k);
            orthCoef.b = pointInter_1.y - pointInter_1.x * orthCoef.k;

            testLine_2.p1 = m_outRoute->back();
            testLine_2.p2 = candidante_2;
            testCoef_2 = calcCoefs(testLine_2);

            pointInter_2.x = ((orthCoef.b - testCoef_2.b) / (testCoef_2.k - orthCoef.k));
            pointInter_2.y = orthCoef.k * pointInter_2.x + orthCoef.b;
            testLine_3.p1 = pointInter_1;
            testLine_3.p2 = pointInter_2;
            testCoef_3 = calcCoefs(testLine_3);

            pointInter_3.x = pointInter_1.x + (pointInter_1.x - pointInter_2.x);
            pointInter_3.y = testCoef_3.k * pointInter_3.x + testCoef_3.b;

            m_outRoute->push_back(pointInter_2);            // Рассчётная вершина обрезается на расстоянии deltaTh от искомой вершины. Получается две точки.
            m_outRoute->push_back(pointInter_3);
        }
        else
        {
            candidante_1.x = bisLen * cos(bisAngle) + pointB.x;
            candidante_1.y = bisLen * sin(bisAngle) + pointB.y;

            candidante_2.x = bisLen * cos(bisAngle - GroupFlight::kPi) + pointB.x;
            candidante_2.y = bisLen * sin(bisAngle - GroupFlight::kPi) + pointB.y;

            if(GroupFlight::isRegionContainsPoint(*m_inRoute, candidante_1))
            {
                resultPoint.x = candidante_1.x;
                resultPoint.y = candidante_1.y;
            }
            else
            {
                resultPoint.x = candidante_2.x;
                resultPoint.y = candidante_2.y;
            }

            m_outRoute->push_back(resultPoint);             // Один из кандидатов является расчётной вершиной.
        }
    }
    else
    {
        double shortLineLength = (sideA < sideC) ? sideA : sideC;
        //double longLineLength = (sideA > sideC) ? sideA : sideC;
        GroupFlight::Coef longLineCoef = calcCoefs((sideA > sideC) ? a_side : c_side);
        //GroupFlight::Coef shortLineCoef = calcCoefs((sideA < sideC) ? a_side : c_side);
        GroupFlight::Coef longOrthLineCoef;
        double triangleHeight = shortLineLength * sin(angleB);
        GroupFlight::Point heightIntersectPoint;
        //double heightDivC_Side = shortLineLength * cos(angleB);

        longOrthLineCoef.k = -(1 / longLineCoef.k);
        if(sideA < sideC)
        {
            longOrthLineCoef.b = pointC.y - pointC.x * longOrthLineCoef.k;
        }
        else
        {
            longOrthLineCoef.b = pointA.y - pointA.x * longOrthLineCoef.k;
        }

        heightIntersectPoint.x = (longOrthLineCoef.b - longLineCoef.b) / (longOrthLineCoef.k - longLineCoef.k);
        heightIntersectPoint.y = longOrthLineCoef.k * heightIntersectPoint.x + longOrthLineCoef.b;

        if(triangleHeight > (2 * deltaTh))
        {
            candidante_1.x = bisLen * cos(bisAngle) + pointB.x;
            candidante_1.y = bisLen * sin(bisAngle) + pointB.y;

            candidante_2.x = bisLen * cos(bisAngle - GroupFlight::kPi) + pointB.x;
            candidante_2.y = bisLen * sin(bisAngle - GroupFlight::kPi) + pointB.y;

            if(GroupFlight::isRegionContainsPoint(*m_inRoute, candidante_1))
            {
                resultPoint.x = candidante_1.x;
                resultPoint.y = candidante_1.y;
            }
            else
            {
                resultPoint.x = candidante_2.x;
                resultPoint.y = candidante_2.y;
            }

            m_outRoute->push_back(resultPoint);             // Один из кандидатов является расчётной вершиной.
        }
        else
        {
            if(angleB < 0.5)
            {
                return Status::ZeroOutsideSharp;
            }
            else
            {
                candidante_1.x = bisLen * cos(bisAngle) + pointB.x;
                candidante_1.y = bisLen * sin(bisAngle) + pointB.y;

                candidante_2.x = bisLen * cos(bisAngle - GroupFlight::kPi) + pointB.x;
                candidante_2.y = bisLen * sin(bisAngle - GroupFlight::kPi) + pointB.y;

                if(GroupFlight::isRegionContainsPoint(*m_inRoute, candidante_1))
                {
                    resultPoint.x = candidante_1.x;
                    resultPoint.y = candidante_1.y;
                }
                else
                {
                    resultPoint.x = candidante_2.x;
                    resultPoint.y = candidante_2.y;
                }

                m_outRoute->push_back(resultPoint);             // Один из кандидатов является расчётной вершиной.
            }
        }
    }
    return Status::Ok;
}

/*!
 *  \brief Основная функция.
 *  Данная функция осуществляет вызов всех других функций, методов и процессов, необходимых для формирования вписываемого полигона.
 *  \param inVals - вектор искомых точек. Объект класса std::vector<GroupFlight::Point>.
 *  \param deltaTh - желаемое сжатие расчётного полигона. Число с плавающей точкой с двойной точностью (double).
 *  \return Вектор расчётных величин. Объект класса std::vector<GroupFlight::Point>.
 */
inline std::vector<GroupFlight::Point> parseRoute(std::vector<GroupFlight::Point> inVals, double deltaTh)
{
    std::vector<GroupFlight::Point> outVals; // Вектор расчётных величин
    Status state;
    if(filterVect(inVals))   // Нет ли пересечений в исходном полигоне
    {
        for(quint32 i = 0; i < inVals.size(); i++)     // Отработка каждой вершины в исходном полигоне
        {
            state = parseAngle(i, deltaTh, &inVals, &outVals);
            if(state == Status::ZeroInsideSharp)
            {
                std::rotate(inVals.begin(), inVals.begin() + 1, inVals.end());
                i--;
            }
            else if(state == Status::ZeroOutsideSharp)
            {
                inVals.erase(inVals.begin() + i);
                outVals.clear();
                i = -1;
            }
        }
        if(!filterVect(outVals))  outVals.clear();
    }
    else
    {
        outVals.clear();
    }
    return outVals;
}

}

#endif // REGIONREDUCE_H
