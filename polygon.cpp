#include "polygon.h"
#include "utility/math/point.h"
#include "utility/math.h"

using namespace util::math;

util::math::Rect util::math::Polygon::boundingRect() const
{
    auto pd = begin();
    const auto pe = end();

    if (empty()) return Rect();

    double minx, maxx, miny, maxy;
    minx = maxx = pd->x;
    miny = maxy = pd->y;

    while (pd != pe) {
        if (pd->x < minx)
            minx = pd->x;
        else if (pd->x > maxx)
            maxx = pd->x;
        if (pd->y < miny)
            miny = pd->y;
        else if (pd->y > maxy)
            maxy = pd->y;
        ++pd;
    }
    return Rect(minx, miny, maxx - minx, maxy - miny);
}


/*!
 * \brief Проверка региона на правильнось
 * \return
 */
bool util::math::Polygon::isRegular(Polygon poly) const
{
    if (size() < 3) return false;                                  // Не менее трёх точек - иначе возврат

    std::vector<math::Line> testVector;
    uint32_t vectSize = poly.end() - poly.begin();
    uint32_t distance = poly.end() - poly.begin();

    for (uint32_t i = 0; i < distance; i++)                           // Выполняется заполнение тест-вектора линиями на базе точек исходного вектора
    {
        if (i == (distance - 1))
        {
            testVector.push_back(math::Line(poly.at(i), poly.at(0)));
        }
        else
        {
            testVector.push_back(math::Line(poly.at(i), poly.at(i + 1)));
        }
    }

    while (testVector.size() > 2)                                // Выполняется проверка на пересечение линии n с линиями n + 2 + i (где i - порядковый номер в векторе)
    {
        if (testVector.size() == vectSize)
            distance--;
        else
            distance = uint32_t(testVector.size());

        for (int i = 2; i < distance; i++)
        {
            if (Line::isIntersects(testVector.at(0), testVector.at(i)))               // Если обнаружено пересечение, то выходной поток сигнализирует ошибку - возвращается False
            {
                return false;
            }
        }
        testVector.erase(testVector.begin());                   // Если в ходе проверки, не обнаржено пересечений ни с одной другой линией в массиве, то элемент удаляется из контейнера
    }
    return true;
}

void util::math::Polygon::isLineInsideIntercept(Polygon *polygon) const
{
    qint32 num = qint32(polygon->size());
    std::vector<Line> testVector;
    Point point1;
    Point point2;
    Point resultPoint;
    Line line1;
    Line line2;

    for(quint32 i = 0; i < num - 1; i++)
    {
        if(i == 0)
        {
            for(quint32 i = 0; i < num; i++)                           // Выполняется заполнение тест-вектора линиями на базе точек исходного вектора
            {
                if(i == (num - 1))
                {
                    testVector.push_back(Line(polygon->at(i), polygon->at(0)));
                }
                else
                {
                    testVector.push_back(Line(polygon->at(i), polygon->at(i + 1)));
                }
            }
        }

        line1 = testVector.at(i);
        if(i != testVector.size() - 2)
        {
            line2 = testVector.at(i + 2);
        }
        else
        {
            line2 = testVector.at(0);
        }

        if(line1.isIntersects(line2))                 // Если обнаружено пересечение...
        {
            if(getIntersectPoint(line1, line2, &resultPoint))
            {
                polygon->emplace(polygon->begin() + i + 1, resultPoint);
                polygon->erase(polygon->begin() + i + 2);             // ...удаляем точки формирующие пересекающие кривые...
                if(i != polygon->size() - 2)
                {
                    polygon->erase(polygon->begin() + i + 2);
                }
                else
                {
                    polygon->erase(polygon->begin());
                }

                i = -1;
                num = qint32(polygon->size());
                testVector.clear();
            }
        }
    }
}

bool util::math::Polygon::getIntersectPoint(Line line1, Line line2, Point *point) const
{
   LineEquation coef_1;
   LineEquation coef_2;
   Point interceptPoint;

   if(((line1.p1.x - line1.p2.x) != 0) && ((line1.p1.y - line1.p2.y) != 0))
   {
       coef_1 = line1.equation();

       if(((line2.p1.x - line2.p2.x) != 0) && ((line2.p1.y - line2.p2.y) != 0))
       {
           coef_2 = line2.equation();

           interceptPoint.x = ((coef_2.b - coef_1.b) / (coef_1.k - coef_2.k));
           interceptPoint.y = coef_2.k * interceptPoint.x + coef_2.b;
       }
       else if((line2.p1.x - line2.p2.x) == 0)
       {
           interceptPoint.x = line2.p1.x;
           interceptPoint.y = coef_1.k * interceptPoint.x + coef_1.b;
       }
       else if((line2.p1.y - line2.p2.y) == 0)
       {
           interceptPoint.y = line2.p1.y;
           interceptPoint.x = (interceptPoint.y - coef_1.b) / coef_1.k;
       }
   }
   else if((line1.p1.x - line1.p2.x) == 0)
   {
       if(((line2.p1.x - line2.p2.x) != 0) && ((line2.p1.y - line2.p2.y) != 0))
       {
           coef_2 = line2.equation();

           interceptPoint.x = line1.p1.x;
           interceptPoint.y = coef_2.k * interceptPoint.x + coef_2.b;
       }
       else if((line2.p1.x - line2.p2.x) == 0)
       {
           return false;
       }
       else if((line2.p1.y - line2.p2.y) == 0)
       {
           interceptPoint.x = line1.p1.x;
           interceptPoint.y = line2.p1.y;
       }
   }
   else if((line1.p1.y - line1.p2.y) == 0)
   {
       if(((line2.p1.x - line2.p2.x) != 0) && ((line2.p1.y - line2.p2.y) != 0))
       {
           coef_2 = line2.equation();

           interceptPoint.y = line1.p1.y;
           interceptPoint.x = (interceptPoint.y - coef_2.b) / coef_2.k;
       }
       else if((line2.p1.x - line2.p2.x) == 0)
       {
           interceptPoint.x = line2.p1.x;
           interceptPoint.y = line1.p1.y;
       }
       else if((line2.p1.y - line2.p2.y) == 0)
       {
           return false;
       }
   }

   return true;
}

bool util::math::Polygon::isContains(const Point &point) const
{
    if (empty()) return false;

    int winding_number = 0;
    Point lastPt = front();
    Point lastStart = lastPt;

    for (size_t i = 1; i < size(); ++i)
    {
        const Point &e = at(i);
        isEctLine(lastPt, e, point, &winding_number);
        lastPt = e;
    }

    if (lastPt != lastStart)
        isEctLine(lastPt, lastStart, point, &winding_number);

    return winding_number != 0;
}

bool util::math::Polygon::isIntersectRegion(const Line &line) const
{
    for (size_t i = 1; i < size(); i++)
    {
        if (line.isIntersects(Line(at(i-1), at(i)))) return true;
    }
    return false;
}

std::vector<Line> util::math::Polygon::getNormalDistances(Polygon inPoly, Point testPoint) const
{
    std::vector<Line> resultLines;
    Point intersectPoint;
    LineEquation inLineCoef;
    LineEquation orthLineCoef;
    Line testLine;

    for(qint32 i = 0; i < inPoly.size(); i++)
    {
        if(i != (inPoly.size() - 1))
        {
            testLine = Line(inPoly.at(i), inPoly.at(i + 1));
        }
        else
        {
            testLine = Line(inPoly.at(i), inPoly.at(0));
        }

        if(((testLine.p1.x - testLine.p2.x) != 0) && ((testLine.p1.y - testLine.p2.y) != 0))
        {
            inLineCoef = testLine.equation();

            orthLineCoef.k = -1 / inLineCoef.k;
            orthLineCoef.b = testPoint.y - testPoint.x * orthLineCoef.k;

            intersectPoint.x = (orthLineCoef.b - inLineCoef.b) / (inLineCoef.k - orthLineCoef.k);
            intersectPoint.y = orthLineCoef.k * intersectPoint.x + orthLineCoef.b;
        }
        else if((testLine.p1.x - testLine.p2.x) == 0)
        {
            intersectPoint.x = testLine.p2.x;
            intersectPoint.y = testPoint.y;
        }
        else if((testLine.p1.y - testLine.p2.y) == 0)
        {
            intersectPoint.y = testLine.p2.y;
            intersectPoint.x = testPoint.x;
        }
        resultLines.push_back((Line(testPoint, intersectPoint)));
    }

    return resultLines;
}

void util::math::Polygon::filterNormalDistances(Polygon inPolygon, Polygon *outPolygon, double deltaTh) const
{
    std::vector<Line> testLine_1Vect;
    Polygon testPoly;
    Line testLine_1;
    Line testLine_2;
    Point testPoint;
    for(qint32 i = 0; i < qint32(outPolygon->size()); i++)
    {
        testPoint = outPolygon->at(i);
        testLine_1Vect = getNormalDistances(inPolygon, testPoint);
        double tempLen = 0;
        for(qint32 t = 0; t < testLine_1Vect.size(); t++)
        {
            testLine_1 = testLine_1Vect.at(t);
            tempLen = testLine_1.length();
            if(tempLen < deltaTh * 0.9)
            {
                testLine_2 = Line(inPolygon.at(t), inPolygon.at((t != (testLine_1Vect.size() - 1) ? (t + 1) : 0)));
                if(testLine_1.isIntersects(testLine_2))
                {
                    outPolygon->erase(outPolygon->begin() + i);
                    i--;
                    break;
                }
                else
                {
                    testLine_1.p2.x += testLine_1.p2.x * 0.001;
                    if(testLine_1.isIntersects(testLine_2))
                    {
                        outPolygon->erase(outPolygon->begin() + i);
                        i--;
                        break;
                    }
                    else
                    {
                        testLine_1.p2.x -= testLine_1.p2.x * 0.001;
                        if(testLine_1.isIntersects(testLine_2))
                        {
                            outPolygon->erase(outPolygon->begin() + i);
                            i--;
                            break;
                        }
                        else
                        {
                            testLine_1.p2.y += testLine_1.p2.y * 0.001;
                            if(testLine_1.isIntersects(testLine_2))
                            {
                                outPolygon->erase(outPolygon->begin() + i);
                                i--;
                                break;
                            }
                            else
                            {
                                testLine_1.p2.y -= testLine_1.p2.y * 0.001;
                                if(testLine_1.isIntersects(testLine_2))
                                {
                                    outPolygon->erase(outPolygon->begin() + i);
                                    i--;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void util::math::Polygon::filterLittlePolyLines(Polygon *outPolygon, double deltaTh) const
{
    Polygon testVals;
    Point testPoint;
    Line testLine;

    for(qint32 i = 0; i < qint32(outPolygon->size()); i++)
    {
        testVals = *outPolygon;
        testPoint = outPolygon->at(i);
        testVals.erase(testVals.begin() + i);
        if(!(testVals.isContains(testPoint)))
        {
            testLine.p1 = outPolygon->at(i);
            if(i == qint32(outPolygon->size() - 1))
            {
                testLine.p2 = outPolygon->at(0);
            }
            else
            {
                testLine.p2 = outPolygon->at(i + 1);
            }

            if(testLine.length() <= deltaTh * 0.5)
            {
                if(i != (outPolygon->size() - 1))
                {
                    outPolygon->erase(outPolygon->begin() + i + 1);
                }
                else
                {
                    outPolygon->erase(outPolygon->begin());
                }
                i--;
            }
        }
    }
}

/*!
 *  \brief Основная функция.
 *  Данная функция осуществляет вызов всех других функций, методов и процессов, необходимых для формирования вписываемого полигона.
 *  \param inVals - вектор искомых точек. Объект класса std::vector<Point>.
 *  \param deltaTh - желаемое сжатие расчётного полигона. Число с плавающей точкой с двойной точностью (double).
 *  \return Вектор расчётных величин. Объект класса std::vector<Point>.
*/
util::math::Polygon util::math::Polygon::adjusted(double deltaTh) const
{
    Polygon inVals = *this;
    Polygon outVals; // Вектор расчётных величин
    Polygon testVals; 
    Polygon testInVals = inVals; // Вектор расчётных величин
    Point testPoint;
    Line testLine;
    Line testLine_1;
    Line testLine_2;
    Status state;

    deltaTh *= -1;

    if(isRegular(inVals))   // Нет ли пересечений в исходном полигоне
    {
        for(quint32 i = 0; i < inVals.size(); i++)     // Отработка каждой вершины в исходном полигоне
        {
            if(size() < 3)   break;    
            
            state = parseAngle(i, deltaTh, &inVals, &outVals);
            
            if(state == Status::ZeroInsideSharp)
            {
                std::rotate(inVals.begin(), inVals.begin() + 1, inVals.end());
                i--;
            }
            else if(state == Status::ZeroOutsideSharp)
            {                
                outVals.clear();
                i = -1;
            }
        }

        if(deltaTh >= 0)
        {
            if(outVals.size() > 3)
            {
                filterLittlePolyLines(&outVals, deltaTh);                   // Данная функция проверяет нет ли во внутреннем полигоне особо малых отрезков
                isLineInsideIntercept(&outVals);                            // Данная функция ищете пересечения линий в полигоне "через одну"
                filterNormalDistances(inVals, &outVals, deltaTh);           // Данная функция проверяет не расположены ли какие-либо точки внутреннего полигона близко к сторонам внешнего полигона
            }

            for(size_t i = 0; i < outVals.size(); i++)      // Проверяем все ли точки внутреннего полигона находятся внутри внешнего полигона
            {
                if(!(testInVals.isContains(outVals.at(i))))
                {
                    outVals.erase(outVals.begin() + i);
                    i--;
                }
            }

            for(size_t i = 0; i < outVals.size(); i++)      // Проверяем не пересекаются ли линии внутреннего и внешнего полигонов
            {
                testLine_1.p1 = outVals.at(i);
                if(i == outVals.size() - 1)
                {
                    testLine_1.p2 = outVals.at(0);
                }
                else
                {
                    testLine_1.p2 = outVals.at(i + 1);
                }

                for(size_t t = 0; t < testInVals.size(); t++)
                {
                    testLine_2.p1 = testInVals.at(t);
                    if(t == testInVals.size() - 1)
                    {
                        testLine_2.p2 = testInVals.at(0);
                    }
                    else
                    {
                        testLine_2.p2 = testInVals.at(t + 1);
                    }

                    if(testLine_1.isIntersects(testLine_2))
                    {
                        outVals.clear();
                        break;
                    }
                }
            }

            if(!isRegular(outVals))  outVals.clear();          // Проверяем нет ли пересечений во внутреннем полигоне
        }
    }
    else
    {
        outVals.clear();
    }

    return outVals;
}

void util::math::Polygon::isEctLine(const Point &p1, const Point &p2, const Point &pos, int *winding) const
{
    double x1 = p1.x;
    double y1 = p1.y;
    double x2 = p2.x;
    double y2 = p2.y;
    double y = pos.y;
    int dir = 1;

    if (fuzzyCompare(y1, y2))
    {
        return;
    }
    else if (y2 < y1)
    {
        double x_tmp = x2; x2 = x1; x1 = x_tmp;
        double y_tmp = y2; y2 = y1; y1 = y_tmp;
        dir = -1;
    }
    if (y >= y1 && y < y2)
    {
        double x = x1 + ((x2 - x1) / (y2 - y1)) * (y - y1);
        if (x <= pos.x) (*winding) += dir;
    }
}

/*!
 *  \brief Данная функция выполняет расчёт угла между биссектрисой угла текущей вершины и осью x.
 *  Выполняется расчёт угла на одну сторону треугольника, затем на вторую.
 *  Результатом является среднее арифметическое этих углов, плюс меньший угол.
 *  Важно! Данная функция не учитывает ситуацию в которой в массиве могут быть дублирующие точки.
 *  Пользователь должен исключить применение данной функции с массивом в котором есть дублирующие точки.
 *  \param a - вершина треугольника с индексом "-1" по отношению к искомой вершине. Обьект Point.
 *  \param b - искомая вершина треуголника. Обьект Point.
 *  \param c - вершина треугольника с индексом "+1" по отношению к искомой вершине. Обьект Point.
 *  \return Угол между биссектрисой и осью Х. Число с плавающей точкой с двойной точностью (double).
*/
double util::math::Polygon::calcAngle(const Point &a, const Point &b, const Point &c) const
{
    double result = 0;
    double angle_1 = 0;
    double angle_2 = 0;

    if((a.x - b.x) == 0)                                    // Определяем лежит ли сторона треугольника на оси у
    {
        if((a.y - b.y) > 0)                                 // Где координата у при этом?
        {
            angle_1 = constants::kHalfPi<double>();
        }
        else
        {
            angle_1 = -constants::kHalfPi<double>();
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
            angle_1 = -constants::kPi<double>();
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
            angle_1 = (constants::kHalfPi<double>() - abs(atan((a.y - b.y) / (a.x - b.x)))) + constants::kThreeSecondPi<double>();
        }
    }
    else
    {
        if((a.y - b.y) > 0)
        {
            angle_1 = (constants::kHalfPi<double>() - abs(atan((a.y - b.y) / (a.x - b.x)))) + constants::kHalfPi<double>();
        }
        else
        {
            angle_1 = atan((a.y - b.y) / (a.x - b.x)) + constants::kPi<double>();
        }
    }

    if((c.x - b.x) == 0)            // То же самое для стороны b
    {
        if((c.y - b.y) > 0)
        {
            angle_2 = constants::kHalfPi<double>();
        }
        else
        {
            angle_2 = -constants::kHalfPi<double>();
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
            angle_2 = -constants::kPi<double>();
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
            angle_2 = (constants::kHalfPi<double>() - abs(atan((c.y - b.y) / (c.x - b.x)))) + constants::kThreeSecondPi<double>();
        }
    }
    else
    {
        if((c.y - b.y) > 0)
        {
            angle_2 = (constants::kHalfPi<double>() - abs(atan((c.y - b.y) / (c.x - b.x)))) + constants::kHalfPi<double>();
        }
        else
        {
            angle_2 = atan((c.y - b.y) / (c.x - b.x)) + constants::kPi<double>();
        }
    }
   
    result = angle_1 + ((angle_2 - angle_1) * 0.5);  // Расчёт - среднее между двумя углами плюс сдвиг на меньший из углов сторон

    if(abs(angle_2 - angle_1) > constants::kPi<double>())
    {
        result += constants::kPi<double>();
    }

    return result;
}

/*!
 *  \brief Функция проверяет находиться ли точка с некоторой областью внутри полигона.
 *  Полигон задан треугольников с вершинами в точках a, b и c.
 *  \param a - Вершина треугольника a. Объект Point.
 *  \param b - Вершина треугольника b. Объект Point.
 *  \param c - Вершина треугольника c. Объект Point.
 *  \param point - Точка для проверки. Объект Point.
 *  \return Булевая величина. True - точка point находиться внутри полигона region. False - точка за границей полигона.
*/
inline bool isRegionContainsPointsReg(Point a, Point b, Point c, Point point)
{
    uint truthCounter = 0;
    Polygon region = std::vector<Point>({a, b, c});

    if(region.isContains(point))   truthCounter++;

    for(uint i = 0; i < 4; i++)
    {
        if(i == 0)
        {
            point.x += point.x * 0.001;
        }
        else if(i == 1)
        {
            point.x -= point.x * 0.001;
        }
        else if(i == 2)
        {
            point.y += point.y * 0.001;
        }
        else if(i == 3)
        {
            point.y -= point.y * 0.001;
        }

        if(region.isContains(point))    truthCounter++;
    }

    if(truthCounter != 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

/*!
 *  \brief Перегруженная функция.
 *  \param region - Исследуемый полигон (регион). Объект класса std::vector<Point>.
 *  \param point - Точка для проверки. Объект Point.
 *  \return Булевая величина. True - точка point находиться внутри полигона region. False - точка за границей полигона.
*/
inline bool isRegionContainsPointsReg(Polygon region, Point point)
{
    uint truthCounter = 0;

    if(region.isContains(point))    truthCounter++;

    for(uint i = 0; i < 4; i++)
    {
        if(i == 0)
        {
            point.x += point.x * 0.001;
        }
        else if(i == 1)
        {
            point.x -= point.x * 0.001;
        }
        else if(i == 2)
        {
            point.y += point.y * 0.001;
        }
        else if(i == 3)
        {
            point.y -= point.y * 0.001;
        }

        if(region.isContains(point))    truthCounter++;
    }

    if(truthCounter != 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}


/*!
 *  \brief Основная функция.
 *  Выполняет расчёт точки для всписываемого полигона.
 *  Расчёт выполняется для вершины треугольника образованного точками входного полигона с индексами:
 *  n - 1, n, n + 1.
 *  \param pointNum - индекс текущей вершины во входном контейнере. Число int.
 *  \param deltaTh - желаемое сжатие расчётного полигона. Число с плавающей точкой с двойной точностью (double).
 *  \param *m_inRoute - вектор точек. Указатель на обьект std::vector<Point>.
 *  \param *m_outRoute - вектор точек. Указатель на обьект std::vector<Point>.
 *  \return Статус отработанной вершины. Обьект Enum-класса Status.
*/
Polygon::Status util::math::Polygon::parseAngle(int pointNum, double deltaTh, Polygon *m_inRoute, Polygon *m_outRoute) const
{
    Point pointA;                                  // Вершины треугольника заданы точками, углами и рёбрами.
    Point pointB(m_inRoute->at(pointNum));
    Point pointC;
    Point resultPoint;
    Point candidante_1;
    Point candidante_2;
    Line a_side;
    Line b_side;
    Line c_side;
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

    if(pointNum == qint32(m_inRoute->size() - 1))
    {
        pointC = m_inRoute->front();
    }
    else
    {
        pointC = m_inRoute->at(pointNum + 1);
    }
    a_side.p2 = pointC;
    b_side.p1 = pointC;

    sideA = a_side.length();
    sideB = b_side.length();
    sideC = c_side.length();

    angleA = ((pow(sideB, 2) + pow(sideC, 2) - pow(sideA, 2))/(2 * sideB * sideC));

    if(angleA > 1)
    {
        angleA = 0;
    }
    else if(angleA < -1)
    {
        angleA = constants::kPi<double>();
    }
    else
    {
        angleA = acos(angleA);
    }

    angleB = ((pow(sideA, 2) + pow(sideC, 2) - pow(sideB, 2))/(2 * sideA * sideC));

    if(angleB > 1)
    {
        angleB = 0;
    }
    else if(angleB < -1)
    {
        angleB = constants::kPi<double>();
    }
    else
    {
        angleB = acos(angleB);
    }

    angleC = constants::kPi<double>() - angleB - angleA;                   // Расчёт параметров треугольника закончен

    bisLen = abs(deltaTh) / (cos(constants::kHalfPi<double>() - constants::kHalf<double>() * angleB));       // Длина биссектрисы
    testLine = 0.1 * abs(deltaTh);   // Длина тест-линии для определения типа угла внутренний/внешний

    bisAngle = calcAngle(pointA, pointB, pointC);   // Угол биссектрисы относительно оси Х

    candidante_1.x = testLine * cos(bisAngle) + pointB.x;   // Две тест-точки. Одна вне полигона, а другая внутри.
    candidante_1.y = testLine * sin(bisAngle) + pointB.y;

    candidante_2.x = testLine * cos(bisAngle - constants::kPi<double>()) + pointB.x;
    candidante_2.y = testLine * sin(bisAngle - constants::kPi<double>()) + pointB.y;

    pointState = m_inRoute->isContains(candidante_1);  // Проверяем - находиться ли одна из точек внутри полигона.

    /******************************** Расчёт ***************************************/

    Polygon testRoute(*m_inRoute);
    testRoute.erase(testRoute.begin() + pointNum);
    bool pointBState = testRoute.isContains(pointB);

    pointBState = deltaTh >= 0 ? pointBState : !pointBState;

    if(pointBState)
    {
        candidante_1.x = bisLen * cos(bisAngle) + pointB.x;
        candidante_1.y = bisLen * sin(bisAngle) + pointB.y;

        candidante_2.x = bisLen * cos(bisAngle - constants::kPi<double>()) + pointB.x;
        candidante_2.y = bisLen * sin(bisAngle - constants::kPi<double>()) + pointB.y;

        if(angleB < 0.5)
        {
            if(pointNum == 0)
            {
                return Status::ZeroInsideSharp;
            }

            if(deltaTh < 0) bisAngle += constants::kPi<double>();

            Point pointInter_1;
            Point pointInter_2;
            Point pointInter_3;
            LineEquation testCoef_1;
            LineEquation testCoef_2;
            LineEquation testCoef_3;
            LineEquation orthCoef;
            Line testLine_1;
            Line testLine_2;
            Line testLine_3;

            pointInter_1.x = deltaTh * cos(bisAngle - constants::kPi<double>()) + pointB.x;
            pointInter_1.y = deltaTh * sin(bisAngle - constants::kPi<double>()) + pointB.y;

            testLine_1.p1 = pointB;
            testLine_1.p2 = pointInter_1;

            testLine_2.p1 = m_outRoute->back();
            testLine_2.p2 = candidante_2;

            if(((testLine_1.p1.x - testLine_1.p2.x) != 0) && ((testLine_1.p1.y - testLine_1.p2.y) != 0))
            {
                testCoef_1 = testLine_1.equation();
                orthCoef.k = - (1 / testCoef_1.k);
                orthCoef.b = pointInter_1.y - pointInter_1.x * orthCoef.k;

                if(((testLine_2.p1.x - testLine_2.p2.x) != 0) && ((testLine_2.p1.y - testLine_2.p2.y) != 0))
                {
                    testCoef_2 = testLine_2.equation();

                    pointInter_2.x = ((orthCoef.b - testCoef_2.b) / (testCoef_2.k - orthCoef.k));
                    pointInter_2.y = orthCoef.k * pointInter_2.x + orthCoef.b;
                }
                else if((testLine_2.p1.x - testLine_2.p2.x) == 0)
                {
                    pointInter_2.x = testLine_2.p1.x;
                    pointInter_2.y = orthCoef.k * pointInter_2.x + orthCoef.b;
                }
                else if((testLine_2.p1.y - testLine_2.p2.y) == 0)
                {
                    testCoef_2 = testLine_2.equation();

                    pointInter_2.y = testLine_2.p1.y;
                    pointInter_2.x = (pointInter_2.y - orthCoef.b) / orthCoef.k;
                }

                testLine_3.p1 = pointInter_1;
                testLine_3.p2 = pointInter_2;
                testCoef_3 = testLine_3.equation();

                pointInter_3.x = pointInter_1.x + (pointInter_1.x - pointInter_2.x);
                pointInter_3.y = testCoef_3.k * pointInter_3.x + testCoef_3.b;
            }
            else if((testLine_1.p1.x - testLine_1.p2.x) == 0)
            {
                testCoef_2 = testLine_2.equation();

                pointInter_2.y = pointInter_1.y;
                pointInter_2.x = (pointInter_2.y - testCoef_2.b) / testCoef_2.k;

                pointInter_3.x = pointInter_1.x + (pointInter_1.x - pointInter_2.x);
                pointInter_3.y = pointInter_1.y;
            }
            else if((testLine_1.p1.y - testLine_1.p2.y) == 0)
            {
                testCoef_2 = testLine_2.equation();

                pointInter_2.x = pointInter_1.x;
                pointInter_2.y = testCoef_2.k * pointInter_2.x + testCoef_2.b;

                pointInter_3.x = pointInter_1.x;
                pointInter_3.y = pointInter_1.y + (pointInter_1.y - pointInter_2.y);
            }

            m_outRoute->push_back(pointInter_2);            // Рассчётная вершина обрезается на расстоянии deltaTh от искомой вершины. Получается две точки.
            m_outRoute->push_back(pointInter_3);
        }
        else
        {
            pointBState = m_inRoute->isContains(candidante_1);

            pointBState = deltaTh >= 0 ? pointBState : !pointBState;
            pointState = deltaTh >= 0 ? pointState : !pointState;

            if((m_inRoute->isContains(candidante_1)) && (pointState))
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
        LineEquation longLineCoef;
        LineEquation longOrthLineCoef;
        double triangleHeight = shortLineLength * sin(angleB);
        Point heightIntersectPoint;
        Line longLine;
        Line shortLine;
        Line nextLine;
        LineEquation nextLineCoef;

        if(deltaTh < 0) bisAngle += constants::kPi<double>();

        candidante_1.x = bisLen * cos(bisAngle) + pointB.x;
        candidante_1.y = bisLen * sin(bisAngle) + pointB.y;

        candidante_2.x = bisLen * cos(bisAngle - constants::kPi<double>()) + pointB.x;
        candidante_2.y = bisLen * sin(bisAngle - constants::kPi<double>()) + pointB.y;

        if(m_inRoute->size() > 3)
        {
            if(sideA > sideC)
            {
                nextLine.p1 = pointA;

                longLine = a_side;
                shortLine = c_side;
            }
            else
            {
                nextLine.p1 = pointC;

                longLine = c_side;
                shortLine = a_side;
            }

            if(((longLine.p1.x - longLine.p2.x) != 0) && ((longLine.p1.y - longLine.p2.y) != 0))
            {
                longLineCoef = longLine.equation();

                nextLineCoef.k = -1 / longLineCoef.k;
                nextLineCoef.b = nextLine.p1.y - nextLine.p1.x * nextLineCoef.k;

                heightIntersectPoint.x = (nextLineCoef.b - longLineCoef.b) / (longLineCoef.k - nextLineCoef.k);
                heightIntersectPoint.y = nextLineCoef.k * heightIntersectPoint.x + nextLineCoef.b;
            }
            else if((longLine.p1.x - longLine.p2.x) == 0)
            {
                heightIntersectPoint.x = longLine.p2.x;
                heightIntersectPoint.y = nextLine.p1.y;
            }
            else if((longLine.p1.y - longLine.p2.y) == 0)
            {
                heightIntersectPoint.y = longLine.p2.y;
                heightIntersectPoint.x = nextLine.p1.x;
            }
        }
        if(triangleHeight > (2 * abs(deltaTh)))
        {
            if((m_inRoute->isContains(candidante_1)) && (pointState))
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
            if(angleB < 1.04)
            {
                if(m_inRoute->size() > 3)
                {
                    if(deltaTh >= 0)
                    {
                        if(sideB >= (1 * deltaTh))
                        {
                            m_inRoute->emplace(m_inRoute->begin() + pointNum, heightIntersectPoint);
                            m_inRoute->erase(m_inRoute->begin() + pointNum + 1);
                        }
                        else
                        {
                            m_inRoute->erase(m_inRoute->begin() + pointNum);
                        }
                    }
                    else
                    {
                        m_inRoute->erase(m_inRoute->begin() + pointNum);
                    }
                }
                else
                {
                    m_inRoute->erase(m_inRoute->begin() + pointNum);
                }

                return Status::ZeroOutsideSharp;
            }
            else
            {
                if((m_inRoute->isContains(candidante_1)) && (pointState))
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
