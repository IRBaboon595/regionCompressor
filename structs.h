#include <cstdint>
#include <vector>
#include <cmath>

//! Файл описвает структруры:
//! Point - точка с двумя координатами
//! Line - отрезок между двумя точками
//! Rect - прямоугольник
//! и алгоритмы с ними,
//! необходимые для работы модуля "GroupFlight"

namespace GroupFlight
{

#ifndef GF_STRUCTS_H
#define GF_STRUCTS_H

    constexpr double kDegPerSecond = 2.7777777777777777777777777777778e-04;           // градусов в секунде
    constexpr double kOneDiv360 = 2.77777777777777777777777777778e-03;                // 1/360
    constexpr double kTwoPi = 6.283185307179586476925286766559005768e+00;             // pi*2
    constexpr double kPi = 3.141592653589793238462643383279502884e+00;                // pi
    constexpr double kHalfPi = 1.570796326794896619231321691639751442e+00;            // pi/2
    constexpr double kQuarterPi = 7.8539816339744830961566084581988e-01;              // pi/4
    constexpr double kOneDivTwoPi = 1.5915494309189533576888376337251e-01;            // 1/(2*pi)
    constexpr double kThreeDivTwoPi = 4.7123889803846898576939650749192543e+00;       // (3*pi)/2
    constexpr double kHalf = 5.000000000000000000000000000000000000e-01;              // 1/2
    constexpr double kRadian = 5.729577951308232087679815481410517033e+01;            // градусов в радиане
    constexpr double kDegree = 1.745329251994329576923690768488612713e-02;            // радиан в градусе (2*pi/360)

    //! \brief Point - точка с двумя координатами
    struct Point
    {
        double x, y;

        Point(double _x, double _y): x(_x), y(_y){}
        Point(): x(0.), y(0.){}
    };

    //! \brief Line - отрезок между двумя точками
    struct Line
    {
        Point p1, p2;

        Line(double _x1, double _y1, double _x2, double _y2):
            p1(_x1, _y1), p2(_x2, _y2){}

        Line(Point _p1, Point _p2): p1(_p1), p2(_p2){}
        Line(): Line(0., 0., 0., 0.){}
    };

    //! \brief Coef - коэффициенты линии
    struct Coef
    {
        double k, b;

        Coef(double _k, double _b): k(_k), b(_b){}
        Coef(Coef const &c): k(c.k), b(c.b){}
        Coef(): k(0.), b(0.){}
    };

    //! \brief Rect - прямоугольник
    //! x, y - координаты верхнего левого угла
    //! w, h - длина и шинина, соответственно
    struct Rect
    {
        double x, y, w, h;

        Rect(double _x, double _y, double _w, double _h):
            x(_x), y(_y), w(_w), h(_h){}

        Rect(): Rect(0., 0., 0., 0.){}
    };

    template<typename T>
    inline constexpr T min(const T &val1, const T &val2)
    { return (val1 < val2) ? val1 : val2;}

    template<typename T>
    inline constexpr T max(const T &val1, const T &val2)
    { return (val1 > val2) ? val1 : val2;}

    template<typename T>
    inline constexpr T bound(const T &low, const T &val, const T &high)
    { return max(low, min(high, val)); }

    static inline bool compareDouble(double p1, double p2){ return (std::abs(p1 - p2) * 1000000000000. <= min(std::abs(p1), std::abs(p2))); }
    static inline bool isNullDouble(double d){ return std::abs(d) <= 0.000000000001; }

    static inline bool compareFloat(float p1, float p2){ return (std::abs(p1 - p2) * 100000.f <= min(std::abs(p1), std::abs(p2))); }
    static inline bool isNullFloat(float d){ return std::abs(d) <= 0.00001f; }

    inline bool operator==(const Point &p1, const Point &p2)
    {
        return ((!p1.x || !p2.x) ? isNullDouble(p1.x - p2.x) : compareDouble(p1.x, p2.x))
                && ((!p1.y || !p2.y) ? isNullDouble(p1.y - p2.y) : compareDouble(p1.y, p2.y));
    }
    inline bool operator!=(const Point &p1, const Point &p2){ return !(p1 == p2); }
    inline const Point operator-(const Point &p1, const Point &p2){ return Point(p1.x - p2.x, p1.y - p2.y); }
    inline const Point operator+(const Point &p1, const Point &p2){ return Point(p1.x + p2.x, p1.y + p2.y); }
    inline const Point operator*(const Point &p, double c){ return Point(p.x * c, p.y * c); }
    inline bool operator==(const Line &l1, const Line &l2){ return (l1.p1 == l2.p1 && l1.p2 == l2.p2); }

    inline Rect boundingRect(const std::vector<Point> &list)
    {
        auto pd = list.begin();
        const auto pe = list.end();

        if (list.empty()) return Rect(0, 0, 0, 0);

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
        return Rect(minx,miny, maxx - minx, maxy - miny);
    }

    inline void sincos(double value, double *s, double *c)
    {
        *s = sin(value);
        *c = cos(value);
    }

    inline void isPolygonEctLine(const Point &p1, const Point &p2, const Point &pos, int *winding)
    {
        double x1 = p1.x;
        double y1 = p1.y;
        double x2 = p2.x;
        double y2 = p2.y;
        double y = pos.y;
        int dir = 1;

        if (compareDouble(y1, y2))
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

    inline bool isRegionContainsPoint(const std::vector<Point> &region, const Point &point)
    {
        if (region.empty()) return false;

        int winding_number = 0;
        Point lastPt = region.front();
        Point lastStart = lastPt;

        for (size_t i = 1; i < region.size(); ++i)
        {
            const Point &e = region.at(i);
            isPolygonEctLine(lastPt, e, point, &winding_number);
            lastPt = e;
        }

        if (lastPt != lastStart)
            isPolygonEctLine(lastPt, lastStart, point, &winding_number);

        return winding_number != 0;
    }

    inline bool isIntersects(const Line &l1, const Line &l2, Point *intPos = nullptr)
    {
        const Point a = l1.p2 - l1.p1;
        const Point b = l2.p1 - l2.p2;
        const Point c = l1.p1 - l2.p1;
        const double denominator = a.y * b.x - a.x * b.y;

        if (denominator == 0 || !std::isfinite(denominator))
            return false;

        const double reciprocal = 1. / denominator;
        const double na = (b.y * c.x - b.x * c.y) * reciprocal;

        if (intPos) *intPos = l1.p1 + a * na;

        if (na < 0.) return false;
        if (na > 1.) return false;

        const double nb = (a.x * c.y - a.y * c.x) * reciprocal;
        if (nb < 0.) return false;
        if (nb > 1.) return false;

        return true;
    }

    inline bool isIntersectRegion(const std::vector<Point>& region, const Line &line)
    {
        for (size_t i = 1; i < region.size(); i++)
            if (isIntersects(line, Line(region.at(i-1), region.at(i))))
                return true;

        return false;
    }

    inline double lineLength(const Line &line)
    {
        double x = line.p2.x - line.p1.x;
        double y = line.p2.y - line.p1.y;
        return sqrt(x*x + y*y);
    }

    inline void setLineLength(Line &line, double length)
    {
        double x = line.p2.x - line.p1.x;
        double y = line.p2.y - line.p1.y;
        double len = sqrt(x*x + y*y);
        Line f(line.p1, Point(line.p1.x + x/len, line.p1.y + y/len));

        line.p2 = Point(line.p1.x + (f.p2.x - f.p1.x) * length,
                        line.p1.y + (f.p2.y - f.p1.y) * length);
    }

    inline double lineAngle(const Line &line)
    {
        const double dx = line.p2.x - line.p1.x;
        const double dy = line.p2.y - line.p1.y;
        const double theta = atan2(-dy, dx) * 360. / kTwoPi;
        const double theta_normalized = theta < 0. ? theta + 360. : theta;

        if (compareDouble(theta_normalized, 360.)) return 0.;
        else return theta_normalized;
    }

    inline void setLineAngle(Line &line, double angle)
    {
        const double angleR = angle * kTwoPi / 360.;
        const double l = lineLength(line);
        const double dx = cos(angleR) * l;
        const double dy = -sin(angleR) * l;
        line.p2.x = line.p1.x + dx;
        line.p2.y = line.p1.y + dy;
    }

    inline Point pointAtLine(const Line &line, double t)
    {
        return Point(line.p1.x + (line.p2.x - line.p1.x) * t, line.p1.y + (line.p2.y - line.p1.y) * t);
    }

    inline void adjustRect(Rect &rect, double xp1, double yp1, double xp2, double yp2)
    {
        rect.x += xp1;
        rect.y += yp1;
        rect.w += xp2 - xp1;
        rect.h += yp2 - yp1;
    }

#endif // GF_STRUCTS_H

} // namespace GroupFlight
