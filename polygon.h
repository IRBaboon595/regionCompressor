#ifndef MATHPOLYGON_H
#define MATHPOLYGON_H

#include <QVector>
#include <vector>
#include "point.h"
#include "rect.h"
#include "line.h"

namespace util
{
    namespace math
    {
        class Polygon : public std::vector<Point>
        {
        public:
            Polygon() {}
            Polygon(int size) : std::vector<Point>(size) {}
            Polygon(const std::vector<Point>& vec) : std::vector<Point>(vec) {}
            Polygon(const QVector<Point>& vec) : std::vector<Point>(vec.begin(), vec.end()) {}
            Polygon(const QVector<QPointF>& polygon)
            {
                reserve(polygon.size());
                for (const QPointF& pos : polygon)
                {
                    push_back({pos.x(), pos.y()});
                }
            }

            Rect boundingRect() const;

            bool isRegular(Polygon poly) const;
            bool isContains(const Point &point) const;
            bool isIntersectRegion(const Line &line) const;

            Polygon adjusted(double deltaTh) const;

            enum class Status : uint8_t
            {
                Ok,
                ZeroInsideSharp,
                ZeroOutsideSharp
            };              

        private:
            void isEctLine(const Point &p1, const Point &p2, const Point &pos, int *winding) const;
            double calcAngle(const Point& a, const Point& b, const Point& c) const;
            Status parseAngle(int pointNum, double deltaTh, Polygon *m_inRoute, Polygon *m_outRoute) const;
        };  // class Polygon
    }   // namespace math
}   // namespace util


#endif // POLYGON_H
