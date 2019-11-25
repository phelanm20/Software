#include "software/new_geom/spline.h"


Spline::Spline(const std::vector<Point>& points, unsigned int order) : knots(points)
{
    makeSegments(points, order);
}

Spline::Spline(const std::initializer_list<Point>& points, unsigned int order) : knots(points)
{
    makeSegments(points, order);
}

Point Spline::valueAt(double val) const
{
    if (val < 0.0 || val > 1.0)
    {
        std::stringstream ss;
        ss << "Tried to evaluate spline at " << val
           << ", which is outside of domain of the spline: [0,1]";
        throw std::invalid_argument(ss.str());
    }

    Point retval;

    if (segments.empty())
    {
        retval = knots.front();
    }
    else
    {
        // Note: this could be more performant with binary search
        auto seg_it = std::find_if(segments.begin(), segments.end(),
                                   [&](const SplineSegment& sseg) {
                                       return (val >= sseg.start && val <= sseg.end);
                                   });

        if (seg_it == segments.end())
        {
            std::stringstream ss;
            ss << "Tried to evaluate spline at " << val
               << ", which was not in any interval of any segment";
            throw std::runtime_error(ss.str());
        }

        retval = Point(seg_it->x.valueAt(val), seg_it->y.valueAt(val));
    }

    return retval;
}

size_t Spline::size(void) const
{
    return knots.size();
}

const std::vector<Point> Spline::getKnots(void) const
{
    return knots;
}

const Point Spline::startPoint(void) const
{
    return knots.front();
}

const Point Spline::endPoint(void) const
{
    return knots.back();
}

/*
Construct a polynomial from coefficients
        * s.t. n = coeffs.size() == the degree of the polynomial
* and of the form coeffs[0]*x^(n-1)
* + coeffs[1]*x^(n-2) + ... + coeffs[n-1]
 *
*/

//Make a segment and add it to segment vector
/*
Polynomial poly_x = Polynomial(std::make_pair(input_start, points[i - 1].x()),
                               std::make_pair(input_end, points[i].x()));
Polynomial poly_y = Polynomial(std::make_pair(input_start, points[i - 1].y()),
                               std::make_pair(input_end, points[i].y()));
segments.push_back(SplineSegment(poly_x, poly_y, input_start, input_end));
*/

//void insert_element (size_type i, size_type j, const_reference t)

void Spline::initSegments(const std::vector<Point>& points, unsigned int order)
{
    using namespace boost::numeric::ublas;
    if (points.size() == 0 || points.size() == 1)
    {
        throw std::runtime_error("Cannot create spline with less than two points.");
    }
    else if (points.size() > 1)
    {
        for (size_t m = 1; m < points.size(); m++)
        {

            fillInterpolationMatrix(points[m], points[m+1], order);

            //Get start and end points to interpolate between
            startX = points[m].x();
            startY = points[m].y();
            endX = points[m + 1].x();
            endY = points[m + 1].y();

            //initialize empty matrix
            matrix<double> m(order + 1, order + 1);

            //initialize first row
            double exp = order;
            for (unsigned int j = 0; j < m.size2(); j++)
            {
                m(0, j) = pow(endX, exp);
                exp--;
            }

            //initialize second row
            double exp = order;
            for (unsigned int j = 0; j < m.size2(); j++)
            {
                m(1, j) = pow(startX, exp);
                exp--;
            }

            //initialize third to final rows
            double exp = order;
            for (unsigned int i = 2; i < m.size1(); i++)
            {
                for(unsigned int j = 0; j < m.size2(); j++)
                {
                    m(i, j) = m(i-1, j)/startX      pow(startX, exp);
                    exp--;
                }
            }


        }
    }
}

void fillInterpolationMatrix(Point start, Point end, unsigned int order){

}
