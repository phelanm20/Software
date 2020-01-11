#include "software/new_geom/spline.h"

Spline::Spline(const std::vector<Point> &points, unsigned int order) : knots(points)
{
    initSegments(points);
}

Spline::Spline(const std::initializer_list<Point> &points, unsigned int order) : knots(points)
{
    initSegments(points);
}

/*
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
*/

int Spline::size(void) const
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

void Spline::initSegments(const std::vector<Point>& points, unsigned int order)
{
    if (points.size() < 2)
    {
        throw std::runtime_error("Cannot create spline with less than two points.");
    }
    else if (points.size() > 1)
    {
        //can be a vector of vectors
        //can be a matrix
        getCoeffs(points, order);

        for (int i = 1; i < points.size(); i++)
        {
            Point first = points[i-1];
            Point last = points[i];

            //HERE
            //you need a function that takes stuff and outputs a vector of coefficients
            //use the polynomial constructor that takes a list of coeffs as doubles
            Polynomial poly = Polynomial !!!;

            segments.push_back(SplineSegment(poly, first, last));
        }
    }
}

//hi this is my function
//Return type cannot be void
//use first/second/third/fourth degree polynomials
void getCoeffs(const std::vector<Point> &points, unsigned int order)
{
    if (points.size() < 2)
    {
        throw std::runtime_error("Cannot create spline with less than two points.");
    }
    else if (points.size() > 1)
    {
        //Initialize coeffs matrix with zeros






        using namespace std;
        using namespace Eigen;
        int numSpline = points.size() - 1;

        typedef Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime, Options> MatrixXd;

        //Initialize interpolation matrix with zeros
        int size = (order + 1)*(points.size() -1);
        MatrixXd splineMatrix = MatrixXd::Zero(8,8);

        //Fill half of interpolation matrix with first two continuity conditions
        //For each polynomial
        for(int i = 0; i < numSpline; i++)
        {
            splineMatrix(0, 4*i) = pow(POINTONE.x(), 3);
            splineMatrix(0, 4*i + 1)
        }
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

