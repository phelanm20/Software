#pragma once
#include "thunderbots/bazel-thunderbots/external/boost/boost/numeric/odeint/external/eigen/eigen.hpp"
#include <math.h>
#include "software/new_geom/polynomial.h"
#include "software/new_geom/point.h"
#include "software/new_geom/vector.h"


/**
 * A Spline is a parameterized piecewise function
 * composed of Polynomials of set order between Points
 */
class Spline {
public:
    Spline() = delete;

    /**
     * Construct a spline by drawing Polynomials of set order between
     * consecutive Points
     *
     * @throws std::invalid_argument if points.size() < 2
     *
     * @param points Points Polynomials are drawn between
     */
    explicit Spline(const std::vector<Point> &points, unsigned int order);

    /**
     * Construct a spline by drawing Polynomials of set order between
     * consecutive Points
     *
     * @throws std::invalid_argument if points.size() < 2
     *
     * @param points Points Polynomials are drawn between
     */
    Spline(const std::initializer_list<Point> &points, unsigned int order);

    /**
     * Calculates the value of spline evaluated at value val
     *
     * @param val value to evaluate spline
     *
     * @throws std::invalid_argument if val is outside of domain
     *
     * @return value of spline evaluated at value val
     * if not defined by a spline then return closest start/end point
     */
    Point valueAt(double val) const;

    /**
     * Gets the number of knots in the spline including start and end points
     *
     * @return size of the spline
     */
    int size(void) const;

    /**
     * Gets knots in the spline including start and end points
     *
     * @return knots in the spline
     */
    const std::vector<Point> getKnots(void) const;

    /**
     * Gets start point of spline
     *
     * @return start point of spline
     */
    const Point startPoint(void) const;

    /**
     * Gets end point of spline
     *
     * @return end point of spline
     */
    const Point endPoint(void) const;

private:
    class SplineSegment{
    public:
        SplineSegment(Polynomial polynomial, Point start, Point end): polynomial(polynomial), start(start), end(end)
        {
        }

        const Polynomial polynomial;
        const Point start;
        const Point end;
    };

    // segments represent the polynomials that interpolate between points
    std::vector<SplineSegment> segments;

    // points that connect segments
    const std::vector<Point> knots;

    /**
     * Initialize (points.size() - 1) segments of order 'order' interpolating the points
     * 'Order' is the value of the largest exponent in the segment
     * Initialize start and end points
     *
     * @throws std::runtime_error if points.size() < 2
     *
     * @param points points to interpolate
     * @param order order of interpolating segments
     */
    void initSegments(const std::vector<Point> &points, unsigned int order);
};
