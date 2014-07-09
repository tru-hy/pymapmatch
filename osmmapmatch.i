%module osmmapmatch
%{
#include "osmmapmatch.hpp"
%}
%include "typemaps.i"
%include "std_vector.i"
%include "std_pair.i"

namespace std {
        %template(Point2dVector) vector< Point2d >;
        %template(EdgePoint2d) pair<Point2d, Point2d>;
        %template(EdgePoint2dVector) vector< pair<Point2d, Point2d> >;
        %template(RealVector) vector< double >;
}

%apply real *OUTPUT { real *x, real *y };
%include "osmmapmatch.hpp"
