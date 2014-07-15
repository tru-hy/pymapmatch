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

%callback("%(upper)s");
WayRole busway_filter(const readosm_way*);
WayRole train_filter(const readosm_way*);
WayRole tram_filter(const readosm_way*);
WayRole subway_filter(const readosm_way*);
%nocallback;
%apply real *OUTPUT { real *x, real *y };
%apply double *OUTPUT { double *latitude, double *longitude };
%include "osmmapmatch.hpp";


