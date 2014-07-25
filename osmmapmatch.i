%module osmmapmatch
%{
#include "osmmapmatch.hpp"
%}
%include "typemaps.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_shared_ptr.i"

%shared_ptr(PositionHypothesis);

namespace std {
        %template(Point2dVector) vector< Point2d >;
        %template(EdgePoint2d) pair<Point2d, Point2d>;
        %template(EdgePoint2dVector) vector< pair<Point2d, Point2d> >;
        %template(RealVector) vector< double >;
        %template(HypothesisVector) vector< shared_ptr<PositionHypothesis> >;
        %template(NodeVector) vector< long long >;
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

%extend PositionHypothesis {
        Point2d position2d() {
                return Point2d($self->position);
        }

        Point2d measurement2d() {
                return Point2d($self->measurement);
        }

        std::vector< node_id_t > subpath_nodes(OsmGraph& graph) {
                auto node_id_prop = bst::get(bst::vertex_name, graph.graph);
                std::vector < node_id_t > result;
                for(auto vertex: $self->subpath) {
                        result.push_back(node_id_prop[vertex]);
                }

                return result;
        }
};
