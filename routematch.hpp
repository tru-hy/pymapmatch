#include <vector>
#include <Eigen/Dense>
#include <spatialindex/SpatialIndex.h>
#include <spatialindex/capi/IdVisitor.h>

using namespace SpatialIndex;
using namespace Eigen;
using std::vector;
using std::cerr;

typedef double real;

template <size_t ndim>
real lineseg_point_projection(real *pr, real *ar, real *br, real &error);


template <size_t ndim>
class LineHitVisitor : public IVisitor {
	public:
	vector<id_type> ids;
	vector<real> distances;
	vector<real> errors;
	vector<LineSegment>& segments;
	Point& needle;
	LineHitVisitor(Point& needle, vector<LineSegment>& segments)
		: needle(needle), segments(segments) {}
	
	void visitNode(const INode& node) {}

	void visitData (const IData &in) {
		auto id = in.getIdentifier();
		LineSegment &shape = segments[id];
		real error;
		auto distance = lineseg_point_projection<ndim>(needle.m_pCoords,
			shape.m_pStartPoint, shape.m_pEndPoint, error);
		ids.push_back(id);
		distances.push_back(distance);
		errors.push_back(error);
	}
	void visitData (std::vector< const IData * > &v) {}
};

template <int ndim>
class Route
{
	private:
	ISpatialIndex *index = NULL;
	IStorageManager *index_storage = NULL;
	vector<LineSegment> segments;

	public:
	Array<real,Dynamic,1> *distances;
	Array<real,Dynamic,ndim> waypoints;
	Route(Array<real,Dynamic,ndim> waypoints);

	void hits_in_range(real *needle, real rng);
	LineHitVisitor<ndim> nearest_hits(real *needle, int n);

	~Route() {
		delete distances;
		delete index_storage;
	}

	
};



extern "C" {
	Route<2>* route2d_new(real *waypoints, size_t n) {
		Map<Array<real,Dynamic,2> > wayarr(waypoints, n, 2);
		return new Route<2>(wayarr);
	}

	void route2d_free(Route<2>* route2d) {
		delete route2d;
	}

	void route2d_hits_in_range(Route<2>* r, real *needle, real rng) {
		r->hits_in_range(needle, rng);
	}
	
	LineHitVisitor<2> route2d_nearest_hits(Route<2>* r, real *needle, int n) {
		r->nearest_hits(needle, n);
	}

	void route2d_naive_match(Route<2>* r, real *needles, size_t n, real *distances) {
		double p[2];
		for(size_t i = 0; i < n; i++) {
			p[0] = needles[i];
			p[1] = needles[n+i];
			auto hits = r->nearest_hits(p, 1);
			auto nodedist = (*(r->distances))(hits.ids[0]);
			distances[i] = hits.distances[0] + nodedist;
		}
	}
}
