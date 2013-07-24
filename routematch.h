// routematch.h
//

#ifndef LZZ_routematch_h
#define LZZ_routematch_h
#include <vector>
#include <cmath>
#include <memory>
#include <Eigen/Dense>
#include <spatialindex/SpatialIndex.h>
#include <spatialindex/capi/IdVisitor.h>

using namespace SpatialIndex;
using namespace Eigen;
using std::vector;
using std::cerr;
using std::shared_ptr;
using std::weak_ptr;
class GaussianRouteModel;
extern GaussianRouteModel DEFAULT_ROUTE_MODEL;
#define LZZ_INLINE inline
typedef double real;
template <size_t ndim>
class LineHitVisitor;
template <size_t ndim>
real lineseg_point_projection (real * pr, real * ar, real * br, real & error);
template <int ndim>
class Route
{
private:
  ISpatialIndex * index;
  IStorageManager * index_storage;
public:
  Array <real,Dynamic,1> * distances;
  Array <real,Dynamic,ndim> waypoints;
  vector <LineSegment> segments;
  Route (Array <real,Dynamic,ndim> waypoints);
  ~ Route ();
  LineHitVisitor <ndim> hits_in_range (real * needle, real rng);
  LineHitVisitor <ndim> nearest_hits (real * needle, int n);
};
template <size_t ndim>
class LineHitVisitor : public IVisitor
{
public:
  vector <id_type> ids;
  vector <real> distances;
  vector <real> errors;
  Route <ndim> & route;
  Point & needle;
  LineHitVisitor (Point & needle, Route <ndim> & route);
  void visitNode (INode const & node);
  void visitData (IData const & in);
  void visitData (std::vector < const IData * > & v);
};
struct RouteState
{
  real timestamp;
  real distance;
  real loglikelihood;
  real cumloglikelihood;
  shared_ptr <RouteState> prev_state;
};
class RouteModel
{
public:
  virtual real measurement_loglik (real error) = 0;
  virtual real transition_loglik (real speed) = 0;
};
real gaussian_logpdf (real var, real x);
class GaussianRouteModel : public RouteModel
{
public:
  real measurement_var;
  real speed_var;
  GaussianRouteModel (real measurement_std, real speed_std);
  real measurement_loglik (real error);
  real transition_loglik (real speed);
};
template <int ndim>
class RouteMatcher
{
  shared_ptr <RouteState> * hypotheses;
  size_t n_hypotheses;
public:
  int path_len;
  Route <ndim> & route;
  RouteModel * route_model;
  real search_range;
  RouteMatcher (Route <ndim> & route, real search_range = 100, RouteModel * route_model = &DEFAULT_ROUTE_MODEL);
  void measurement (real ts, real (point) [ndim]);
  void get_path (real * ts, real * dist);
private:
  void find_hypothesis (RouteState & state);
};
extern "C"
{
  Route <2> * route2d_new (real * waypoints, size_t n);
}
extern "C"
{
  void route2d_free (Route <2> * route2d);
}
extern "C"
{
  void route2d_distances (Route <2> * r, real * distances);
}
extern "C"
{
  void route2d_naive_match (Route <2> * r, real * needles, size_t n, real * distances);
}
extern "C"
{
  size_t route2d_hmm_match (Route <2> * r, real * ts, real * pos, size_t n, real * outts, real * outdist);
}
template <size_t ndim>
real lineseg_point_projection (real * pr, real * ar, real * br, real & error)
{
	Map<Matrix<real,1,ndim> > p(pr, 1, ndim);
	Map<Matrix<real,1,ndim> > a(ar, 1, ndim);
	Map<Matrix<real,1,ndim> > b(br, 1, ndim);
	auto segd = b - a;
	auto seglen = segd.norm();
	auto normstart = p - a;
	auto t = normstart.dot(segd)/(seglen*seglen);
	if(t > 1) {
		error = (p - b).norm();
		return seglen;
	}

	if(t < 0) {
		error = normstart.norm();
		return 0.0;
	}
	
	auto proj = a + t*segd;
	error = (p - proj).norm();
	return t*seglen;

}
template <int ndim>
Route <ndim>::Route (Array <real,Dynamic,ndim> waypoints)
  : waypoints (waypoints)
        {
		distances = new Array<real,Dynamic,1>(waypoints.rows(), 1);
		auto& distref = *distances;
		auto wpm = waypoints.matrix();
		
		distref(0) = 0.0;
		size_t n = distref.rows();
		for(size_t i = 1; i < n; ++i) {
			distref(i) = (wpm.row(i) - wpm.row(i-1)).norm() + distref(i-1);
		}
		
		// TODO: Bulk load
		index_storage = StorageManager::createNewMemoryStorageManager();
		id_type indexIdentifier;
		index = RTree::createNewRTree(*index_storage,
			0.7, 10, 10, ndim, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
		
		real startpoint[ndim];
		real endpoint[ndim];
		for(size_t i = 0; i < n-1; ++i) {
			for(int d = 0; d < ndim; ++d) {
				startpoint[d] = waypoints(i,d);
				endpoint[d] = waypoints(i+1,d);
			}
			LineSegment seg(startpoint, endpoint, ndim);
			index->insertData(0, NULL, seg, i);
			segments.push_back(seg);
		}
	}
template <int ndim>
Route <ndim>::~ Route ()
                 {
		delete distances;
		delete index;
		delete index_storage;
	}
template <int ndim>
LineHitVisitor <ndim> Route <ndim>::hits_in_range (real * needle, real rng)
        {
		Point np(needle, ndim);
		double start[ndim];
		double end[ndim];
		for(int d = 0; d < ndim; ++d) {
			start[d] = needle[d] - rng;
			end[d] = needle[d] + rng;
		}

		Region bbox(start, end, ndim);
		LineHitVisitor<ndim> result(np, *this);
		index->intersectsWithQuery(bbox, result);
		return result;

	}
template <int ndim>
LineHitVisitor <ndim> Route <ndim>::nearest_hits (real * needle, int n)
        {
		Point np(needle, ndim);
		LineHitVisitor<ndim> result(np, *this);
		index->nearestNeighborQuery(n, np, result);
		return result;
	}
template <size_t ndim>
LineHitVisitor <ndim>::LineHitVisitor (Point & needle, Route <ndim> & route)
  : needle (needle), route (route)
                                               {}
template <size_t ndim>
void LineHitVisitor <ndim>::visitNode (INode const & node)
                                          {}
template <size_t ndim>
void LineHitVisitor <ndim>::visitData (IData const & in)
                                         {
		auto id = in.getIdentifier();
		LineSegment &shape = route.segments[id];
		real error;
		auto distance = lineseg_point_projection<ndim>(needle.m_pCoords,
			shape.m_pStartPoint, shape.m_pEndPoint, error);
		ids.push_back(id);
		distances.push_back(distance + (*route.distances)(id));
		errors.push_back(error);
	}
template <size_t ndim>
void LineHitVisitor <ndim>::visitData (std::vector < const IData * > & v)
                                                         {}
template <int ndim>
RouteMatcher <ndim>::RouteMatcher (Route <ndim> & route, real search_range, RouteModel * route_model)
  : route (route), search_range (search_range), route_model (route_model), n_hypotheses (0), path_len (0)
        {
		
	}
template <int ndim>
void RouteMatcher <ndim>::measurement (real ts, real (point) [ndim])
        {
		auto hits = route.hits_in_range(point, search_range);
		auto n_hits = hits.ids.size();
		if(n_hits == 0) return;
		path_len++;

		auto states = new shared_ptr<RouteState>[hits.ids.size()];
		for(size_t i = 0; i < n_hits; ++i) {
			auto state = shared_ptr<RouteState>(new RouteState);
			state->timestamp = ts;
			state->distance = hits.distances[i];
			state->loglikelihood = route_model->measurement_loglik(hits.errors[i]);
			state->cumloglikelihood = state->loglikelihood;
			states[i] = state;
		}
		
		if(!hypotheses) {
			hypotheses = states;
			n_hypotheses = n_hits;
			return;
		}
		
		for(size_t i = 0; i < n_hits; ++i) find_hypothesis(*states[i]);

		delete [] hypotheses;
		hypotheses = states;
		n_hypotheses = n_hits;
	}
template <int ndim>
void RouteMatcher <ndim>::get_path (real * ts, real * dist)
                                            {
		real max_lik = -INFINITY;
		shared_ptr<RouteState> state;
		
		for (size_t i = 0; i < n_hypotheses; ++i) {
			auto lik = hypotheses[i]->cumloglikelihood;
			if (lik > max_lik) {
				state = hypotheses[i];
				max_lik = lik;
			}
		}
		
		vector<real> tsv;
		vector<real> distv;
		size_t i = path_len-1;
		while(state) {
			ts[i] = state->timestamp;
			dist[i] = state->distance;
			state = state->prev_state;
			--i;
		}
	}
template <int ndim>
void RouteMatcher <ndim>::find_hypothesis (RouteState & state)
                                                {
		real max_lik = -INFINITY;
		shared_ptr<RouteState> winner;

		for (size_t i = 0; i < n_hypotheses; ++i) {
			auto dt = state.timestamp - hypotheses[i]->timestamp;
			auto dd = state.distance - hypotheses[i]->distance;
			auto speed = dd/dt;
			auto trans_lik = route_model->transition_loglik(speed);
			auto lik = hypotheses[i]->cumloglikelihood + trans_lik;
			if (lik > max_lik) {
				winner = hypotheses[i];
				max_lik = lik;
			}
		}
		
		state.prev_state = winner;
		state.cumloglikelihood += max_lik;
	}
#undef LZZ_INLINE
#endif
