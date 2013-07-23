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

typedef double real;

template <size_t ndim>
real lineseg_point_projection(real *pr, real *ar, real *br, real &error);

template <size_t ndim>
class LineHitVisitor;

template <int ndim>
class Route
{
	private:
	ISpatialIndex *index = NULL;
	IStorageManager *index_storage = NULL;

	public:
	Array<real,Dynamic,1> *distances;
	Array<real,Dynamic,ndim> waypoints;
	Route(Array<real,Dynamic,ndim> waypoints);
	vector<LineSegment> segments;
	
	// TODO: Verify that these returns don't copy nor
	// leak!
	LineHitVisitor<ndim> hits_in_range(real *needle, real rng);

	LineHitVisitor<ndim> nearest_hits(real *needle, int n);

	~Route() {
		delete distances;
		delete index_storage;
	}
};

template <size_t ndim>
class LineHitVisitor : public IVisitor {
	public:
	vector<id_type> ids;
	vector<real> distances;
	vector<real> errors;
	Route<ndim>& route;
	Point& needle;
	LineHitVisitor(Point& needle, Route<ndim>& route)
		: needle(needle), route(route) {}
	
	void visitNode(const INode& node) {}

	void visitData (const IData &in) {
		auto id = in.getIdentifier();
		LineSegment &shape = route.segments[id];
		real error;
		auto distance = lineseg_point_projection<ndim>(needle.m_pCoords,
			shape.m_pStartPoint, shape.m_pEndPoint, error);
		ids.push_back(id);
		distances.push_back(distance + (*route.distances)(id));
		errors.push_back(error);
	}
	void visitData (std::vector< const IData * > &v) {}
};


struct RouteState {
	real timestamp;
	real distance;
	real loglikelihood;
	real cumloglikelihood;
	shared_ptr<RouteState> prev_state;
};

class RouteModel {
	public:
	virtual real measurement_loglik(real error) = 0;
	virtual real transition_loglik(real speed) = 0;
};

real gaussian_logpdf(real var, real x)
{
	auto normer = std::log(1.0/std::sqrt(2.0*M_PI*var));
	return normer - x*x/(2*var);
}

class GaussianRouteModel : public RouteModel {
	public:
	real measurement_var;
	real speed_var;
	
	GaussianRouteModel(real measurement_std, real speed_std)
	{
		measurement_var = measurement_std*measurement_std;
		speed_var = speed_std*speed_std;
	}

	real measurement_loglik(real error)
	{
		return gaussian_logpdf(measurement_var, error);
	}
	
	real transition_loglik(real speed)
	{
		return gaussian_logpdf(speed_var, speed);
	}
};

GaussianRouteModel DEFAULT_ROUTE_MODEL(30, 30);

template <int ndim>
class RouteMatcher {
	shared_ptr<RouteState>* hypotheses = NULL;
	size_t n_hypotheses = 0;
	public:
	int path_len = 0;
	Route<ndim>& route;
	RouteModel* route_model;
	real search_range;

	RouteMatcher(Route<ndim>& route, real search_range=100,
			RouteModel* route_model=&DEFAULT_ROUTE_MODEL)
		: route(route), search_range(search_range), route_model(route_model)
	{
		
	}

	void measurement(real ts, real point[ndim])
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

	void get_path(real *ts, real *dist) {
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

	private:
	void find_hypothesis(RouteState& state) {
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
};

extern "C" {
	Route<2>* route2d_new(real *waypoints, size_t n) {
		Map<Array<real,Dynamic,2> > wayarr(waypoints, n, 2);
		return new Route<2>(wayarr);
	}

	void route2d_free(Route<2>* route2d) {
		delete route2d;
	}
	
	void route2d_naive_match(Route<2>* r, real *needles, size_t n, real *distances) {
		double p[2];
		for(size_t i = 0; i < n; i++) {
			p[0] = needles[i];
			p[1] = needles[n+i];
			auto hits = r->nearest_hits(p, 1);
			distances[i] = hits.distances[0];
		}
	}

	size_t route2d_hmm_match(Route<2>* r, real *ts, real *pos, size_t n,
			real *outts, real *outdist) {
		RouteMatcher<2> matcher(*r);
		double p[2];
		for(size_t i = 0; i < n; i++) {
			p[0] = pos[i];
			p[1] = pos[n+i];
			matcher.measurement(ts[i], p);
		}
		
		matcher.get_path(outts, outdist);
		return matcher.path_len;
	}
}

