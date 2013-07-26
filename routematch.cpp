// routematch.cpp
//

#include "routematch.h"
GaussianRouteModel DEFAULT_ROUTE_MODEL(30, 30);
#define LZZ_INLINE inline
real gaussian_logpdf (real var, real x)
{
	auto normer = std::log(1.0/std::sqrt(2.0*M_PI*var));
	return normer - x*x/(2*var);
}
GaussianRouteModel::GaussianRouteModel (real measurement_std, real speed_std)
        {
		measurement_var = measurement_std*measurement_std;
		speed_var = speed_std*speed_std;
	}
real GaussianRouteModel::measurement_loglik (real error)
        {
		return gaussian_logpdf(measurement_var, error);
	}
real GaussianRouteModel::transition_loglik (real speed)
        {
		return gaussian_logpdf(speed_var, speed);
	}
extern "C"
{
  Route <2> * route2d_new (real * waypoints, size_t n)
                                                         {
		Map<Array<real,Dynamic,2> > wayarr(waypoints, n, 2);
		return new Route<2>(wayarr);
	}
}
extern "C"
{
  void route2d_free (Route <2> * route2d)
                                             {
		delete route2d;
	}
}
extern "C"
{
  void route2d_distances (Route <2> * r, real * distances)
                                                             {
		auto& distref = *(r->distances);
		size_t n = distref.size();

		for(size_t i = 0; i < n; ++i) {
			distances[i] = distref(i);
		}
	}
}
extern "C"
{
  void route2d_naive_match (Route <2> * r, real * needles, size_t n, real * distances)
                                                                                        {
		double p[2];
		for(size_t i = 0; i < n; i++) {
			p[0] = needles[i];
			p[1] = needles[n+i];
			auto hits = r->nearest_hits(p, 1);
			distances[i] = hits.distances[0];
		}
	}
}
extern "C"
{
  size_t route2d_hmm_match (Route <2> * r, real * ts, real * pos, size_t n, real measurement_std, real transition_std, real * outts, real * outdist)
                                                    {
		GaussianRouteModel model(measurement_std, transition_std);
		RouteMatcher<2> matcher(*r, &model);
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
#undef LZZ_INLINE
