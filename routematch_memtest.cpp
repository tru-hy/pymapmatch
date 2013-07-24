#include <routematch.h>

int main() {
	const size_t ndim = 2;
	auto route_n = 100;
	auto waypoints = Array<real, Dynamic, ndim>::Random(route_n, ndim);
	Route<ndim> route(waypoints);

	RouteMatcher<ndim> matcher(route);
	auto point_n = 100;
	auto ts = Array<real, Dynamic, 1>::LinSpaced(point_n, 0, point_n);
	auto points = Array<real, Dynamic, ndim, RowMajor>::Random(point_n, ndim);
	double point[ndim] = {0}; 
	for(auto i = 0; i < point_n; ++i) {
		Map<Array<real, 1, ndim, RowMajor> >(point, 1, ndim) = points.row(i);
		matcher.measurement(ts(i), point);
	}
}
