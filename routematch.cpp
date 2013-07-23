#include "routematch.hpp"

#include <iostream>
#include <vector>

using std::cerr;
using std::vector;


template <int ndim>
Route<ndim>::Route(Array<real,Dynamic,ndim> waypoints)
	: waypoints(waypoints)
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

template <size_t ndim>
real lineseg_point_projection(real *pr, real *ar, real *br, real &error)
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
LineHitVisitor<ndim> Route<ndim>::hits_in_range(real *needle, real rng)
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
LineHitVisitor<ndim> Route<ndim>::nearest_hits(real *needle, int n) {
	Point np(needle, ndim);
	LineHitVisitor<ndim> result(np, *this);
	index->nearestNeighborQuery(n, np, result);
	return result;
}

