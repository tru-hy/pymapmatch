from itertools import chain

import numpy as np
from .ext import rtree

def polyline_dists(points):
	diffs = np.diff(points, axis=1)
	dists = np.empty(len(points))
	dists[1:] = np.sqrt(np.sum(diffs**2, axis=0))
	dists[0] = 0.0
	dists = np.cumsum(dists)
	return dists

def polyline_bbox(points):
	for i in range(len(points)-1):
		dims = zip(points[i], points[i+1])
		bbox = chain(*((min(dim), max(dim)) for dim in dims))
		yield tuple(bbox)

def lineseg_point_projection(seg, point):
	segd = seg[1] - seg[0]
	seglen = np.linalg.norm(segd)
	p = np.array(point)
	t = np.dot(p - seg[0], segd)/seglen**2
	if t > seglen: return seglen, np.linalg.norm(p - seg[1])
	if t < 0: return 0, np.linalg.norm(p - seg[0])
	
	proj = seg[0] + t * segd
	return t*seglen, np.linalg.norm(p - proj)
	

class Route:
	def __init__(self, waypoints, distances=None):
		self.waypoints = waypoints
		if distances is None:
			distances = polyline_dists(waypoints)
		self.distances = distances
		
		indexitems = ((i, bbox, None)
			for (i, bbox) in enumerate(polyline_bbox(self.waypoints)))
		self.index = rtree.index.Index(indexitems, interleaved=False)
		# Use numpy arrays to speed up some calculations later
		self.segments = [np.array([self.waypoints[i], self.waypoints[i+1]])
			for i in range(len(waypoints)-1)]
	
	def hits_in_range(self, point, rng):
		bbox = tuple(chain(*((d-rng, d+rng) for d in point)))
		candidates = self.index.intersection(bbox, objects=False)
		for i in candidates:
			(distance, error) = lineseg_point_projection(
				self.segments[i], point)
			if error > rng: continue
			yield i, self.distances[i]+distance, error
