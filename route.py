from itertools import chain

import numpy as np
from scipy.interpolate import interp1d
from .ext import rtree

def polyline_dists(points):
	diffs = np.diff(points, axis=0)
	dists = np.empty(len(points))
	dists[1:] = np.sqrt(np.sum(diffs**2, axis=1))
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
	if t > 1: return seglen, np.linalg.norm(p - seg[1])
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

		self.position_at = interp1d(distances, waypoints, axis=0)
	
	def hits_in_range(self, point, rng):
		bbox = tuple(chain(*((d-rng, d+rng) for d in point)))
		candidates = self.index.intersection(bbox, objects=False)
		for i in candidates:
			(distance, error) = lineseg_point_projection(
				self.segments[i], point)
			if error > rng: continue
			yield i, self.distances[i]+distance, error
	
	def closest_hits(self, point, n=1):
		bbox = tuple(chain(*((d, d) for d in point)))
		hits = self.index.nearest(bbox, n)
		for c, i in enumerate(hits):
			if c > n: break
			(distance, error) = lineseg_point_projection(
				self.segments[i], point)

			yield i, self.distances[i]+distance, error

def posnorm_logpdf(std):
	var = std**2
	normer = np.log(2.0/np.sqrt(2*np.pi*var))
	def logpdf(d):
		lik = normer - d**2/(2*var)
		return lik*np.NINF**(d < 0)
	return logpdf


def norm_logpdf(std):
	var = std**2
	normer = np.log(1.0/np.sqrt(2*np.pi*var))
	def logpdf(d):
		return normer - d**2/(2*var)
	return logpdf

def speed_logpdf(dist):
	def logpdf(dt, dd):
		return dist(dd/dt)
	return pdf

class State:
	def __init__(self, ts, dist, lik):
		self.ts = ts
		self.dist = dist
		self.lik = lik
		self.prev_state = None

class NaiveRouteMatcher:
	def __init__(self, route):
		if not hasattr(route, 'hits_in_range'):
			route = Route(route)
		self.route = route
		self.path = []
	
	def __call__(self, ts, point):
		self.path.append([ts, self.route.closest_hits(point).next()[1]])
	
	def get_path(self):
		return self.path
	

class RouteMatcher:
	def __init__(self, route,
			measurement_model=norm_logpdf(30),
			transition_model=norm_logpdf(25),
			search_range=3*30):
		if not hasattr(route, 'hits_in_range'):
			route = Route(route)
		self.route = route
		self.hypotheses = None
		self.measurement_model = measurement_model
		self.transition_model = transition_model
		self.search_range = search_range

	def __call__(self, ts, point):
		hits = self.route.hits_in_range(point, self.search_range)
		states = []
		for _, dist, err in hits:
			states.append(State(ts, dist, self.measurement_model(err)))
		
		if len(states) == 0:
			return

		if self.hypotheses is None:
			self.hypotheses = states
			return
		
		prev_dists = np.array([h.dist for h in self.hypotheses])
		prev_ts = np.array([h.ts for h in self.hypotheses])
		prev_lik = np.array([h.lik for h in self.hypotheses])

		for state in states:
			dd = state.dist - prev_dists
			dt = state.ts - prev_ts
			speeds = dd/dt
			liks = self.transition_model(speeds)
			
			liks += prev_lik
			prev = np.argmax(liks)

			state.prev_state = self.hypotheses[prev]

			state.lik += liks[prev]

		self.hypotheses = states
	
	def get_path(self):
		winner = np.argmax([h.lik for h in self.hypotheses])
		state = self.hypotheses[winner]

		path = []
		while state:
			path.append([state.ts, state.dist])
			state = state.prev_state
		return path[::-1]

