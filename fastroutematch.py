import ctypes
import os

import numpy as np
from numpy.ctypeslib import ndpointer


dllpath = os.path.join(os.path.dirname(__file__), 'routematch.so')
dllpath = os.path.abspath(dllpath)
rmatch = ctypes.CDLL(dllpath)


route2d = rmatch.route2d_new
route2d.restype = ctypes.c_void_p
route2d.argtypes = [ndpointer(ctypes.c_double), ctypes.c_size_t]
rmatch.route2d_naive_match.argtypes = [
	ctypes.c_void_p,
	ndpointer(ctypes.c_double),
	ctypes.c_int,
	ndpointer(ctypes.c_double)
	]

c_double_p = ctypes.POINTER(ctypes.c_double)
rmatch.route2d_hmm_match.argtypes = [
	ctypes.c_void_p,
	ndpointer(ctypes.c_double),
	ndpointer(ctypes.c_double),
	ctypes.c_int,
	ctypes.c_double,
	ctypes.c_double,
	ndpointer(ctypes.c_double),
	ndpointer(ctypes.c_double),
	]
rmatch.route2d_hmm_match.restype = ctypes.c_size_t

rmatch.route2d_distances.argtypes = [
	ctypes.c_void_p, ndpointer(ctypes.c_double)]

def _naive_match(r, coords, out=None):
	coords = np.asfortranarray(coords)
	if out is None:
		out = np.asfortranarray(np.empty(len(coords)))
	rmatch.route2d_naive_match(r, coords, len(coords), out)
	return out

def _route_match(r, ts, coords, measurement_std, transition_std):
	ts = np.asfortranarray(ts)
	coords = np.asfortranarray(coords)
	new_ts = np.empty(len(ts))
	new_dist = np.empty(len(ts))
	n = rmatch.route2d_hmm_match(r, ts, coords, len(ts),
		measurement_std, transition_std,
		new_ts, new_dist)
	new_ts = new_ts[:n]
	new_dist = new_dist[:n]
	return new_ts, new_dist

class RouteMatcher2d:
	def __init__(self, waypoints,
			measurement_std=30.0,
			transition_std=30.0):
		self.measurement_std = measurement_std
		self.transition_std = transition_std
		self._route = rmatch.route2d_new(np.asfortranarray(waypoints), len(waypoints))
		self.distances = np.empty(len(waypoints))
		rmatch.route2d_distances(self._route, self.distances)
	
	def __call__(self, ts, points):
		return _route_match(self._route, ts, points,
			self.measurement_std, self.transition_std)
	
	def __del__(self):
		rmatch.route2d_free(self._route)
		
