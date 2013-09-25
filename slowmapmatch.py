from itertools import izip
import sys

import numpy as np
import rtree
import networkx


def lineseg_point_projection(p, a, b):
	(p, a, b) = map(np.array, (p, a, b))
	segd = b - a
	seglen = np.linalg.norm(segd)
	normstart = p - a
	t = np.dot(normstart, segd)/(seglen**2)
	if t > 1:
		error = np.linalg.norm(p - b)
		t = 1.0
		return t, error, b

	if t < 0:
		error = np.linalg.norm(normstart)
		t = 0.0
		return t, error, a
	
	proj = a + t*segd
	error = np.linalg.norm(p - proj)
	return t, error, proj

class _State(object):
	def __init__(self, ts, position, point, parent):
		self.ts = ts
		self.position = position
		self.point = point
		self.parent = parent
		self.cumlik = 0.0
		self.path = []
		self.transition_lik = np.nan

def are_same_node((se, st), (e, t)):
	if st not in (0.0, 1.0):
		return
	if t not in (0.0, 1.0):
		return
	return se[int(st)] == e[int(t)]

class _Matcher(object):
	def __init__(self, model):
		self.model = model
		self.graph = self.model.graph
		self.states = []
		self._prev_pos = None
	
	def __call__(self, ts, p):
		r = self.model.search_radius
		region = (
			p[0]-r, p[1]-r,
			p[0]+r, p[1]+r
			)
		
		# Find edges in the search radius and
		# project the current measurement to them
		candidates = self.model.segindex.intersection(region, objects='raw')
		positions = []
		for cand in candidates:
			a = self.model.node_coords[cand[0]]
			b = self.model.node_coords[cand[1]]
			t, error, point = lineseg_point_projection(p, a, b)
			if error > r:
				continue
			positions.append((cand, t, error, point))
		
		# Prune duplicates that hit exactly on
		# the nodes.
		# TODO: For some reason this pruning affects
		# results. Probably some logic flaw in the programmer.
		# Which is a shame as this does remove a lot of hypotheses.
		"""
		inter_node_hits = []
		exact_node_hits = {}
		for (e, t, error, point) in positions:
			node = None
			if t == 0.0:
				node = e[0]
			if t == 1.0:
				node = e[1]
			if node is None:
				inter_node_hits.append((e, t, error, point))
				continue
			if node in exact_node_hits: continue
			exact_node_hits[node] = (e, t, error, point)

		hits = inter_node_hits + exact_node_hits.values()
		"""
		hits = positions
		
		# Do nothing if no hits found. Not catastrophical
		# if this is due to an outlier.
		# TODO: We could track "out of map" trajectories
		if len(hits) == 0:
			print >>sys.stderr, "NO HITS"
			return

		# No previous states, start with current hits
		if len(self.states) == 0:
			for e, t, error, point in positions:
				s = _State(ts, (e, t), point, None)
				s.measurement_lik = self.model.measurement_logpdf(error)
				s.cumlik += s.measurement_lik
				self.states.append(s)
			self._prev_pos = p
			return

		# Find most likely predecessor for each
		# hit.
		new_states = []
		for (e, t, error, point) in hits:
			# TODO: No need to add temporary nodes for
			# exact hits
			# Add temporary target
			dist = self.model.edge_costs[e]
			self.graph.add_edge(e[0], "tmptarget", weight=t*dist)
			self.graph.add_edge("tmptarget", e[1], weight=(1.0-t)*dist)
			self.model.node_coords['tmptarget'] = point

			max_lik = -np.inf
			max_hit = None
			for prev in self.states:
				# Add temporary target node
				se, st = prev.position
				dist = self.model.edge_costs[se]
				if e == se and st <= t:
					# The source is on the same edge than the
					# start and before, so it can jump straight
					# there. 
					self.graph.add_edge("tmpsource", "tmptarget",
						weight=dist*(t-st))
				
				# Allow self-transition if the nodes are exactly the same.
				# TODO: Get rid of this by getting rid of the whole temporary
				#	node mess.
				if are_same_node((se, st), (e, t)):
					self.graph.add_edge("tmpsource", "tmptarget",
						weight=0.0)
				
				self.graph.add_edge(se[0], "tmpsource", weight=st*dist)
				self.graph.add_edge("tmpsource", se[1], weight=(1.0-st)*dist)
				self.model.node_coords['tmpsource'] = prev.point
				# Find shortest path
				try:
					distance, path = networkx.bidirectional_dijkstra(
						self.graph,
						"tmpsource", "tmptarget",
						weight='weight')
					"""
					def euclidean_dist(a, b):
						return np.linalg.norm(np.subtract(
							self.model.node_coords[a],
							self.model.node_coords[b]))
					path = networkx.astar_path(self.graph, "tmpsource", "tmptarget",
						euclidean_dist, weight='weight')
					distance = sum(self.graph[path[i]][path[i+1]]['weight']
						for i in range(len(path)-1))
					"""
				except networkx.exception.NetworkXNoPath:
					continue
				finally:
					# Remove temporary target node
					self.graph.remove_node("tmpsource")
				
				
				dt = ts - prev.ts

				path_coords = ([prev.point] +
					[self.model.node_coords[n] for n in path[1:-1]] +
					[point])
				totlik = prev.cumlik + self.model.transition_logpdf(distance, dt,
						(self._prev_pos, p), path_coords)
				if totlik > max_lik:
					max_lik = totlik
					max_hit = (prev, path)
			
			# Remove temporary source node
			self.graph.remove_node("tmptarget")

			if max_hit is None:
				# Non-reachable target
				continue
			
			s = _State(ts, (e, t), point, max_hit[0])
			s.path = max_hit[1][1:-1]
			s.cumlik += max_lik
			s.cumlik += self.model.measurement_logpdf(error)
			new_states.append(s)

		if len(new_states) == 0:
			# No valid transitions, ignore current measurement
			print >>sys.stderr, "NO TRANSITIONS"
			return
		
		self._prev_pos = p
		self.states = new_states
	
	def get_path(self):
		state = self.states[np.argmax([s.cumlik for s in self.states])]
		path = []
		while state:
			path.extend([state.position] + state.path[::-1])
			state = state.parent
		return path[::-1]
	
	def get_winner_state_path(self):
		state = self.states[np.argmax([s.cumlik for s in self.states])]
		states = []
		while state:
			states.append(state)
			state = state.parent
		return states[::-1]

	def get_map_path(self):
		# Path with "non-map-nodes" except for start and end
		# skipped. Gives same geometry with less points.
		states = self.get_winner_state_path()
		coords = []
		coords.append(states[0].position)
		for state in states[1:]:
			coords.extend(n for n in state.path)
		coords.append(states[-1].position)
		return coords
		
	
	def get_map_coordinates(self):
		states = self.get_winner_state_path()
		coords = []
		coords.append(states[0].point)
		for state in states[1:]:
			coords.extend(self.model.node_coords[n] for n in state.path)
		coords.append(states[-1].point)
		return coords
		
	
	def get_coordinates(self):
		state = self.states[np.argmax([s.cumlik for s in self.states])]
		path = []
		while state:
			#path.extend([state.point] +
			#	[self.model.node_coords[n] for n in state.path[::-1]])
			path.append(state.point)
			state = state.parent
		return path[::-1]

def gaussian_logpdf(std):
	var = std**2
	normer = np.log(1.0/(np.sqrt(2.0*np.pi)*std))
	return lambda x: normer - x**2/(2*var)

def speed_gaussian_logpdf(std):
	logpdf = gaussian_logpdf(std)
	return lambda dist, dt, points, path: logpdf(dist/dt)

class MapMatcher2d(object):
	def __init__(self, edge_costs, node_coordinates,
			search_radius=50.0,
			measurement_logpdf=gaussian_logpdf(30),
			transition_logpdf=speed_gaussian_logpdf(30)
			):
		self.edge_costs = dict(edge_costs)
		self.node_coords = node_coordinates
		self.graph = networkx.DiGraph()

		self.search_radius = search_radius
		self.measurement_logpdf=measurement_logpdf
		self.transition_logpdf=transition_logpdf
		
		def iter_segments():
			for i, (e, cost) in enumerate(self.edge_costs.iteritems()):
				try:
					n1 = self.node_coords[e[0]]
					n2 = self.node_coords[e[1]]
				except KeyError:
					continue
				
				self.graph.add_edge(*e, weight=cost)
				bbox = (
					min(n1[0], n2[0]),
					min(n1[1], n2[1]),
					max(n1[0], n2[0]),
					max(n1[1], n2[1]),
					)
				yield (i, bbox, e)
		self.segindex = rtree.index.Index(iter_segments())

	def __call__(self, ts, points):
		matcher = _Matcher(self)
		vals = []
		for i, (t, p) in enumerate(izip(ts, points)):
			matcher(t, p)
		return matcher

#if __name__ == '__main__':
	
