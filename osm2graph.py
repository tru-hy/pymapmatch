import sys

import numpy as np
from imposm.parser import OSMParser

def is_oneway(tags):
	if tags.get('oneway') == 'yes': return True
	if tags.get('highway') == 'motorway': return True
	if tags.get('junction') == 'roundabout': return True

def get_graph(filename):
	node_coords = {}
	edges = []
	etags = {}

	def node(id, *coords):
		node_coords[id] = coords

	def way(id, tags, refs):
		is_way = False
		if 'highway' in tags:
			is_way = True
		if 'busway' in tags:
			is_way = True
		if not is_way: return
		if tags['highway'] == 'footway': return
		if tags['highway'] == 'cycleway': return
		if tags['highway'] == 'steps': return
		
		for i in range(len(refs)-1):
			edge = refs[i], refs[i+1]
			etags[edge] = tags
			edges.append(edge)
			

		if not is_oneway(tags):
			tags['oneway'] = 'yes'
			way(None, tags, refs[::-1])
		
	def nodes(nodes):
		for args in nodes:
			node(*args)

	def ways(ways):
		for args in ways:
			way(*args)
	

	parser = OSMParser(
		ways_callback=ways,
		coords_callback=nodes,
		)
	parser.parse(filename)

	return node_coords, edges, etags

def euclidean_edge_costs(node_coords, edges):
	for e in edges:
		try:
			n1 = node_coords[e[0]]
			n2 = node_coords[e[1]]
			dist = np.sqrt(np.sum(np.subtract(n1, n2)**2))
		except KeyError:
			continue
		yield e, dist

def fastlines(segments, *args, **kwargs):
	import matplotlib.pyplot as plt
	invert = kwargs.pop('invert_dims', False)
	x = []
	y = []
	for a, b in segments:
		x.append(a[0])
		x.append(b[0])
		x.append(None)
		y.append(a[1])
		y.append(b[1])
		y.append(None)
	if invert:
		x, y = y, x
	return plt.plot(x, y, *args, **kwargs)

def plot_graph(nodes, edges, *args, **kwargs):
	def segments():
		for e in edges:
			try:
				yield nodes[e[0]], nodes[e[1]]
			except KeyError:
				continue
	return fastlines(segments(), *args, **kwargs)

def fasttest(f):
	print >>sys.stderr, "Building graph"
	nodes, edges, tags = get_graph(f)
	print >>sys.stderr, "%i nodes, %i edges"%(len(nodes), len(edges))
	
	from graph_tool.all import Graph

	g = Graph()
	node_map = {}
	
	print "Adding vertices"
	for osm_id in nodes.iterkeys():
		node_map[osm_id] = g.add_vertex()
	
	print "Adding edges"
	for a, b in edges:
		try:
			g.add_edge(node_map[a], node_map[b])
		except KeyError:
			continue
	
	
		

def plottest(f):
	from matplotlib.collections import LineCollection
	import matplotlib.pyplot as plt
	import numpy as np
	import random
	#from scipy.sparse import lil_matrix
	#from scipy.sparse.csgraph import dijkstra
	import networkx
	import time
	
	print >>sys.stderr, "Building graph"
	nodes, edges, tags = get_graph(f)
	print >>sys.stderr, "%i nodes, %i edges"%(len(nodes), len(edges))
	return
	#distances = lil_matrix((len(node_list), len(node_list)))
	durations = networkx.DiGraph()
	used_nodes = set()

	lines = []
	print >>sys.stderr, "Calculating durations"
	for edge in edges:
		try:
			a = nodes[edge[0]]
			b = nodes[edge[1]]
		except KeyError:
			#print >>sys.stderr, "Node not found!"
			continue
		# WOW! Doing cartesian!
		distance = np.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)
		#maxspeed = float(tags[edge].get('maxspeed', 10))
		duration = distance
		#durations[node_id_map[edge[0]],node_id_map[edge[1]]] = distance
		durations.add_edge(edge[0],edge[1], distance=distance, duration=duration)
		#plt.plot(*zip(a, b), color='black')
	

	#durations = networkx.strongly_connected_component_subgraphs(durations)[0]
	used_nodes = durations.nodes()

	lines = []
	for edge in durations.edges():
		a = nodes[edge[0]]
		b = nodes[edge[1]]
		lines.append((a, b))

	start, end = random.sample(used_nodes, 2)


	print >>sys.stderr, "Dijkstra...",
	t = time.time()
	path = networkx.shortest_path(durations, start, end, weight='duration')
	#line = [nodes[n] for n in path]
	dijt = time.time() - t
	print >>sys.stderr, dijt

	print >>sys.stderr, "Done"
	return
	"""
	print >>sys.stderr, "Astar...",
	t = time.time()
	apath = networkx.astar_path(distances, start, end,
		heuristic=lambda a, b: np.sqrt((nodes[a][0] - nodes[b][0])**2 + (nodes[a][1] - nodes[b][1])**2)
		)
	at = time.time() - t
	print >>sys.stderr, at
	print >>sys.stderr, "Astar times faster", dijt/at

	if path != apath:
		print path, apath
		print >>sys.stderr, "HEURISTICS FAILURE!"
	"""
	
	print >>sys.stderr, "Plotting"
	fastlines(lines, '.-k', alpha=0.1)
	#lines = LineCollection(lines)
	#lines.set_color('black')
	#lines.set_alpha(0.4)
	#ax = plt.subplot(1,1,1)
	#ax.add_collection(lines)
	#ax.autoscale_view()

	plt.scatter(*zip(nodes[start], nodes[end]))
	line = [nodes[n] for n in path]
	plt.plot(*zip(*line), color='red', linewidth=3, alpha=0.5)
	
	path = networkx.shortest_path(durations, start, end, weight='distance')
	line = [nodes[n] for n in path]
	plt.plot(*zip(*line), color='green', linewidth=3, alpha=0.5)

	
	#path = networkx.shortest_path(distances, start, end)
	#plt.plot(*zip(*line), color='green', linewidth=3, alpha=0.5)

	plt.show()

if __name__ == '__main__':
	#plottest(sys.argv[1])
	fasttest(sys.argv[1])

