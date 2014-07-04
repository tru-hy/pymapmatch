// TODO: This is currently double as float causes
// headaches with PROJ.4.
typedef double real;

#ifndef SWIG
/* The code is horrible. I'm deeply sorry. C++, and especially
 * Boost, is hell. */
#include <readosm.h>
#include <malloc.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <ctime>
#include <queue>
#include <set>
#include <list>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/random.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <proj_api.h>

using std::unique_ptr;
using std::shared_ptr;
using std::vector;
using std::string;
namespace bst = boost;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef long long node_id_t;


typedef bst::adjacency_list<bst::vecS, bst::vecS, bst::bidirectionalS,
	bst::property<bst::vertex_name_t, node_id_t>,
	bst::property<bst::edge_weight_t, real>
	> Graph;
typedef bst::graph_traits < Graph >::vertex_descriptor Vertex;
typedef bst::graph_traits < Graph >::edge_descriptor Edge;

real rad_to_deg(real rad) {
	return rad / (M_PI/180.0);
}

real deg_to_rad(real deg) {
	return deg * (M_PI/180.0);
}

class OsmReaderError : public std::runtime_error {
	using std::runtime_error::runtime_error;
};

enum WayRole {
	IgnoreWay = 0,
	OneWayWay,
	TwoWayWay
};

typedef bg::model::point<real, 2, bg::cs::cartesian> Point;
typedef bg::model::segment<Point> LineSegment;
typedef bg::model::box<Point> Bbox;


// OMFG! I have to implement this manually!?!
// Also giving up on generality, as, you know, Boost.
inline real vect_norm(Point p) {
	real c = p.get<0>();
	real result = c*c;
	c = p.get<1>();
	result += c*c;

	return std::sqrt(result);
}

// There's probably something like this in the depths
// of Boost.Geometry, but couldn't find it, so here we go again.
std::pair<real, real> lineseg_point_projection(LineSegment seg, Point p)
{
	auto a = get<0>(seg);
	auto b = get<1>(seg);

	// I guess Boost tries to make the APIs as difficut
	// as possible. This should be: auto segd = b - a;
	// but boost fucking wants this shit:
	auto segd = b; bg::subtract_point(segd, a);
	auto seglen = vect_norm(segd);
	auto normstart = p; bg::subtract_point(normstart, a);
	auto t = bg::dot_product(normstart, segd)/(seglen*seglen);

	real error;

	if(t > 1.0) {
		t = 1.0;
	}

	if(t < 0.0) {
		t = 0.0;
	}
	
	// Reusing stuff because the boost api is a fucking disaster.
	// With Eigen this was: auto proj = a + t*segd;
	bg::multiply_value(segd, t);
	auto proj = a;
	bg::add_point(proj, segd);
	
	bg::subtract_point(p, proj);
	error = vect_norm(p);
	return std::make_pair(error, t*seglen);
}


class SearchTargetFound {};
class SearchTargetNotFound {};

/*
class _SingleTargetVisitor : public bst::default_astar_visitor {
	const Vertex& target;
	public:
	_SingleTargetVisitor(Vertex& t) : target(t) {

	}

	template <class Graph>
	void examine_vertex(Vertex v, Graph& g) {
		if(v == target) {
			throw SearchTargetFound();
		}
	}
}*/

template <class DistanceMap>
class _SingleTargetVisitor : public bst::default_dijkstra_visitor {
	const Vertex& target;
	DistanceMap& dist_map;
	public:
	
	_SingleTargetVisitor(Vertex& t, DistanceMap& m)
		: target(t), dist_map(m) {

	}

	template <class Graph>
	void examine_vertex(Vertex v, Graph& g) {
		/*if(dist_map[v] > 100.0) {
			std::cout << "Bailing out" << std::endl;
			throw SearchTargetNotFound();
		}*/
		
		if(v == target) {
			throw SearchTargetFound();
		}
	}
};

class distance_heuristic : public bst::astar_heuristic<Graph, real>
{
	public:
	real operator()(Vertex u) {
		return 0.0;
	}
};

template <class Key, class Value>
class default_map
{	
	private:
	std::unordered_map<Key, Value> m_map;
	Value deflt;

	public:
	typedef std::pair<Key, Value> value_type;
	typedef Key key_type;

	default_map(const Value def): deflt(def) {

	}

	Value& operator[](const Key& k) {
		if(m_map.count(k) == 0) {
			m_map[k] = deflt;
		}
		
		return m_map[k];
	}

	size_t count(const Key& k) {
		return m_map.count(k);
	}
};

template <class Key, class Value>
class more_hacky_shit_maps_for_the_crappy_boost_api
{	
	private:
	std::unordered_map<Key, Value> m_map;
	Value deflt;

	public:
	typedef std::pair<Key, Value> value_type;
	typedef Key key_type;

	Value& operator[](const Key& k) {
		if(m_map.count(k) == 0) {
			m_map[k] = k;
		}
		
		return m_map[k];
	}
};

real single_target_shortest_path(Vertex source, Vertex target, const Graph& graph, real cutoff=1.0/0.0) {
	default_map<Vertex, real> distances(1.0/0.0);
	default_map<Vertex, real> cost(1.0/0.0);
	std::unordered_map<Vertex, Vertex> predecessor;
	auto weights = bst::get(bst::edge_weight, graph);

	typedef std::pair<real, Vertex> Node;
	std::priority_queue<Node, std::vector<Node>, std::greater<Node>> queue;
	queue.push(std::make_pair(0.0, source));
	predecessor[source] = source;
	
	while(!queue.empty()) {
		auto current = queue.top();
		queue.pop();
		if(current.first > cutoff) {
			break;
		}
		distances[current.second] = current.first;
		if(current.first == target) {
			break;
		}

		auto children = bst::out_edges(current.second, graph);
		for(auto edge=children.first; edge!=children.second; edge++) {
			auto alt = distances[current.second] + weights[*edge];
			auto child = bst::target(*edge, graph);
			
			if(distances.count(child)) {
				// Sanity check that we don't get a new distance
				// for an already finished node
				assert(distances[child] <= alt);
				// We've already finished this one
				continue;
			}
			
			if(alt < cost[child]) {
				cost[child] = alt;
				queue.push(std::make_pair(alt, child));
				predecessor[child] = current.first;
			}
		}
	}

	return distances[target];
}

template <class Callback>
void reverse_shortest_path(Vertex source, const Graph& graph, Callback callback) {
	default_map<Vertex, real> distances(1.0/0.0);
	default_map<Vertex, real> cost(1.0/0.0);
	std::unordered_map<Vertex, Vertex> successor;
	auto weights = bst::get(bst::edge_weight, graph);

	typedef std::pair<real, Vertex> Node;
	std::priority_queue<Node, std::vector<Node>, std::greater<Node>> queue;
	queue.push(std::make_pair(0.0, source));
	
	while(!queue.empty()) {
		auto top = queue.top();
		auto current = top.second;
		auto current_dist = top.first;
		queue.pop();
		distances[current] = current_dist;
		if(!callback(current, current_dist, successor)) {
			break;
		}

		auto children = bst::in_edges(current, graph);
		for(auto edge=children.first; edge!=children.second; edge++) {
			auto alt = current_dist + weights[*edge];
			auto child = bst::source(*edge, graph);
			
			if(distances.count(child)) {
				// Sanity check that we don't get a new distance
				// for an already finished node
				assert(distances[child] <= alt);
				// We've already finished this one
				continue;
			}
			
			if(alt < cost[child]) {
				cost[child] = alt;
				queue.push(std::make_pair(alt, child));
				successor[child] = current;
			}
		}
	}
}

real single_target_shortest_path_slow(Vertex source, Vertex target, Graph graph) {
	default_map<Vertex, real> distances(1.0/0.0);
	default_map<Vertex, real> costs(1.0/0.0);
	distances[target] = 0.0;
	
	//std::map<Vertex, Vertex> parents;
	more_hacky_shit_maps_for_the_crappy_boost_api<Vertex, Vertex> parents;

	auto parent_map = bst::associative_property_map<decltype(parents)>(parents);
	auto dist_map = bst::associative_property_map<decltype(distances)>(distances);
	auto costs_map = bst::associative_property_map<decltype(costs)>(costs);
	
	/*auto colorshit = default_map<Vertex, bst::default_color_type>(bst::default_color_type::white);
	auto color_map = bst::associative_property_map<decltype(colorshit)>(colorshit);*/
	
	auto found = false;
	try {
	bst::astar_search(graph, source, distance_heuristic(),
		bst::predecessor_map(parent_map).
		visitor(_SingleTargetVisitor<decltype(distances)>(target, distances)).
		distance_map(dist_map).
		rank_map(costs_map)
		);
	} catch(SearchTargetFound e) {
		found = true;
	} catch(SearchTargetNotFound e) {

	}

	
	typedef typename bst::property_traits<decltype(dist_map)>::value_type D;
	D inf = std::numeric_limits<D>::max();
	

	/*
	try {
	bst::dijkstra_shortest_paths(
		graph,
		source,
		bst::associative_property_map<decltype(parents)>(parents),
		dist_map,
		bst::get(bst::edge_weight, graph), 
		bst::get(bst::vertex_index, graph), 
		std::less<D>(), // compare
		bst::closed_plus<D>(inf), // combine
		D(), // zero
		_SingleTargetDijkstraVisitor<decltype(distances)>(target, distances) // visitor
		);
	} catch(SearchTargetFound e) {
		found = true;
	}
	*/

	if(!found) return 1.0/0.0;
	return distances[target];
	//std::cout << "Found " << found << std::endl;
	//std::cout << "Distance " << distances[target] << std::endl;
}

#endif // SWIG

// A hack to return stuff to python as the
// Boost.Geometry Point is too bizarre
struct Point2d {
	real x; real y;
	Point2d(){}
	Point2d(Point& p) {
		x = p.get<0>();
		y = p.get<1>();
	}

	Point2d(real x, real y)
		:x(x), y(y){}
};

class CoordinateProjector {
	private:
	projPJ projector;
	projPJ wgs;

	public:
	CoordinateProjector(const char* proj_string)
			:projector(NULL), wgs(NULL) {
		projector = pj_init_plus(proj_string);
		if(!projector) {
			throw OsmReaderError("Couldn't initialize coordinate projection");
		}

		projPJ tmp_wgs = pj_init_plus("+init=epsg:4326");
		if(!tmp_wgs) {
			throw OsmReaderError("Couldn't initialize WGS coordinate system");
		}

		if(!(wgs = pj_latlong_from_proj(tmp_wgs))) {
			pj_free(tmp_wgs);
			throw OsmReaderError("Couldn't initialize WGS latlong coordinate system");
		}
		pj_free(tmp_wgs);
	}

	void operator()(double latitude, double longitude, real* x, real* y) {
		*y = deg_to_rad(latitude);
		*x = deg_to_rad(longitude);
		
		if(pj_transform(wgs, projector, 1, 1, x, y, NULL)) {
			throw OsmReaderError("Failed to project a coordinate");
		}
	}

	/*void inverse(real* xy, Wgs84Point& wgs) {
		
	}*/

	~CoordinateProjector() {
		if(projector) pj_free(projector);
		if(wgs) pj_free(wgs);
	}
};



class OsmGraph {
	private:
	CoordinateProjector& coord_proj;

	public:
	Graph graph;
	std::unordered_map<node_id_t, Point> node_coordinates;
	std::unordered_map<node_id_t, Vertex> id_to_vertex;
	bgi::rtree< std::pair<Bbox, std::pair<Edge, LineSegment> >, bgi::quadratic<16> > edge_index;
	
	private:
	static int handle_osm_node(const void *usr_data, const readosm_node *node) {
		auto self = (OsmGraph*)usr_data;

		real x, y;
		self->coord_proj(node->latitude, node->longitude, &x, &y);
		self->node_coordinates[node->id] = Point(x, y);
		
				
		return READOSM_OK;
	}

	// TODO: Make configurable
	static WayRole get_way_role(const readosm_way *way) {
		const char *highway = NULL;
		const char *busway = NULL;
		const char *oneway = NULL;
		const char *junction = NULL;
		
		for(size_t i=0; i < way->tag_count; ++i) {
			if(string("highway").compare(way->tags[i].key) == 0)
				highway = way->tags[i].value;
			if(string("busway").compare(way->tags[i].key) == 0)
				busway = way->tags[i].value;
			if(string("oneway").compare(way->tags[i].key) == 0)
				oneway = way->tags[i].value;
			if(string("junction").compare(way->tags[i].key) == 0)
				junction = way->tags[i].value;
		}
		
		if(!(busway || highway)) return IgnoreWay;
		if(!highway) return TwoWayWay;

		if(string("footway").compare(highway) == 0) return IgnoreWay;
		if(string("cycleway").compare(highway) == 0) return IgnoreWay;
		if(string("steps").compare(highway) == 0) return IgnoreWay;
		
		if(string(highway).compare("motorway") == 0) return OneWayWay;
		if(oneway && string(oneway).compare("yes") == 0) return OneWayWay;
		if(junction && string(junction).compare("roundabout") == 0) return OneWayWay;
		
		return TwoWayWay;
	}

	void ensure_node(node_id_t node_id) {
		if(id_to_vertex.count(node_id)) return;
		if(node_coordinates.count(node_id) == 0) return;

		auto node_id_prop = bst::get(bst::vertex_name, graph);
		Vertex v = bst::add_vertex(graph);
		id_to_vertex[node_id] = v;
		node_id_prop[v] = node_id;
	}
	
	void add_new_edge(node_id_t src, node_id_t dst) {
		ensure_node(src); ensure_node(dst);
		if(!(id_to_vertex.count(src) && id_to_vertex.count(dst))) {
			return;
		}
		auto& s = node_coordinates[src];
		auto& e = node_coordinates[dst];
		auto seg = LineSegment(s, e);
		
		auto new_edge = bst::add_edge(
			id_to_vertex[src],
			id_to_vertex[dst],
			graph).first;
		
		auto length = bg::length(seg);
		auto edge_id = num_edges(graph);
		
		get(bst::edge_weight, graph)[new_edge] = length;
		Bbox bbox;
		bg::envelope(seg, bbox);
		auto entry = std::make_pair(bbox, std::make_pair(new_edge, seg));
		edge_index.insert(entry);

	}

	static int handle_osm_way(const void *usr_data, const readosm_way *way) {
		// TODO: The graph search could be optimized quite a bit by
		//	just inserting one edge per way and calculating the
		//	edge length from the segments.
		//	Also duplicating both ways in the RTree causes almost
		//	double the storage/computation. The implementation is
		//	simpler this way though.
		auto self = (OsmGraph*)usr_data;
		auto& i2v = self->id_to_vertex;
		WayRole role = get_way_role(way);
		if(role == IgnoreWay) {
			return READOSM_OK;
		}

		for(size_t i=0; i < way->node_ref_count - 1; ++i) {
			self->add_new_edge(way->node_refs[i], way->node_refs[i+1]);
		}

		if(role == OneWayWay) return READOSM_OK;
		
		for(size_t i=0; i < way->node_ref_count - 1; ++i) {
			self->add_new_edge(way->node_refs[i+1], way->node_refs[i]);
		}


		return READOSM_OK;
	}
	
	public:
	OsmGraph(const char *filename, CoordinateProjector& proj)
			:coord_proj(proj) {
		const void *osm_handle;
		auto status = readosm_open(filename, &osm_handle);
		if(status != READOSM_OK) {
			throw new OsmReaderError("Failed to open input file");
		}
		// NOTE: Assumes currently that all nodes are read before
		//	the ways. Seems to be so, but not really assured
		//	anywhere.
		readosm_parse(osm_handle, this, handle_osm_node, handle_osm_way, NULL);
		readosm_close(osm_handle);
	}

	Point& get_vertex_point(Vertex vertex) {
		auto node_id = get(bst::vertex_name, graph)[vertex];
		return node_coordinates[node_id];
	}
	
	std::vector< std::pair<Point2d, Point2d> > get_edge_coordinates() {
		vector< std::pair<Point2d, Point2d> > result;
		auto iters = bst::edges(graph);
		for(auto edge=iters.first; edge != iters.second; edge++) {
			auto src = Point2d(get_vertex_point(bst::source(*edge, graph)));
			auto dst = Point2d(get_vertex_point(bst::target(*edge, graph)));
			result.push_back(std::make_pair(src, dst));
		}

		return result;
	}

	~OsmGraph() {
	}
};



struct PositionHypothesis {
	shared_ptr<PositionHypothesis> parent;
	real timestamp;
	real total_likelihood;
	real measurement_likelihood;
	real transition_likelihood;
	real measurement_error;
	Edge edge;
	real edge_offset;
	vector<Vertex> subpath; 
};

real gaussian_logpdf(real x, real m, real s) {
	real normer = std::log(1.0/(s*std::sqrt(2*M_PI)));
	return normer - (x-m)*(x-m)/(2*s*s);
}

template <class SuccessorMap>
vector<Vertex> build_path_from_successors(Vertex start, Graph &g, SuccessorMap& suc) {
	vector<Vertex> result;

	auto end = suc.end();
	auto current = suc.find(start);
	while(current != end) {
		result.push_back(current->first);
		current = suc.find(current->second);
	}

	return result;
}

class MapMatcher2d {
	// TODO: Probably doing a lot of large object copies. Let's hope
	//	for copy elision, and of course optimize later.
	OsmGraph& graph;
	real search_radius;
	vector< shared_ptr<PositionHypothesis> >* hypotheses = NULL;
	Point previous_measurement;

	public:
	MapMatcher2d(OsmGraph& g, real search_radius=100.0)
		: graph(g), search_radius(search_radius)  {

	}

	~MapMatcher2d() {
		if(hypotheses) delete hypotheses;
	}

	// TODO: Take as template parameter
	// TODO: Give measurement and target
	real measurement_likelihood(real error) {
		// TODO: Distances aren't distributed like this
		return gaussian_logpdf(error, 0.0, 30.0);
	}
	real transition_likelihood(real dist) {
		return gaussian_logpdf(dist, 0.0, 10.0);
	}

	void measurement(real ts, real x, real y) {
		Bbox search_area(
			Point(x - search_radius, y - search_radius),
			Point(x + search_radius, y + search_radius));

		auto query = bgi::intersects(search_area);
		std::vector< std::pair<Bbox, std::pair<Edge, LineSegment> > > results;
		graph.edge_index.query(query, std::back_inserter(results));
		
		auto point = Point(x, y);
		if(!hypotheses) {
			previous_measurement = point;
		}
		auto measured_dist = bg::length(LineSegment(point, previous_measurement));
		previous_measurement = point;

		auto new_hypotheses = new vector< shared_ptr<PositionHypothesis> >;

		for(auto result: results) {
			auto edge = result.second.first;
			auto seg = result.second.second;
			auto proj = lineseg_point_projection(seg, Point(x, y));
			auto error = proj.first;
			auto t = proj.second;
			if(error > search_radius) continue;
			
			auto hypo = shared_ptr<PositionHypothesis>(new PositionHypothesis {
				.timestamp=ts,
				.edge=edge,
				.edge_offset=t,
				.measurement_error=error,
				.measurement_likelihood=measurement_likelihood(error),
				.parent=NULL
				});
			hypo->total_likelihood = hypo->measurement_likelihood;
			if(!hypotheses) {
				// First round, so don't search for
				// parents
				new_hypotheses->push_back(hypo);
				continue;
			}
			
			std::unordered_map<Vertex, vector<shared_ptr<PositionHypothesis>> > targets;
			for(auto prev: *hypotheses) {
				targets[bst::target(prev->edge, graph.graph)].push_back(prev);
			}
			
			shared_ptr<PositionHypothesis> best_parent;
			auto best_likelihood = -1.0/0.0;
			auto best_transition = -1.0/0.0;
			vector<Vertex> best_path;
			
			auto visit = [&] (Vertex current, real dist, std::unordered_map<Vertex, Vertex>& successors) {

				// TODO: Make configurable!
				// TODO: With certain transition functions we could probably
				//	infer when no other hypothesis can win the current best,
				//	allowing for an "optimal" return and probably a lot earlier.
				if(dist > (10.0+measured_dist)*5.0) {
					return false;
				}

				auto it = targets.find(current);
				if(it == targets.end()) return true;
				auto path = build_path_from_successors(
						current, graph.graph, successors);

				for(auto parent: it->second) {
					auto parent_likelihood = parent->total_likelihood;
					// TODO! Take account the mid-edge hits!!
					auto transition = transition_likelihood(dist);
					auto likelihood = parent_likelihood + transition;
					if(likelihood > best_likelihood) {
						best_parent = parent;
						best_likelihood = likelihood;
						best_transition = transition;
						best_path = path;
					}
				
				}
				
				targets.erase(it);
				if(!targets.size()) {
					// All hypotheses found
					return false;
				}
					
				return true;
			};
		
			reverse_shortest_path(
				bst::source(hypo->edge, graph.graph), graph.graph,
				visit
				);
			if(!best_parent) {
				// Not reachable from parents, skipidiskip
				continue;
			}

			hypo->parent = best_parent;
			hypo->transition_likelihood = best_transition;
			hypo->total_likelihood += best_likelihood + best_transition;
			hypo->subpath = best_path;
			new_hypotheses->push_back(hypo);
		}

		if(new_hypotheses->size() == 0) {
			// Found no hypotheses! Let's hope it was an outlier
			// and ignore this round.
			delete new_hypotheses;
			return;
		}
		
		delete hypotheses;
		hypotheses = new_hypotheses;
	}
	
	void measurement(real ts, Point p) {
		measurement(ts, p.get<0>(), p.get<1>());
	}

	std::vector< Point2d > best_match_coordinates() {
		vector< Point2d > path;
		shared_ptr<PositionHypothesis> current;
		auto max_likelihood = -1.0/0.0;
		// Find the best current state
		for(auto hypo: *hypotheses) {
			if(hypo->total_likelihood < max_likelihood) continue;
			max_likelihood = hypo->total_likelihood;
			current = hypo;
		}
		
		std::list< shared_ptr<PositionHypothesis> > states;

		while(current) {
			states.push_front(current);
			current = current->parent;
		}
		
		for(auto current: states){
			for(auto vertex: current->subpath) {
				auto point = graph.get_vertex_point(vertex);
				Point2d coords(point);
				path.push_back(coords);
			}
			auto point = graph.get_vertex_point(
				bst::source(current->edge, graph.graph));
			path.push_back(Point2d(point));
			current = current->parent;
		}

		return path;
	}
};

