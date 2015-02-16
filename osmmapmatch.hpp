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

#endif
typedef long long node_id_t;
#ifndef SWIG

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

struct LinesegProjectionResult {
	Point projected;
	real error;
	real t;
};
// There's probably something like this in the depths
// of Boost.Geometry, but couldn't find it, so here we go again.
LinesegProjectionResult lineseg_point_projection(LineSegment seg, Point p)
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
	LinesegProjectionResult result = {
		.projected=proj,
		.error=error,
		.t=t*seglen};
	return result;
}

template <class Key, class Value>
class default_map
{	
	public:
	typedef std::pair<const Key, Value> value_type;
	typedef Key key_type;

	private:
	std::unordered_map<Key, Value, std::hash<Key>, std::equal_to<Key>> m_map;
	Value deflt;

	public:
	
	default_map(const Value def): deflt(def),
		m_map() {

	}

	Value& operator[](const Key& k) {
		auto it = m_map.find(k);
		if(it == m_map.end()) {
			return m_map[k] = deflt;
		}
		
		return it->second;
	}

	size_t count(const Key& k) {
		return m_map.count(k);
	}
};

template <class Callback>
void reverse_shortest_path(Vertex target, Graph& graph, Callback callback) {
	default_map<Vertex, real> distances(1.0/0.0);
	default_map<Vertex, real> cost(1.0/0.0);
	std::unordered_map<Vertex, Vertex> successor;
	auto weights = bst::get(bst::edge_weight, graph);

	typedef std::pair<real, Vertex> Node;
	std::priority_queue<Node, std::vector<Node>, std::greater<Node>> queue;
	queue.push(std::make_pair(0.0, target));
	
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

template <class SuccessorMap>
vector<Vertex> build_path_from_successors(Vertex start, Graph &g, SuccessorMap& suc) {
	vector<Vertex> result;

	auto end = suc.end();
	auto current = suc.find(start);
	result.push_back(start);
	while(current != end) {
		result.push_back(current->second);
		current = suc.find(current->second);
	}

	return result;
}

auto single_target_shortest_path = [](Vertex source, Vertex target, Graph& graph) {
	vector<Vertex> path;
	auto visitor = [&](Vertex current, real distance, 
		std::unordered_map<Vertex, Vertex>& successors) {
		if(current != source) return true;
		path = build_path_from_successors(source, graph, successors);
		return false;
	};

	reverse_shortest_path(target, graph, visitor);

	return path;
};


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
			throw OsmReaderError("Failed to project a coordinate to XY " + std::to_string(latitude) + " " + std::to_string(longitude));
		}
	}
	
	void inverse(real x, real y, double *latitude, double *longitude) {
		*latitude = y;
		*longitude = x;
		
		if(pj_transform(projector, wgs, 1, 1, longitude, latitude, NULL)) {
			throw OsmReaderError("Failed to project a coordinate to LatLng "  + std::to_string(x) + " " + std::to_string(y));
		}

		*latitude = rad_to_deg(*latitude);
		*longitude = rad_to_deg(*longitude);
	}

	~CoordinateProjector() {
		if(projector) pj_free(projector);
		if(wgs) pj_free(wgs);
	}
};

enum WayRole {
	WayRoleIgnore = 0,
	WayRoleOneWay,
	WayRoleTwoWay
};

WayRole busway_filter(const readosm_way *way) {
	const char *highway = NULL;
	const char *busway = NULL;
	const char *oneway = NULL;
	const char *junction = NULL;
	const char *ferry = NULL;

	for(size_t i=0; i < way->tag_count; ++i) {
		if(string("highway").compare(way->tags[i].key) == 0)
			highway = way->tags[i].value;
		if(string("busway").compare(way->tags[i].key) == 0)
			busway = way->tags[i].value;
		if(string("oneway").compare(way->tags[i].key) == 0)
			oneway = way->tags[i].value;
		if(string("junction").compare(way->tags[i].key) == 0)
			junction = way->tags[i].value;
		if(string("ferry").compare(way->tags[i].key) == 0)
			ferry = way->tags[i].value;
	}
		
	if(!(busway || highway || ferry)) return WayRoleIgnore;
	if(!highway) return WayRoleTwoWay;

	if(string("footway").compare(highway) == 0) return WayRoleIgnore;
	if(string("cycleway").compare(highway) == 0) return WayRoleIgnore;
	if(string("steps").compare(highway) == 0) return WayRoleIgnore;
	if(string("path").compare(highway) == 0) return WayRoleIgnore;
	if(string("construction").compare(highway) == 0) return WayRoleIgnore;
	if(string("proposed").compare(highway) == 0) return WayRoleIgnore;
	if(string("bridleway").compare(highway) == 0) return WayRoleIgnore;

	
	if(string(highway).compare("motorway") == 0) return WayRoleOneWay;
	if(oneway && string(oneway).compare("yes") == 0) return WayRoleOneWay;
	if(junction && string(junction).compare("roundabout") == 0) return WayRoleOneWay;
		
	return WayRoleTwoWay;
}

WayRole train_filter(const readosm_way *way) {
	const char *railway = NULL;
	const char *oneway = NULL;
	for(size_t i=0; i < way->tag_count; ++i) {
		if(string("railway").compare(way->tags[i].key) == 0)
			railway = way->tags[i].value;
		if(string("oneway").compare(way->tags[i].key) == 0)
			oneway = way->tags[i].value;
	};

	if(!railway) return WayRoleIgnore;
	if(string("rail").compare(railway) != 0) return WayRoleIgnore;
	if(oneway && string(oneway).compare("yes") == 0) return WayRoleOneWay;
	return WayRoleTwoWay;
}

WayRole tram_filter(const readosm_way *way) {
	const char *railway = NULL;
	const char *oneway = NULL;
	for(size_t i=0; i < way->tag_count; ++i) {
		if(string("railway").compare(way->tags[i].key) == 0)
			railway = way->tags[i].value;
		if(string("oneway").compare(way->tags[i].key) == 0)
			oneway = way->tags[i].value;
	};

	if(!railway) return WayRoleIgnore;
	if(string("tram").compare(railway) != 0) return WayRoleIgnore;
	if(oneway && string(oneway).compare("yes") == 0) return WayRoleOneWay;
	return WayRoleTwoWay;
}

WayRole subway_filter(const readosm_way *way) {
	const char *railway = NULL;
	for(size_t i=0; i < way->tag_count; ++i) {
		if(string("railway").compare(way->tags[i].key) == 0)
			railway = way->tags[i].value;
	};

	if(!railway) return WayRoleIgnore;
	if(string("subway").compare(railway) != 0) return WayRoleIgnore;
	return WayRoleTwoWay;
}

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
		auto self = (OsmGraph*)usr_data;
		return self->do_handle_osm_way(way);
	}

	int do_handle_osm_way(const readosm_way *way) {
		// TODO: The graph search could be optimized quite a bit by
		//	just inserting one edge per way and calculating the
		//	edge length from the segments.
		//	Also duplicating both ways in the RTree causes almost
		//	double the storage/computation. The implementation is
		//	simpler this way though.
		auto& i2v = id_to_vertex;
		WayRole role = get_way_role(way);
		if(role == WayRoleIgnore) {
			return READOSM_OK;
		}

		for(size_t i=0; i < way->node_ref_count - 1; ++i) {
			add_new_edge(way->node_refs[i], way->node_refs[i+1]);
		}

		if(role == WayRoleOneWay) return READOSM_OK;
		
		for(size_t i=0; i < way->node_ref_count - 1; ++i) {
			add_new_edge(way->node_refs[i+1], way->node_refs[i]);
		}


		return READOSM_OK;
	}
	
	WayRole(&get_way_role)(const readosm_way *);

	public:
	OsmGraph(const char *filename, CoordinateProjector& proj, WayRole(*get_way_role)(const readosm_way *))
			:coord_proj(proj), get_way_role(*get_way_role) {
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

	node_id_t get_vertex_id(Vertex vertex) {
		return get(bst::vertex_name, graph)[vertex];
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
	std::shared_ptr<PositionHypothesis> parent;
	real timestamp;
	real total_likelihood;
	real measurement_likelihood;
	real transition_likelihood;
 	real measurement_error;
	Edge edge;
	real edge_offset;
	std::vector<Vertex> subpath;
	Point position;
	Point measurement;
};

real gaussian_logpdf(real x, real m, real s) {
	real normer = std::log(1.0/(s*std::sqrt(2*M_PI)));
	return normer - (x-m)*(x-m)/(2*s*s);
}

class StateLikelihoodModel {
	public:
	virtual real measurement(const Point& m, const Point& p) = 0;
	virtual real transition(const PositionHypothesis& parent, const PositionHypothesis& next, const vector<Vertex>& path, real path_length) = 0;
	virtual ~StateLikelihoodModel() {};

	virtual real best_transition_still_possible(real measured_dist, real distance) {
		return 0.0;
	}
};

inline real vector_angle(Point& a, Point &b) {
	auto cosangle = bg::dot_product(a, b)/(vect_norm(a)*vect_norm(b));
	// Avoid acos returning NaNs due to rounding errors
	if(cosangle < -1.0) {
		cosangle = -1.0;
	} else if(cosangle > 1.0) {
		cosangle = 1.0;
	}
	auto angle = std::acos(cosangle);
	return angle;
}

#if defined __FAST_MATH__
#error This will not currently work with fast math due to needed NaN checks
#endif
class DrawnGaussianStateModel : public StateLikelihoodModel {
	real measurement_std;
	real length_error_std;
	OsmGraph& graph;
	public:

	DrawnGaussianStateModel(real measurement_std, real length_error_std, OsmGraph& graph)
		:measurement_std(measurement_std), graph(graph), length_error_std(length_error_std) {}

	real measurement(const Point& m, const Point& p) {
		auto dx = m.get<0>() - p.get<0>();
		auto dy = m.get<1>() - p.get<1>();
		return gaussian_logpdf(dx, 0.0, measurement_std) +
			gaussian_logpdf(dy, 0.0, measurement_std);
	}

	real transition(const PositionHypothesis& parent, const PositionHypothesis& next, const vector<Vertex> &path, real path_length) {
		// The idea here is to calculate how big a share of the
		// path is in different direction than the measurement.
		// TODO: I have a hunch that this could be reduced to a very
		//	simple form of just some distance ratios.
		auto orig_direction = next.measurement;
		bg::subtract_point(orig_direction, parent.measurement);
		
		auto total_angle = 0.0;
		auto total_length = 0.0;
		int n_spans = 0;
		
		auto evaluate_span = [&](Point& prev_point, Point& next_point) {
			auto direction = next_point;
			bg::subtract_point(direction, prev_point);
			auto length = vect_norm(direction);
			prev_point = next_point;
			auto angle = vector_angle(orig_direction, direction);
			if(std::isnan(angle)) {
				return;
			}
			total_angle += angle/M_PI*length;
			total_length += length;
			n_spans++;

		};

		auto prev_point = parent.position;
		for(auto i=0; i < path.size(); i++) {
			auto next_point = graph.get_vertex_point(path[i]);
			evaluate_span(prev_point, next_point);
			prev_point = next_point;
		}

		if(path.size() > 0) {
			auto next_point = graph.get_vertex_point(bst::source(next.edge, graph.graph));
			evaluate_span(prev_point, next_point);
			prev_point = next_point;
		}
		
		auto next_point = next.position;
		evaluate_span(prev_point, next_point);
		
		if(total_length < 1e-6) {
			// Give some penalty for hanging around in the same
			// node. Otherwise truncates start and endpoints.
			// TODO: Could probably be handled more elegantly.
			total_angle = length_error_std;
		} else {
			total_angle /= total_length;
		}
		
		return gaussian_logpdf(total_angle, 0.0, length_error_std);
		
		//auto measured_length = bg::length(LineSegment(parent.measurement, next.measurement));
		//auto relative_dist = (path_length+length_error_std)/(measured_length+length_error_std);
		//return gaussian_logpdf(relative_dist-1.0, 0, 0.1);
		//return gaussian_logpdf(path_length, 0.0, length_error_std);
	}

	real best_transition_still_possible(real measured_length, real path_length) {
		return gaussian_logpdf(0.0, 0.0, length_error_std);
		//auto diff = path_length - measured_length;
		//if(diff < 0.0) diff = 0.0;
		//return gaussian_logpdf(path_length/(measured_length+1.0), 0, length_error_std);
		//return gaussian_logpdf(path_length, 0.0, length_error_std);
	}
};

class MapMatcher2d {
	// TODO: Probably doing a lot of large object copies. Let's hope
	//	for copy elision, and of course optimize later.
	// TODO: There's a huge problem in the approach! By examining just
	//	a single point on the edge, if the noise causes the measurement
	//	to "backtrack", the result is very bad! Should consider the whole
	//	edge instead.
	OsmGraph& graph;
	real search_radius;
	vector< shared_ptr<PositionHypothesis> >* hypotheses = NULL;
	Point previous_measurement;
	StateLikelihoodModel& state_model;

	public:
	int n_outliers = 0;

	MapMatcher2d(OsmGraph& g, StateLikelihoodModel& state_model, real search_radius=100.0)
		: graph(g), search_radius(search_radius), state_model(state_model)  {

	}

	~MapMatcher2d() {
		if(hypotheses) delete hypotheses;
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
		auto best_total_likelihood = -1.0/0.0;
		
		#pragma omp parallel for
		for(auto i=results.begin(); i < results.end(); i++) {
			auto& result = *i;
			auto edge = result.second.first;
			auto seg = result.second.second;
			auto proj = lineseg_point_projection(seg, Point(x, y));
			auto error = proj.error;
			auto t = proj.t;
			if(error > search_radius) {
				continue;
			}
			
			auto hypo = shared_ptr<PositionHypothesis>(new PositionHypothesis {
				.parent=NULL,
				.timestamp=ts,
				.total_likelihood=0.0,
				.measurement_likelihood=state_model.measurement(point, proj.projected),
				.transition_likelihood=0.0,
				.measurement_error=error,
				.edge=edge,
				.edge_offset=t,
				.subpath=vector<Vertex>(),
				.position=proj.projected,
				.measurement=point
				});
			hypo->total_likelihood = hypo->measurement_likelihood;
			if(!hypotheses) {
				// First round, so don't search for
				// parents
				#pragma omp critical
				{
				new_hypotheses->push_back(hypo);
				}
				continue;
			}
			
						
			shared_ptr<PositionHypothesis> best_parent;
			auto best_likelihood = -1.0/0.0;
			auto best_transition = -1.0/0.0;
			vector<Vertex> best_path;
			
			/* No polymorphic lambdas in C++11 yet :( */
			auto consider_parent = [&] (shared_ptr<PositionHypothesis>& parent, real dist,
					vector<Vertex>& path ) {
				auto parent_likelihood = parent->total_likelihood;
				auto transition = state_model.transition(*parent, *hypo, path, dist);
				auto likelihood = parent_likelihood + transition;
				if(likelihood > best_likelihood) {
					best_parent = parent;
					best_likelihood = likelihood;
					best_transition = transition;
					best_path = path;
				}
			};
			
			std::unordered_map<Vertex, vector<shared_ptr<PositionHypothesis>> > targets;
			for(auto& prev: *hypotheses) {
				if(prev->edge == hypo->edge) {
					// The graph search doesn't work in the special case
					// where the parent and hypo are the same edge.
					auto real_dist = hypo->edge_offset - prev->edge_offset;
					if(real_dist < 0.0) {
						// Won't go into reverse direction. This is a prime example
						// of the problem that considering only single point in the
						// edge causes!
						continue;
					}
					// No path as it's the same edge
					vector<Vertex> path;
					consider_parent(prev, real_dist, path);
					continue;
				}
				targets[bst::target(prev->edge, graph.graph)].push_back(prev);
			}

			auto visit = [&] (Vertex current, real dist, std::unordered_map<Vertex, Vertex>& successors) {
				dist += hypo->edge_offset;
				
				
				// TODO: Make configurable!
				// TODO: With certain transition likelihood functions we could probably
				//	infer when no other hypothesis can win the current best,
				//	allowing for an "optimal" return and probably a lot earlier.
				if(dist > (10.0+measured_dist)*5.0) {
					//std::cout << "Hit the limit :(" << std::endl;
					return false;
				}

				auto it = targets.find(current);
				if(it == targets.end()) return true;
				auto hypo = it->second;
				auto path = build_path_from_successors(
						current, graph.graph, successors);

				for(auto& parent: it->second) {
					auto parent_left = get(bst::edge_weight, graph.graph)[parent->edge];
					parent_left -= parent->edge_offset;
					auto real_dist = dist + parent_left;
					consider_parent(parent, real_dist, path);
				}
				
				targets.erase(it);

				
				if(!targets.size()) {
					// All hypotheses found
					return false;
				}
				
				// Prune parents that can't win
				// TODO: This would prune a lot more if done
				//	for each found node and track the current
				//	worst.
				auto itr = targets.begin();
				while(itr != targets.end()) {
					for(auto& parent: itr->second) {
						auto parent_left = get(bst::edge_weight, graph.graph)[parent->edge];
						parent_left -= parent->edge_offset;
						auto real_dist = dist + parent_left;
						auto best_possible = state_model.best_transition_still_possible(
							measured_dist, real_dist);
						best_possible += parent->total_likelihood;
						if(best_possible >= best_likelihood) {
							goto dontremove;
						}
					}
					itr = targets.erase(itr);
					if(targets.size() == 0) {
						return false;
					}
					continue;

					dontremove:
					++itr;
				}
					
				return true;
			};
		
			reverse_shortest_path(
				bst::source(hypo->edge, graph.graph), graph.graph,
				visit
				);
			if(!best_parent) {
				// Not reachable from parents, skipidiskip
				// TODO: Getting here is very expensive. A nicer
				// way would be to simultaneously do multi-target
				// multi-source search (or just parallelise the different searches)
				// and give up when the search starts to look too bad.
				continue;
			}

			hypo->parent = best_parent;
			hypo->transition_likelihood = best_transition;
			hypo->total_likelihood += best_likelihood + best_transition;
			hypo->subpath = best_path;
			#pragma omp critical
			{
			new_hypotheses->push_back(hypo);
			if(hypo->total_likelihood > best_total_likelihood) {
				best_total_likelihood = 0.0;
			}
			}
		}
		
		if(new_hypotheses->size() == 0) {
			// Found no hypotheses! Let's hope it was an outlier
			// and ignore this round.
			delete new_hypotheses;
			n_outliers += 1;
			return;
		}
		
		delete hypotheses;
		hypotheses = new_hypotheses;
	}
	
	void measurement(real ts, Point p) {
		measurement(ts, p.get<0>(), p.get<1>());
	}
	
	
	void measurements(const std::vector<real>& ts, const std::vector<Point2d>& points) {
		auto s = ts.size();
		for(int i = 0; i < s; i++) {
			measurement(ts[i], points[i].x, points[i].y);
		}
	}

	std::list<Vertex> route_vertex_path(std::list<shared_ptr<PositionHypothesis> >& path) {
		std::list<Vertex> vertices;
		
		Vertex prev_inserted;
		auto insert_unique = [&](Vertex vertex) {
			if(vertices.size() and prev_inserted == vertex) return;
			vertices.push_back(vertex);
			prev_inserted = vertex;
		};

		for(auto current: path) {
			if(current->subpath.size() == 0) {
				// If the subpath length is zero, we are
				// on the same edge as the parent, so
				// don't add the source.
				continue;
			}
			for(auto vertex: current->subpath) {
				insert_unique(vertex);
			}
			insert_unique(bst::source(current->edge, graph.graph));
		}

		return vertices;
	}

	std::list<shared_ptr<PositionHypothesis> > get_hypothesis_path(shared_ptr<PositionHypothesis> current) {
		std::list< shared_ptr<PositionHypothesis> > states;
		while(current) {
			states.push_front(current);
			current = current->parent;
		}

		return states;
	}

	std::vector< Point2d > best_match_coordinates() {
		vector< Point2d > path;
		auto current = best_current_hypothesis();		
		auto states = get_hypothesis_path(current);
		path.push_back(states.front()->position);
		states.pop_front();

		for(auto vertex: route_vertex_path(states)) {
			auto point = graph.get_vertex_point(vertex);
			path.push_back(point);
		}
		path.push_back(states.back()->position);

		return path;
	}

	std::vector< node_id_t > best_match_node_ids() {
		vector< node_id_t > path;
		auto current = best_current_hypothesis();
		auto states = get_hypothesis_path(current);
		// OSM id's have to be >= 1, so hacking
		// zero to be missing.
		path.push_back(0);
		states.pop_front();

		for(auto vertex: route_vertex_path(states)) {
			auto point = graph.get_vertex_id(vertex);
			path.push_back(point);
		}
		path.push_back(0);

		return path;
	}

	std::shared_ptr<PositionHypothesis> best_current_hypothesis() {
		shared_ptr<PositionHypothesis> current;
		auto max_likelihood = -1.0/0.0;
		// Find the best current state
		for(auto hypo: *hypotheses) {
			if(hypo->total_likelihood < max_likelihood) continue;
			max_likelihood = hypo->total_likelihood;
			current = hypo;
		}

		return current;
	}
	
	std::vector<std::shared_ptr<PositionHypothesis> > current_hypotheses() {
		return *hypotheses;
	}
	
};

template<class Rndgen>
std::vector<Point2d> get_random_path_custom(OsmGraph& graph, int n_waypoints, Rndgen& gen) {

	std::vector<Point2d> result;

	auto src = bst::random_vertex(graph.graph, gen);
	for(int i=0; i < n_waypoints; ++i) {
		auto dst = bst::random_vertex(graph.graph, gen);
		auto path = single_target_shortest_path(src, dst, graph.graph);
		if(!path.size()) {
			// Recurse until a path is found
			return get_random_path_custom(graph, n_waypoints, gen);
		}

		for(auto node: path) {
			result.push_back(graph.get_vertex_point(node));
		}

		src = dst;
	}

	return result;
}

std::vector<Point2d> get_random_path(OsmGraph& graph, int n_waypoints=1) {
	std::random_device rd;
	std::mt19937 gen(rd());
	return get_random_path_custom(graph, n_waypoints, gen);
}

std::vector<Point2d> get_shortest_node_path(OsmGraph& graph, node_id_t start, node_id_t end) {
	std::vector<Point2d> result;
	auto s = graph.id_to_vertex[start];
	auto e = graph.id_to_vertex[end];
	for(auto node: single_target_shortest_path(s, e, graph.graph)) {
		result.push_back(graph.get_vertex_point(node));
	}

	return result;
};
