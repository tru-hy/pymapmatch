#include "osmmapmatch.hpp"

int main(int argc, char **argv) {
	CoordinateProjector proj("+init=epsg:3067");
	OsmGraph graph(argv[1], proj, busway_filter);
	
	//std::random_device rd;
	std::mt19937 gen(42);

	int n_iters = 5;
	for(int i = 0; i < n_iters; i++) {
		DrawnGaussianStateModel model(30, 0.1, graph);
		MapMatcher2d matcher(graph, model);
		auto path = get_random_path_custom(graph, 5, gen);
		for(auto coord: path) {
			// TODO: Sparseness and random noise
			matcher.measurement(0, coord.x, coord.y);
		}
		std::cout << i << std::endl;
	}
	
	/*
	auto bounds = graph.edge_index.bounds();
	std::uniform_real_distribution<> xgen(
		bounds.min_corner().get<0>(),
		bounds.max_corner().get<0>()
		);
	
	std::uniform_real_distribution<> ygen(
		bounds.min_corner().get<1>(),
		bounds.max_corner().get<1>()
		);
	
	auto n = 10;
	auto t = std::clock();
	for(int i=0; i < n; ++i) {
		real x = xgen(gen);
		real y = ygen(gen);
		
		std::cout << "Running measurement" << std::endl;
		matcher.measurement(0, x, y);
		std::cout << "Done" << std::endl;
	}

	std::cout << double(std::clock() - t)/double(CLOCKS_PER_SEC)/double(n)*1000.0 << " ms per measurement" << std::endl;
	
	//auto nodes = get(std::vertex_name, graph.graph);
	//std::cout << nodes.size() << std::endl;
	t = std::clock();
	n = 1000;
	for(int i=0; i < n; ++i) {
		auto source = bst::random_vertex(graph.graph, gen);
		auto target = bst::random_vertex(graph.graph, gen);
		//auto dist_true = single_target_shortest_path_slow(source, target, graph.graph);
		auto dist_guess = single_target_shortest_path(source, target, graph.graph, 1000.0);
		//std::cout << dist_true << " " << ((dist_true == dist_guess)?"win":"fail") << std::endl;
	}
	std::cout << double(std::clock() - t)/double(CLOCKS_PER_SEC)/double(n)*1000.0 << " ms per dijkstra query" << std::endl;
	//std::cout << random_vertex(graph.graph, gen) << std::endl;
	//string tmp; std::getline(std::cin, tmp);
	*/
}
