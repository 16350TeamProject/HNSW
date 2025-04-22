/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include <random>
#include <vector>
#include <array>
#include <algorithm>
#include <utility> 
#include <queue>

#include <tuple>
#include <string>
#include <stdexcept>
#include <regex> // For regex and split logic
#include <iostream> // cout, endl
#include <fstream> // For reading/writing files
#include <assert.h> 
#include <unordered_set>
#include <cfloat>
#include <ctime>

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define PRM         0
#define RRT         1
#define PRM_HNSW    2
#define RRT_HNSW    3
#define PRM_KDTree  4
#define RRT_KDTree  5

/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

//the length of each link in the arm
#define LINKLENGTH_CELLS 10

// Some potentially helpful imports
using std::vector;
using std::array;
using std::string;
using std::runtime_error;
using std::tuple;
using std::make_tuple;
using std::tie;
using std::cout;
using std::endl;
using std::swap;

#include "hnswlib/hnswlib/hnswlib.h"
#include "nanoflann/include/nanoflann.hpp"
using namespace std;
using namespace nanoflann;

double PI=3.141592654;
/// @brief 
/// @param filepath 
/// @return map, x_size, y_size
tuple<double*, int, int> loadMap(string filepath) {
	std::FILE *f = fopen(filepath.c_str(), "r");
	if (f) {
	}
	else {
		printf("Opening file failed! \n");
		throw runtime_error("Opening map file failed!");
	}
	int height, width;
	if (fscanf(f, "height %d\nwidth %d\n", &height, &width) != 2) {
		throw runtime_error("Invalid loadMap parsing map metadata");
	}
	
	////// Go through file and add to m_occupancy
	double* map = new double[height*width];

	double cx, cy, cz;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			char c;
			do {
				if (fscanf(f, "%c", &c) != 1) {
					throw runtime_error("Invalid parsing individual map data");
				}
			} while (isspace(c));
			if (!(c == '0')) { 
				map[y+x*width] = 1; // Note transposed from visual
			} else {
				map[y+x*width] = 0;
			}
		}
	}
	fclose(f);
	return make_tuple(map, width, height);
}

// Splits string based on deliminator
vector<string> split(const string& str, const string& delim) {   
		// https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c/64886763#64886763
		const std::regex ws_re(delim);
		return { std::sregex_token_iterator(str.begin(), str.end(), ws_re, -1), std::sregex_token_iterator() };
}


double* doubleArrayFromString(string str) {
	vector<string> vals = split(str, ",");
	double* ans = new double[vals.size()];
	for (int i = 0; i < vals.size(); ++i) {
		ans[i] = std::stod(vals[i]);
	}
	return ans;
}

bool equalDoubleArrays(double* v1, double *v2, int size) {
    for (int i = 0; i < size; ++i) {
        if (abs(v1[i]-v2[i]) > 1e-3) {
            cout << endl;
            return false;
        }
    }
    return true;
}

typedef struct {
	int X1, Y1;
	int X2, Y2;
	int Increment;
	int UsingYIndex;
	int DeltaX, DeltaY;
	int DTerm;
	int IncrE, IncrNE;
	int XIndex, YIndex;
	int Flipped;
} bresenham_param_t;


void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size) {
	double cellsize = 1.0;
	//take the nearest cell
	*pX = (int)(x/(double)(cellsize));
	if( x < 0) *pX = 0;
	if( *pX >= x_size) *pX = x_size-1;

	*pY = (int)(y/(double)(cellsize));
	if( y < 0) *pY = 0;
	if( *pY >= y_size) *pY = y_size-1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params) {
	params->UsingYIndex = 0;

	if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
		(params->UsingYIndex)++;

	if (params->UsingYIndex)
		{
			params->Y1=p1x;
			params->X1=p1y;
			params->Y2=p2x;
			params->X2=p2y;
		}
	else
		{
			params->X1=p1x;
			params->Y1=p1y;
			params->X2=p2x;
			params->Y2=p2y;
		}

	 if ((p2x - p1x) * (p2y - p1y) < 0)
		{
			params->Flipped = 1;
			params->Y1 = -params->Y1;
			params->Y2 = -params->Y2;
		}
	else
		params->Flipped = 0;

	if (params->X2 > params->X1)
		params->Increment = 1;
	else
		params->Increment = -1;

	params->DeltaX=params->X2-params->X1;
	params->DeltaY=params->Y2-params->Y1;

	params->IncrE=2*params->DeltaY*params->Increment;
	params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
	params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

	params->XIndex = params->X1;
	params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y) {
	if (params->UsingYIndex) {
        *y = params->XIndex;
        *x = params->YIndex;
        if (params->Flipped)
            *x = -*x;
    }
	else {
        *x = params->XIndex;
        *y = params->YIndex;
        if (params->Flipped)
            *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params) {
	if (params->XIndex == params->X2) {
        return 0;
    }
	params->XIndex += params->Increment;
	if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
		params->DTerm += params->IncrE;
	else {
        params->DTerm += params->IncrNE;
        params->YIndex += params->Increment;
	}
	return 1;
}



int IsValidLineSegment(double x0, double y0, double x1, double y1, double*	map,
			 int x_size, int y_size) {
	bresenham_param_t params;
	int nX, nY; 
	short unsigned int nX0, nY0, nX1, nY1;

	//printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
		
	//make sure the line segment is inside the environment
	if(x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

	//printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	//iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do {
		get_current_point(&params, &nX, &nY);
		if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
			return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double*	map,
			 int x_size, int y_size) {
    double x0,y0,x1,y1;
    int i;
		
	 //iterate through all the links starting with the base
	x1 = ((double)x_size)/2.0;
	y1 = 0;
	for(i = 0; i < numofDOFs; i++){
		//compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
		y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

		//check the validity of the corresponding line segment
		if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
			return 0;
	}    
	return 1;
}

// Compute the Euclidean distance between two configurations
double computeDistance(const vector<double>& a, const vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        sum += pow(a[i] - b[i], 2);
    }
    return sqrt(sum);
}

// PRM Planner
static void PRMPlanner(
	double* map,
	int x_size,
	int y_size,
	double* armstart_anglesV_rad,
	double* armgoal_anglesV_rad,
	int numofDOFs,
	double*** plan,
	int* planlength)
{
	const int NUM_SAMPLES = 10000;   // sample size of random configurations
	const int K_NEAREST = 13;      // number of nearest neighbors to attempt to connect

    struct Node {
        vector<double> angles;
        vector<Node*> neighbors;
        Node* parent = nullptr;
        double cost = 0.0;
    };

	vector<Node*> roadmap;

	// check if a configuration is valid
	auto isValid = [&](const vector<double>& config) {
		return IsValidArmConfiguration(const_cast<double*>(config.data()), numofDOFs, map, x_size, y_size);
	};

	// generate random valid samples, so we search for a path in the selected space
	for (int i = 0; i < NUM_SAMPLES; i++) {
		vector<double> sample(numofDOFs);
		for (int j = 0; j < numofDOFs; j++) {
			sample[j] = ((double)rand() / RAND_MAX) * 2 * PI;  // random angle [0, 2π]
		}

		if (isValid(sample)) {
			Node* new_node = new Node{sample};
			roadmap.push_back(new_node);
		}
	}

	// add start and goal configurations to the roadmap
	Node* start_node = new Node{{armstart_anglesV_rad, armstart_anglesV_rad + numofDOFs}};
	Node* goal_node = new Node{{armgoal_anglesV_rad, armgoal_anglesV_rad + numofDOFs}};

	if (!isValid(start_node->angles) || !isValid(goal_node->angles)) {
		printf("Start or Goal is in collision!\n");
		return;
	}
	roadmap.push_back(start_node);
	roadmap.push_back(goal_node);

	// trying to connect to nearest neighbors
	for (auto& node : roadmap) {
		std::vector<std::pair<double, Node*>> neighbors;
		
		for (auto& other : roadmap) {
			if (node == other) continue;
			double dist = computeDistance(node->angles, other->angles);
			neighbors.emplace_back(dist, other);
		}

		// sort by distance and connect up to K_NEAREST valid paths
		sort(neighbors.begin(), neighbors.end());
		for (int i = 0; i < std::min(K_NEAREST, (int)neighbors.size()); i++) {
			if (isValid(neighbors[i].second->angles)) {
				node->neighbors.push_back(neighbors[i].second);
			}
		}
	}

	auto compare = [](const std::pair<double, Node*>& a, const std::pair<double, Node*>& b) {
		return a.first > b.first;
	};	

    std::priority_queue<
        std::pair<double, Node*>,
        vector<std::pair<double, Node*>>,
        decltype(compare)
    > pq(compare);

	std::unordered_set<Node*> visited;

	start_node->cost = 0.0;
	pq.emplace(0.0, start_node);

	bool found = false;

	// Dijkstra graph search to find the shortest path
	while (!pq.empty()) {
		Node* current = pq.top().second;
		pq.pop();

		// skip visited nodes
		if (visited.find(current) != visited.end())
        continue; 

		visited.insert(current);

		// find the goal node
		if (current == goal_node) {
			found = true;
			break;
		}

		for (Node* neighbor : current->neighbors) {
			double new_cost = current->cost + computeDistance(current->angles, neighbor->angles);
			if (visited.find(neighbor) == visited.end() && 
				(neighbor->parent == nullptr || new_cost < neighbor->cost)) {
				
				neighbor->cost = new_cost;
				neighbor->parent = current;
				pq.emplace(new_cost, neighbor);
			}
		}
	}

	if (!found) {
		printf("No valid path found!\n");
		return;
	}

	vector<vector<double>> path;

	// ensure the goal node is valid
	if (goal_node == nullptr) {
		printf("ERROR DETECTED: Goal node is NULL. No valid path found!\n");
		return;
	}
	if (goal_node->parent == nullptr) {
		printf("ERROR DETECTED: Goal node is not connected to the graph. No valid path found!\n");
		return;
	}
	
	// backtrack to extract the path
	int MAX_ITER = 10000;
	int iteration = 0;
	for (Node* node = goal_node; node != nullptr; node = node->parent) {
		if (iteration++ > MAX_ITER) {
			printf("ERROR DETECTED: Infinite loop detected while backtracking!\n");
			return;
		}
		path.push_back(node->angles);
		// printf("  Step %d: Node at %p\n", iteration, node);
	}
	
	reverse(path.begin(), path.end());
	// printf("DEBUG: Extracted path length = %d\n", (int)path.size());
	

	// store the path
	*planlength = path.size();

	//printf("DEBUG: Extracted path length = %d\n", *planlength);

	*plan = (double**) malloc(*planlength * sizeof(double*));

	for (int i = 0; i < *planlength; i++) {
		(*plan)[i] = (double*) malloc(numofDOFs * sizeof(double));
		memcpy((*plan)[i], path[i].data(), numofDOFs * sizeof(double));
	}

	printf("PRM successfully found a path with %d waypoints.\n", *planlength);

	for (Node* node : roadmap) {
		delete node;
	}
}


struct Node {
	std::vector<double> angles;
	std::vector<Node*> neighbors;
	Node* parent = nullptr;
	double cost = 0.0;
};

static void debugSummary(const std::vector<Node*>& roadmap, Node* start_node, Node* goal_node) {
	printf("[DEBUG] Total nodes in roadmap: %lu\n", roadmap.size());
	printf("[DEBUG] Start node neighbors: %lu\n", start_node->neighbors.size());
	printf("[DEBUG] Goal node neighbors: %lu\n", goal_node->neighbors.size());
	double dist = computeDistance(start_node->angles, goal_node->angles);
	printf("[DEBUG] Start-to-Goal Euclidean distance: %.6f\n", dist);

	int total_edges = 0;
	for (auto node : roadmap) total_edges += node->neighbors.size();
	printf("[DEBUG] Total roadmap edges: %d (avg %.2f per node)\n", total_edges, total_edges / (double)roadmap.size());
}

static void PRMHNSWPlanner(
	double* map,
	int x_size,
	int y_size,
	double* armstart_anglesV_rad,
	double* armgoal_anglesV_rad,
	int numofDOFs,
	double*** plan,
	int* planlength)
{
	const int NUM_SAMPLES = 10000;
	const int K_NEAREST = 500;

	std::vector<Node*> roadmap;
	auto isValid = [&](const std::vector<double>& config) {
		return IsValidArmConfiguration(const_cast<double*>(config.data()), numofDOFs, map, x_size, y_size);
	};

	hnswlib::L2Space space(numofDOFs);
	hnswlib::HierarchicalNSW<float> hnsw(&space, NUM_SAMPLES + 2);
	std::vector<std::vector<float>> hnsw_data_storage;
	std::unordered_map<int, Node*> label_to_node;

	int idx = 0;
	for (int i = 0; i < NUM_SAMPLES; i++) {
		std::vector<double> sample(numofDOFs);
		for (int j = 0; j < numofDOFs; j++) {
			sample[j] = ((double)rand() / RAND_MAX) * 2 * PI;
		}
		if (isValid(sample)) {
			Node* node = new Node{sample};
			roadmap.push_back(node);
			hnsw_data_storage.emplace_back(sample.begin(), sample.end());
			hnsw.addPoint(hnsw_data_storage.back().data(), idx);
			label_to_node[idx] = node;
			idx++;
		}
	}

	Node* start_node = new Node{{armstart_anglesV_rad, armstart_anglesV_rad + numofDOFs}};
	Node* goal_node = new Node{{armgoal_anglesV_rad, armgoal_anglesV_rad + numofDOFs}};

	if (!isValid(start_node->angles) || !isValid(goal_node->angles)) {
		printf("Start or Goal is in collision!\n");
		*planlength = 0;
		*plan = nullptr;
		return;
	}

	hnsw_data_storage.emplace_back(start_node->angles.begin(), start_node->angles.end());
	hnsw.addPoint(hnsw_data_storage.back().data(), idx);
	label_to_node[idx] = start_node;
	roadmap.push_back(start_node);
	int start_idx = idx++;

	hnsw_data_storage.emplace_back(goal_node->angles.begin(), goal_node->angles.end());
	hnsw.addPoint(hnsw_data_storage.back().data(), idx);
	label_to_node[idx] = goal_node;
	roadmap.push_back(goal_node);
	int goal_idx = idx++;

	for (int i = 0; i < roadmap.size(); ++i) {
		Node* node = roadmap[i];
		std::priority_queue<std::pair<float, hnswlib::labeltype>> result = hnsw.searchKnn(node->angles.data(), K_NEAREST + 1);
		while (!result.empty()) {
			auto [dist, label] = result.top(); result.pop();
			if (label_to_node.find(label) == label_to_node.end()) continue;
			if (label_to_node[label] == node) continue;
			Node* neighbor = label_to_node[label];
			if (isValid(neighbor->angles)) {
				node->neighbors.push_back(neighbor);
			}
		}
	}

	// debugSummary(roadmap, start_node, goal_node);

	auto compare = [](const std::pair<double, Node*>& a, const std::pair<double, Node*>& b) {
		return a.first > b.first;
	};

	std::priority_queue<std::pair<double, Node*>, std::vector<std::pair<double, Node*>>, decltype(compare)> pq(compare);
	std::unordered_set<Node*> visited;
	start_node->cost = 0.0;
	pq.emplace(0.0, start_node);

	// printf("[DEBUG] Starting Dijkstra from start_node...\n");

	bool found = false;
	while (!pq.empty()) {
		Node* current = pq.top().second;
		pq.pop();
		if (visited.find(current) != visited.end()) continue;
		visited.insert(current);

		// printf("  Visiting node with cost %.3f, neighbors: %lu\n", current->cost, current->neighbors.size());
		double dist_to_goal = computeDistance(current->angles, goal_node->angles);
		// printf("    --> Distance to goal: %.6f\n", dist_to_goal);

		if (current == goal_node) {
			found = true;
			break;
		}
		for (Node* neighbor : current->neighbors) {
			double new_cost = current->cost + computeDistance(current->angles, neighbor->angles);
			if (visited.find(neighbor) == visited.end() && (neighbor->parent == nullptr || new_cost < neighbor->cost)) {
				neighbor->cost = new_cost;
				neighbor->parent = current;
				pq.emplace(new_cost, neighbor);
			}
		}
	}

	if (goal_node->parent == nullptr) {
		printf("[DEBUG] Goal node was never reached in search.\n");
	}

	if (!found) {
		printf("No valid path found!\n");
		*planlength = 0;
		*plan = nullptr;
		return;
	}

	std::vector<std::vector<double>> path;
	int MAX_ITER = 10000, iteration = 0;
	for (Node* node = goal_node; node != nullptr; node = node->parent) {
		if (iteration++ > MAX_ITER) {
			printf("ERROR DETECTED: Infinite loop while backtracking!\n");
			return;
		}
		path.push_back(node->angles);
	}
	reverse(path.begin(), path.end());
	*planlength = path.size();
	*plan = (double**) malloc(*planlength * sizeof(double*));
	for (int i = 0; i < *planlength; i++) {
		(*plan)[i] = (double*) malloc(numofDOFs * sizeof(double));
		memcpy((*plan)[i], path[i].data(), numofDOFs * sizeof(double));
	}

	printf("PRM-HNSW found a path with %d waypoints.\n", *planlength);
	for (Node* node : roadmap) delete node;
}

// Standard RRT Planner (Unidirectional-Extending)
static void RRTPlanner(
    double* map,
    int x_size,
    int y_size,
    double* armstart_anglesV_rad,
    double* armgoal_anglesV_rad,
    int numofDOFs,
    double*** plan,
    int* planlength)
{
    const int MAX_ITER = 7000;       // maximum iterations
    const double INITIAL_STEP_SIZE = 0.05; // dynamic step size
    const double MIN_STEP_SIZE = 0.005; // smaller step size near obstacles
    const double GOAL_BIAS_INITIAL = 0.05;   // probability of sampling goal directly (5%)
    const double GOAL_THRESHOLD = 0.02; // threshold connection to goal

    struct Node {
        vector<double> angles;
        Node* parent;
        Node(vector<double> ang, Node* p = nullptr) : angles(ang), parent(p) {}
    };

    vector<Node*> tree;
    double goal_bias = GOAL_BIAS_INITIAL;

    auto isValid = [&](const vector<double>& config) {
        return IsValidArmConfiguration(const_cast<double*>(config.data()), numofDOFs, map, x_size, y_size);
    };

    // find the nearest neighbor node in the tree to a given target configuration
    auto nearestNeighbor = [](const vector<Node*>& tree, const vector<double>& target) {
        Node* best = nullptr;
        double min_dist = DBL_MAX;
        for (Node* node : tree) {
            double dist = computeDistance(node->angles, target);
            if (dist < min_dist) {
                min_dist = dist;
                best = node;
            }
        }
        return best;
    };

    auto extend = [&](Node* nearest, const std::vector<double>& target, double step_size) {
        std::vector<double> new_ang = nearest->angles;
        double dist = computeDistance(nearest->angles, target);
        double alpha = std::min(step_size / dist, 1.0);
        for (int i = 0; i < numofDOFs; i++) {
            new_ang[i] = nearest->angles[i] + alpha * (target[i] - nearest->angles[i]);
        }
        if (isValid(new_ang)) {
            return new Node(new_ang, nearest);
        }
        return static_cast<Node*>(nullptr);
    };

    tree.push_back(new Node({armstart_anglesV_rad, armstart_anglesV_rad + numofDOFs}));
    vector<double> goal_config(armgoal_anglesV_rad, armgoal_anglesV_rad + numofDOFs);

    if (!isValid(tree.back()->angles) || !isValid(goal_config)) {
        printf("Start or Goal is in collision!\n");
        return;
    }

    // main RRT extending
    for (int iter = 0; iter < MAX_ITER; iter++) {
        // goal bias
        vector<double> random_config(numofDOFs);
        if ((double)rand() / RAND_MAX < goal_bias) {
            random_config = goal_config;
        } else {
            for (int j = 0; j < numofDOFs; j++) {
                random_config[j] = ((double)rand() / RAND_MAX) * 2 * M_PI;
            }
        }

        Node* nearest = nearestNeighbor(tree, random_config);
        double step_size = std::max(MIN_STEP_SIZE, INITIAL_STEP_SIZE * (1.0 - (double)iter / MAX_ITER));
        Node* new_node = extend(nearest, random_config, step_size);
        if (!new_node) continue; // trapped
        tree.push_back(new_node);

        // check the distance to the goal
        double dist_to_goal = computeDistance(new_node->angles, goal_config);
        if (dist_to_goal < GOAL_THRESHOLD) {
            Node* final_node = extend(new_node, goal_config, MIN_STEP_SIZE);
            if (final_node) {
                tree.push_back(final_node);
                break;
            }
        }

        goal_bias = std::min(0.3, goal_bias + 0.0005); // dynamically increase goal bias over iterations
    }

    // backtrack the path
    vector<vector<double>> path;
    Node* node = tree.back();
    while (node) {
        path.push_back(node->angles);
        node = node->parent;
    }
    reverse(path.begin(), path.end());

    // ensure the path starts at start position and ends at goal position
    memcpy(path.front().data(), armstart_anglesV_rad, numofDOFs * sizeof(double));
    memcpy(path.back().data(), armgoal_anglesV_rad, numofDOFs * sizeof(double));

    *planlength = path.size();
    *plan = (double**) malloc(*planlength * sizeof(double*));
    for (int i = 0; i < *planlength; i++) {
        (*plan)[i] = (double*) malloc(numofDOFs * sizeof(double));
        memcpy((*plan)[i], path[i].data(), numofDOFs * sizeof(double));
    }

    for (Node* node : tree) delete node;

    printf("RRT successfully found a path with %d waypoints.\n", *planlength);
}


// RRT-HNSW Planner
static void RRTHNSWPlanner(
    double* map,
    int x_size,
    int y_size,
    double* armstart_anglesV_rad,
    double* armgoal_anglesV_rad,
    int numofDOFs,
    double*** plan,
    int* planlength)
{
    const int MAX_ITER = 20000;
    const double INITIAL_STEP_SIZE = 0.05;
    const double MIN_STEP_SIZE = 0.005;
    const double GOAL_BIAS_INITIAL = 0.1;
    const double GOAL_THRESHOLD = 0.02;

    struct Node {
        std::vector<double> angles;
        Node* parent;
        Node(std::vector<double> ang, Node* p = nullptr) : angles(ang), parent(p) {}
    };

    std::vector<Node*> tree;
    hnswlib::L2Space space(numofDOFs);
    hnswlib::HierarchicalNSW<float> hnsw(&space, MAX_ITER + 2);
    std::vector<std::vector<float>> hnsw_data;
    std::unordered_map<int, Node*> label_to_node;

    double goal_bias = GOAL_BIAS_INITIAL;
    auto isValid = [&](const std::vector<double>& config) {
        return IsValidArmConfiguration(const_cast<double*>(config.data()), numofDOFs, map, x_size, y_size);
    };

    auto extend = [&](Node* nearest, const std::vector<double>& target, double step_size) {
        std::vector<double> new_ang = nearest->angles;
        double dist = computeDistance(nearest->angles, target);
        double alpha = std::min(step_size / dist, 1.0);
        for (int i = 0; i < numofDOFs; i++) {
            new_ang[i] = nearest->angles[i] + alpha * (target[i] - nearest->angles[i]);
        }
        if (isValid(new_ang)) return new Node(new_ang, nearest);
        return static_cast<Node*>(nullptr);
    };

    Node* start_node = new Node({armstart_anglesV_rad, armstart_anglesV_rad + numofDOFs});
    std::vector<double> goal_config(armgoal_anglesV_rad, armgoal_anglesV_rad + numofDOFs);
    if (!isValid(start_node->angles) || !isValid(goal_config)) {
        printf("Start or Goal is in collision!\n");
        return;
    }

    tree.push_back(start_node);
    hnsw_data.emplace_back(start_node->angles.begin(), start_node->angles.end());
    hnsw.addPoint(hnsw_data.back().data(), 0);
    label_to_node[0] = start_node;
    int hnsw_id = 1;

    Node* final_node = nullptr;
    for (int iter = 0; iter < MAX_ITER; iter++) {
        std::vector<double> random_config(numofDOFs);
        if ((double)rand() / RAND_MAX < goal_bias) {
            random_config = goal_config;
        } else {
            for (int j = 0; j < numofDOFs; j++) {
                random_config[j] = ((double)rand() / RAND_MAX) * 2 * M_PI;
            }
        }

        auto result = hnsw.searchKnn(random_config.data(), 1);
        Node* nearest = label_to_node[result.top().second];

        double step_size = std::max(MIN_STEP_SIZE, INITIAL_STEP_SIZE * (1.0 - (double)iter / MAX_ITER));
        Node* new_node = extend(nearest, random_config, step_size);
        if (!new_node) continue;

        tree.push_back(new_node);
        hnsw_data.emplace_back(new_node->angles.begin(), new_node->angles.end());
        hnsw.addPoint(hnsw_data.back().data(), hnsw_id);
        label_to_node[hnsw_id] = new_node;
        hnsw_id++;

        double dist_to_goal = computeDistance(new_node->angles, goal_config);
        if (dist_to_goal < GOAL_THRESHOLD) {
            Node* connect_goal = extend(new_node, goal_config, MIN_STEP_SIZE);
            if (connect_goal) {
                final_node = connect_goal;
                final_node->parent = new_node;
                tree.push_back(final_node);
                break;
            }
        }

        goal_bias = std::min(0.3, goal_bias + 0.0005);
    }

    if (!final_node) {
        printf("RRT-HNSW failed to find a path.\n");
        *planlength = 0;
        *plan = nullptr;
        return;
    }

    std::vector<std::vector<double>> path;
    for (Node* node = final_node; node != nullptr; node = node->parent) {
        path.push_back(node->angles);
    }
    std::reverse(path.begin(), path.end());

    memcpy(path.front().data(), armstart_anglesV_rad, numofDOFs * sizeof(double));
    memcpy(path.back().data(), armgoal_anglesV_rad, numofDOFs * sizeof(double));

    *planlength = path.size();
    *plan = (double**) malloc(*planlength * sizeof(double*));
    for (int i = 0; i < *planlength; i++) {
        (*plan)[i] = (double*) malloc(numofDOFs * sizeof(double));
        memcpy((*plan)[i], path[i].data(), numofDOFs * sizeof(double));
    }

    for (Node* node : tree) delete node;
    printf("RRT-HNSW successfully found a path with %d waypoints.\n", *planlength);
}

struct NodeCloud {
    vector<Node*> nodes;

    // Required for nanoflann
    inline size_t kdtree_get_point_count() const { return nodes.size(); }
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return nodes[idx]->angles[dim];
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<double, NodeCloud>,
    NodeCloud,
    -1, // Dynamic dimensions
    size_t
> KDTree;

// PRM-KDTree Planner
static void PRMKDTreePlanner(
    double* map, int x_size, int y_size,
    double* armstart_anglesV_rad, double* armgoal_anglesV_rad, int numofDOFs,
    double*** plan, int* planlength)
{
    const int NUM_SAMPLES = 10000;
    const int K_NEAREST = 13;

    vector<Node*> roadmap;
    NodeCloud cloud;

    auto isValid = [&](const vector<double>& config) {
        return IsValidArmConfiguration(const_cast<double*>(config.data()), numofDOFs, map, x_size, y_size);
    };

    // Generate valid samples
    for (int i = 0; i < NUM_SAMPLES; i++) {
        vector<double> sample(numofDOFs);
        for (int j = 0; j < numofDOFs; j++)
            sample[j] = ((double)rand() / RAND_MAX) * 2 * PI;

        if (isValid(sample)) {
            Node* new_node = new Node{sample};
            roadmap.push_back(new_node);
            cloud.nodes.push_back(new_node);
        }
    }

    Node* start_node = new Node{{armstart_anglesV_rad, armstart_anglesV_rad + numofDOFs}};
    Node* goal_node = new Node{{armgoal_anglesV_rad, armgoal_anglesV_rad + numofDOFs}};
    roadmap.push_back(start_node);
    roadmap.push_back(goal_node);
    cloud.nodes.push_back(start_node);
    cloud.nodes.push_back(goal_node);

    if (!isValid(start_node->angles) || !isValid(goal_node->angles)) {
        printf("Start or Goal is in collision!\n");
        return;
    }

    // Build KDTree
    KDTree index(numofDOFs, cloud, KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    // Connect neighbors using KDTree
    for (Node* node : roadmap) {
        vector<size_t> ret_indices(K_NEAREST + 1);
        vector<double> out_distances_sqr(K_NEAREST + 1);

        KNNResultSet<double> resultSet(K_NEAREST + 1);
        resultSet.init(&ret_indices[0], &out_distances_sqr[0]);
        index.findNeighbors(resultSet, node->angles.data(), nanoflann::SearchParameters(10));

        for (size_t i = 1; i < resultSet.size(); i++) { // Skip first as it's self
            Node* neighbor = cloud.nodes[ret_indices[i]];
            if (isValid(neighbor->angles))
                node->neighbors.push_back(neighbor);
        }
    }

    // Use existing Dijkstra algorithm for searching
    auto compare = [](const pair<double, Node*>& a, const pair<double, Node*>& b) {
        return a.first > b.first;
    };

    priority_queue<pair<double, Node*>, vector<pair<double, Node*>>, decltype(compare)> pq(compare);
    unordered_set<Node*> visited;

    start_node->cost = 0.0;
    pq.emplace(0.0, start_node);

    bool found = false;
    while (!pq.empty()) {
        Node* current = pq.top().second;
        pq.pop();

        if (visited.count(current)) continue;
        visited.insert(current);

        if (current == goal_node) {
            found = true;
            break;
        }

        for (Node* neighbor : current->neighbors) {
            double new_cost = current->cost + computeDistance(current->angles, neighbor->angles);
            if (!visited.count(neighbor) && (neighbor->parent == nullptr || new_cost < neighbor->cost)) {
                neighbor->cost = new_cost;
                neighbor->parent = current;
                pq.emplace(new_cost, neighbor);
            }
        }
    }

    if (!found) {
        printf("No valid path found!\n");
        return;
    }

    vector<vector<double>> path;
    for (Node* node = goal_node; node != nullptr; node = node->parent)
        path.push_back(node->angles);

    reverse(path.begin(), path.end());
    *planlength = path.size();
    *plan = (double**)malloc(*planlength * sizeof(double*));

    for (int i = 0; i < *planlength; i++) {
        (*plan)[i] = (double*)malloc(numofDOFs * sizeof(double));
        memcpy((*plan)[i], path[i].data(), numofDOFs * sizeof(double));
    }

    printf("PRM-KDTree successfully found a path with %d waypoints.\n", *planlength);
    for (Node* node : roadmap) delete node;
}

// RRT-KDTree Planner
static void RRTKDTreePlanner(
    double* map, int x_size, int y_size,
    double* armstart_anglesV_rad, double* armgoal_anglesV_rad, int numofDOFs,
    double*** plan, int* planlength)
{
    const int MAX_ITER = 10000;
    const double STEP_SIZE = 0.05;
    const double GOAL_BIAS = 0.1;
    const double GOAL_THRESHOLD = 0.02;

    vector<Node*> tree;
    NodeCloud cloud;

    auto isValid = [&](const vector<double>& config) {
        return IsValidArmConfiguration(const_cast<double*>(config.data()), numofDOFs, map, x_size, y_size);
    };

    Node* start_node = new Node{{armstart_anglesV_rad, armstart_anglesV_rad + numofDOFs}};
    vector<double> goal_config(armgoal_anglesV_rad, armgoal_anglesV_rad + numofDOFs);

    if (!isValid(start_node->angles) || !isValid(goal_config)) {
        printf("Start or Goal is in collision!\n");
        return;
    }

    tree.push_back(start_node);
    cloud.nodes.push_back(start_node);

    KDTree index(numofDOFs, cloud, KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    Node* final_node = nullptr;

    for (int iter = 0; iter < MAX_ITER; iter++) {
        vector<double> sample(numofDOFs);
        if ((double)rand() / RAND_MAX < GOAL_BIAS) {
            sample = goal_config;
        } else {
            for (int j = 0; j < numofDOFs; j++)
                sample[j] = ((double)rand() / RAND_MAX) * 2 * PI;
        }

        size_t nearest_idx;
        double dist_sq;
        KNNResultSet<double> resultSet(1);
        resultSet.init(&nearest_idx, &dist_sq);
        index.findNeighbors(resultSet, sample.data(), nanoflann::SearchParameters(10));
        Node* nearest = cloud.nodes[nearest_idx];

        vector<double> direction(numofDOFs);
        double dist = sqrt(dist_sq);
        for (int j = 0; j < numofDOFs; j++)
            direction[j] = nearest->angles[j] + STEP_SIZE * (sample[j] - nearest->angles[j]) / dist;

        if (!isValid(direction)) continue;

		Node* new_node = new Node{direction};
		new_node->parent = nearest;

        tree.push_back(new_node);
        cloud.nodes.push_back(new_node);
		
		if (cloud.nodes.size() % 500 == 0) {
			index.buildIndex();
		}
        if (computeDistance(new_node->angles, goal_config) < GOAL_THRESHOLD) {
			final_node = new Node{goal_config};
			final_node->parent = new_node;
            tree.push_back(final_node);
            break;
        }
    }

    if (!final_node) {
        printf("RRT-KDTree failed to find a path.\n");
        return;
    }

    vector<vector<double>> path;
    for (Node* node = final_node; node != nullptr; node = node->parent)
        path.push_back(node->angles);

    reverse(path.begin(), path.end());
    *planlength = path.size();
    *plan = (double**)malloc(*planlength * sizeof(double*));

    for (int i = 0; i < *planlength; i++) {
        (*plan)[i] = (double*)malloc(numofDOFs * sizeof(double));
        memcpy((*plan)[i], path[i].data(), numofDOFs * sizeof(double));
    }

    printf("RRT-KDTree successfully found a path with %d waypoints.\n", *planlength);
    for (Node* node : tree) delete node;
}

/** Your final solution will be graded by an grading script which will
 * send the default 6 arguments:
 *    map, numOfDOFs, commaSeparatedStartPos, commaSeparatedGoalPos, 
 *    whichPlanner, outputFilePath
 * An example run after compiling and getting the planner.out executable
 * >> ./planner.out map1.txt 5 1.57,0.78,1.57,0.78,1.57 0.392,2.35,3.14,2.82,4.71 0 output.txt
 * See the hw handout for full information.
 * If you modify this for testing (e.g. to try out different hyper-parameters),
 * make sure it can run with the original 6 commands.
 * Programs that do not will automatically get a 0.
 * */
int main(int argc, char** argv) {
	double* map;
	int x_size, y_size;

	tie(map, x_size, y_size) = loadMap(argv[1]);
	const int numOfDOFs = std::stoi(argv[2]);
	double* startPos = doubleArrayFromString(argv[3]);
	double* goalPos = doubleArrayFromString(argv[4]);
	int whichPlanner = std::stoi(argv[5]);
	string outputFile = argv[6];

	if(!IsValidArmConfiguration(startPos, numOfDOFs, map, x_size, y_size)||
			!IsValidArmConfiguration(goalPos, numOfDOFs, map, x_size, y_size)) {
		throw runtime_error("Invalid start or goal configuration!\n");
	}

	///////////////////////////////////////
	//// Feel free to modify anything below. Be careful modifying anything above.

	double** plan = NULL;
	int planlength = 0;
    
	if (whichPlanner == PRM) {
		PRMPlanner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	} else if (whichPlanner == RRT) {
		RRTPlanner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	} else if (whichPlanner == PRM_HNSW) {
		PRMHNSWPlanner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	} else if (whichPlanner == RRT_HNSW) {
		RRTHNSWPlanner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	} else if (whichPlanner == PRM_KDTree) {
		PRMKDTreePlanner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	} else if (whichPlanner == RRT_KDTree) {
		RRTKDTreePlanner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	}

	//// Feel free to modify anything above.
	//// If you modify something below, please change it back afterwards as the 
	//// grading script will not work.
	///////////////////////////////////////

    // Your solution's path should start with startPos and end with goalPos
    if (!equalDoubleArrays(plan[0], startPos, numOfDOFs) || 
    	!equalDoubleArrays(plan[planlength-1], goalPos, numOfDOFs)) {
		throw std::runtime_error("Start or goal position not matching");
	}

	/** Saves the solution to output file
	 * Do not modify the output log file output format as it is required for visualization
	 * and for grading.
	 */
	std::ofstream m_log_fstream;
	m_log_fstream.open(outputFile, std::ios::trunc); // Creates new or replaces existing file
	if (!m_log_fstream.is_open()) {
		throw std::runtime_error("Cannot open file");
	}
	m_log_fstream << argv[1] << endl; // Write out map name first
	/// Then write out all the joint angles in the plan sequentially
	for (int i = 0; i < planlength; ++i) {
		for (int k = 0; k < numOfDOFs; ++k) {
			m_log_fstream << plan[i][k] << ",";
		}
		m_log_fstream << endl;
	}
}
