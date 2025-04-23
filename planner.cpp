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
#include <queue>

#include <tuple>
#include <string>
#include <stdexcept>
#include <regex> // For regex and split logic
#include <iostream> // cout, endl
#include <fstream> // For reading/writing files
#include <assert.h> 

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define PRM         0
#define PRM_HNSW    1
#define RRT         2


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

#define PI 3.141592654

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

#ifdef DEBUG
  #define DEBUG_PRINT(x) std::cout << x
#else
  #define DEBUG_PRINT(x)
#endif

/** @brief 
 * @param filepath 
 * @return map, x_size, y_size
 */
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
std::vector<string> split(const string& str, const string& delim) {   
		// https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c/64886763#64886763
		const std::regex ws_re(delim);
		return { std::sregex_token_iterator(str.begin(), str.end(), ws_re, -1), std::sregex_token_iterator() };
}

double* doubleArrayFromString(string str) {
	std::vector<string> vals = split(str, ",");
	double* ans = new double[vals.size()];
	for (int i = 0; i < vals.size(); ++i) {
		ans[i] = std::stod(vals[i]);
	}
	return ans;
}

bool equalDoubleArrays(double* v1, double *v2, int size) {
    for (int i = 0; i < size; ++i) {
        if (abs(v1[i]-v2[i]) > 1e-3) {
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

static void planner(
			double* map,
			int x_size,
			int y_size,
			double* armstart_anglesV_rad,
			double* armgoal_anglesV_rad,
            int numofDOFs,
            double*** plan,
            int* planlength) 
{
	//no plan by default
	*plan = NULL;
	*planlength = 0;
		
    //for now just do straight interpolation between start and goal checking for the validity of samples

    double distance = 0;
    int i,j;
    for (j = 0; j < numofDOFs; j++){
        if(distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
            distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
    }
    int numofsamples = (int)(distance/(PI/20));
    if(numofsamples < 2){
        printf("The arm is already at the goal\n");
        return;
    }
	int countNumInvalid = 0;
    *plan = (double**) malloc(numofsamples*sizeof(double*));
    for (i = 0; i < numofsamples; i++){
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = armstart_anglesV_rad[j] + ((double)(i)/(numofsamples-1))*(armgoal_anglesV_rad[j] - armstart_anglesV_rad[j]);
        }
        if(!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size)) {
			++countNumInvalid;
        }
    }
	printf("Linear interpolation collided at %d instances across the path\n", countNumInvalid);
    *planlength = numofsamples;
    
    return;
}

double distance(std::vector<double>& q1, std::vector<double>& q2) {
    double dist = 0.0;
    for (size_t i = 0; i < q1.size(); i++) {
        double diff = q1[i] - q2[i];
        dist += diff * diff;
    }
    return sqrt(dist);
}

struct TreeNode {
    std::vector<double> angles;
    TreeNode* parent;
    double cost;

    TreeNode(const std::vector<double>& a, TreeNode* p = nullptr, double c = 0.0) 
        : angles(a), parent(p), cost(c) {}
};

std::vector<double> randomConfig(int numofDOFs, double* map, int x_size, int y_size) {
    std::vector<double> q(numofDOFs);
    // Continue sampling until a valid configuration is found.
    while (!IsValidArmConfiguration(q.data(), numofDOFs, map, x_size, y_size)) {
        for (int i = 0; i < numofDOFs; i++) {
            q[i] = ((double)rand() / RAND_MAX) * 2 * PI;
        }
    }
    return q;
}

TreeNode* nearestNeighbor(std::vector<TreeNode*>& tree, std::vector<double>& q_rand) {
    TreeNode* nearest = nullptr;
    double min_dist = INFINITY;
    
    for (TreeNode* node : tree) {
        double dist = 0;
        for (int i = 0; i < q_rand.size(); i++) {
			double diff = node->angles[i] - q_rand[i];
            dist += diff * diff;
        }
        if (dist < min_dist) {
            min_dist = dist;
            nearest = node;
        }
    }
    return nearest;
}

std::vector<double> steer(std::vector<double>& q_from, std::vector<double>& q_to, double epsilon, int numofDOFs)
{
    std::vector<double> q_new(numofDOFs);
    double dist = distance(q_from, q_to);

    // If the points are effectively the same, return q_from
    if (dist < 1e-3) {
        return q_from;
    }

    double step_size = std::min(epsilon, dist);
    double factor = step_size / dist;

    for (int i = 0; i < numofDOFs; i++) {
        q_new[i] = q_from[i] + factor * (q_to[i] - q_from[i]);
    }
    return q_new;
}

std::vector<double> extend(std::vector<TreeNode*>& T,
                      std::vector<double>& q_target,
                      double* map,
                      int x_size,
                      int y_size,
                      double epsilon,
                      int numofDOFs)
{
    // 1) Find the nearest node in T to q_target
    TreeNode* q_near = nearestNeighbor(T, q_target);
    if (!q_near) {
        DEBUG_PRINT("[extend] ERROR: No nearest node found! Tree size: " << T.size() << std::endl);
        return {}; // Tree is empty or something went wrong
    }

    DEBUG_PRINT("[extend] q_near: ");
    for (double v : q_near->angles) {
        DEBUG_PRINT(v << " ");
    }
    DEBUG_PRINT(std::endl);

    // 2) Compute q_new by stepping from q_near->angles toward q_target
    std::vector<double> q_new = steer(q_near->angles, q_target, epsilon, numofDOFs);

    DEBUG_PRINT("[extend] q_target: ");
    for (double v : q_target) {
        DEBUG_PRINT(v << " ");
    }
    DEBUG_PRINT(std::endl);

    DEBUG_PRINT("[extend] q_new: ");
    for (double v : q_new) {
        DEBUG_PRINT(v << " ");
    }
    DEBUG_PRINT(std::endl);

    // If steer returned the same as q_near->angles => no progress
    if (equalDoubleArrays(q_new.data(), q_near->angles.data(), numofDOFs)) {
        DEBUG_PRINT("[extend] WARNING: q_new equals q_near (no progress made)." << std::endl);
        return {};  // "Trapped" => empty
    }

    // 3) Check if q_new is valid (collision-free)
    if (!IsValidArmConfiguration(q_new.data(), numofDOFs, map, x_size, y_size)) {
        DEBUG_PRINT("[extend] WARNING: q_new is invalid (collision detected)." << std::endl);
        return {};  // "Trapped" => empty
    }

    // 4) Add q_new as a new TreeNode in T
    TreeNode* new_node = new TreeNode(q_new, q_near);
    T.push_back(new_node);

    // 5) Return the newly created configuration
    return q_new;
}

bool connect(std::vector<TreeNode*>& T,
             std::vector<double>& q_target,
             double* map,
             int x_size,
             int y_size,
             double epsilon,
             int numofDOFs)
{
    while (true) {
        // Attempt to extend the tree toward q_target
        std::vector<double> q_new = extend(T, q_target, map, x_size, y_size, epsilon, numofDOFs);

        // If extend(...) returned an empty std::vector, we're trapped (no progress).
        if (q_new.empty()) {
            return false; // TRAPPED
        }

        // If q_new == q_target, we've reached the target exactly.
        if (equalDoubleArrays(q_new.data(), q_target.data(), numofDOFs)) {
            return true; // REACHED
        }

        // Otherwise, we ADVANCED but haven't reached q_target yet.
        // We keep looping to extend further from the newly added node.
    }
}

void getRRTConnectPath(std::vector<TreeNode*>& T_a, std::vector<TreeNode*>& T_b, double*** plan, int* planlength) {
    // Extract path from T_a (from start to connection)
    std::vector<std::vector<double>> path_a;
    TreeNode* node = T_a.back();
    DEBUG_PRINT("[getRRTConnectPath] T_a connection node: ");
    for (double v : node->angles) {
        DEBUG_PRINT(v << " ");
    }
    DEBUG_PRINT(std::endl);
    while (node) {
        path_a.push_back(node->angles);
        node = node->parent;
    }
    reverse(path_a.begin(), path_a.end()); // From start -> connection

    // Extract path from T_b (from connection to goal)
    std::vector<std::vector<double>> path_b;
    node = T_b.back();
    DEBUG_PRINT("[getRRTConnectPath] T_b connection node: ");
    for (double v : node->angles) {
        DEBUG_PRINT(v << " ");
    }
    DEBUG_PRINT(std::endl);
    while (node) {
        path_b.push_back(node->angles);
        node = node->parent;
    }

    // Check if the connection nodes match.
    // If they do, the first element of path_b should equal the last element of path_a.
    if (!path_a.empty() && !path_b.empty() &&
        equalDoubleArrays(path_a.back().data(), path_b.front().data(), path_a.back().size())) {
        DEBUG_PRINT("[getRRTConnectPath] Duplicate connection node detected, removing duplicate." << std::endl);
        path_b.erase(path_b.begin());
    } else {
        DEBUG_PRINT("[getRRTConnectPath] Warning: Connection nodes do not match!" << std::endl);
    }

    // Concatenate the two paths.
    std::vector<std::vector<double>> full_path = path_a;
    full_path.insert(full_path.end(), path_b.begin(), path_b.end());

    // Print the full path for debugging.
    DEBUG_PRINT("[getRRTConnectPath] Full path:" << std::endl);
    for (size_t i = 0; i < full_path.size(); i++) {
        DEBUG_PRINT("  Node " << i << ": ");
        for (double val : full_path[i]) {
            DEBUG_PRINT(val << " ");
        }
        DEBUG_PRINT(std::endl);
    }

    // Allocate memory for the plan output.
    *planlength = full_path.size();
    *plan = new double*[*planlength];
    for (size_t i = 0; i < full_path.size(); i++) {
        (*plan)[i] = new double[full_path[i].size()];
        copy(full_path[i].begin(), full_path[i].end(), (*plan)[i]);
    }
}


struct PRMNode {
    std::vector<double> angles;
    std::vector<int> neighbors;
};

static void PRM_planner(
    double* map,
    int x_size,
    int y_size,
    double* armstart_anglesV_rad,
    double* armgoal_anglesV_rad,
    int numofDOFs,
    double*** plan,
    int* planlength)
{
    *plan = NULL;
    *planlength = 0;
    
    int numSamples = 20000;         // Number of random samples to add.
    int k = 500;                   // Number of nearest neighbors to connect.
    
    // Build the roadmap.
    // First, add the start and goal nodes.
    vector<PRMNode> roadmap;
    {
        PRMNode startNode;
        startNode.angles = vector<double>(armstart_anglesV_rad, armstart_anglesV_rad + numofDOFs);
        roadmap.push_back(startNode);
        DEBUG_PRINT("[PRM] Start node added." << endl);

        PRMNode goalNode;
        goalNode.angles = vector<double>(armgoal_anglesV_rad, armgoal_anglesV_rad + numofDOFs);
        roadmap.push_back(goalNode);
        DEBUG_PRINT("[PRM] Goal node added." << endl);
    }
    
    // Sample random, collision-free configurations.
    for (int i = 0; i < numSamples; i++) {
        vector<double> q = randomConfig(numofDOFs, map, x_size, y_size);
        PRMNode node;
        node.angles = q;
        roadmap.push_back(node);
    }
    DEBUG_PRINT("[PRM] Roadmap sampled: " << roadmap.size() << " nodes." << endl);

    // Build roadmap edges using k nearest neighbors.
    int N = roadmap.size();
    for (int i = 0; i < N; i++) {
        // Build a list of (distance, index) pairs for node i.
        vector<std::pair<double, int>> dists;
        for (int j = 0; j < N; j++) {
            if (j == i) continue;
            double d = distance(roadmap[i].angles, roadmap[j].angles);
            if (d > 1e-3) {
                dists.push_back(std::make_pair(d, j));
            }
        }
        // Sort the pairs by distance.
        std::sort(dists.begin(), dists.end(), [](auto &a, auto &b) { return a.first < b.first; });
        int num_neighbors = std::min(k, (int)dists.size());
        for (int n = 0; n < num_neighbors; n++) {
            int neighbor_index = dists[n].second;
            roadmap[i].neighbors.push_back(neighbor_index);
            // Add reverse edge if not already present.
            roadmap[neighbor_index].neighbors.push_back(i);
            DEBUG_PRINT("[PRM] Edge added between nodes " << i << " and " << neighbor_index 
                        << " (d = " << dists[n].first << ")." << endl);
        }
    }
    DEBUG_PRINT("[PRM] Roadmap edges built." << endl);

    // Compute the shortest path from start (index 0) to goal (index 1) using Dijkstra
    std::vector<double> cost(N, std::numeric_limits<double>::infinity());
    std::vector<int> prev(N, -1);
    typedef std::pair<double, int> P;
    std::priority_queue<P, vector<P>, std::greater<P>> pq;
    cost[0] = 0.0;
    pq.push(P(0.0, 0));
    DEBUG_PRINT("[PRM] Dijkstra: Starting from node 0." << endl);

    while (!pq.empty()) {
        auto cur_entry = pq.top();
        pq.pop();
        double cur_cost = cur_entry.first;
        int cur = cur_entry.second;
        DEBUG_PRINT("[PRM] Dijkstra: Processing node " << cur << " with cost " << cur_cost << "." << endl);
        if (cur_cost > cost[cur]) continue; // Outdated entry.
        if (cur == 1) {
            DEBUG_PRINT("[PRM] Dijkstra: Goal node reached with cost " << cur_cost << "." << endl);
            break; // Goal reached.
        }
        for (int neighbor : roadmap[cur].neighbors) {
            double weight = distance(roadmap[cur].angles, roadmap[neighbor].angles);
            if (cost[cur] + weight < cost[neighbor]) {
                cost[neighbor] = cost[cur] + weight;
                prev[neighbor] = cur;
                pq.push(P(cost[neighbor], neighbor));
                DEBUG_PRINT("[PRM] Dijkstra: Updating node " << neighbor 
                            << " with new cost " << cost[neighbor] 
                            << " via node " << cur << "." << endl);
            }
        }
    }

    if (cost[1] == std::numeric_limits<double>::infinity()) {
        DEBUG_PRINT("[PRM] Failed to connect start and goal." << endl);
        return;
    }
    
    // Reconstruct the path from goal (index 1) back to start (index 0).
    std::vector<int> pathIndices;
    for (int cur = 1; cur != -1; cur = prev[cur]) {
        pathIndices.push_back(cur);
    }
    reverse(pathIndices.begin(), pathIndices.end());
    
    DEBUG_PRINT("[PRM] Path found with " << pathIndices.size() << " nodes." << endl);
    for (size_t i = 0; i < pathIndices.size(); i++) {
        DEBUG_PRINT("  Node " << i << ": ");
        for (double v : roadmap[pathIndices[i]].angles)
            DEBUG_PRINT(v << " ");
        DEBUG_PRINT(std::endl);
    }
    
    // Build the output plan from the roadmap.
    vector<vector<double>> full_path;
    for (int idx : pathIndices) {
        full_path.push_back(roadmap[idx].angles);
    }
    
    // Allocate the plan output.
    *planlength = full_path.size();
    *plan = new double*[*planlength];
    for (size_t i = 0; i < full_path.size(); i++) {
        (*plan)[i] = new double[full_path[i].size()];
        std::copy(full_path[i].begin(), full_path[i].end(), (*plan)[i]);
    }
}

static void RRT_planner(
    double* map,
    int x_size,
    int y_size,
    double* armstart_anglesV_rad,
    double* armgoal_anglesV_rad,
    int numofDOFs,
    double*** plan,
    int* planlength)
{
    *plan = NULL;
    *planlength = 0;

    int K = 5000;
    double epsilon = 1;

    if (numofDOFs > 5) {
        K = 50000;      // More iterations for a larger space.
        epsilon = 2;  
    }

    // Initialize the tree with the start configuration.
    std::vector<TreeNode*> T;
    T.push_back(new TreeNode(vector<double>(armstart_anglesV_rad, armstart_anglesV_rad + numofDOFs)));

    DEBUG_PRINT("RRT Planner Started" << std::endl);
    DEBUG_PRINT("Start Configuration: ");
    for (double angle : T.front()->angles)
        DEBUG_PRINT(angle << " ");
    DEBUG_PRINT(std::endl);

    std::vector<double> goalConfig(armgoal_anglesV_rad, armgoal_anglesV_rad + numofDOFs);
    DEBUG_PRINT("Goal Configuration: ");
    for (double angle : goalConfig)
        DEBUG_PRINT(angle << " ");
    DEBUG_PRINT(std::endl);

    TreeNode* goal_node = nullptr;
    for (int i = 0; i < K; i++) {
        DEBUG_PRINT("Iteration " << i + 1 << std::endl);
        
        vector<double> q_rand = randomConfig(numofDOFs, map, x_size, y_size);

        DEBUG_PRINT("  q_rand: ");
        for (double val : q_rand)
            DEBUG_PRINT(val << " ");
        DEBUG_PRINT(std::endl);

        vector<double> q_new = extend(T, q_rand, map, x_size, y_size, epsilon, numofDOFs);
        if (!q_new.empty()) {
            DEBUG_PRINT("  q_new added to tree: ");
            for (double val : q_new)
                DEBUG_PRINT(val << " ");
            DEBUG_PRINT(std::endl);

            // Check if the new configuration matches the goal (within tolerance).
            if (distance(q_new, goalConfig) <= epsilon) {
                T.push_back(new TreeNode(goalConfig, T.back()));
                DEBUG_PRINT("  Goal reached!" << std::endl);
                goal_node = T.back();
                break;
            }
        } else {
            DEBUG_PRINT("  q_new is empty (extension failed)" << std::endl);
        }
    }

    if (goal_node == nullptr) {
        DEBUG_PRINT("Failed to find a path within " << K << " iterations." << std::endl);
        return;
    }

    // Extract the path from the goal node back to the start.
    std::vector<std::vector<double>> path;
    TreeNode* node = goal_node;
    while (node) {
        path.push_back(node->angles);
        node = node->parent;
    }
    reverse(path.begin(), path.end());
    DEBUG_PRINT("Path found with length: " << path.size() << std::endl);
    for (size_t i = 0; i < path.size(); i++) {
        DEBUG_PRINT("  Node " << i << ": ");
        for (double val : path[i])
            DEBUG_PRINT(val << " ");
        DEBUG_PRINT(std::endl);
    }

    // Allocate the plan output.
    *planlength = path.size();
    *plan = new double*[*planlength];
    for (size_t i = 0; i < path.size(); i++) {
        (*plan)[i] = new double[path[i].size()];
        std::copy(path[i].begin(), path[i].end(), (*plan)[i]);
    }
}

std::vector<TreeNode*> getNeighborsWithinRadius(
    std::vector<TreeNode*>& tree,
    std::vector<double>& q_new,
    double neighborRadius,
    int numofDOFs)
{
    std::vector<TreeNode*> neighbors;
    for (TreeNode* node : tree) {
        double d = distance(node->angles, q_new);
        if (d < neighborRadius) {
            neighbors.push_back(node);
        }
    }
    return neighbors;
}


#include "hnswlib/hnswlib/hnswlib.h" 

// PRM_planner with HNSW
static void PRM_HNSW_planner(
    double* map,
    int x_size,
    int y_size,
    double* armstart_anglesV_rad,
    double* armgoal_anglesV_rad,
    int numofDOFs,
    double*** plan,
    int* planlength)
{
    *plan = NULL;
    *planlength = 0;

    int numSamples = 20000;
    int k = 500;
    int dim = numofDOFs;

    using namespace hnswlib;

    // Create L2 space and index
    L2Space l2space(dim);
    HierarchicalNSW<float> hnsw_index(&l2space, numSamples + 2); // +2 for start & goal

    std::vector<PRMNode> roadmap;
    std::vector<std::vector<float>> node_data; // Keep data for HNSW

    // Add start node
    PRMNode startNode;
    startNode.angles = std::vector<double>(armstart_anglesV_rad, armstart_anglesV_rad + numofDOFs);
    roadmap.push_back(startNode);
    node_data.push_back(std::vector<float>(startNode.angles.begin(), startNode.angles.end()));
    hnsw_index.addPoint(node_data.back().data(), 0);

    // Add goal node
    PRMNode goalNode;
    goalNode.angles = std::vector<double>(armgoal_anglesV_rad, armgoal_anglesV_rad + numofDOFs);
    roadmap.push_back(goalNode);
    node_data.push_back(std::vector<float>(goalNode.angles.begin(), goalNode.angles.end()));
    hnsw_index.addPoint(node_data.back().data(), 1);

    // Sample valid configurations
    for (int i = 0; i < numSamples; i++) {
        std::vector<double> q = randomConfig(numofDOFs, map, x_size, y_size);
        PRMNode node;
        node.angles = q;
        roadmap.push_back(node);
        node_data.push_back(std::vector<float>(q.begin(), q.end()));
        hnsw_index.addPoint(node_data.back().data(), i + 2);
    }

    // Connect using HNSW k-NN
    for (int i = 0; i < roadmap.size(); ++i) {
		auto result = hnsw_index.searchKnn(node_data[i].data(), k + 1);

        while (!result.empty()) {
            auto [dist, j] = result.top();
            result.pop();
            if (i == j) continue; // skip self
            if (IsValidArmConfiguration(roadmap[i].angles.data(), numofDOFs, map, x_size, y_size) &&
                IsValidArmConfiguration(roadmap[j].angles.data(), numofDOFs, map, x_size, y_size) &&
                IsValidLineSegment(roadmap[i].angles[0], roadmap[i].angles[1],
                                   roadmap[j].angles[0], roadmap[j].angles[1],
                                   map, x_size, y_size))
            {
                roadmap[i].neighbors.push_back(j);
                roadmap[j].neighbors.push_back(i);
            }
        }
    }

    // Dijkstra from start(0) to goal(1)
    int N = roadmap.size();
    std::vector<double> cost(N, std::numeric_limits<double>::infinity());
    std::vector<int> prev(N, -1);
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> pq;
    cost[0] = 0.0;
    pq.emplace(0.0, 0);

    while (!pq.empty()) {
        auto [c, u] = pq.top(); pq.pop();
        if (c > cost[u]) continue;
        if (u == 1) break;
        for (int v : roadmap[u].neighbors) {
            double weight = distance(roadmap[u].angles, roadmap[v].angles);
            if (cost[u] + weight < cost[v]) {
                cost[v] = cost[u] + weight;
                prev[v] = u;
                pq.emplace(cost[v], v);
            }
        }
    }

    if (cost[1] == std::numeric_limits<double>::infinity()) return;

    std::vector<int> pathIndices;
    for (int cur = 1; cur != -1; cur = prev[cur]) pathIndices.push_back(cur);
    std::reverse(pathIndices.begin(), pathIndices.end());

    *planlength = pathIndices.size();
    *plan = new double*[*planlength];
    for (int i = 0; i < *planlength; ++i) {
        (*plan)[i] = new double[numofDOFs];
        std::copy(roadmap[pathIndices[i]].angles.begin(), roadmap[pathIndices[i]].angles.end(), (*plan)[i]);
    }
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

	if (whichPlanner == RRT) RRT_planner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	else if (whichPlanner == PRM) PRM_planner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	else if (whichPlanner == PRM_HNSW) PRM_HNSW_planner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	else throw std::runtime_error("Invalid planner number!\n");

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
