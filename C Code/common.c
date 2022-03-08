#include <stdio.h>
#include <string.h>
#include "common.h"
#include "ziggurat.h"

/**
 * Check to see if the graph represented by the adjacency matric M
 * is connected or not. 
 * @param n size is nxn
 * @param M is the matric 
 * @return 0 for false, 1 for true.
 * 
 */
int is_connected(int n, int M[n][n]){
    // convert M to upper triangular matric
    // for (int j = 0; j < n; j++){
    //     for (int k = 0; k <= j; k++){
    //         M[j][k] = 0;
    //     }
    // }
    int Map[30] = {0}; // To mark vertices. 

    int c = 1; //count
    int current_traverse[30] = {0};

    current_traverse[0] = 1;

    int done = 0;

    while (done == 0){
        int temp[30] = {0};
        int count = 0;
        //mark all in temp
        for (int i = 0; i < n; i++){
            if (current_traverse[i] > 0){
                //mark
                if (Map[i] != 1){
                    Map[i] = 1;
                    temp[count] = i;
                    count++;
                }
            }
        }
        // check count and terminate.
        if (count == 0){
            break;
        }
        //reset current_traverse
        for (int i = 0; i < n; i++){
            current_traverse[i] = 0;
        }
        //recomupte temp using OR
        for (int i = 0; i < count; i++){
            for (int j = 0; j < n; j++){
                current_traverse[j] = current_traverse[j] + M[temp[i]][j];
            }
        }
    }
    // Check and return 
    double sum = 0;
    for (int i = 0; i < n; i++){
        sum += Map[i];
    }
    if (sum == n){
        return 1;
    }
    else{
        return 0;
    }
}

/**
 * Print a 2d matric
 * @param n size is nxn
 * @param M is the matric 
 */
void print_matric(int n, int M[n][n]){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            printf("%d ", M[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/**
 * Proximity net of the flock encoded as a graph.
 * The graph is encoded as adjacency matrix.
 * @param s The state of the system.
 * @param M The adjacency matric.
 * @param rs The radius of sensing.
 */
void state_to_adjacency_matrix( int n, double sqDistances[n][n], int M[n][n], double rs){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (sqDistances[i][j] < rs * rs){
                M[i][j] = 1;
            }
            else{
                M[i][j] = 0;
            }
        }
        M[i][i] = 0;
    }
}


/**
 * The braking distance between agents i and j.
 * @param s The state of the system.
 * @param i The agent number.
 * @param j The agent number.
 * @param c The constraints.
 * @return The distance required to come to a stop appying maximum acc.
 */
double braking_distance(State *s, int i, int j, Constraint *c)
{  
    double P_ij[2] = {s->x[i] - s->x[j], s->y[i] - s->y[j]}; 
    double v_ij[2] = {s->vx[i] - s->vx[j], s->vy[i] - s->vy[j]};

    double distance = norm(P_ij[0], P_ij[1]);
    double v_normal[2] = { P_ij[0] * v_ij[0] / distance, P_ij[1] * v_ij[1] / distance }; 

    double v_mag = -(P_ij[0] * v_ij[0] + P_ij[1] * v_ij[1] ) / distance;
    double d_brake = 0;
    int steps = 0;

    double rel_v = v_mag;
    double rel_a = c->amax;
    if (v_mag > 0){
        steps =  ceil( rel_v / (rel_a * c->dt) );
        d_brake = c->dt * (steps * rel_v - 0.5 * (steps - 1) * steps * rel_a * c->dt);
    }
    // printf("i:%d,%d, dist: %f, v_mag: %f, steps: %d, d_br: %f\n",i, j, distance, v_mag, steps, d_brake);
    return d_brake;
}

/**
 * The braking distance between agents i and j.
 * @param s The state of the system.
 * @param i The agent number.
 * @param j The agent number.
 * @param pos The Worst case postion
 * @param vel The Worst case Velocity
 * @param c The constraints.
 * @return The distance required to come to a stop appying maximum acc.
 */
double braking_distance_worst(State *s, int i, int j, double pos[2], double vel[2], Constraint *c)
{  
    double P_ij[2] = {s->x[i] - pos[0], s->y[i] - pos[1]}; 
    double v_ij[2] = {s->vx[i] - vel[0], s->vy[i] - vel[1]};

    double distance = norm(P_ij[0], P_ij[1]);

    double v_mag = -(P_ij[0] * v_ij[0] + P_ij[1] * v_ij[1] ) / distance;
    double d_brake = 0;
    int steps = 0;

    double rel_v = v_mag;
    double rel_a = c->amax;
    if (v_mag > 0){
        steps = ceil( rel_v / (rel_a * c->dt) );
        d_brake = c->dt * (steps * rel_v - 0.5 * (steps - 1) * steps * rel_a * c->dt);
    }
    // printf("i:%d,%d, dist: %f, v_mag: %f, steps: %d, d_br: %f\n",i, j, distance, v_mag, steps, d_brake);
    return d_brake;
}

/**
 * The braking distance between agents i and j.
 * @param s The state of the system.
 * @param i The agent number.
 * @param j The agent number.
 * @param pos The Worst case postion
 * @param vel The Worst case Velocity
 * @param c The constraints.
 * @return The distance required to come to a stop as proptional to the
 * square of the projection of the relative velocity in the direction from one 
 * agent to the other.
 */
double braking_distance_propotional(State *s, int i, int j, double pos[2], double vel[2], Constraint *c)
{  
    double P_ij[2] = {s->x[i] - pos[0], s->y[i] - pos[1]}; 
    double v_ij[2] = {s->vx[i] - vel[0], s->vy[i] - vel[1]};

    double distance = norm(P_ij[0], P_ij[1]);
    double d_brake = 0;

    double v_mag = -(P_ij[0] * v_ij[0] + P_ij[1] * v_ij[1] ) / distance;
    if (v_mag > 0){
        d_brake = c->beta * (pow(v_mag, 2) / (4 * c->amax));
    }
    return d_brake;
}

void appendToFile(char * name, double val){
    FILE *fp = fopen(name, "a");
    fprintf(fp, "%.30f\n", val);
    fclose(fp);
}

double pointToRectangle(Rectangle *R, double x, double y){
    double cx = R->cx;
    double cy = R->cy;
    double theta = R->theta; //in 
    double l = R->l;
    double w = R->w;
    //Translation by centre.
    x -= cx;
    y -= cy;
    // printf("x:\t%f y:\t%f\n", x, y);
    //Rotation
    double nx = cos(-theta) * x - sin(-theta) * y;
    double ny = sin(-theta) * x + cos(-theta) * y;
    // printf("nx:\t%f ny:\t%f\n", nx, ny);
    //Distance
    double dx = MAX(fabs(nx) - w / 2, 0);
    double dy = MAX(fabs(ny) - l / 2, 0);
    // printf("dx:\t%f dy:\t%f\n",dx , dy);
    return sqrt (dx * dx + dy * dy);
}

/**
 * Deep clones a state.
 * @param s The state to clone.
 * @return Return a deep clone of the specified state.
 * @remark Caller is responsible for freeing the memory of the clone.
 */
State cloneState(State *s) 
{
    int n = s->numBirds;
    State nextS;
    nextS.numBirds = n;
    nextS.x = (double *)malloc(n * sizeof(double));
    nextS.y = (double *)malloc(n * sizeof(double));
    nextS.vx = (double *)malloc(n * sizeof(double));
    nextS.vy = (double *)malloc(n * sizeof(double));
    nextS.policy = (int *)malloc(n * sizeof(int));

    // BAIG WAS HERE
    malloc_count += 5;

    memcpy(nextS.x, s->x, n * sizeof(double));
    memcpy(nextS.y, s->y, n * sizeof(double));
    memcpy(nextS.vx, s->vx, n * sizeof(double));
    memcpy(nextS.vy, s->vy, n * sizeof(double));

    memcpy(nextS.policy, s->policy, n * sizeof(int));

    return nextS;
}


void countRandom(double *r, int len, double sigma)
{
    int nBin = 22;
    double ticks[nBin];
    int hist[nBin];
    memset(hist, 0, nBin * sizeof(int));
    int k = 0;
    for (int i = -10; i <= 10; i++) {
        ticks[k++] = i*sigma;
    }
    ticks[nBin-1] = DBL_MAX;
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < nBin; j++) {
            if (ticks[j] > r[i]) {
                hist[j]++;
                break;
            }
        }    
    }
    printf("Distribution: \n");
    for (int i = 0; i < nBin; i+=11) {
        for (int j = 0; j < 11; j++) {
            printf("B%d: %d, %.3f, ", (i+j+1), hist[i+j], r[i]);
        }
        printf("\n");
    }
}

/**
 * Adds multiplicative noise to an array.
 * @param x The array to add noise to.
 * @param len The number of elements in the array, starting from x[0], to add noise to.
 * @param noiseParams The parameters of noise.
 */
void addNoiseToVariable(double *x, int len, NoiseParameters *noiseParams, uint32_t *seed, uint32_t kn[128], float fn[128], float wn[128])
{
    if (!noiseParams->enabled) {
        return;
    }
    
    if (!noiseParams->independent) {
        double noise = 0;
        if (noiseParams->uniform) {
            noise = randd(noiseParams->lower, noiseParams->upper);
        } else {
            noise = r4_nor(seed, kn, fn, wn) * noiseParams->sigma + noiseParams->mu;
        }
        // if (noise < noiseParams->lower) {
            // noise = noiseParams->lower;
        // }
        // if (noise > noiseParams->upper) {
            // noise = noiseParams->upper;
        // }
        for (int i = 0; i < len; i++) {
            x[i] += noise;
        }
    } else {
        double r[len];
        if (noiseParams->uniform) {
            randdVector(noiseParams->lower, noiseParams->upper, len, r);
        } else {
            for (int i = 0; i < len; i++) {
                r[i] = r4_nor(seed, kn, fn, wn) * noiseParams->sigma + noiseParams->mu;
            }
            // normrndVector(noiseParams->mu, noiseParams->sigma, seed, len, r);
            // printf("mu: %.3f, sigma: %.3f\n", noiseParams->mu, noiseParams->sigma);
            //countRandom(r, len, noiseParams->sigma);
        }
        for (int i = 0; i < len; i++) {
            // if (r[i] < noiseParams->lower) {
                // r[i] = noiseParams->lower;
            // }
            // if (r[i] > noiseParams->upper) {
                // r[i] = noiseParams->upper;
            // }
            x[i] += r[i];
        }
    }
}

/**
 * Adds noise to a specified state.
 * @param s The state to add noise to.
 * @param noiseCfg The noise configuration.
 */
void addNoiseToState(State *s, NoiseConfig *noiseCfg, uint32_t *seed, uint32_t kn[128], float fn[128], float wn[128])
{
    addNoiseToVariable(s->x, s->numBirds, &noiseCfg->positionSensorNoise, seed, kn, fn, wn);
    addNoiseToVariable(s->y, s->numBirds, &noiseCfg->positionSensorNoise, seed, kn, fn, wn);
    addNoiseToVariable(s->vx, s->numBirds, &noiseCfg->velocitySensorNoise, seed, kn, fn, wn);
    addNoiseToVariable(s->vy, s->numBirds, &noiseCfg->velocitySensorNoise, seed, kn, fn, wn); 
}

/**
 * Computes pair-wise squared distances between every two birds at a specified state.
 * @param s The state at which to compute the pair-wise squared distances.
 * @param sqDistances The adjacency matrix to store the pair-wise squared distances.
 *        sqDistances[i][j] is the squared distance between bird i and bird j
 *        sqDistances is a symmetric matrix, i.e., sqDistances[i][j] == sqDistances[j][i]
 *        sqDistances[i][i] == 0
 */
void computeSquaredDistances(State *s, double sqDistances[s->numBirds][s->numBirds])
{
    double dx = 0;
    double dy = 0;
    for (int i = 0; i < s->numBirds; i++){
        sqDistances[i][i] = 0;
        for (int j = i + 1; j < s->numBirds; j++){
            dx = s->x[i] - s->x[j];
            dy = s->y[i] - s->y[j];
            sqDistances[i][j] = dx * dx + dy * dy;
            sqDistances[j][i] = sqDistances[i][j]; 
        }
    }  
}

/**
 * Computes pair-wise squared distances between every two birds in a neighborhood of a specified agent.
 * @param s The state at which to compute the pair-wise squared distances.
 * @param sqDistances The adjacency matrix to store the pair-wise squared distances.
 *        sqDistances[i][j] is the squared distance between bird i and bird j
 *        sqDistances is a symmetric matrix, i.e., sqDistances[i][j] == sqDistances[j][i]
 * @param birdIdx The index of the specified agent.
 * @param neighborIds The array containing the indices of the neighbors of the specified agent.
 * @param numNeighbors The number of neighbors of the specified agent.
 */
void computeSquaredDistancesInNeighborhood(State *s, double sqDistances[s->numBirds][s->numBirds], int birdIdx, int *neighborIds, int numNeighbors)
{
    int ni = 0;
    int nj = 0;
    double dx = 0;
    double dy = 0;
    
    for (int i = 0; i < numNeighbors; i++) {
        ni = neighborIds[i]; // index of i-th neighbor of the specified agent
        // Computes distances between the specified agent and its neighbors
        dx = s->x[birdIdx] - s->x[ni];
        dy = s->y[birdIdx] - s->y[ni];
        sqDistances[birdIdx][ni] = dx * dx + dy * dy;
        sqDistances[ni][birdIdx] = sqDistances[birdIdx][ni]; 
        // Computes the distances between neighbors themselves
        for (int j = i + 1; j < numNeighbors; j++) {
            nj = neighborIds[j];
            dx = s->x[nj] - s->x[ni];
            dy = s->y[nj] - s->y[ni];
            sqDistances[nj][ni] = dx * dx + dy * dy;
            sqDistances[ni][nj] = sqDistances[nj][ni]; 
        }
    }  
}

/**
 * Computes pair-wise distances between every two birds at a specified state.
 * @param s The state at which to compute the pair-wise distances.
 * @param distances The adjacency matrix to store the pair-wise distances.
 *        distances[i][j] is the distance between bird i and bird j
 *        distances is a symmetric matrix, i.e., distances[i][j] == distances[j][i]
 *        distances[i][i] == 0
 */
void computeDistances(State *s, double distances[s->numBirds][s->numBirds])
{
    for (int i = 0; i < s->numBirds; i++){
        distances[i][i] = 0;
        for (int j = i + 1; j < s->numBirds; j++){
            distances[i][j] = norm(s->x[i] - s->x[j], s->y[i] - s->y[j]);
            distances[j][i] = distances[i][j]; 
        }
    }        
}

/**
 * Computes the number of neighbors and the neighbor indices for all agents given a squared distances matrix and a radius.
 * @param s The state at which to compute the neighbor information.
 * @param neighbors The N-by-N matrix to store the neighbor indices, where N is the number of agents.
 *        neighbors[i][j] is the index of j-th neighbor of agent i.
 * @param numNeighbors The array of size N to store the number of neighbors for each agent, where N is the number of agents.
 *        numNeighbors[i] is the number of neighbors of agent i.
 * @param sqDistances The adjacency matrix containing the pair-wise squared distances.
 *        sqDistances[i][j] is the squared distance between bird i and bird j
 *        sqDistances is a symmetric matrix, i.e., sqDistances[i][j] == sqDistances[j][i]
 *        sqDistances[i][i] can be undefined.
 * @param radius The radius to determine neighborhood. It can be thought of as the sensing range of each agent.
 */
void computeNeighbors(State *s, int neighbors[s->numBirds][s->numBirds], int numNeighbors[s->numBirds], double sqDistances[s->numBirds][s->numBirds], double radius)
{
    double sqRadius = radius * radius;
    int n = s->numBirds;
    memset(numNeighbors, 0, n * sizeof(int));
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (sqDistances[i][j] <= sqRadius) {
                neighbors[i][numNeighbors[i]] = j;
                neighbors[j][numNeighbors[j]] = i;
                numNeighbors[i]++;
                numNeighbors[j]++;
            }
        }
    }
}

int compare (const void * a, const void * b)
{
  if (*(double*)a > *(double*)b)
    return 1;
  else if (*(double*)a < *(double*)b)
    return -1;
  else
    return 0;  
}

/**
 * Computes the neighbor indices for all agents for n-nearest neighbors.
 * @param s The state at which to compute the neighbor information.
 * @param neighbors The N-by-N matrix to store the neighbor indices, where N is the number of agents.
 *        neighbors[i][j] is the index of j-th neighbor of agent i.
 * @param numNeighbors The array of size N to store the number of neighbors for each agent, where N is the number of agents.
 *        numNeighbors[i] is the number of neighbors of agent i.
 * @param sqDistances The adjacency matrix containing the pair-wise squared distances.
 *        sqDistances[i][j] is the squared distance between bird i and bird j
 *        sqDistances is a symmetric matrix, i.e., sqDistances[i][j] == sqDistances[j][i]
 *        sqDistances[i][i] can be undefined.
 * @param radius The radius to determine neighborhood. It can be thought of as the sensing range of each agent.
 */
void getKNN(State *s, int neighbors[s->numBirds][s->numBirds], double sqDistances[s->numBirds][s->numBirds], FitnessParameter *param)
{
    int n = s->numBirds;
    double array[n][2]; 
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++){
            array[j][0] = sqDistances[i][j];
            array[j][1] = j;
        }
        qsort(array, n, sizeof array[0], compare);
        for (int k = 0; k < param->knn ; k++){
            neighbors[i][k] = (int) array[k+1][1];
        }
    }
}

/**
 * Applies a specified control action to a given state to advance the system to a new state.
 * @param s The current state.
 * @param a The control action. It is assumed that action satisfies the constraints.
 * @param c The constraints.
 * @return The new state obtained by applying action a to state s.
 */
State dynamics(State *s, Action *a, Constraint *c)
{
    int n = s->numBirds;
    State nextS;
    nextS.numBirds = n;
    nextS.x = (double *)malloc(n * sizeof(double));
    nextS.y = (double *)malloc(n * sizeof(double));
    nextS.vx = (double *)malloc(n * sizeof(double));
    nextS.vy = (double *)malloc(n * sizeof(double));
    nextS.policy = (int *)malloc(n * sizeof(int));
    // initState(&nextS, n);

    // BAIG WAS HERE
    malloc_count += 5;

    for (int i = 0; i < n; i++) {
        nextS.policy[i] = s->policy[i];
        if (c->predator && i == n - 1){
            nextS.vx[i] = s->vx[i] + c->dt * a->ax[i];
            nextS.vy[i] = s->vy[i] + c->dt * a->ay[i];
            trim(&nextS.vx[i], &nextS.vy[i], c->pFactor * c->vmax);
            nextS.x[i] = s->x[i] + c->dt * nextS.vx[i];
            nextS.y[i] = s->y[i] + c->dt * nextS.vy[i];
            continue;
        }
        nextS.vx[i] = s->vx[i] + c->dt * a->ax[i];
        nextS.vy[i] = s->vy[i] + c->dt * a->ay[i];
        trim(&nextS.vx[i], &nextS.vy[i], c->vmax);
        nextS.x[i] = s->x[i] + c->dt * nextS.vx[i];
        nextS.y[i] = s->y[i] + c->dt * nextS.vy[i];
    }
    
    return nextS;
}

/**
 * Applies a specified control action to a given state to advance the system to a new state.
 * @param s The current state.
 * @param a The control action. It is assumed that action satisfies the constraints.
 * @param c The constraints.
 * @return The new state obtained by applying action a to state s.
 */
State dynamics_max_braking(State *s, Action *a, Constraint *c)
{
    int n = s->numBirds;
    State nextS;
    nextS.numBirds = n;
    nextS.x = (double *)malloc(n * sizeof(double));
    nextS.y = (double *)malloc(n * sizeof(double));
    nextS.vx = (double *)malloc(n * sizeof(double));
    nextS.vy = (double *)malloc(n * sizeof(double));
    nextS.policy = (int *)malloc(n * sizeof(int));
    // initState(&nextS, n);

    // BAIG WAS HERE
    malloc_count += 5;

    for (int i = 0; i < n; i++) {
        nextS.policy[i] = s->policy[i];
        // Ensure acceleration doesn't exceed amax
        // Remove this step as the optimizer already trimmed the acceleration
        // trim(&a->ax[i], &a->ay[i], c->amax);
        if(s->policy[i] == 1){
            nextS.vx[i] = s->vx[i] + c->dt * a->ax[i];
            nextS.vy[i] = s->vy[i] + c->dt * a->ay[i];            
            if (s->vx[i] > 0){
                nextS.vx[i] = MAX(0, nextS.vx[i]);
            }else if (s->vx[i] < 0){
                nextS.vx[i] = MIN(0, nextS.vx[i]);
            }
            else{
                nextS.vx[i] = 0;
            }

            if (s->vy[i] > 0){
                nextS.vy[i] = MAX(0, nextS.vy[i]);
            }else if (s->vy[i] < 0){
                nextS.vy[i] = MIN(0, nextS.vy[i]);
            }
            else{
                nextS.vy[i] = 0;
            }
        }
        else{
            nextS.vx[i] = s->vx[i] + c->dt * a->ax[i];
            nextS.vy[i] = s->vy[i] + c->dt * a->ay[i];            
        }
        // Ensure velocity doesn't exceed vmax
        double normV = norm(nextS.vx[i], nextS.vy[i]);
        if (normV > c->vmax) {
            nextS.vx[i] *= c->vmax / normV;
            nextS.vy[i] *= c->vmax / normV;
        }
        //trim(&nextS.vx[i], &nextS.vy[i], c->vmax);
        nextS.x[i] = s->x[i] + c->dt * nextS.vx[i];
        nextS.y[i] = s->y[i] + c->dt * nextS.vy[i];
    }
    
    return nextS;
}

/**
 * Calculates the cost of a specified control action. The cost is defined as the 
 * sum of squares of control action magnitudes.
 * @param a The control action to calculate the cost.
 * @return Return the control action cost, which is the sum of squares of control action magnitudes.
 */
double controlCost(Action *a) {
    double cost = 0.0;
    for (int i = 0; i < a->numBirds; i++) {
        cost += a->ax[i] * a->ax[i] + a->ay[i] * a->ay[i];
    }
    return cost;
}


/**
 * Writes an array of length n to a line in a file stream. 
 * The written line ends with a space and a new line character, i.e., it ends with " \n".
 * @param fp The file stream to write to.
 * @param var The array to write to file.
 * @param n The length of the array.
 */
void logStateVar(FILE *fp, double *var, int n)
{
    for (int i = 0; i < n; i++) {
        fprintf(fp, "%f ", var[i]);
    }
    fprintf(fp, "\n");
}

void logStateVarInt(FILE *fp, int *var, int n)
{    
    for (int i = 0; i < n; i++) {
        fprintf(fp, "%d ", var[i]);
    }
    fprintf(fp, "\n");

}
/**
 * Writes evolutions of state variables to corresponding files specified in log config.
 * @param states The array of states to write to files.
 * @param steps The number of states, i.e., length of the array that states points to.
 * @param logCfg The configuration, e.g., file names.
 */
void logStates(State *states, int steps, LogConfig *logCfg)
{
    FILE *fpx = fopen(logCfg->logFileX, "w+");
    FILE *fpy = fopen(logCfg->logFileY, "w+");
    FILE *fpvx = fopen(logCfg->logFileVx, "w+");
    FILE *fpvy = fopen(logCfg->logFileVy, "w+");
    FILE *fppolicy = fopen(logCfg->logFilePolicy, "w+");
    
    if (fpx == NULL) {
        perror("fopen");
        printf("ERROR opening file x %s\n", logCfg->logFileX);
    }
    if (fpy == NULL) {
        perror("fopen");
        printf("ERROR opening file y %s\n", logCfg->logFileY);
    }
    
    if (fpvx == NULL) {
        perror("fopen");
        printf("ERROR opening file vx %s\n", logCfg->logFileVx);
    }
    
    if (fpvy == NULL) {
        perror("fopen");
        printf("ERROR opening file  vy %s\n", logCfg->logFileVy);
    }

    if (fppolicy == NULL) {
        perror("fopen");
        printf("ERROR opening file p %s\n", logCfg->logFilePolicy);
    }
    
    if (fpx == NULL || fpy == NULL || fpvx == NULL || fpvy == NULL || fppolicy == NULL) {
        if (fpx){
            fclose(fpx);
        }
        if (fpy){
            fclose(fpy);
        }
        if (fpvx){
            fclose(fpvx);
        }
        if (fpvy){
            fclose(fpvy);
        }
        if (fppolicy){
            fclose(fppolicy);
        }
        return;
    }

    State *s;
    int n = states->numBirds;
    for (int i = 0; i < steps; i++) {
        s = &states[i];
        logStateVar(fpx, s->x, n);
        logStateVar(fpy, s->y, n);
        logStateVar(fpvx, s->vx, n);
        logStateVar(fpvy, s->vy, n);
        logStateVarInt(fppolicy, s->policy, n);
    }

    fclose(fpx);
    fclose(fpy);
    fclose(fpvx);
    fclose(fpvy);
    fclose(fppolicy);
}

/**
 * Writes the sequence of control actions to files specified in log config.
 * @param actions The sequence of control actions to write to files.
 * @param steps The number of control actions, i.e., the length of the actions array.
 * @param logCfg The configuration, e.g., file names.
 */
void logActions(Action *actions, int steps, LogConfig *logCfg)
{
    FILE *fpax = fopen(logCfg->logFileAx, "w+");
    FILE *fpay = fopen(logCfg->logFileAy, "w+");
    if (fpax == NULL) {
        perror("fopen");        
        printf("ERROR opening file %s\n", logCfg->logFileAx);
    }
    
    if (fpay == NULL) {
        perror("fopen");
        printf("ERROR opening file %s\n", logCfg->logFileAy);
    }
    
    if (fpax == NULL || fpay == NULL) {
        if (fpax){
            fclose(fpax);
        }        
        if (fpay){
            fclose(fpay);
        }
        return;
    }
    
    Action *a;
    int n = actions->numBirds;
    for (int i = 0; i < steps; i++) {
        a = &actions[i];
        logStateVar(fpax, a->ax, n);
        logStateVar(fpay, a->ay, n);
    }
    fclose(fpax);
    fclose(fpay);
}

/**
 * Writes the evolution of fitness values to a file sepcified in log config.
 * @param fitnesses The array of fitness values.
 * @param steps The number of fitness values, i.e., the length of fitnesses array.
 * @param logCfg The configuration, e.g., the file name.
 */
void logFitnesses(double *fitnesses, int steps, LogConfig *logCfg)
{
    FILE *fp = fopen(logCfg->logFileFitness, "w+");
        
    if (fp == NULL) {
        printf("ERROR opening file %s\n", logCfg->logFileFitness);
        return;
    }
    
    for (int i = 0; i < steps; i++)
    {
        fprintf(fp, "%f\n", fitnesses[i]);
    }
    fclose(fp);
}

/**
 * Writes the evolution of lie-derivative values to a file sepcified in log config.
 * @param LD The array of fitness values.
 * @param steps The number of fitness values, i.e., the length of fitnesses array.
 * @param logCfg The configuration, e.g., the file name.
 */
void logLD(double *LD, int steps, LogConfig *logCfg)
{
    FILE *fp = fopen(logCfg->logFileLieDerivative, "w+");
        
    if (fp == NULL) {
        printf("ERROR opening file %s\n", logCfg->logFileFitness);
        return;
    }
    
    for (int i = 0; i < steps; i++)
    {
        fprintf(fp, "%f\n", LD[i]);
    }
    fclose(fp);
}

int getNumberOfCores() {
#ifdef WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#elif MACOS
    int nm[2];
    size_t len = 4;
    uint32_t count;
 
    nm[0] = CTL_HW; nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);
 
    if(count < 1) {
    nm[1] = HW_NCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);
    if(count < 1) { count = 1; }
    }
    return count;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}
