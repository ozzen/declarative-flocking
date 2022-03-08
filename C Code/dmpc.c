/**
 * dmpc.c
 *
 * Purpose: Distributed MPC-GD declarative flocking controller.
 *
 * How to compile:
 *    gcc dmpc.c conf.c common.c ziggurat.c -O3 -pthread -o dmpc
 *
 * How to use: run
 *    ./dmpc
 * to print usage.
 */

#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "common.h"
#include "conf.h"
#include "ziggurat.h"
#include "reynolds.h"

// BAIG WAS HERE
long free_count = 0;
long malloc_count = 0;

/**
 * Defines parameters needed to run a simulation.
 */
typedef struct Simulation {
    int simNum;
    char *initFileTemplate;
    char *targetFileTemplate;
    char *resultFileTemplate;
    int horizon;
    int steps;
    Constraint *constraints;
    FitnessParameter *fitnessParams;
    GradientDescentOptions *gdOptions;
    NoiseConfig *noiseCfg;
    double runTimeSeconds;
} Simulation;

double separationConstrained(State *s, double sqDistances[s->numBirds][s->numBirds], int birdIdx, int* neighborIds, FitnessParameter *param, Constraint *c);
double obstacleCost(State *s, Constraint *c, FitnessParameter *param, int birdId);
double targetCost(State *s, int birdIdx, double target[2]);
double predatorAvoidance(State *s, double sqDistances[s->numBirds][s->numBirds], int birdIdx, int * neighborIds, int numNeighbors, int numPreds, FitnessParameter *param, Constraint *c);
double controller_predator(State *s, Action *a, Constraint *c);
double averageSquaredDistance(State *s, double sqDistances[s->numBirds][s->numBirds], FitnessParameter *param, int birdIdx, int *neighborIds, int numNeighbors);


/*
* Calculates the rule-based control action for the predator.
* The predator seeks the flock centre.
*/
double controller_predator(State *s, Action *a, Constraint *c){
    //Compute the flock centre.
    int n = s->numBirds;
    double cx = 0;
    double cy = 0;
    for(int i = 0; i < n - 1; i++){
        cx += s->x[i];
        cy += s->y[i];
    }
    cx /= n - 1;
    cy /= n - 1;
    //Compute acceleration.
    a[0].ax[n - 1] = cx - s->x[n - 1];
    a[0].ay[n - 1] = cy - s->y[n - 1];
    normalize(&a[0].ax[n - 1], &a[0].ay[n - 1]);
    //Trim acceleration
    a[0].ax[n - 1] *= c->pFactor * c->amax;
    a[0].ay[n - 1] *= c->pFactor * c->amax;
    // trim(&a[0].ax[n - 1], &a[0].ay[n - 1], c->pFactor * c->amax);
    return norm(cx - s->x[n - 1], cy - s->y[n - 1]);
}

/**
Calculates the combined fitness for predator - prey flock.
**/
double predatorAvoidance(State *s, double sqDistances[s->numBirds][s->numBirds], int birdIdx, int * neighborIds, int numNeighbors, int numPreds, FitnessParameter *param, Constraint *c){
    int n = s->numBirds;
    double centreX = 0;
    double centreY = 0;
    if (numPreds >= n){
        printf("Error in predatorAvoidance()");
        return 0;
    }

    s->numBirds = n - 1;
    double sqD[s->numBirds][s->numBirds];
    computeSquaredDistances(s, sqD);
    double asd = averageSquaredDistance(s, sqD, param, birdIdx, neighborIds, numNeighbors);
    double spC = separationConstrained(s, sqD, birdIdx, neighborIds, param, c);
    s->numBirds = n;

    double pa = pow(MAX(c->predAvD - sqrt(sqDistances[birdIdx][n]), 0), 2);

    double constraint_violations = sqrt(pow(spC, 2) + pa);
    return  param->wd * asd +  param->wspc * constraint_violations;
    // return  param->wd * asd +  param->wspc * spC;

}

double separationConstrained(State *s, double sqDistances[s->numBirds][s->numBirds], int birdIdx, int* neighborIds, FitnessParameter *param, Constraint *c){
    int n = s->numBirds;
    // printf("n in spC: %d\n", n);
    double sp = 0;
    for (int i = 0; i < param->knn ; i++){
        sp += pow(MAX(c->dmin - sqrt(sqDistances[birdIdx][neighborIds[i]]), 0), 2);
    }
    // count = n * (n-1) / 2;
    return sqrt(sp);// / (double)count;
}

/**
 * Calculates sum of minimum distance constraint violations of the flock and rectangular obstacles.
 * @param s The state at which to compute average distance.
 * @param R The rectangular obstacle.
 * @param param The fitness parameters.
 * @return Returns the average distance of the whole flock at the Rectangle.
 */
double obstacleCost(State *s, Constraint *c, FitnessParameter *param, int birdId)
{
    Obstacles *O = c->obs;
    int n = s->numBirds;
    int numO = O->n;
    int count = 0;
    double dist = 0;
    // double sqRadius = param->rc * param->rc;
    for (int k = 0; k < numO; k++){
            double p2R = pointToRectangle(&(O->R[k]), s->x[birdId], s->y[birdId]);
            if (p2R >= 0){
                dist += pow(MAX(0, c->dmin - p2R), 2);
            }
    }
    return sqrt(dist);
}


 /* Calculates the separation of agents in the neighborhood of a specified agent.
 * @param s The state to calculate fitness.
 * @param sqDistances The squared distances between any two agents.
 *        sqDistances[i][j] means the squared distance between agent i and agent j.
 *        sqDistances[i][i] can be undefined.
 * @param param The parameters for fitness calculation.
 * @param birdIdx The index of the specified agent to calculate the separation.
 * @param neighborIds The array of neighbor indices of the specified agent.
 * @param numNeighbors The number of neighbors of the specified agent.
 * @return Returns the separation of agents in the neighborhood of the specified agent.
 */
double separation(State *s, double sqDistances[s->numBirds][s->numBirds], FitnessParameter *param, int birdIdx, int* neighborIds, int numNeighbors)
{
    int count = 0;
    double ret = 0;
    if (numNeighbors <= 0) {
        return 0;
    }
    double sp = 0.0f;
    for (int i = 0; i < numNeighbors; i++) {
        sp += 1.0f / sqDistances[birdIdx][neighborIds[i]];
        count++;
    }
    //sp /= numNeighbors;
    // if (count != 0){
    //     sp /= count;
    // }
    // sp /= numNeighbors;
    return sp;
}

/*
 * Calculates the average squared distance of agents in the neighborhood of a specified agent.
 * @param s The state to calculate fitness.
 * @param sqDistances The squared distances between any two agents.
 *        sqDistances[i][j] means the squared distance between agent i and agent j.
 *        sqDistances[i][i] can be undefined.
 * @param param The parameters for fitness calculation.
 * @param birdIdx The index of the specified agent to calculate the average squared distance.
 * @param neighborIds The array of neighbor indices of the specified agent.
 * @param numNeighbors The number of neighbors of the specified agent.
 * @return Returns the average squared distance of agents in the neighborhood of the specified agent.
 */
double averageSquaredDistance(State *s, double sqDistances[s->numBirds][s->numBirds], FitnessParameter *param, int birdIdx, int *neighborIds, int numNeighbors)
{
    // printf("knn, %d\n", numNeighbors);
    if (numNeighbors <= 0) {
        return 0;
    }
    double avgSqDist = 0.0f;
    for (int i = 0; i < numNeighbors; i++) {
        avgSqDist += sqDistances[birdIdx][neighborIds[i]];
    }
    avgSqDist /= numNeighbors;
    return avgSqDist;
}

/**
 * Calculates the squared diameter of the neighborhood of a specified agent.
 * @param s The state to calculate fitness.
 * @param sqDistances The squared distances between any two agents.
 *        sqDistances[i][j] means the squared distance between agent i and agent j.
 *        sqDistances[i][i] can be undefined.
 * @param param The parameters for fitness calculation.
 * @param birdIdx The index of the specified agent whose neighborhood's diameter will be computed.
 * @param neighborIds The array of neighbor indices of the specified agent.
 * @param numNeighbors The number of neighbors of the specified agent.
 * @return Returns the squared diameter of the neighborhood of the specified agent.
 */
double squaredNeighborhoodDiameter(State *s, double sqDistances[s->numBirds][s->numBirds], FitnessParameter *param, int birdIdx, int *neighborIds, int numNeighbors)
{
    if (numNeighbors <= 0) {
        return 0;
    }
    double sqDiameter = -DBL_MAX;
    // Check distances between agent i and its neighbors
    for (int i = 0; i < numNeighbors; i++) {
        if (sqDiameter < sqDistances[birdIdx][neighborIds[i]]) {
            sqDiameter = sqDistances[birdIdx][neighborIds[i]];
        }
    }
    // Check distances among neighbors
    for (int i = 0; i < numNeighbors; i++) {
        for (int j = i + 1; j < numNeighbors; j++) {
            int ni = neighborIds[i];
            int nj = neighborIds[j];
            if (sqDiameter < sqDistances[ni][nj]) {
                sqDiameter = sqDistances[ni][nj];
            }
        }
    }
    return sqDiameter;
}

/**
 * Calculates the distance from the birdIdx to the target.
 * @param s The state to calculate fitness.
 * @param birdIdx The index of the specified agent to calculate the fitness of its neighborhood.
 * @param target The loaction of the goal.
 * @return The squared distance from the bird to the targer.
 */
double targetCost(State *s, int birdIdx, double target[2])
{
    double ret = pow(s->x[birdIdx] - target[0], 2) + pow(s->y[birdIdx] - target[1], 2);
    return ret;
}

/**
 * Calculates the fitness value of a neighborhood of a specified agent at a given state.
 * @param s The state to calculate fitness.
 * @param sqDistances The squared distances between any two agents.
 *        sqDistances[i][j] means the squared distance between agent i and agent j.
 *        sqDistances[i][i] can be undefined.
 * @param param The parameters for fitness calculation.
 * @param birdIdx The index of the specified agent to calculate the fitness of its neighborhood.
 * @param neighborIds The array of  neighbor indices of the specified agent.
 * @param numNeighbors The number of neighbors of the specified agent.
 * @return The fitness of the neighborhood of the specified agent at the specified state.
 */
double fitness(State *s, double sqDistances[s->numBirds][s->numBirds], FitnessParameter *param, Constraint *c, int birdIdx, int *neighborIds, int numNeighbors)
{
    /// Fitness function based on average squared distance in the neighborhood and separation
    // double target[2] = {60, 50};
    double asd = averageSquaredDistance(s, sqDistances, param, birdIdx, neighborIds, numNeighbors);
    double sp = separation(s, sqDistances, param, birdIdx, neighborIds, numNeighbors);
    // double spC = separationConstrained(s, sqDistances, birdIdx, neighborIds, param, c);
    // double obst = obstacleCost(s, c, param, birdIdx);
    // double tgt = targetCost(s, birdIdx, target);
    // double penalty = sqrt( pow(spC, 2 ) + pow(obst, 2) );

    double f = param->wd * asd + param->ws * sp;
    // double f = param->wd * asd + param->wspc * spC;
    
    // double f = param->wd * asd + param->wspc * penalty + param->wt * tgt;

    // double f = predatorAvoidance(s, sqDistances, birdIdx, neighborIds, numNeighbors, 1, param, c);

    return f;

    /// Fitness function based on squared neighborhood diameter and separation
    // double sqDiameter = squaredNeighborhoodDiameter(s, sqDistances, param, birdIdx, neighborIds, numNeighbors);
    // double sp = separation(s, sqDistances, param, birdIdx, neighborIds, numNeighbors);
    // double f = sqDiameter + sp;
    // return f;
}

/**
 * Calculates the cost of applying a sequence of control actions to a specified bird.
 * The neighbors choose optimal action sequence w.r.t the central agent.
 * @param s The current system state containing states of all agents.
 * @param a The sequence of control actions of the specified bird and its neighbors.
 *          The length of the sequence (size of array) is the horizon.
 * @param h The prediction horizon.
 * @param c The system constraints.
 * @param sqDistances The squared distances between any two agents.
 *        sqDistances[i][j] means the squared distance between agent i and agent j.
 *        sqDistances[i][i] can be undefined.
 * @param param The parameters for fitness calculation.
 * @param birdIdx The index of the specified bird to calculate the cost.
 * @param neighborIds The array of neighbor indices of the specified bird.
 * @param numNeighbors The number of neighbors of the specified bird.
 * @return Returns the sum of the fitness values over the horizon.
 */
double costSumOptim(State *s, Action *a, int h, Constraint *c, FitnessParameter *param, int birdIdx, int *neighborIds, int numNeighbors)
{
    // Make a copy of the current state to play with
    State cs = cloneState(s);
    int n = s->numBirds;
    double sqDistances[n][n];
    computeSquaredDistances(s, sqDistances);
    double cost = 0;
    for (int i = 0; i < h; i++) {
        // Apply acceleration to the specified bird
        cs.vx[birdIdx] += a[i].ax[0];
        cs.vy[birdIdx] += a[i].ay[0];
        // Make sure velocity does not exceed vmax
        trim(&cs.vx[birdIdx], &cs.vy[birdIdx], c->vmax);
        cs.x[birdIdx] += cs.vx[birdIdx];
        cs.y[birdIdx] += cs.vy[birdIdx];
        // Assume neighbors keep the same velocities during the horizon
        for (int j = 0; j < numNeighbors; j++) {
            int nIdx = neighborIds[j];
            cs.vx[nIdx] += a[i].ax[j+1];
            cs.vy[nIdx] += a[i].ay[j+1];
            trim(&cs.vx[nIdx], &cs.vy[nIdx], c->vmax);

            cs.x[nIdx] += cs.vx[nIdx];
            cs.y[nIdx] += cs.vy[nIdx];
        }
        // Calculate fitness at the predicted state and add it to the cost
        computeSquaredDistancesInNeighborhood(&cs, sqDistances, birdIdx, neighborIds, numNeighbors);
        cost += fitness(&cs, sqDistances, param, c, birdIdx, neighborIds, numNeighbors);
        cost += param->wu * (a[i].ax[0] * a[i].ax[0] + a[i].ay[0] * a[i].ay[0]); // control input term
    }
    freeState(&cs);
    return cost;
}

/**
 * Calculates the cost of applying a sequence of control actions to a specified bird.
 * Neighbors move with linear velocities.
 * @param s The current system state containing states of all agents.
 * @param a The sequence of control actions of the specified bird.
 *          The length of the sequence (size of array) is the horizon.
 * @param h The prediction horizon.
 * @param c The system constraints.
 * @param sqDistances The squared distances between any two agents.
 *        sqDistances[i][j] means the squared distance between agent i and agent j.
 *        sqDistances[i][i] can be undefined.
 * @param param The parameters for fitness calculation.
 * @param birdIdx The index of the specified bird to calculate the cost.
 * @param neighborIds The array of neighbor indices of the specified bird.
 * @param numNeighbors The number of neighbors of the specified bird.
 * @return Returns the sum of the fitness values over the horizon.
 */
double costSum(State *s, Action *a, int h, Constraint *c, FitnessParameter *param, int birdIdx, int *neighborIds, int numNeighbors)
{
    // printf("%d, cs:1\n", birdIdx);
    // Make a copy of the current state to play with
    State cs = cloneState(s);
    int n = s->numBirds;
    double sqDistances[n][n];
    computeSquaredDistances(s, sqDistances);
    double cost = 0;
    int controlSteps = c->ct/c->dt;
    // printf("%d, cs:2\n", birdIdx);

    for (int i = 0; i < h; i++) {
        // Apply acceleration to the specified bird
        for (int k = 0; k<controlSteps; k++){
            cs.vx[birdIdx] += c->dt * a[i].ax[0];
            cs.vy[birdIdx] += c->dt * a[i].ay[0];
            // Make sure velocity does not exceed vmax
            // printf("%d, cs:2.1\n", birdIdx);

            trim(&cs.vx[birdIdx], &cs.vy[birdIdx], c->vmax);
            cs.x[birdIdx] += c->dt * cs.vx[birdIdx];
            cs.y[birdIdx] += c->dt * cs.vy[birdIdx];
            // printf("%d, cs:2.11\n", birdIdx);
            // printf("num_n: %d, cs:2.12\n", numNeighbors);

            // Assume neighbors keep the same velocities during the horizon
            for (int j = 0; j < numNeighbors; j++) {
                int nIdx = neighborIds[j];
                // printf("n_id: %d, cs:2.13\n", nIdx);
                cs.x[nIdx] += c->dt * cs.vx[nIdx];
                cs.y[nIdx] += c->dt * cs.vy[nIdx];
            }
            // printf("%d, cs:2.2\n", birdIdx);
    
            // Model the behavior of the predator, if predator is present
            if (c->predator){
                trim(&cs.vx[n-1], &cs.vy[n-1], c->pFactor * c->vmax);
                cs.x[n-1] += c->dt * cs.vx[n-1];
                cs.y[n-1] += c->dt * cs.vy[n-1];                
            }
        }
        // printf("%d, cs:3\n", birdIdx);
    
        // Calculate fitness at the predicted state and add it to the cost
        computeSquaredDistancesInNeighborhood(&cs, sqDistances, birdIdx, neighborIds, numNeighbors);
        cost += fitness(&cs, sqDistances, param, c, birdIdx, neighborIds, numNeighbors);
        cost += param->wu * (a[i].ax[0] * a[i].ax[0] + a[i].ay[0] * a[i].ay[0]); // control input term
    }
    freeState(&cs);
    return cost;
}


/**
 * Finds the best sequence of control actions that yields the minimum cost when applied to a specified state.
 * @param s The state.
 * @param a The best sequence of actions found.
 * @param h The horizon --- number of control actions in a.
 * @param c The system constraints.
 * @param param The parameters for fitness calculation.
 * @param options The parameters for gradient descent algorithm.
 * @return Returns true if a sequence of actions is found. Returns false otherwise.
 * @remark The optimizer may fail to find a solution that reduces the fitness at the specified state.
 */
bool optimize(State *s, Action *a, int h, Constraint *c, FitnessParameter *param, GradientDescentOptions *options, int birdIdx, int *neighborIds, int numNeighbors)
{
    // printf("i:%d, n:%d\n", birdIdx, s->numBirds);

    // Initialize a
    for (int i = 0; i < h; i++) {
        memset(a[i].ax, 0, sizeof(double));
        memset(a[i].ay, 0, sizeof(double));
    }
    double dfax[h];  // Store derivative of cost function w.r.t ax
    double dfay[h];  // Store derivative of cost function w.r.t ay
    double currCost;
    double newCost;
    double gradMag = 0;
    // printf("%d, op:2\n", birdIdx);

    for (int i = 0; i < options->maxIters; i++) {
        // Calculate the cost for the current accelerations
        currCost = costSum(s, a, h, c, param, birdIdx, neighborIds, numNeighbors);
        // printf("%d, op:2.1\n", birdIdx);

        gradMag = 0;
        for (int j = 0; j < h; j++) {
            // Calculate derivative of cost function w.r.t ax and ay of the specified agent
            a[j].ax[0] += options->delta;
            newCost = costSum(s, a, h, c, param, birdIdx, neighborIds, numNeighbors);
            a[j].ax[0] -= options->delta;
            dfax[j] = (newCost - currCost) / options->delta;

            a[j].ay[0] += options->delta;
            newCost = costSum(s, a, h, c, param, birdIdx, neighborIds, numNeighbors);
            a[j].ay[0] -= options->delta;
            dfay[j] = (newCost - currCost) / options->delta;

            gradMag += pow(dfax[j], 2) + pow(dfay[j], 2);
        }
    // printf("%d, op:3\n", birdIdx);

        // printf("%d: %f\n", i, sqrt(gradMag));
        if (gradMag <= options->precision) {
            //printf("Bird %d GD terminates after %d iters\n", birdIdx, (i + 1));
            break;
        }

        // Update accelerations
        for (int j = 0; j < h; j++) {
            a[j].ax[0] -= options->learningRate * dfax[j];
            a[j].ay[0] -= options->learningRate * dfay[j];
            // Ensure acceleration doesn't exceed amax
            trim(&a[j].ax[0], &a[j].ay[0], c->amax);
        }
    }

    return true;
}

bool optimize_n(State *s, Action *a, int h, Constraint *c, FitnessParameter *param, GradientDescentOptions *options, int birdIdx, int *neighborIds, int numNeighbors)
{
    // Initialize a, Index 0 is my acc, 1 to end is for neighbors.
    for (int i = 0; i < h; i++) {
        memset(a[i].ax, 0, (numNeighbors + 1) * sizeof(double));
        memset(a[i].ay, 0, (numNeighbors + 1) * sizeof(double));
    }
    // printf("%d\n", numNeighbors);
    // printf("1\n");

    double dfax[h][numNeighbors + 1];  // Store derivative of cost function w.r.t ax
    double dfay[h][numNeighbors + 1];  // Store derivative of cost function w.r.t ay
    double currCost;
    double newCost;
    double gradMag = 0;

    for (int i = 0; i < options->maxIters; i++) {
        // Calculate the cost for the current accelerations
        currCost = costSumOptim(s, a, h, c, param, birdIdx, neighborIds, numNeighbors);
        gradMag = 0;
        for (int k = 0; k < numNeighbors+1; k++) {
            for (int j = 0; j < h; j++) {
                // Calculate derivative of cost function w.r.t ax and ay of the specified agent
                a[j].ax[k] += options->delta;
                newCost = costSumOptim(s, a, h, c, param, birdIdx, neighborIds, numNeighbors);
                a[j].ax[k] -= options->delta;
                dfax[j][k] = (newCost - currCost) / options->delta;

                a[j].ay[k] += options->delta;
                newCost = costSumOptim(s, a, h, c, param, birdIdx, neighborIds, numNeighbors);
                a[j].ay[k] -= options->delta;
                dfay[j][k] = (newCost - currCost) / options->delta;

                gradMag += fabs(dfax[j][k]) + fabs(dfay[j][k]);
            }
        }
        // printf("2\n");

        if (gradMag <= options->precision) {
            //printf("Bird %d GD terminates after %d iters\n", birdIdx, (i + 1));
            break;
        }

        // Update accelerations
        for (int k = 0; k < numNeighbors+1; k++) {
            for (int j = 0; j < h; j++) {
                a[j].ax[k] -= options->learningRate * dfax[j][k];
                a[j].ay[k] -= options->learningRate * dfay[j][k];
                // Ensure acceleration doesn't exceed amax
                trim(&a[j].ax[k], &a[j].ay[k], c->amax);
            }
        }
        // printf("3\n");

    }
    // printf("4");

    return true;
}

/**
 * Computed the worst case action.
 * @param pos the position of the agent
 * @param vel the velocity of the agent
 * @param c the constraints
 * @param param the parameters.
 */

void worst_case(State *s, int i, int j, double pos[2], double vel[2], Constraint *c)
{
    double Centre[2] = { s->x[j] + c->ct * s->vx[j], s->y[j] + c->ct * s->vy[j] };
    double rad = 0.5 * c->amax * c->ct * c->ct;
    double unit[2] = {s->x[i] - Centre[0], s->y[i] - Centre[1]};
    double u[2] = {rad * unit[0] / norm(unit[0], unit[1]),  rad * unit[1] / norm(unit[0], unit[1]) };
    pos[0] = Centre[0] + u[0];
    pos[1] = Centre[1] + u[1];

    double acc[2] = {pos[0] - Centre[0], pos[1] - Centre[1]};
    set_magnitude(&acc[0], &acc[1], c->amax);

    vel[0] = s->vx[j];
    vel[1] = s->vy[j];
    vel[0] += c->dt * acc[0];
    vel[1] += c->dt * acc[1];

    //Constrain velocity
    // trim(&vel[0], &vel[1], c->vmax);

}
/**
 * The forward switching logic for Distributed simplex architecture.
 * @param s The state of the system.
 * @param i The agent number for which FSL is called.
 * @param h The prediction horizon.
 * @param param The parameters.
 * @param c The constraints.
 */
bool FSL(State *s, int i, int h, FitnessParameter *param, ReynoldsParameters *ReynParams, Constraint *c, GradientDescentOptions *options, int *neighborIds, int numNeighbors)
{
    int n = s->numBirds;
    Action a[h];
    for (int i = 0; i < h; i++) {
        initAction(&a[i], n);
    }
    //Simulate Advance Controller.
    // optimize(s, a, h, c, param, options, i, neighborIds, numNeighbors); // Find the best sequence of acceleration for agent i
    Action tempAction;
    initAction(&tempAction, s->numBirds);
    steer_cohesion_reynolds(s, &tempAction, i, param, c);
    a->ax[0] = tempAction.ax[0];
    a->ay[0] = tempAction.ay[0];

    State nextS = cloneState(s);
    int controlSteps = c->ct / c->dt;
    for (int m = 0; m < controlSteps; m++){
        nextS.vx[i] += c->dt * a->ax[0];
        nextS.vy[i] += c->dt * a->ay[0];
        // Ensure velocity doesn't exceed vmax
        trim(&nextS.vx[i], &nextS.vy[i], c->vmax);
        nextS.x[i] += c->dt * nextS.vx[i];
        nextS.y[i] += c->dt * nextS.vy[i];
    }

    bool ret = false;
    //Switching Condition on d_braking
    for (int j = 0; j < n; j++)
    {
        //Check if neighbor
        double dist = norm(s->x[i] - s->x[j], s->y[i] - s->y[j]);
        if (dist < param->rs && i != j ){
            //Worst Action for neighbor
            double pos[2];
            double vel[2];
            worst_case(&nextS, i, j, pos, vel, c);

            double d_ij_wc = norm(s->x[j] - pos[0], s->y[j] - pos[1]);
            double d_ij_ac = norm(s->x[i] - nextS.x[i], s->x[i] - nextS.x[i]);
            double d_ij_o = norm(s->x[i] - s->x[j], s->y[i] - s->y[j]);
            double d_ij = norm(nextS.x[i] - pos[0], nextS.y[i] - pos[1]);
            // double d_ij = norm(NextS.x[i] - NextS.x[j], NextS.y[i] - NextS.y[j]);
            // double d_brake = braking_distance(s, i, j, c);
            double d_brake = braking_distance_propotional(&nextS, i, j, pos, vel, c);
            // if (i == 9){
            //     printf("d_ij_ac: %f, d_ij_wc = %f, d_ij_o: %f, d_ij: %f, d_brake: %f\n", d_ij_ac, d_ij_wc, d_ij_o, d_ij, d_brake);
            // }
            if (d_ij <= c->dmin + d_brake){
                ret = true;
                break;
            }
        }
    }

    //Clean Up
    freeState(&nextS);
    for (int i = 0; i < h; i++)
    {
        freeAction(&a[i]);
    }
    freeAction(&tempAction);
    return ret;
}

/**
 * The reverse switching logic for Distributed simplex architecture.
 * @param s The state of the system.
 * @param i The agent number for which RSL is called.
 * @param h The prediction horizon.
 * @param param The parameters.
 * @param c The constraints.
 */
bool RSL(State *s, int i, int h, FitnessParameter *param, Constraint *c, GradientDescentOptions *options, int *neighborIds, int numNeighbors)
{
    int n = s->numBirds;
    double minDistance = 10000;
    double dist = 0;
    for (int j = 0; j < n; j++)
    {
        dist = norm(s->x[i] - s->x[j], s->y[i] - s->y[j]);
        if (i != j && dist < minDistance)
        {
            minDistance = dist;
        }
    }
    for (int j = 0; j < n; j++)
    {
        if(i != j){
            // double d_brake = braking_distance(s, i, j, c);
            // double d_ij = norm(s->x[i] - s->x[j], s->y[i] - s->y[j]);
            double rel_v = 2 * c->vmax;
            double rel_a = c->amax;
            int steps = ceil( (rel_v) / (rel_a * c->dt) );
            // double d_br_max = c->dt * (steps * rel_v - 0.5 * (steps - 1) * steps * rel_a * c->dt);
            double d_br_max = c->beta * pow(rel_v, 2)/ (4 * c->amax);
            if (minDistance >= c->dmin + d_br_max){
                return true;
            }
        }
    }
    return false;
}

void RandomControl(Action * a, Constraint *c)
{
    memset(a[0].ax, 0, sizeof(double));
    memset(a[0].ay, 0, sizeof(double));

    double mag = randd(0, c->amax);
    double theta = randd(-PI, PI);

    a[0].ax[0] = mag * cos(theta);
    a[0].ay[0] = mag * sin(theta);
}
/**
 * Uses MPC with a specified horizon over a specified number of time steps to
 * control the system starting from a given state.
 * @param s The initial state.
 * @param h The prediction horizon.
 * @param steps The number of simulation time steps.
 * @param param The parameters for fitness calculation.
 * @param options The parameters for gradient descent algorithm.
 * @param logCfg The configurations for logging data to files.
 */
void distributedControl(State *s, int h, int steps, Constraint *c, FitnessParameter *param, GradientDescentOptions *options, LogConfig *logCfg, NoiseConfig *noiseCfg)
{
    char targetFile[500];
    sprintf(targetFile, "experiments_circle/src/target_conf_%d.txt", 1);

    pointObstacles targets;
    initPointObstacles(&targets, 30);
    readConfTarget(targetFile, &targets);
    c->pointObs = &targets;

    /// Random generator initialization
    uint32_t seed = initRandomSeed();
    uint32_t kn[128];
    float fn[128];
    float wn[128];
    r4_nor_setup(kn, fn, wn);
    int numBs = s->numBirds;

    State states[steps]; // Stores the evolution of state.
    Action actions[steps];
    ReynoldsParameters ReynParams;

    char *confFile = "C Code/reynolds.conf";
    int reynSteps;
    Constraint dummy;
    readReynoldsConfig(confFile, &dummy, &ReynParams, &reynSteps);

    int n = s->numBirds;
    states[0] = *s; // initial state
    // There are no accelerations before the initial state
    actions[0].numBirds = n;
    actions[0].ax = (double *)malloc(n*sizeof(double));
    actions[0].ay = (double *)malloc(n*sizeof(double));

    //Baig was here
    malloc_count = malloc_count + 2;

    memset(actions[0].ax, 0, n*8);
    memset(actions[0].ay, 0, n*sizeof(double));

    // Initialize a sequence of accelerations of one(or more if neighbors optimize for the central bird) bird for optimizer to update at each time step
    Action a[h];

    // Action tempAction;
    // initAction(&tempAction, numBs);

    for (int i = 0; i < h; i++) {
        initAction(&a[i], numBs);
        // initAction(&a[i], 1);
    }

    int neighbors[n][n];    
    int numNeighbors[n];
    double sqDistances[n][n];
    int controlStep = c->ct/c->dt;
    clock_t timeStart = 0;
    clock_t timeEnd = 0;
    double total_time = 0;

    for (int t = 1; t < steps; t++) {
        initAction(&actions[t], n);
        State *s = &states[t-1];  // Current state
        if ((t - 1)%controlStep == 0){
            for (int i = 0; i < n; i++) {
                State s_sensed = cloneState(s);
                // addNoiseToState(&s_sensed, noiseCfg, &seed, kn, fn, wn);
                computeSquaredDistances(&s_sensed, sqDistances);
                // At each time step, assume then neighborhoods don't change in the horizon
                computeNeighbors(&s_sensed, neighbors, numNeighbors, sqDistances, param->rs);
                if (c->predator){
                    if (i == n - 1){
                        double dist = 0;
                        dist = controller_predator(&s_sensed, a, c);
                        if (dist < 10){
                            a[0].ax[n-1] = actions[t-1].ax[n-1];
                            a[0].ay[n-1] = actions[t-1].ay[n-1];
                        }

                        actions[t].ax[i] = a[0].ax[n-1];
                        actions[t].ay[i] = a[0].ay[n-1]; 
                        continue;
                    }
                    // s_sensed.numBirds = n - 1;
                    // double sqD[s_sensed.numBirds][s_sensed.numBirds];
                    // computeSquaredDistances(&s_sensed, sqD);
                    getKNN(&s_sensed, neighbors, sqDistances, param);
                    // for (int idx = 0; idx < param->knn; idx++){
                    //     printf("%d, ", neighbors[i][idx]);
                    // }
                    // printf("\n");

                    //MPC controller
                    // printf("%d, dc:4.1\n", i);
                    // s_sensed.numBirds = n;

                    optimize(&s_sensed, a, h, c, param, options, i, &neighbors[i][0], param->knn); // Find the best sequence of acceleration for agent i
                    // optimize(&s_sensed, a, h, c, param, options, i, &neighbors[i][0], numNeighbors[i]); // Find the best sequence of acceleration for agent i
                    actions[t].ax[i] = a[0].ax[0];
                    actions[t].ay[i] = a[0].ay[0]; 

                }
                else{
                    // getKNN(&s_sensed, neighbors, sqDistances, param);  
                    //MPC controller
                    timeStart = clock();
                    // optimize(&s_sensed, a, h, c, param, options, i, &neighbors[i][0], param->knn); // Find the best sequence of acceleration for agent i
                    optimize(&s_sensed, a, h, c, param, options, i, &neighbors[i][0], numNeighbors[i]); // Find the best sequence of acceleration for agent i
                    timeEnd = clock();
                    total_time += (double) (timeEnd - timeStart) / CLOCKS_PER_SEC;
                    actions[t].ax[i] = a[0].ax[0];
                    actions[t].ay[i] = a[0].ay[0]; 
                }


                freeState(&s_sensed);
            }
            // Advance the system
            states[t] = dynamics(s, &actions[t], c);
        }
        else{
            states[t] = dynamics(s, &actions[t-1], c);
            // actions[t].ax = (double *)malloc(n*sizeof(double));
            // actions[t].ay = (double *)malloc(n*sizeof(double));
            memcpy(actions[t].ax, actions[t-1].ax, n*sizeof(double));
            memcpy(actions[t].ay, actions[t-1].ay, n*sizeof(double));
            //Baig was here
            // malloc_count = malloc_count + 2;
        }

    }

    for (int i = 0; i < h; i++) {
        freeAction(&a[i]);
    }

    // freeAction(&tempAction);
    logStates(states, steps, logCfg);
    logActions(actions, steps, logCfg);

    // Free memory
    for (int t = 0; t < steps; t++) {
        if (t > 0) { // Do not free the initial state as it was created outside of this function
            freeState(&states[t]);
        }
        freeAction(&actions[t]);
    }
    printf("avg controller time: %f\n", (double) (total_time / steps) * (c->ct / c->dt)  );
    // freePointObstacles(c->pointObs);
}

/**
 * Runs a simulation with the specified parameters.
 * This method serves as the entry point to a pthread.
 * @param simPtr Pointer to a Simulation struct holding simulation parameters.
 * @return Function returns NULL. However, the runTimeSeconds in the specified
 *         simulation struct will be updated with the running time of the simulation.
 */
void *runSimulation(void *simPtr)
{
    double timeBegin = clock();
    Simulation *sim = (Simulation *)simPtr;
    // Construct file names
    char initFile[500];
    // char targetFile[500];
    char xFile[500];
    char yFile[500];
    char vxFile[500];
    char vyFile[500];
    char axFile[500];
    char ayFile[500];
    char policyFile[500];
    char fitnessFile[500];

    int simNum = sim->simNum;

    sprintf(initFile, sim->initFileTemplate, simNum);
    sprintf(xFile, sim->resultFileTemplate, simNum, "x");
    sprintf(yFile, sim->resultFileTemplate, simNum, "y");
    sprintf(vxFile, sim->resultFileTemplate, simNum, "vx");
    sprintf(vyFile, sim->resultFileTemplate, simNum, "vy");
    sprintf(axFile, sim->resultFileTemplate, simNum, "ax");
    sprintf(ayFile, sim->resultFileTemplate, simNum, "ay");
    sprintf(fitnessFile, sim->resultFileTemplate, simNum, "fitness");
    sprintf(policyFile, sim->resultFileTemplate, simNum, "policy");


    LogConfig logCfg;
    logCfg.logFileX = xFile;
    logCfg.logFileY = yFile;
    logCfg.logFileVx = vxFile;
    logCfg.logFileVy = vyFile;
    logCfg.logFileAx = axFile;
    logCfg.logFileAy = ayFile;
    logCfg.logFileFitness = fitnessFile;
    logCfg.logFilePolicy = policyFile;

    // Read initial state
    int numBirds = readNumBirds(initFile);
    int n = numBirds;
    // double x[n];
    // double y[n];
    // double vx[n];
    // double vy[n];
    // int p[n] = {0};

    State is;
    initState(&is, n);
    readConfState(initFile, &is);

    // Run global MPC controller
    distributedControl(&is, sim->horizon, sim->steps, sim->constraints, sim->fitnessParams, sim->gdOptions, &logCfg, sim->noiseCfg);
    double timeEnd = clock();
    sim->runTimeSeconds = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;
    freeState(&is);

    return NULL;
}

////////////////////////////////////////////////////////////////////////////
// Implementation of function prototypes                                  //
////////////////////////////////////////////////////////////////////////////


/**
 * Calculates the cohesion metric at a specified state.
 * @param s The state at which to compute cohesion metric.
 * @param distances The pair-wise distances matrix.
 * @param param The fitness parameters for calculating cohesion metric.
 * @return Returns the cohesion metric at the specified state.
 */
double cohesion(State *s, double distances[s->numBirds][s->numBirds], FitnessParameter *param)
{
    double c = 0;
    int count = 0;
    int n = s->numBirds;
    double radius = param->rc;
    double minDistance = DBL_MAX;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (distances[i][j] < radius){
                c -= pow(param->mc, -distances[i][j] / radius) - 1;
                count++; //The number of piars that are at most radius apart.
            }
            else if (distances[i][j] < minDistance) {
                minDistance = distances[i][j];
            }
        }
    }
    if (count == 0){
        return minDistance;
    }
    c /= (double)count; //Mean
    return c;
}


/**
 * Calculates the velocity matching metric at a specified state.
 * @param s The state at which to compute velocity matching metric.
 * @param distances The pair-wise distances matrix.
 * @param param The fitness parameters for calculating velocity matching metric.
 * @return Returns the velocity matching metric at the specified state.
 */
double velocityMatching(State *s, double distances[s->numBirds][s->numBirds], FitnessParameter *param)
{
    int n = s->numBirds;
    double sum = 0.0;
    int count = 0;
    double diff = 0;
    double radius = param->rv;
    double minDistance = DBL_MAX;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (distances[i][j] < radius){
                diff = norm(s->vx[i] - s->vx[j], s->vy[i] - s->vy[j])
                     / (norm(s->vx[i], s->vy[i]) + norm(s->vx[j], s->vy[j]));
                sum += diff;
                count++;
            }
            else {
                minDistance = distances[i][j] < minDistance? distances[i][j]: minDistance;
            }
        }
    }
    if (count == 0){
        return minDistance;
    }
    return sum / count;
}

/**
 * Calculates the average distance of the whole flock.
 * @param s The state at which to compute average distance.
 * @param distances The pair-wise distances matrix.
 * @param param The fitness parameters.
 * @return Returns the average distance of the whole flock at the specified state.
 */
double averageDistance(State *s, double distances[s->numBirds][s->numBirds], FitnessParameter *param)
{
    int n = s->numBirds;
    int count = 0;
    double sd = 0;
    // double sqRadius = param->rc * param->rc;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            sd += distances[i][j];
        }
    }
    count = n * (n-1) / 2;
    sd /= (double)count;
    return sd;
}

/**
 * Finds the maximal value in a specified symmetric matrix, ignoring the diagonal elements.
 * @param n The dimension of the matrix.
 * @param matrix The symmetric matrix.
 * @return Returns the maximal value in the specified matrix.
 */
double maxValue(int n, double matrix[n][n])
{
    double maxVal = -DBL_MAX;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (matrix[i][j] > maxVal) {
                maxVal = matrix[i][j];
            }
        }
    }
    return maxVal;
}

void printUsage(char *programName)
{
    char *usage = "%s -b[-sc] init_file_template output_file_template\n\n"
                  "    -b numSims    Number of simulations in batch mode. Required.\n"
                  "    -s startNum   Starting simulation number. Optional. Default is 1.\n"
                  "    -c conf_file  Specify config file to read parameters.\n"
                  "                  Default config file is 'dmpc.conf'.\n"
                  "    init_file_template    Template to look for initial configuration files. Required.\n"
                  "    output_file_template  Template to create file names to write results. Required.\n"
                  "\n"
                  "Example: run 3 simulations\n"
                  "    %s -b 3 -s 6 initstate_%%d.txt output_%%d_%%s.txt\n"
                  "\nThe program will read initial states for 3 simulations from 'initstate_6.txt', "
                  "'initstate_7.txt', and 'initstate_8.txt'.\n"
                  "It will write:\n"
                  "    evolution of x to 'output_6_x.txt', 'output_7_x.txt', and 'output_8_x.txt'\n"
                  "    evolution of y to 'output_6_y.txt', 'output_7_y.txt', and 'output_8_y.txt'\n"
                  "    evolution of vx to 'output_6_vx.txt', 'output_7_vx.txt', and 'output_8_vx.txt'\n"
                  "    evolution of vy to 'output_6_vy.txt', 'output_7_vy.txt', and 'output_8_vy.txt'\n"
                  "    accelerations to 'output_6_ax.txt', 'output_7_ax.txt', 'output_8_ax.txt', "
                  "'output_6_ay.txt', 'output_7_ay.txt', and 'output_8_ay.txt' \n"
                  "    evolution of fitness values to 'output_6_fitness.txt', 'output_7_fitness.txt', and 'output_8_fitness.txt'\n";
    printf(usage, programName, programName);
}

int main(int argc, char *argv[])
{
    printf("----------DMPC\n");
    Obstacles O;
    char* filename1 = "rectangles.txt";
    readRectObstacles(filename1, &O);

    double timeBegin = clock();
    int numSims = 0;
    int startNum = 1;
    char *confFile = "C Code/dmpc.conf";
    char *initFileTemplate = NULL;
    char *resultFileTemplate = NULL;
    bool quietMode = false;

    int opt;
    while ((opt = getopt(argc, argv, "b:s:c:q:h")) != -1) {
        switch (opt) {
        case 'b': numSims = atoi(optarg); break;
        case 's': startNum = atoi(optarg); break;
        case 'c': confFile = optarg; break;
        case 'q': quietMode = true; break;
        case 'h': printUsage(argv[0]); exit(EXIT_SUCCESS);
        default:
            printUsage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (optind < argc - 1) {
        initFileTemplate = argv[optind];
        resultFileTemplate = argv[optind + 1];
    }
    else {
        printUsage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if (numSims < 1) {
        printf("ERROR: Number of simulations (-b) is required and must be greater than zero.\n");
        exit(EXIT_FAILURE);
    }

    // Parse MPC parameters from config file
    Constraint c;
    FitnessParameter param;
    GradientDescentOptions gdOptions;
    NoiseConfig noiseCfg;
    int h;
    int steps;
    int numCpuCores = getNumberOfCores();
    int maxThreads = numCpuCores * 2; // Impose this threshold to conserve resources

    readGlobalMpcConfig(confFile, &c, &param, &gdOptions, &h, &steps);
    c.obs = &O;

    // Manually set noise config
    noiseCfg.positionSensorNoise.enabled = false;
    noiseCfg.positionSensorNoise.uniform = false;
    noiseCfg.positionSensorNoise.independent = true;
    noiseCfg.positionSensorNoise.mu = 0;
    noiseCfg.positionSensorNoise.sigma = 2.0;
    noiseCfg.positionSensorNoise.lower = -2;
    noiseCfg.positionSensorNoise.upper = 2;

    noiseCfg.velocitySensorNoise.enabled = false;
    noiseCfg.velocitySensorNoise.uniform = false;
    noiseCfg.velocitySensorNoise.independent = true;
    noiseCfg.velocitySensorNoise.mu = 0;
    noiseCfg.velocitySensorNoise.sigma = 1.0;
    noiseCfg.velocitySensorNoise.lower = -1;
    noiseCfg.velocitySensorNoise.upper = 1;

    noiseCfg.actuatorNoise.enabled = false;
    noiseCfg.actuatorNoise.uniform = false;
    noiseCfg.actuatorNoise.independent = true;
    noiseCfg.actuatorNoise.mu = 0;
    noiseCfg.actuatorNoise.sigma = 0.02;

    if (!quietMode) {
        printf("horizon = %d \t time steps = %d\n", h, steps);
        printf("amax = %.3f \t vmax = %.3f\n", c.amax, c.vmax);
        printf("rc = %.3f \t wc = %.3f \t mc = %.3f\n", param.rc, param.wc, param.mc);
        printf("rs = %.3f \t ws = %.3f \t ms = %.3f\n", param.rs, param.ws, param.ms);
        printf("rv = %.3f \t wv = %.3f\n", param.rv, param.wv);
        printf("learningRate = %.6f \t delta = %.9f \t maxIters = %d\n", gdOptions.learningRate, gdOptions.delta, gdOptions.maxIters);
    }

    int noiseExps = 1;

    char noiseOutFileTemplate[500];
    for (int nIdx = 1; nIdx <= noiseExps; nIdx++) {
        noiseCfg.positionSensorNoise.sigma = 0.2 * nIdx;
        noiseCfg.velocitySensorNoise.sigma = 0.1 * nIdx;
        if (!quietMode) {
            printf("Noise Experiment #%d: sigma_x = %.3f, sigma_v = %.3f\n", nIdx,
               noiseCfg.positionSensorNoise.sigma, noiseCfg.velocitySensorNoise.sigma);
        }
        sprintf(noiseOutFileTemplate, "experiments/exp3/ndmpc_sn_x%02d_v%02d_%%d_%%s.txt", 2 * nIdx, 1 * nIdx);

        pthread_t simThreads[numSims];
        Simulation simParams[numSims];

        int error_code = 0;
        int finishedThreadIndex = 0;

        for (int i = 0; i < numSims; i++) {
            int simNum = i + startNum;
            if (!quietMode) {
                printf("Creating a pthread for simulation #%d...\n", simNum);
            }
            simParams[i].simNum = simNum;
            simParams[i].initFileTemplate = initFileTemplate;
            simParams[i].resultFileTemplate = resultFileTemplate;//noiseOutFileTemplate; //resultFileTemplate;
            simParams[i].horizon = h;
            simParams[i].steps = steps;
            simParams[i].constraints = &c;
            simParams[i].fitnessParams = &param;
            simParams[i].gdOptions = &gdOptions;
            simParams[i].noiseCfg = &noiseCfg;
            simParams[i].runTimeSeconds = -1;

            error_code = pthread_create(&simThreads[i], NULL, runSimulation, (void *)&simParams[i]);
            if (error_code == EAGAIN || ((i + 1) % maxThreads) == 0) { // system lacks resources, try again
                printf("Waiting for some threads to finish before creating new ones...\n");
                // wait for some threads to finish before creating new ones
                for (; finishedThreadIndex < i - numCpuCores; finishedThreadIndex++) {
                    error_code = pthread_join(simThreads[finishedThreadIndex], NULL);
                    if (error_code) {
                        fprintf(stderr, "Error joining thread for simulation #%d. Error code: %d\n", simNum, error_code);
                        exit(EXIT_FAILURE);
                    }
                }
            }
            else if (error_code) {
                fprintf(stderr, "Error creating thread: %d\n", error_code);
                exit(EXIT_FAILURE);
            }
        }

        // Wait for all threads to finish
        for (int i = finishedThreadIndex; i < numSims; i++) {
            int simNum = i + startNum;
            error_code = pthread_join(simThreads[i], NULL);
            if (error_code) {
                fprintf(stderr, "Error joining thread for simulation #%d. Error code: %d\n", simNum, error_code);
                exit(EXIT_FAILURE);
            }
        }

        if (!quietMode) {
            for (int i = 0; i < numSims; i++) {
               printf("Simulation #%d finished in %f seconds\n", simParams[i].simNum, simParams[i].runTimeSeconds);
            }
        }
    }

    double timeEnd = clock();
    double timeSeconds = (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC;
    if (!quietMode) {
        printf("Total running time = %f seconds\n", timeSeconds);
    }
    printf("malloc: %ld, free: %ld\n", malloc_count, free_count);
}
