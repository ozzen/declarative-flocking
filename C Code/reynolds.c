/**
 * gmpc.c
 * 
 * Purpose: Reynold's rule-based flocking controller.
 *
 * How to compile:
 *    gcc reynolds.c conf.c common.c normal.c -O3 -pthread -o reynolds
 *
 * On Windows, use MinGW-64 to compile:
 *    gcc reynolds.c conf.c common.c -O3 -pthread -o reynolds.exe
 *
 * How to use: run
 *    ./reynolds
 * to print usage.
 */
 
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <getopt.h> 
#include "common.h"
#include "conf.h"
#include "reynolds.h"
#include "ziggurat.h"

// BAIG WAS HERE
// long free_count = 0;
// long malloc_count = 0;

/**
 * Defines parameters needed to run a simulation.
 */
typedef struct Simulation {
    int simNum;
    char *initFileTemplate;
    char *resultFileTemplate;
    int steps;
    Constraint *constraints;
    ReynoldsParameters *reynoldsParams;
    NoiseConfig *noiseCfg;
    double runTimeSeconds;
} Simulation;

/**
 * Calculates the weighted control action (accelerations, or steering forces) using Reynold's cohesion rule
 * to move birds closer together.
 * @param s The current state.
 * @param param The parameters for Reynold's boid model.
 * @return Returns the accelerations for all birds.
 */
Action steer_cohesion(State *s, double distances[s->numBirds][s->numBirds], ReynoldsParameters *param)
{
    int n = s->numBirds;
    Action a;
    initAction(&a, n);
    memset(a.ax, 0, n*sizeof(double));
    memset(a.ay, 0, n*sizeof(double));
    int neighbors;
    // Calculate action for each bird i
    for (int i = 0; i < n; i++) {
        neighbors = 0;
        // Loop through all birds to find neighbors of bird i
        for (int j = 0; j < n; j++) {
            // Only count those in bird i's neighborhood
            if ((i != j) && (distances[i][j] <= param->rc)) {
                // Accumulate sum of neighbor's positions
                a.ax[i] += s->x[j];
                a.ay[i] += s->y[j];
                neighbors++;
            }
        }

        if (neighbors > 0) {
            // Take average positions of the neighbors (centroid)
            a.ax[i] /= (double)neighbors;
            a.ay[i] /= (double)neighbors;
            // Calculate the direction from bird i to the centroid
            a.ax[i] -= s->x[i];
            a.ay[i] -= s->y[i];
            // Normalize to get unit-vector
            // normalize(&a.ax[i], &a.ay[i]);
            // Weigh
            a.ax[i] *= param->wc;
            a.ay[i] *= param->wc;
        }
        
        //printf("Cohesion: Bird %d's # of neighbors: %d, ax = %.4f, ay = %.4f\n", i, neighbors, a.ax[i], a.ay[i]);
    }
    return a;
}
/**
 * Calculates the weighted control action (accelerations, or steering forces) using Reynold's alignment rule.
 * @param s The current state.
 * @param param The parameters for Reynold's boid model.
 * @return Returns the accelerations for all birds.
 */
Action steer_alignment(State *s, double distances[s->numBirds][s->numBirds], ReynoldsParameters *param)
{
    int n = s->numBirds;
    Action a;
    initAction(&a, n);
    memset(a.ax, 0, n*sizeof(double));
    memset(a.ay, 0, n*sizeof(double));
    int neighbors;
    // Calculate action for each bird i
    for (int i = 0; i < n; i++) {
        neighbors = 0;
        // Loop through all birds to find neighbors of bird i
        for (int j = 0; j < n; j++) {
            // Only count those in bird i's neighborhood
            if ((i != j) && (distances[i][j] <= param->rv)) {
                // Accumulate sum of neighbor's heading direction
                a.ax[i] += s->vx[j];
                a.ay[i] += s->vy[j];
                neighbors++;
            }
        }

        if (neighbors > 0) {
            // Take average heading direction
            a.ax[i] /= (double)neighbors;
            a.ay[i] /= (double)neighbors;
            // Calculate the offset between bird i and the average heading direction
            a.ax[i] -= s->vx[i];
            a.ay[i] -= s->vy[i];
            // Normalize to get unit-vector
            // normalize(&a.ax[i], &a.ay[i]);
            // Weigh
            a.ax[i] *= param->wv;
            a.ay[i] *= param->wv;
        }
        //printf("Alignment: Bird %d's # of neighbors: %d, ax = %.4f, ay = %.4f\n", i, neighbors, a.ax[i], a.ay[i]);
    }
    return a;    
}

/**
 * Calculates the weighted control action (accelerations, or steering forces) using Reynold's separation rule.
 * @param s The current state.
 * @param param The parameters for Reynold's boid model.
 * @return Returns the accelerations for all birds.
 */
Action steer_separation(State *s, double distances[s->numBirds][s->numBirds], ReynoldsParameters *param)
{
    int n = s->numBirds;
    Action a;
    initAction(&a, n);
    memset(a.ax, 0, n*sizeof(double));
    memset(a.ay, 0, n*sizeof(double));
    int neighbors;
    // Calculate action for each bird i
    for (int i = 0; i < n; i++) {
        neighbors = 0;
        // Loop through all birds to find neighbors of bird i
        for (int j = 0; j < n; j++) {
            // Only count those in bird i's neighborhood
            if ((i != j) && (distances[i][j] <= param->rs)) {
                a.ax[i] += (s->x[i] - s->x[j]) / distances[i][j] / distances[i][j];
                a.ay[i] += (s->y[i] - s->y[j]) / distances[i][j] / distances[i][j];
                neighbors++;
            }
        }
        // Normalize
        // normalize(&a.ax[i], &a.ay[i]);
        // Weigh
        a.ax[i] *= param->ws;
        a.ay[i] *= param->ws;
        //printf("Separation: Bird %d's # of neighbors: %d, ax = %.4f, ay = %.4f\n", i, neighbors, a.ax[i], a.ay[i]);
    }
    return a;
}

/**
 * Calculates the control action (accelerations, or steering forces) using Reynold's rules.
 * @param s The current state.
 * @param param The parameters for Reynold's boid model.
 * @param c The system constraints.
 * @return Returns the accelerations for all birds.
 */
Action steer_total(State *s, ReynoldsParameters *param, Constraint *c)
{
    Action total;
    initAction(&total, s->numBirds);
    double distances[s->numBirds][s->numBirds];
    computeDistances(s, distances);
    Action ac = steer_cohesion(s, distances, param);
    Action as = steer_separation(s, distances, param);
    Action av = steer_alignment(s, distances, param);
    addActions(&ac, &as, &total);
    addActions(&av, &total, &total);
    freeAction(&ac);
    freeAction(&as);
    freeAction(&av);
    // Ensure acceleration doesn't exceed amax
    for (int i = 0; i < total.numBirds; i++) {
        double normA = norm(total.ax[i], total.ay[i]);
        if (normA > c->amax) {
            total.ax[i] *= c->amax / normA;
            total.ay[i] *= c->amax / normA;
        }
        //trim(&total.ax[i], &total.ay[i], c->amax);
    }
    return total;
}

Action steer_cohesion_distributed(State *s, int i, double distances[s->numBirds][s->numBirds], FitnessParameter *param, Constraint *c)
{
    int n = s->numBirds;
    Action a;
    initAction(&a, n);
    memset(a.ax, 0, n*sizeof(double));
    memset(a.ay, 0, n*sizeof(double));
    int neighbors;
    // Calculate action for each bird i
    neighbors = 0;
    // printf("rs: %f\n");
    // Loop through all birds to find neighbors of bird i
    for (int j = 0; j < n; j++) {
        // Only count those in bird i's neighborhood
        if ((i != j) && (distances[i][j] <= param->rs)) {
            // Accumulate sum of neighbor's positions
            a.ax[0] += s->x[j];
            a.ay[0] += s->y[j];
            neighbors++;
        }
    }
    if (neighbors > 0) {
        // Take average positions of the neighbors (centroid)
        a.ax[0] /= (double)neighbors;
        a.ay[0] /= (double)neighbors;

        // Calculate the direction from bird i to the centroid
        a.ax[0] -= s->x[i];
        a.ay[0] -= s->y[i];
        // Normalize to get unit-vector
        normalize(&a.ax[0], &a.ay[0]);
        // Weigh
        a.ax[0] *= c->amax;
        a.ay[0] *= c->amax;
    }
    //printf("Cohesion: Bird %d's # of neighbors: %d, ax = %.4f, ay = %.4f\n", i, neighbors, a.ax[i], a.ay[i]);
    return a;
}

Action steer_target(State *s, int i, double distances[s->numBirds][s->numBirds], FitnessParameter *param, Constraint *c)
{
    int n = s->numBirds;
    Action a;
    initAction(&a, n);
    memset(a.ax, 0, n*sizeof(double));
    memset(a.ay, 0, n*sizeof(double));

    // Calculate the direction from bird i to the target
    a.ax[0] = c->pointObs->px[i] - s->x[i];
    a.ay[0] = c->pointObs->py[i] - s->y[i];

    // Normalize to get unit-vector
    normalize(&a.ax[0], &a.ay[0]);

    // Weigh
    a.ax[0] *= c->amax;
    a.ay[0] *= c->amax;
    return a;
}

Action steer_separation_distributed(State *s, int i, double distances[s->numBirds][s->numBirds], FitnessParameter *param, Constraint *c)
{
    int n = s->numBirds;
    Action a;
    initAction(&a, n);
    memset(a.ax, 0, n*sizeof(double));
    memset(a.ay, 0, n*sizeof(double));
    int neighbors;
    // Calculate action for each bird i
    neighbors = 0;
    // Loop through all birds to find neighbors of bird i
    for (int j = 0; j < n; j++) {
        // Only count those in bird i's neighborhood
        if ((i != j) && (distances[i][j] <= param->rs)) {
            double d_brake = braking_distance(s, i, j, c);
            // double scale =(distances[i][j] - c->dmin - d_brake) * (distances[i][j] - c->dmin - d_brake);
            double scale = MAX(0.000001, (distances[i][j] - c->dmin - d_brake) * (distances[i][j] - c->dmin - d_brake));
            // double scale = (distances[i][j]) * (distances[i][j]);
            a.ax[0] += (s->x[i] - s->x[j]) / scale;
            a.ay[0] += (s->y[i] - s->y[j]) / scale;
            // a.ax[0] += (s->x[i] - s->x[j]) / distances[i][j] / distances[i][j];
            // a.ay[0] += (s->y[i] - s->y[j]) / distances[i][j] / distances[i][j];
            neighbors++;
        }
    }
    // if (neighbors != 0){
    //     a.ax[0] /= neighbors;
    //     a.ay[0] /= neighbors;
    // }
    // Normalize
    normalize(&a.ax[0], &a.ay[0]);
    // Weigh
    a.ax[0] *= c->amax;
    a.ay[0] *= c->amax;    

    //printf("Separation: Bird %d's # of neighbors: %d, ax = %.4f, ay = %.4f\n", i, neighbors, a.ax[i], a.ay[i]);
    return a;
}

double shouviks_weight(double r, double dmin){
    double y = 1 - pow(M_E, -pow(r - dmin, 2));
    double ret = 1/(1 + y);
    return ret;
}

Action steer_alignment_distributed(State *s, int i, double distances[s->numBirds][s->numBirds], FitnessParameter *param)
{
    int n = s->numBirds;
    Action a;
    initAction(&a, n);
    memset(a.ax, 0, n*sizeof(double));
    memset(a.ay, 0, n*sizeof(double));
    int neighbors;
    // Calculate action for each bird i
    neighbors = 0;
    // Loop through all birds to find neighbors of bird i
    for (int j = 0; j < n; j++) {
        // Only count those in bird i's neighborhood
        if ((i != j) && (distances[i][j] <= param->rs)) {
            // Accumulate sum of neighbor's heading direction
            a.ax[0] += s->vx[j];
            a.ay[0] += s->vy[j];
            neighbors++;
        }
    }

    if (neighbors > 0) {
        // Take average heading direction
        a.ax[0] /= (double)neighbors;
        a.ay[0] /= (double)neighbors;
        // Calculate the offset between bird i and the average heading direction
        a.ax[0] -= s->vx[i];
        a.ay[0] -= s->vy[i];
        // Normalize to get unit-vector
        // normalize(&a.ax[i], &a.ay[i]);
        // Weigh
        a.ax[0] *= param->wv;
        a.ay[0] *= param->wv;
    }
    //printf("Alignment: Bird %d's # of neighbors: %d, ax = %.4f, ay = %.4f\n", i, neighbors, a.ax[i], a.ay[i]);
    return a;    
}

bool steer_total_distributed(State *s, Action *total, int i, FitnessParameter *param, Constraint *c)
{
    //Initialize a
    memset(total->ax, 0, sizeof(double));
    memset(total->ay, 0, sizeof(double));

    // Action total;
    // initAction(&total, s->numBirds);
    double distances[s->numBirds][s->numBirds];
    computeDistances(s, distances);
    Action ac = steer_cohesion_distributed(s, i, distances, param, c);
    Action as = steer_separation_distributed(s, i, distances, param, c);
    Action av = steer_alignment_distributed(s, i, distances, param);

    total->ax[0] = ac.ax[0] + as.ax[0] + av.ax[0];
    total->ay[0] = ac.ay[0] + as.ay[0] + av.ay[0];

    freeAction(&ac);
    freeAction(&as);
    freeAction(&av);
    // Ensure acceleration doesn't exceed amax
    double normA = norm(total->ax[0], total->ay[0]);
    if (normA > c->amax) {
        total->ax[0] *= c->amax / normA;
        total->ay[0] *= c->amax / normA;
    }
    //trim(&total.ax[i], &total.ay[i], c->amax);
    return true;
}


Action steer_separation_closest_distributed(State *s, int i, double distances[s->numBirds][s->numBirds], ReynoldsParameters *param, Constraint *c)
{
    int n = s->numBirds;
    Action a;
    initAction(&a, n);
    memset(a.ax, 0, n*sizeof(double));
    memset(a.ay, 0, n*sizeof(double));
    int neighbors;
// Calculate action for each bird i
    neighbors = 0;
    // Loop through all birds to find neighbors of bird i
    // Find the closest neighbor.
    int closestN = 0;
    double minDist = 100000;
    for (int k = 0; k < n; k++) {
        if (distances[k][i] < param->rs && k != i && distances[k][i] < minDist){
            double d_brake = braking_distance(s, i, k, c);
            // minDist = distances[k][i] - d_brake;
            minDist = distances[k][i];
            closestN = k;
        }
    }
    
    a.ax[0] = (s->x[i] - s->x[closestN]) / distances[i][closestN] / distances[i][closestN];
    a.ay[0] = (s->y[i] - s->y[closestN]) / distances[i][closestN] / distances[i][closestN];
    // Normalize
    // normalize(&a.ax[i], &a.ay[i]);
    // Weigh
    a.ax[0] *= param->ws;
    a.ay[0] *= param->ws;
    //printf("Separation: Bird %d's # of neighbors: %d, ax = %.4f, ay = %.4f\n", i, neighbors, a.ax[i], a.ay[i]);
    return a;
}


bool steer_max_avoidance(State *s, Action *total, int i, ReynoldsParameters *param, Constraint *c)
{
    //Initialize a
    memset(total->ax, 0, sizeof(double));
    memset(total->ay, 0, sizeof(double));

    // Action total;
    // initAction(&total, s->numBirds);
    double distances[s->numBirds][s->numBirds];
    computeDistances(s, distances);
    Action as2 = steer_separation_closest_distributed(s, i, distances, param, c);

    total->ax[0] = as2.ax[0];
    total->ay[0] = as2.ay[0];

    // Ensure acceleration doesn't exceed amax
    double normA = norm(total->ax[0], total->ay[0]);
    if (normA > c->amax) {
        total->ax[0] *= c->amax / normA;
        total->ay[0] *= c->amax / normA;
    }
    freeAction(&as2);
    return true;
}

bool steer_separation_reynolds(State *s, Action *total, int i, FitnessParameter *param, Constraint *c)
{
    //Initialize a
    memset(total->ax, 0, sizeof(double));
    memset(total->ay, 0, sizeof(double));

    // Action total;
    // initAction(&total, s->numBirds);
    double distances[s->numBirds][s->numBirds];
    computeDistances(s, distances);
    Action as2 = steer_separation_distributed(s, i, distances, param, c);

    total->ax[0] = as2.ax[0];
    total->ay[0] = as2.ay[0];

    // Ensure acceleration doesn't exceed amax
    double normA = norm(total->ax[0], total->ay[0]);
    if (normA > c->amax) {
        total->ax[0] *= c->amax / normA;
        total->ay[0] *= c->amax / normA;
    }
    freeAction(&as2);
    return true;
}

bool steer_cohesion_reynolds(State *s, Action *total, int i, FitnessParameter *param, Constraint *c)
{
    //Initialize a
    memset(total->ax, 0, sizeof(double));
    memset(total->ay, 0, sizeof(double));

    // Action total;
    // initAction(&total, s->numBirds);
    double distances[s->numBirds][s->numBirds];
    computeDistances(s, distances);
    // Action as2 = steer_cohesion_distributed(s, i, distances, param, c);
    Action as2 = steer_target(s, i, distances, param, c);

    total->ax[0] = as2.ax[0];
    total->ay[0] = as2.ay[0];

    // Ensure acceleration doesn't exceed amax
    double normA = norm(total->ax[0], total->ay[0]);
    if (normA > c->amax) {
        total->ax[0] *= c->amax / normA;
        total->ay[0] *= c->amax / normA;
    }
    freeAction(&as2);
    return true;
}

bool steer_max_braking(State *s, Action *total, int i, Constraint *c)
{
        //Initialize a
    memset(total->ax, 0, sizeof(double));
    memset(total->ay, 0, sizeof(double));

    double vel[2] = {s->vx[i], s->vy[i]};
    double normVel = norm(vel[0], vel[1]); 

    total->ax[0] = -vel[0] * (c->amax / normVel);
    total->ay[0] = -vel[1] * (c->amax / normVel);

    return true;
}

/**
 * Run Reynold's rule-based controller over a specified number of time steps to 
 * control the system starting from a given state.
 * @param s The initial state.
 * @param steps The number of simulation time steps.
 * @param param The parameters for Reynold's boid model.
 * @param c The system constraints.
 * @param logCfg The configurations for logging data to files.
 * @param noiseCfg The configuration for generating noise.
 */
void reynolds_control(State *s, int steps, ReynoldsParameters *param, Constraint *c, LogConfig *logCfg, NoiseConfig *noiseCfg)
{
    /// Random generator initialization
    uint32_t seed = initRandomSeed();
    uint32_t kn[128];
    float fn[128];
    float wn[128];
    r4_nor_setup(kn, fn, wn);
    
    State states[steps]; // Stores the evolution of state.
    Action actions[steps];
    
    int n = s->numBirds;
    states[0] = *s; // initial state
    // There are no accelerations before the initial state
    initAction(&actions[0], n);
    memset(actions[0].ax, 0, n*sizeof(double));
    memset(actions[0].ay, 0, n*sizeof(double));
    
    printf("here\n");
    for (int t = 1; t < steps; t++) {
        //printf("Time step: %d\n", t);
        printf("inside1\n");
        State *s = &states[t-1];  // Current state
        printf("inside2\n");
        State s_sensed = cloneState(s);
        printf("inside3\n");
        // addNoiseToState(&s_sensed, noiseCfg, &seed, kn, fn, wn);
        actions[t] = steer_total(&s_sensed, param, c);
        printf("inside\n");
        states[t] = dynamics(s, &actions[t], c);
        freeState(&s_sensed);
    }
    
    logStates(states, steps, logCfg);
    logActions(actions, steps, logCfg);
    
    // Free memory
    for (int t = 0; t < steps; t++) {
        if (t > 0) { // Do not free the initial state as it was created outside of this function
            freeState(&states[t]);
        }
        freeAction(&actions[t]);
    }
}

/**
 * Runs a simulation with the specified parameters.
 * This method serves as the entry point to a pthread.
 * @param simPtr Pointer to a Simulation struct holding simulation parameters.
 * @return Function returns NULL. However, the runTimeSeconds in the specified 
 *         simulation struct will be updated with the running time of the simulation.
 */
void *runSimulationR(void *simPtr)
{
    double timeBegin = clock();
    Simulation *sim = (Simulation *)simPtr;
    // Construct file names
    char initFile[500];
    char xFile[500];
    char yFile[500];
    char vxFile[500];
    char vyFile[500];
    char axFile[500];
    char ayFile[500];
    
    int simNum = sim->simNum;

    sprintf(initFile, sim->initFileTemplate, simNum);
    sprintf(xFile, sim->resultFileTemplate, simNum, "x");
    sprintf(yFile, sim->resultFileTemplate, simNum, "y");
    sprintf(vxFile, sim->resultFileTemplate, simNum, "vx");
    sprintf(vyFile, sim->resultFileTemplate, simNum, "vy");
    sprintf(axFile, sim->resultFileTemplate, simNum, "ax");
    sprintf(ayFile, sim->resultFileTemplate, simNum, "ay");
    
    LogConfig logCfg;
    logCfg.logFileX = xFile;
    logCfg.logFileY = yFile;
    logCfg.logFileVx = vxFile;
    logCfg.logFileVy = vyFile;
    logCfg.logFileAx = axFile;
    logCfg.logFileAy = ayFile;
    
    // Read initial state
    int numBirds = readNumBirds(initFile);
    int n = numBirds;
    // double x[n];
    // double y[n];
    // double vx[n];
    // double vy[n];
    // readConf(initFile, x, y, vx, vy, n);
    // Initial state
    State is;
    is.numBirds = n;
    initState(&is, n);
    readConfState(initFile, &is);
    // is.x = x;
    // is.y = y;
    // is.vx = vx;
    // is.vy = vy;
    
    // Run global MPC controller
    reynolds_control(&is, sim->steps, sim->reynoldsParams, sim->constraints, &logCfg, sim->noiseCfg);
    double timeEnd = clock();
    sim->runTimeSeconds = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;
    // freeState(&is);
    return NULL;
}

// void printUsage(char *programName)
// {
//     char *usage = "%s -b[-sc] init_file_template output_file_template\n\n"
//                   "    -b numSims    Number of simulations in batch mode. Required.\n"
//                   "    -s startNum   Starting simulation number. Optional. Default is 1.\n"
//                   "    -c conf_file  Specify config file to read parameters.\n"
//                   "                  Default config file is 'reynolds.conf'.\n"
//                   "    init_file_template    Template to look for initial configuration files. Required.\n"
//                   "    output_file_template  Template to create file names to write results. Required.\n"
//                   "\n"
//                   "Example: run 3 simulations\n"
//                   "    %s -b 3 -s 6 initstate_%%d.txt output_%%d_%%s.txt\n"
//                   "\nThe program will read initial states for 3 simulations from 'initstate_6.txt', "
//                   "'initstate_7.txt', and 'initstate_8.txt'.\n"
//                   "It will write:\n"
//                   "    evolution of x to 'output_6_x.txt', 'output_7_x.txt', and 'output_8_x.txt'\n"
//                   "    evolution of y to 'output_6_y.txt', 'output_7_y.txt', and 'output_8_y.txt'\n"
//                   "    evolution of vx to 'output_6_vx.txt', 'output_7_vx.txt', and 'output_8_vx.txt'\n"
//                   "    evolution of vy to 'output_6_vy.txt', 'output_7_vy.txt', and 'output_8_vy.txt'\n"
//                   "    evolution of ax to 'output_6_ax.txt', 'output_7_ax.txt', 'output_8_ax.txt'\n"
//                   "    evolution of ay to 'output_6_ay.txt', 'output_7_ay.txt', and 'output_8_ay.txt'\n";
//     printf(usage, programName, programName);
// }

// int main(int argc, char *argv[])
// {
//     double timeBegin = clock();
//     int numSims = 0;
//     int startNum = 1;
//     char *confFile = "C Code/reynolds.conf";
//     char *initFileTemplate = NULL;
//     char *resultFileTemplate = NULL;
//     bool quietMode = false;
    
//     int opt;
//     while ((opt = getopt(argc, argv, "b:s:c:qh")) != -1) {
//         switch (opt) {
//         case 'b': numSims = atoi(optarg); break;
//         case 's': startNum = atoi(optarg); break;
//         case 'c': confFile = optarg; break;
//         case 'q': quietMode = true; break;
//         case 'h': printUsage(argv[0]); exit(EXIT_SUCCESS);
//         default:
//             printUsage(argv[0]);
//             exit(EXIT_FAILURE);
//         }
//     }
    
//     if (optind < argc - 1) {
//         initFileTemplate = argv[optind];
//         resultFileTemplate = argv[optind + 1];
//     }
//     else {
//         printUsage(argv[0]);
//         exit(EXIT_FAILURE);
//     }
    
//     if (numSims < 1) {
//         printf("ERROR: Number of simulations (-b) is required and must be greater than zero.\n");
//         exit(EXIT_FAILURE);
//     }
    
//     // Parse Reynold's parameters from config file
//     Constraint c;
//     ReynoldsParameters param;
//     NoiseConfig noiseCfg;
//     int steps;
    
//     readReynoldsConfig(confFile, &c, &param, &steps);
    
//     // Manually set noise config
//     noiseCfg.positionSensorNoise.enabled = true;
//     noiseCfg.positionSensorNoise.uniform = false;
//     noiseCfg.positionSensorNoise.independent = true;
//     noiseCfg.positionSensorNoise.mu = 0;
//     noiseCfg.positionSensorNoise.sigma = 0.2;
//     noiseCfg.positionSensorNoise.lower = -2;
//     noiseCfg.positionSensorNoise.upper = 2;
    
//     noiseCfg.velocitySensorNoise.enabled = true;
//     noiseCfg.velocitySensorNoise.uniform = false;
//     noiseCfg.velocitySensorNoise.independent = true;
//     noiseCfg.velocitySensorNoise.mu = 0;
//     noiseCfg.velocitySensorNoise.sigma = 0.1;
//     noiseCfg.velocitySensorNoise.lower = -1;
//     noiseCfg.velocitySensorNoise.upper = 1;
    
//     noiseCfg.actuatorNoise.enabled = false;
//     noiseCfg.actuatorNoise.uniform = false;
//     noiseCfg.actuatorNoise.independent = true;
//     noiseCfg.actuatorNoise.mu = 0;
//     noiseCfg.actuatorNoise.sigma = 0.02;
    
//     if (!quietMode) {
//         printf("# of time steps = %d \t time step = %.3f\n", steps, c.dt);
//         printf("amax = %.3f \t vmax = %.3f\n", c.amax, c.vmax);
//         printf("rc = %.3f \t wc = %.3f\n", param.rc, param.wc);
//         printf("rs = %.3f \t ws = %.3f\n", param.rs, param.ws);
//         printf("rv = %.3f \t wv = %.3f\n", param.rv, param.wv);
//     }
    
//     int noiseExps = 1; // Number of noise experiments
//     char noiseOutFileTemplate[500];
//     for (int nIdx = 1; nIdx <= noiseExps; nIdx++) {
//         noiseCfg.positionSensorNoise.sigma = 0.2 * nIdx;
//         noiseCfg.velocitySensorNoise.sigma = 0.1 * nIdx;
//         sprintf(noiseOutFileTemplate, "experiments\\sensor_noise\\reynolds_sn_x%02d_v%02d_%%d_%%s.txt", 2 * nIdx, 1 * nIdx);
        
//         pthread_t simThreads[numSims];
//         Simulation simParams[numSims];
        
//         for (int i = 0; i < numSims; i++) {
//             int simNum = i + startNum;
//             if (!quietMode) {
//                 printf("Creating a pthread for simulation #%d...\n", simNum);
//             }
//             simParams[i].simNum = simNum;
//             simParams[i].initFileTemplate = initFileTemplate;
//             simParams[i].resultFileTemplate = resultFileTemplate;
//             simParams[i].steps = steps;
//             simParams[i].constraints = &c;
//             simParams[i].reynoldsParams = &param;
//             simParams[i].noiseCfg = &noiseCfg;
//             simParams[i].runTimeSeconds = -1;
            
//             if (pthread_create(&simThreads[i], NULL, runSimulationR, (void *)&simParams[i])) {
//                 fprintf(stderr, "Error creating thread\n");
//             }
//         }
        
//         // Wait for all threads to finish
//         for (int i = 0; i < numSims; i++) {
//             int simNum = i + startNum;
//             if (pthread_join(simThreads[i], NULL)) {
//                 fprintf(stderr, "Error joining thread for simulation #%d\n", simNum);
//             }
//         }
        
//         if (!quietMode) {
//             for (int i = 0; i < numSims; i++) {
//                printf("Simulation #%d finished in %f seconds\n", simParams[i].simNum, simParams[i].runTimeSeconds);
//             }
//         }
//     }
//     double timeEnd = clock();
//     double timeSeconds = (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC;
//     if (!quietMode) {
//         printf("Total running time = %f seconds\n", timeSeconds);
//     }
// }