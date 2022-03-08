/**
 * gmpc.c
 * 
 * Purpose: Global MPC flocking controller.
 *
 * How to compile:
 *    gcc gmpc.c conf.c common.c normal.c -O3 -pthread -o gmpc
 *
 * How to use: run
 *    ./gmpc
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
#include "common.h"
#include "conf.h"
#include "ziggurat.h"

// BAIG WAS HERE
long free_count = 0;
long malloc_count = 0;
/**
 * Defines parameters needed to run a simulation.
 */
typedef struct Simulation {
    int simNum;
    char *initFileTemplate;
    char *resultFileTemplate;
    int horizon;
    int steps;
    Constraint *constraints;
    FitnessParameter *fitnessParams;
    GradientDescentOptions *gdOptions;
    NoiseConfig *noiseCfg;
    double runTimeSeconds;
} Simulation;

double cohesion(State *s, double distances[s->numBirds][s->numBirds], FitnessParameter *param);
double separationExponential(State *s, double distances[s->numBirds][s->numBirds], FitnessParameter *param);
double separationContrained(State *s, double sqDistances[s->numBirds][s->numBirds], FitnessParameter *param, Constraint *c);
double velocityMatching(State *s, double distances[s->numBirds][s->numBirds], FitnessParameter *param);
double averageSquaredDistance(State *s, double sqDistances[s->numBirds][s->numBirds], FitnessParameter *param);
double seperationPolynomial(State *s, double sqDistances[s->numBirds][s->numBirds], FitnessParameter *param);
double averageDistance(State *s, double distances[s->numBirds][s->numBirds], FitnessParameter *param);
double maxValue(int n, double matrix[n][n]);
double obstacleCost1(State *s, Obstacles *O, FitnessParameter *param);
double obstacleCost2(State *s, Rectangle *R, FitnessParameter *param);
double pointObstacleCost(State *s, pointObstacles *P, FitnessParameter *param);
double target(State *s,  double target[2] );
double predatorAvoidance(State *s, double sqDistances[s->numBirds][s->numBirds], int numPreds, FitnessParameter *param, Constraint *c);
double caBarrierFunction_ij(double * p1, double * p2, double * v1, double * v2, Constraint *c);
double CompositeBarrierFunction(State *s, Constraint *c);
double barrierFunctionConsraint(State *s, Action *a, Constraint *c);
double alpha(double val, Constraint * c);
double Lie_Derivative(State *s, Action *a, Constraint *c);
double controller_predator(State *s, Action *a, Constraint *c);

/**
 * Calculates the fitness value for a given state.
 * @param s The state to calculate fitness.
 * @param param The parameters for fitness calculation.
 * @return The fitness value for the specified state.
 */
double fitness(State *s, Action *a, FitnessParameter *param, Constraint *c)
{
    double sqDistances[s->numBirds][s->numBirds];
    computeSquaredDistances(s, sqDistances);

    // double sp = seperationPolynomial(s, sqDistances, param);
    double spC = separationContrained(s, sqDistances, param, c);
    double asd = averageSquaredDistance(s, sqDistances, param);
    // double obst = obstacleCost1(s, c->obs, param);
    // double P2Obst = pointObstacleCost(s, c->pointObs, param);
    // double cbf = CompositeBarrierFunction(s, c);
    // double constrainedBf = barrierFunctionConsraint(s, a, c);
    // double f = predatorAvoidance(s, sqDistances, 1, param, c);

    // double destination[2] = {60, 50};
    // double tgt = target(s, destination);

    //double f = param->wd * pow(asd,2) + param->ws * pow(sp,2) + pow(tgt,2) + 800000*pow(P2Obst,2) ; 
    // double f = param->wd * pow(asd,1) + param->ws * pow(sp,1);// + 20000* pow(spC,2) + 100000*pow(obst,2) + 30*pow(tgt,1); 
    // printf("BF: %f\n", barrierFunctionConsraint(s, a, c));
    // double f = param->wd * asd + param->ws * sp;
    double f = param->wd * asd + param->wspc * spC;
    // Obstacle_avoidance_final:
    // double f = param->wd * asd + 20000 * spC + 200000 * obst + 1000 * tgt;
    // double constraint_violation = pow(spC, 2) + pow(obst, 2);
    // double f = param->wd * asd + param->wspc * sqrt(constraint_violation) + param->wt * tgt;
    // double pred = predatorAvoidance(s, sqDistances, 1, param);
    return f;    
}

/**
Calculates the combined fitness for predator - prey flock. 
**/
double predatorAvoidance(State *s, double sqDistances[s->numBirds][s->numBirds], int numPreds, FitnessParameter *param, Constraint *c){
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
    double asd = averageSquaredDistance(s, sqD, param);
    double spC = separationContrained(s, sqD, param, c);
    s->numBirds = n;


    // for (int j = 0; j < f.numBirds; j++){
    //     centreX += f.x[j];
    //     centreY += f.y[j];
    // }

    // centreX /= f.numBirds;
    // centreY /= f.numBirds;

    double pa = 0;
    // double pAttack = 0;

    for(int k = n-numPreds; k < n; k++){
        for(int i = 0; i < n-numPreds; i++){
                pa += pow(MAX(c->predAvD - sqrt(sqDistances[i][k]), 0), 2);
                // pa += MAX(c->predAvD*c->predAvD  - sqDistances[i][k], 0);
        }
    }

    // for(int k = n-numPreds; k < n; k++){
    //     pAttack += pow(s->x[k] - centreX, 2) + pow(s->y[k] - centreY, 2);
    // }
    double constraint_violations = sqrt(pow(spC, 2) + pa);
    return  param->wd * asd +  param->wspc * constraint_violations;

}
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

/*
 * alpha is a class kappa function.
 */
double alpha(double val, Constraint * c){
    double gamma = 0.25e-05;
    // return val/20;
    // return -val / (val - c->alpha);
    return gamma * pow(val, 2);

}

/**
 * Calculates the composite barrier function constraint (as in eqn 21 of Wang et. al.)
 * @param s The current state.
 * @param a The sequence of actions.
 * @param c The constraints on the system. 
 * @return The value for the barrier function. 
 */
double barrierFunctionConsraint(State *s, Action *a, Constraint *c){
    double LD = Lie_Derivative(s, a, c);
    // if (LD < 0){
    //     printf("LD:%f\n", LD);
    // }
    // printf("LD:%f\nbf: %f\n", Lie_Derivative, CompositeBarrierFunction(s, c));
    double alphaVal = alpha( CompositeBarrierFunction(s, c) , c);
    // alphaVal = 0.5;
    // printf("Alpha%f\n", alphaVal);
    // return MAX(-(Lie_Derivative + alphaVal), 0);
    return MAX((-LD-alphaVal), 0);
    // return ;
}

/**
 * Calculates the Lie derivative of composite barrier function constraint (as in eqn 21 of Wang et. al.)
 * @param s The current state.
 * @param a The sequence of actions.
 * @param c The constraints on the system. 
 * @return The value for the barrier function. 
 */
double Lie_Derivative(State *s, Action *a, Constraint *c){
    Direction d;
    d.numBirds = s->numBirds;
    int n = s->numBirds;
    double delta = 0.0000000001;

    d.d_vx = (double *)malloc(n*sizeof(double));
    d.d_vy = (double *)malloc(n*sizeof(double));
    d.d_ux = (double *)malloc(n*sizeof(double));
    d.d_uy = (double *)malloc(n*sizeof(double));
    // BAIG WAS HERE
    malloc_count += 4;
    for (int i = 0; i < n; i++){
         double next_vx = s->vx[i]  + c->dt * a->ax[i];
         double next_vy = s->vy[i]  + c->dt * a->ay[i];
         d.d_vx[i] = c->dt * next_vx;
         d.d_vy[i] = c->dt * next_vy;
         d.d_ux[i] = c->dt * a->ax[i];
         d.d_uy[i] = c->dt * a->ay[i];
    }
    State cs = cloneState(s);
    for (int i = 0; i < n; i++){
        cs.x[i] -= delta * d.d_vx[i];
        cs.y[i] -= delta * d.d_vy[i];
        cs.vx[i] -= delta * d.d_ux[i];
        cs.vy[i] -= delta * d.d_uy[i];
    }
    double Lie_Derivative = (CompositeBarrierFunction(&cs, c) - CompositeBarrierFunction(s, c)) / delta;

    freeState(&cs);
    freeDirection(&d);
    return -Lie_Derivative;
}


/**
 * Calculates the composite barrier function (as in Wang et. al.)
 * @param s The current state.
 * @param c The constraints on the system. 
 * @return The value for the barrier function. 
 */
double CompositeBarrierFunction(State *s, Constraint *c){
    int n = s->numBirds;
    double pi[2];
    double pj[2];
    double vi[2];
    double vj[2];
    double ret = 1;
    int count = 0;

    for (int i = 0; i < n; i++){
        for (int j = i+1; j < n; j++){
            count++;
            pi[0] = s->x[i];
            pi[1] = s->y[i];
            pj[0] = s->x[j];
            pj[1] = s->y[j];            

            vi[0] = s->vx[i];
            vi[1] = s->vy[i];
            vj[0] = s->vx[j];
            vj[1] = s->vy[j];
            // printf("pw_hij: %f\n", caBarrierFunction_ij(pi, pj, vi, vj, c));
            ret *= caBarrierFunction_ij(pi, pj, vi, vj, c);
        }
    }
    return ret;
}

/**
 * Calculates the Collision Avoidance Barrier function for a pair of agents. (as in Wang et. al.)
 * @param p1, p2 are the positions.
 * @param v1, v2 are the velocities. 
 * @param dmin 
 * @param maxAcc 
 * @return The value for the barrier function. 
 */
double caBarrierFunction_ij(double * pi, double * pj, double * vi, double * vj, Constraint *c){
    double p_ij[2]; 
    double v_ij[2]; 
    double maxAcc = c->amax;
    double dmin = c->dmin;

    p_ij[0] = pi[0] - pj[0];
    p_ij[1] = pi[1] - pj[1];
    double mag = norm(p_ij[0], p_ij[1]);  
    if (mag > 4*dmin){
        return 1;
    }

    v_ij[0] = vi[0] - vj[0];
    v_ij[1] = vi[1] - vj[1];

    double ret = dot(p_ij[0]/mag, p_ij[1]/mag, v_ij[0], v_ij[1]) + sqrt( 4 * maxAcc * c->dt * ( norm(p_ij[0], p_ij[1]) - dmin) );
    return MAX(ret, 0);
}

/**
 * Calculates the cost of applying a sequence of control actions to a given state.
 * @param s The state.
 * @param a The array of actions to apply in sequence.
 * @param h The horizon --- number of control actions in a to apply.
 * @param c The system constraints.
 * @param param The parameters for fitness calculation.
 * @return The cost of applying a to s. It is the fitness value for the final state in the horizon.
 */
double costFinalState(State *s, Action *a, int h, Constraint *c, FitnessParameter *param)
{
    State currentS;
    State nextS = dynamics(s, &a[0], c);
    for (int i = 1; i < h; i++) {
        if (i > 1) {
            freeState(&currentS);
        }
        currentS = nextS;
        nextS = dynamics(&currentS, &a[i], c);
    }
    
    double cost = fitness(&nextS, &a[h], param, c); //Probably smth wrong.
    freeState(&nextS);
    if (h > 1) {
        freeState(&currentS);
    }
    return cost;
}

/**
 * Calculates the cost of applying a sequence of control actions to a given state.
 * @param s The state.
 * @param a The array of actions to apply in sequence.
 * @param h The horizon --- number of control actions in a to apply.
 * @param c The system constraints.
 * @param param The parameters for fitness calculation.
 * @return The cost of applying a to s. It is the minimum fitness value achieved during the horizon.
 */
// double costMinimum(State *s, Action *a, int h, Constraint *c, FitnessParameter *param)
// {
//     State currentS;
//     State nextS = dynamics(s, &a[0], c);
//     double cost = fitness(&nextS, &a[0], param, c);
//     double f;
//     for (int i = 1; i < h; i++) {
//         if (i > 1) {
//             freeState(&currentS);
//         }
//         currentS = nextS;
//         nextS = dynamics(&currentS, &a[i], c);
//         f = fitness(&nextS, &a[h], param, c);
//         if (f < cost) {
//             cost = f;
//         }
//     }
//     freeState(&nextS);
//     if (h > 1) {
//         freeState(&currentS);
//     }
//     return cost;
// }

/**
 * Calculates the cost of applying a sequence of control actions to a given state.
 * @param s The state.
 * @param a The array of actions to apply in sequence.
 * @param h The horizon --- number of control actions in a to apply.
 * @param c The system constraints.
 * @param param The parameters for fitness calculation.
 * @return The cost of applying a to s. It is the sum of fitness values at states during the horizon 
 *         and the cost of the control action sequence.
 */
// double costSum(State *s, Action *a, int h, Constraint *c, FitnessParameter *param)
// {
//     State currentS;
//     State nextS = dynamics(s, &a[0], c);
//     double cost = fitness(&nextS, param) + param->wu * controlCost(&a[0]);
//     double f;
//     for (int i = 1; i < h; i++) {
//         if (i > 1) {
//             freeState(&currentS);
//         }
//         currentS = nextS;
//         nextS = dynamics(&currentS, &a[i], c);
//         cost += fitness(&nextS, param) + param->wu * controlCost(&a[i]);
//     }
//     freeState(&nextS);
//     if (h > 1) {
//         freeState(&currentS);
//     }
//     return cost;
// }


/**
 * Calculates the cost of applying a sequence of control actions to a given state.
 * @param s The state.
 * @param a The array of actions to apply in sequence.
 * @param h The horizon --- number of control actions in a to apply.
 * @param c The system constraints.
 * @param param The parameters for fitness calculation.
 * @return The cost of applying a to s. It is the sum of fitness values at states during the horizon 
 *         and the cost of the control action sequence.
 */
double costSum(State *s, Action *a, int h, Constraint *c, FitnessParameter *param){
    State cs = cloneState(s);
    int n = s->numBirds;
    double sqDistances[n][n];
    computeSquaredDistances(s, sqDistances);
    double cost = 0;
    double actionCost = 0;
    int controlSteps = c->ct/c->dt;
    for(int j = 0; j<h; j++){
        for(int i = 0; i<n; i++){
            for (int k = 0; k<controlSteps; k++){
                if (c->predator && i == n - 1){
                    trim(&cs.vx[i], &cs.vy[i], c->pFactor * c->vmax);
                    cs.x[i] += c->dt * cs.vx[i];
                    cs.y[i] += c->dt * cs.vy[i]; 
                    continue;                  
                }
                cs.vx[i] += c->dt * a[j].ax[i];
                cs.vy[i] += c->dt * a[j].ay[i];
                trim(&cs.vx[i], &cs.vy[i], c->vmax);
                cs.x[i] += c->dt * cs.vx[i];
                cs.y[i] += c->dt * cs.vy[i];
            }
            actionCost += a[j].ax[i]*a[j].ax[i] + a[j].ay[i]*a[j].ay[i];
        }
        cost += fitness(&cs, &a[j], param, c);

    }
    cost += fitness(&cs, &a[h-1], param, c);
    freeState(&cs); 
    return (double) cost/h + actionCost/h;
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
bool optimize_adapt(State *s, Action *a, int h, Constraint *c, FitnessParameter *param, GradientDescentOptions *options)
{
    int n =  s->numBirds;
    double learningRate = options->learningRate;
    int count = 0;
    // Initialize a
    // for (int i = 0; i < h; i++) {
    //     memset(a[i].ax, 0, n*sizeof(double));
    //     memset(a[i].ay, 0, n*sizeof(double));
    // }
    
    double dfax[h][n];  // Store derivative of cost function w.r.t ax
    double dfay[h][n];  // Store derivative of cost function w.r.t ay

    double currCost;
    double newCost;
    char name[500] = "fitOverHorizon.txt";
    char lambdaFile[500] = "lambda.txt";

    Action aCopy[h];
    for (int i = 0; i < h; i++) {
        aCopy[i].numBirds = n;
        aCopy[i].ax = (double *)malloc(n*sizeof(double));
        aCopy[i].ay = (double *)malloc(n*sizeof(double));

    // BAIG WAS HERE
    malloc_count += 2;
    }

    appendToFile(name, costSum(s, a, h, c, param));
    appendToFile(lambdaFile, learningRate);
    for (int i = 0; i < options->maxIters; i++) {
        count++;
        if (count > 200){
            break;
        }
        // Calculate the cost for the current accelerations
        currCost = costSum(s, a, h, c, param);

        double gradMagnitude = 0;
        for (int j = 0; j < h; j++) {
            for (int k = 0; k < n; k++) {
                // Calculate derivative of cost function w.r.t ax and ay of each bird
                a[j].ax[k] += options->delta;
                newCost = costSum(s, a, h, c, param);
                a[j].ax[k] -= options->delta;
                dfax[j][k] = (newCost - currCost) / options->delta;
                
                a[j].ay[k] += options->delta;
                newCost = costSum(s, a, h, c, param);
                a[j].ay[k] -= options->delta;
                dfay[j][k] = (newCost - currCost) / options->delta;
                
                gradMagnitude += fabs(dfax[j][k]) + fabs(dfay[j][k]);
            }
        }
        // Update accelerations
        for (int j = 0; j < h; j++) {
            // double temp = 0;
            for (int k = 0; k < n; k++) {
                aCopy[j].ax[k] = a[j].ax[k];
                aCopy[j].ay[k] = a[j].ay[k];

                a[j].ax[k] -= learningRate * dfax[j][k];
                a[j].ay[k] -= learningRate * dfay[j][k];
            }
        } 
        // trimActionSeq(a, h, c->amax); 
        trimActionSeqNaive(a, h, c->amax);

        double newCost =  costSum(s, a, h, c, param);
        if (newCost - currCost <= 0.00001){
            learningRate = 1.05 * learningRate;
            appendToFile(name, newCost);
            appendToFile(lambdaFile, learningRate);


        }
        else{
            i -= 1;
            //Undo Change:
            for (int j = 0; j < h; j++) {
                for (int k = 0; k < n; k++) {
                    a[j].ax[k] = aCopy[j].ax[k];
                    a[j].ay[k] = aCopy[j].ay[k];
                }
            }
            //Decrease learning-rate
            learningRate = 0.50 * learningRate;
            // appendToFile(lambdaFile, learningRate);
        }

    }
    for (int i = 0; i < h; i++){
        freeAction(&aCopy[i]);
    }
    return true;
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
bool optimize(State *s, Action *a, int h, Constraint *c, FitnessParameter *param, GradientDescentOptions *options)
{
    int n =  s->numBirds;
    // Initialize a
    for (int i = 0; i < h; i++) {
        memset(a[i].ax, 0, n*sizeof(double));
        memset(a[i].ay, 0, n*sizeof(double));
    }
    
    double dfax[h][n];  // Store derivative of cost function w.r.t ax
    double dfay[h][n];  // Store derivative of cost function w.r.t ay

    double currCost;
    double newCost;
    char name[500] = "fitOverHorizon.txt";
    
    for (int i = 0; i < options->maxIters; i++) {
        // Calculate the cost for the current accelerations
        currCost = costSum(s, a, h, c, param);
        appendToFile(name, currCost);

        double gradMagnitude = 0;
        for (int j = 0; j < h; j++) {
            for (int k = 0; k < n; k++) {
                // Calculate derivative of cost function w.r.t ax and ay of each bird
                a[j].ax[k] += options->delta;
                newCost = costSum(s, a, h, c, param);
                a[j].ax[k] -= options->delta;
                dfax[j][k] = (newCost - currCost) / options->delta;
                
                a[j].ay[k] += options->delta;
                newCost = costSum(s, a, h, c, param);
                a[j].ay[k] -= options->delta;
                dfay[j][k] = (newCost - currCost) / options->delta;
                
                gradMagnitude += pow(dfax[j][k],2) + pow(dfay[j][k],2);
            }
        }
        gradMagnitude = sqrt(gradMagnitude);
        // printf("k: %d, Grad: %f\n", i, gradMagnitude);
        // Update accelerations
        for (int j = 0; j < h; j++) {
            // double temp = 0;
            for (int k = 0; k < n; k++) {
                a[j].ax[k] -= options->learningRate * dfax[j][k];
                a[j].ay[k] -= options->learningRate * dfay[j][k];
                if (c->predator && i == n - 1){
                    trim(&a[j].ax[k], &a[j].ay[k], c->pFactor * c->amax);
                }
                else{
                    trim(&a[j].ax[k], &a[j].ay[k], c->amax);
                }
            }
        } 
        // trimActionSeq(a, h, c->amax); 
        // trimActionSeqNaive(a, h, c->amax);
    }
    appendToFile(name, costSum(s, a, h, c, param));
    return true;
}

bool optimize_momentum(State *s, Action *a, int h, Constraint *c, FitnessParameter *param, GradientDescentOptions *options)
{
    int n =  s->numBirds;
    // Initialize a
    for (int i = 0; i < h; i++) {
        memset(a[i].ax, 0, n*sizeof(double));
        memset(a[i].ay, 0, n*sizeof(double));
    }
    
    double dfax[h][n];  // Store derivative of cost function w.r.t ax
    double dfay[h][n];  // Store derivative of cost function w.r.t ay

    double veeteeX[h][n];
    double veeteeY[h][n];
    for (int j = 0; j < h; j++) {
        for (int k = 0; k < n; k++) {
            veeteeX[j][k] = 0;
            veeteeY[j][k] = 0;
        }
    }

    double currCost;
    double newCost;
    double gradMagnitude = 0;
    double gamma = 0.25;
    
    
    for (int i = 0; i < options->maxIters; i++) {
        // Calculate the cost for the current accelerations
        currCost = costSum(s, a, h, c, param);
        gradMagnitude = 0;
        for (int j = 0; j < h; j++) {
            for (int k = 0; k < n; k++) {
                // Calculate derivative of cost function w.r.t ax and ay of each bird
                a[j].ax[k] += options->delta;
                newCost = costSum(s, a, h, c, param);
                a[j].ax[k] -= options->delta;
                dfax[j][k] = (newCost - currCost) / options->delta;
                
                a[j].ay[k] += options->delta;
                newCost = costSum(s, a, h, c, param);
                a[j].ay[k] -= options->delta;
                dfay[j][k] = (newCost - currCost) / options->delta;
                
                gradMagnitude += fabs(dfax[j][k]) + fabs(dfay[j][k]);
            }
        }
        
        double * maxAgentAcc;   
        maxAgentAcc = (double *)malloc(h*sizeof(double));
        //BAIG WAS HERE
        malloc_count += 1;
        memset(maxAgentAcc, 0, h*sizeof(double));

        
        // Update accelerations
        for (int j = 0; j < h; j++) {
            double temp = 0;
            for (int k = 0; k < n; k++) {
                veeteeX[j][k] = gamma * veeteeX[j][k] + options->learningRate * dfax[j][k];
                a[j].ax[k] -= veeteeX[j][k];
                veeteeY[j][k] = gamma * veeteeY[j][k] + options->learningRate * dfay[j][k];
                a[j].ay[k] -= veeteeY[j][k];

                temp = norm(a[j].ax[k], a[j].ay[k]);
                if (maxAgentAcc[j] > temp){
                    maxAgentAcc[j] = temp;
                }
                // Ensure acceleration doesn't exceed amax
                // trim(&a[j].ax[k], &a[j].ay[k], c->amax);
                // double normA = norm(a[j].ax[k], a[j].ay[k]);
                // if (normA > c->amax) {
                //     a[j].ax[k] *= c->amax / normA;
                //     a[j].ay[k] *= c->amax / normA;
                // }
            }
        } 
        // Alternate, scale accelerations with the maximum acc for any agent at any time step.
        for (int j = 0; j < h; j++) {
            if (maxAgentAcc[j] > c->amax){
                for (int k = 0; k < n; k++) { 
                    a[j].ax[k] *= c->amax / maxAgentAcc[j];
                    a[j].ay[k] *= c->amax / maxAgentAcc[j];   
                }
            }
        }     
    }
    // for (int j = 0; j < h; j++) {
    //         for (int k = 0; k < n; k++) {
    //             a[j].ax[k] /= c->ct/c->dt;
    //             a[j].ay[k] /= c->ct/c->dt;
    //     }
    // }

    return true;
}

bool optimize_Nesterov(State *s, Action *a, int h, Constraint *c, FitnessParameter *param, GradientDescentOptions *options)
{
    int n =  s->numBirds;
    // Initialize a
    for (int i = 0; i < h; i++) {
        memset(a[i].ax, 0, n*sizeof(double));
        memset(a[i].ay, 0, n*sizeof(double));
    }

    Action a_f[h];
    for (int i = 0; i < h; i++) {
        a_f[i].numBirds = n;
        a_f[i].ax = (double *)malloc(n*sizeof(double));
        a_f[i].ay = (double *)malloc(n*sizeof(double));

        memset(a_f[i].ax, 0, n*sizeof(double));
        memset(a_f[i].ay, 0, n*sizeof(double));

    // BAIG WAS HERE
    malloc_count += 2;
    }


    
    double dfax[h][n];  // Store derivative of cost function w.r.t ax
    double dfay[h][n];  // Store derivative of cost function w.r.t ay

    double veeteeX[h][n];
    double veeteeY[h][n];
    for (int j = 0; j < h; j++) {
        for (int k = 0; k < n; k++) {
            veeteeX[j][k] = 0;
            veeteeY[j][k] = 0;
        }
    }

    double currCost;
    double newCost;
    double gradMagnitude = 0;
    double gamma = 0.90;
    char name[500] = "fitOverHorizon.txt";
    
    for (int i = 0; i < options->maxIters; i++) {
        //set a-f
        for (int j = 0; j < h; j++) {
            for (int k = 0; k < n; k++) {
                a_f[j].ax[k] = a[j].ax[k] - gamma * veeteeX[j][k];
                a_f[j].ay[k] = a[j].ay[k] - gamma * veeteeY[j][k];
            }
        }
        // Calculate the cost for the current accelerations
        currCost = costSum(s, a_f, h, c, param);
        appendToFile(name, currCost);

        gradMagnitude = 0;
        for (int j = 0; j < h; j++) {
            for (int k = 0; k < n; k++) {
                // Calculate derivative of cost function w.r.t ax and ay of each bird
                a_f[j].ax[k] += options->delta;
                newCost = costSum(s, a_f, h, c, param);
                a_f[j].ax[k] -= options->delta;
                dfax[j][k] = (newCost - currCost) / options->delta;
                
                a_f[j].ay[k] += options->delta;
                newCost = costSum(s, a_f, h, c, param);
                a_f[j].ay[k] -= options->delta;
                dfay[j][k] = (newCost - currCost) / options->delta;
                
                gradMagnitude += fabs(dfax[j][k]) + fabs(dfay[j][k]);
            }
        }
        

        
        // Update accelerations
        for (int j = 0; j < h; j++) {
            for (int k = 0; k < n; k++) {
                veeteeX[j][k] = gamma * veeteeX[j][k] + options->learningRate * dfax[j][k];
                a[j].ax[k] -= veeteeX[j][k];
                veeteeY[j][k] = gamma * veeteeY[j][k] + options->learningRate * dfay[j][k];
                a[j].ay[k] -= veeteeY[j][k];
                // Ensure acceleration doesn't exceed amax
                // trim(&a[j].ax[k], &a[j].ay[k], c->amax);
                // double normA = norm(a[j].ax[k], a[j].ay[k]);
                // if (normA > c->amax) {
                //     a[j].ax[k] *= c->amax / normA;
                //     a[j].ay[k] *= c->amax / normA;
                // }
            }
        } 
        trimActionSeq(a, h, c->amax); 
    }
    appendToFile(name, costSum(s, a, h, c, param));

    //Free a_future
    for(int i = 0; i < h; i++){
        freeAction(&a_f[i]);
    }

    return true;
}

bool optimize_adadelta(State *s, Action *a, int h, Constraint *c, FitnessParameter *param, GradientDescentOptions *options)
{
    int n =  s->numBirds;
    // Initialize a
    for (int i = 0; i < h; i++) {
        memset(a[i].ax, 0, n*sizeof(double));
        memset(a[i].ay, 0, n*sizeof(double));
    }
    
    double dfax[h][n];  // Store derivative of cost function w.r.t ax
    double dfay[h][n];  // Store derivative of cost function w.r.t ay

    double avgSqGradX[h][n]; 
    double avgSqGradY[h][n]; 
    //Zero initialize avgSqGrad
    memset( avgSqGradX, 0, h*n*sizeof(double) );
    memset( avgSqGradY, 0, h*n*sizeof(double) );

    double avgSqStepX[h][n];
    double avgSqStepY[h][n];

    //Zero initialize avgSqGrad
    memset( avgSqStepX, 0, h*n*sizeof(double) );
    memset( avgSqStepY, 0, h*n*sizeof(double) );

    double currCost;
    double newCost;
    double gradMagnitude = 0;
    double gamma = 0.90;
    double epsilon = 0.00000001;
    char name[500] = "fitOverHorizon.txt";

    
    
    for (int i = 0; i < options->maxIters; i++) {
        // Calculate the cost for the current accelerations
        currCost = costSum(s, a, h, c, param);
        appendToFile(name, currCost);

        gradMagnitude = 0;
        for (int j = 0; j < h; j++) {
            for (int k = 0; k < n; k++) {
                // Calculate derivative of cost function w.r.t ax and ay of each bird
                a[j].ax[k] += options->delta;
                newCost = costSum(s, a, h, c, param);
                a[j].ax[k] -= options->delta;
                dfax[j][k] = (newCost - currCost) / options->delta;
                
                a[j].ay[k] += options->delta;
                newCost = costSum(s, a, h, c, param);
                a[j].ay[k] -= options->delta;
                dfay[j][k] = (newCost - currCost) / options->delta;

                avgSqGradX[j][k] = gamma * avgSqGradX[j][k] + (1 - gamma) * pow(dfax[j][k], 2);
                avgSqGradY[j][k] = gamma * avgSqGradY[j][k] + (1 - gamma) * pow(dfay[j][k], 2);
                
                gradMagnitude += fabs(dfax[j][k]) + fabs(dfay[j][k]);
            }
        }
        


        // Update accelerations
        for (int j = 0; j < h; j++) {
            for (int k = 0; k < n; k++) {
                double stepX =  - 2*(sqrt(avgSqStepX[j][k] + epsilon) / sqrt(avgSqGradX[j][k] + epsilon)) * dfax[j][k];
                double stepY =  - 2*(sqrt(avgSqStepY[j][k] + epsilon) / sqrt(avgSqGradY[j][k] + epsilon)) * dfay[j][k];
                avgSqStepX[j][k] = gamma * avgSqStepX[j][k] + (1 - gamma) * pow(stepX,2);
                avgSqStepY[j][k] = gamma * avgSqStepY[j][k] + (1 - gamma) * pow(stepY,2);

                a[j].ax[k] += stepX;
                a[j].ay[k] += stepY;

                // if (j == 1 && k == 1){
                //     printf("step(%d):%.5f, acc:%f\n", i, sqrt( avgSqStepX[1][1])/sqrt(avgSqGradX[1][1]), a[1].ax[1] );
                // }
                // Ensure acceleration doesn't exceed amax
                // trim(&a[j].ax[k], &a[j].ay[k], c->amax);
                // double normA = norm(a[j].ax[k], a[j].ay[k]);
                // if (normA > c->amax) {
                //     a[j].ax[k] *= c->amax / normA;
                //     a[j].ay[k] *= c->amax / normA;
                // }
            }
        } 
        trimActionSeq(a, h, c->amax);
    }
    appendToFile(name, costSum(s, a, h, c, param));

    // for (int j = 0; j < h; j++) {
    //         for (int k = 0; k < n; k++) {
    //             a[j].ax[k] /= c->ct/c->dt;
    //             a[j].ay[k] /= c->ct/c->dt;
    //     }
    // }

    return true;
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
void globalControl(State *s, int h, int steps, Constraint *c, FitnessParameter *param, GradientDescentOptions *options, LogConfig *logCfg, NoiseConfig *noiseCfg)
{
    /// Random generator initialization
    uint32_t seed = initRandomSeed();
    uint32_t kn[128];
    float fn[128];
    float wn[128];
    r4_nor_setup(kn, fn, wn);
    
    State states[steps];     // Stores the evolution of state.
    Action actions[steps];   // Stores the evolution of acceleration.
    double fitnesses[steps]; // Stores the evolution of fitness value.
    
    int n = s->numBirds;
    states[0] = *s; // initial state
    // There are no accelerations before the initial state
    actions[0].numBirds = n;
    actions[0].ax = (double *)malloc(n*sizeof(double));
    actions[0].ay = (double *)malloc(n*sizeof(double));

    // BAIG WAS HERE
    malloc_count += 2;

    memset(actions[0].ax, 0, n*sizeof(double));
    memset(actions[0].ay, 0, n*sizeof(double));
    
    // Initialize a sequence of accelerations for optimizer to update at each time step
    Action a[h];
    for (int i = 0; i < h; i++) {
        a[i].numBirds = n;
        a[i].ax = (double *)malloc(n*sizeof(double));
        a[i].ay = (double *)malloc(n*sizeof(double));

	// BAIG WAS HERE
	malloc_count += 2;
    }

    double LD[steps];
    // Fitness at initial state
    fitnesses[0] = fitness(s, &a[0], param, c);
    int controlStep = c->ct/c->dt;
    LD[0] = 0;
    clock_t timeStart = 0;
    clock_t timeEnd = 0;
    double total_time = 0;
    // printf("%d\n", controlStep);
    for (int t = 1; t < steps; t++) {
        if ( t % 100 == 0){
            printf("Time: %d\n", t);
        }
        State *s = &states[t-1];  // Current state

        State s_sensed = cloneState(s);
        
        // addNoiseToState(&s_sensed, noiseCfg, &seed, kn, fn, wn); // Add sensing noise
        // if ( (t - 1 ) % controlStep == 0){
        //     // printf("%d\n", t);
        //     optimize(&s_sensed, a, h, c, param, options); // Find the best sequence of acceleration
        //     states[t] = dynamics(s, &a[0], c);
        //     actions[t] = a[0];   
        //     // Since we keep a[0], we should allocate new memory for the next time step
        //     a[0].ax = (double *)malloc(n*sizeof(double));
        //     a[0].ay = (double *)malloc(n*sizeof(double));   

        //     // BAIG WAS HERE
        //     malloc_count += 2;      
        // }
        if ( (t - 1 ) % controlStep == 0){
            // printf("%d\n", t);
            timeStart = clock();
            optimize(&s_sensed, a, h, c, param, options); // Find the best sequence of acceleration
            timeEnd = clock();
            total_time += (double) (timeEnd - timeStart) / CLOCKS_PER_SEC;
            if (c->predator){
                double dist = 0;
                dist = controller_predator(s, a, c);
                if (dist < 10){
                    a[0].ax[n-1] = actions[t-1].ax[n-1];
                    a[0].ay[n-1] = actions[t-1].ay[n-1];
                }
            }
            states[t] = dynamics(s, &a[0], c);
            actions[t] = a[0];   
            // Since we keep a[0], we should allocate new memory for the next time step
            a[0].ax = (double *)malloc(n*sizeof(double));
            a[0].ay = (double *)malloc(n*sizeof(double));

            for (int l = 0; l < h-1; l++){
                memcpy(a[l].ax, a[l+1].ax, n*sizeof(double));
                memcpy(a[l].ay, a[l+1].ay, n*sizeof(double));
            }   
            memset(a[h-1].ax, 0, sizeof(double));
            memset(a[h-1].ay, 0, sizeof(double));
            // BAIG WAS HERE
            malloc_count += 2;      
        }
        else{
            states[t] = dynamics(s, &actions[t-1], c);
            actions[t].ax = (double *)malloc(n*sizeof(double));
            actions[t].ay = (double *)malloc(n*sizeof(double));
            memcpy(actions[t].ax, actions[t-1].ax, n*sizeof(double));
            memcpy(actions[t].ay, actions[t-1].ay, n*sizeof(double));

            //BAIG WAS HERE
            malloc_count += 2; 
            // actions[t] = actions[t-1];            
        }
        // Store and apply the first acceleration

        // Calculate and store fitness value
        fitnesses[t] = fitness(&states[t-1], &actions[t], param, c);
        LD[t] = Lie_Derivative(&states[t-1], &actions[t], c);

        freeState(&s_sensed);
    }
    
    logStates(states, steps, logCfg);
    logActions(actions, steps, logCfg);
    // logFitnesses(fitnesses, steps, logCfg);
    // logLD(LD, steps, logCfg);

    
    // Free memory
    for (int t = 0; t < steps; t++) {
        if (t > 0) { // Do not free the initial state as it was created outside of this function
            freeState(&states[t]);
        }
        // printf("Action: %d\n", t);
        freeAction(&actions[t]);
    }

    for (int i = 1; i < h; i++) {
        // printf("a: %d\n", i);
        freeAction(&a[i]);
    }
    printf("avg controller time: %f\n", (double)(total_time / steps) * (c->ct / c->dt)  );
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
    char xFile[500];
    char yFile[500];
    char vxFile[500];
    char vyFile[500];
    char axFile[500];
    char ayFile[500];
    char policy[500];
    char fitnessFile[500];
    char lieDerivativeFile[500];
    
    int simNum = sim->simNum;

    sprintf(initFile, sim->initFileTemplate, simNum);
    sprintf(xFile, sim->resultFileTemplate, simNum, "x");
    sprintf(yFile, sim->resultFileTemplate, simNum, "y");
    sprintf(vxFile, sim->resultFileTemplate, simNum, "vx");
    sprintf(vyFile, sim->resultFileTemplate, simNum, "vy");
    sprintf(axFile, sim->resultFileTemplate, simNum, "ax");
    sprintf(ayFile, sim->resultFileTemplate, simNum, "ay");
    sprintf(policy, sim->resultFileTemplate, simNum, "policy");
    sprintf(fitnessFile, sim->resultFileTemplate, simNum, "fitness");
    sprintf(lieDerivativeFile, sim->resultFileTemplate, simNum, "LD");
    
    LogConfig logCfg;
    logCfg.logFileX = xFile;
    logCfg.logFileY = yFile;
    logCfg.logFileVx = vxFile;
    logCfg.logFileVy = vyFile;
    logCfg.logFileAx = axFile;
    logCfg.logFileAy = ayFile;
    logCfg.logFileFitness = fitnessFile;
    logCfg.logFileLieDerivative = lieDerivativeFile;
    logCfg.logFilePolicy = policy;
    
    // Read initial state
    int numBirds = readNumBirds(initFile);
    int n = numBirds;

    State is;
    initState(&is, n);
    readConfState(initFile, &is);
    
    //Where is the allocation for positions and vel.
    // Run global MPC controller
    globalControl(&is, sim->horizon, sim->steps, sim->constraints, sim->fitnessParams, sim->gdOptions, &logCfg, sim->noiseCfg);

    double timeEnd = clock();
    sim->runTimeSeconds = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;
    // freeState(&is);
    return NULL;
}

void genTraningData(int h, Constraint *c, FitnessParameter *param, GradientDescentOptions *options){
    // Load data from file to array 
    int n = 2;
    char *infilename = "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/Training_data/halfSamples_2_2618276.txt";
    char *outfilename = "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/Training_data/fullSamples_2_2618276.txt";
    FILE *fp = fopen(infilename, "r");
    FILE *fpO = fopen(outfilename, "w+");
    float t1,t2,t3,t4,t5,t6,t7,t8;
    State s;
    initState(&s, n);
    s.numBirds = n;
    int i = 0;

    // Initialize a sequence of accelerations for optimizer to update at each time step
    Action a[h];
    for (int i = 0; i < h; i++) {
        a[i].numBirds = n;
        a[i].ax = (double *)malloc(n*sizeof(double));
        a[i].ay = (double *)malloc(n*sizeof(double));

	// BAIG WAS HERE
	malloc_count += 2;
    }
    fprintf(fpO, "px1,px2,py1,py2,qx1,qx2,qy1,qy2,ax1,ax2,ay1,ay2\n");
    while (fscanf(fp, "%f %f %f %f %f %f %f %f ", &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8) != EOF) {
        if ((i%100) == 0 ){
            printf("Trainig examples done: %d\n", i);
        }
        s.x[0] = t1;
        s.y[0] = t2;
        s.x[1] = t3;
        s.y[1] = t4;

        s.vx[0] = t5;
        s.vy[0] = t6;
        s.vx[1] = t7;
        s.vy[1] = t8;

        // Call controller on the file
        optimize(&s, a, h, c, param, options); // Find the best sequence of acceleration
        
        //Write to out file.
        fprintf(fpO, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", t1, t3, t2, t4, t5, t7, t6, t8, a[0].ax[0], a[0].ax[1], a[0].ay[0], a[0].ay[1]);
        i++;
    }
    fclose(fp);
    fclose(fpO);
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
 * Calculates the separation metric at a specified state.
 * @param s The state at which to compute separation metric.
 * @param distances The pair-wise distances matrix.
 * @param param The fitness parameters for calculating separation metric.
 * @return Returns the separation metric at the specified state.
 */
double separationExponential(State *s, double distances[s->numBirds][s->numBirds], FitnessParameter *param)
{
    int n = s->numBirds;
    int count = 0;
    double sp = 0;
    double radius = param->rs;
    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++){
            if ( i != j && distances[i][j] < radius){
                sp += pow(param->ms, -3 * distances[i][j] / radius) ;
                count++; 
            }
        }
    }
    if (count == 0){
        return 0;
    }
    return sp / (double)count;
}

/**
 * Calculates the separation metric at a specified state.
 * @param s The state at which to compute separation metric.
 * @param sqDistances The pair-wise squared distances matrix.
 * @param param The fitness parameters for calculating separation metric.
 * @return Returns the separation metric at the specified state.
 */
double seperationPolynomial(State *s, double sqDistances[s->numBirds][s->numBirds], FitnessParameter *param)
{
    int n = s->numBirds;
    int count = 0;
    double sp = 0;
    double sqRadius = param->rs * param->rs;
    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++) {
            if (sqDistances[i][j] < sqRadius) {
               sp += 1.0f / sqDistances[i][j];
               count++; 
            }
        }
    }
    if (count == 0) {
         return 0;
    }
    // count = n * (n-1) / 2;
    return sp; //(double)count;
}

double separationContrained(State *s, double sqDistances[s->numBirds][s->numBirds], FitnessParameter *param, Constraint *c){

    int n = s->numBirds;
    // printf("n in spC: %d\n", n);
    int count = 0;
    double sp = 0;
    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++) {
            sp += pow(MAX(c->dmin - sqrt(sqDistances[i][j]), 0), 2);
        }
    }
    // count = n * (n-1) / 2;
    return sqrt(sp);// / (double)count;
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
 * Calculates the average squared distance metric at a specified state.
 * @param s The state at which to compute average squared distance metric.
 * @param sqDistances The pair-wise squared distances matrix.
 * @param param The fitness parameters for calculating average squared distance metric.
 * @return Returns the average squared distance metric at the specified state.
 */
double averageSquaredDistance(State *s, double sqDistances[s->numBirds][s->numBirds], FitnessParameter *param)
{
    int n = s->numBirds;
    // printf("n in asd: %d\n", n);
    int count = 0;
    double sd = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if ( i != j) { 
                sd += sqDistances[i][j];
            }
        }
    }
    count = n * (n-1) / 2;
    return sd / (double)count;
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
 * Calculates sum of minimum distance constraint violations of the flock and point-like obstacles.
 * @param s The state at which to compute average distance.
 * @param R The rectangular obstacle.
 * @param param The fitness parameters.
 * @return Returns the average distance of the whole flock at the Rectangle.
 */
double pointObstacleCost(State *s, pointObstacles *P, FitnessParameter *param)
{
    double dMin = 2.0;

    int n = s->numBirds;
    int numO = P->n - 1;
    int count = 0;
    double dist = 0;
    // double sqRadius = param->rc * param->rc;
    for (int k = 0; k < numO; k++){
        for (int i = 0; i < n; i++) {
            double A2P = norm(s->x[i] - P->px[k], s->y[i] - P->py[k]);
            if (A2P > 0){
            //     dist += 1/p2R;
                dist += MAX(0, dMin - A2P);
            }
        }
    }
    // count = numO * n * (n-1) / 2; //Change this
    // dist /= (double)count;
    return dist;
}

/**
 * Calculates sum of minimum distance constraint violations of the flock and rectangular obstacles.
 * @param s The state at which to compute average distance.
 * @param R The rectangular obstacle.
 * @param param The fitness parameters.
 * @return Returns the average distance of the whole flock at the Rectangle.
 */
double obstacleCost1(State *s, Obstacles *O, FitnessParameter *param)
{
    double dMin = 2.0;

    int n = s->numBirds;
    int numO = O->n;
    int count = 0;
    double dist = 0;
    // double sqRadius = param->rc * param->rc;
    for (int k = 0; k < numO; k++){
        for (int i = 0; i < n; i++) {
            double p2R = pointToRectangle(&(O->R[k]), s->x[i], s->y[i]);
            if (p2R >= 0){
            //     dist += 1/p2R;
                dist += pow(MAX(0, dMin - p2R), 2);
            }
        }
    }
    // count = numO * n * (n-1) / 2; //Change this
    // dist /= (double)count;
    return sqrt(dist);
}

double obstacleCost2(State *s, Rectangle *R, FitnessParameter *param)
{
    int n = s->numBirds;
    int count = 0;
    double dist = 0;
    // double sqRadius = param->rc * param->rc;
    for (int i = 0; i < n; i++) {
        double p2R = pointToRectangle(R, s->x[i], s->y[i]);
        if (p2R > 0){
            dist += 1/p2R;
        }
    }
    count = n * (n-1) / 2;
    dist /= (double)count;
    return dist;
}

/**
 * Calculates the average um of squares of distances of the flock to a target.
 * @param s The state at which to compute average distance.
 * @param target[2] is the target location
 * @return Returns the average squared distance of the whole flock at the Rectangle.
 */
double target(State *s,  double target[2]){
    int n = s->numBirds;
    double ret = 0;
    for(int i = 0; i < n-1; i++){
        ret += norm(s->x[i] - target[0], s->y[i] - target[1]);
    }
    return ret / (double)n;
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
                  "                  Default config file is 'gmpc.conf'.\n"
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
    printf("----------GMPC\n");
    int printNum = 0;

    Obstacles O;
    char* filename1 = "rectangles.txt";
    readRectObstacles(filename1, &O);

    // pointObstacles P;
    // char* filename2 = "pointObstacles25.txt";
    // readPointObstacles(filename2, &P);

    for (int i = 0; i<O.n; i++){
        printf("Rectangle obstacles\n");
        printf("%f\t%f\t%f\t%f\t%f\t\n", O.R[i].cx, O.R[i].cy, O.R[i].w, O.R[i].l, O.R[i].theta);
    }

    // for (int i = 0; i<P.n; i++){
    //     printf("Point obstacles\n");
    //     printf("%f\t%f\n", P.px[i], P.py[i]);
    // }
    double timeBegin = clock();
    int numSims = 1;
    int startNum = 1;

    char *confFile = "C Code/gmpc.conf";
    char *initFileTemplate = NULL;
    char *resultFileTemplate = NULL;
    bool quietMode = false;
    
    int opt;
    while ((opt = getopt(argc, argv, "b:s:c:qh")) != -1) {
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
    // int maxThreads = 1;
    readGlobalMpcConfig(confFile, &c, &param, &gdOptions, &h, &steps);

    //Generate and save synthetic training data.
    genTraningData(h, &c, &param, &gdOptions);

    c.obs = &O; 
    // c.pointObs = &P;
    // Manually set noise config
    noiseCfg.positionSensorNoise.enabled = true;
    noiseCfg.positionSensorNoise.uniform = false;
    noiseCfg.positionSensorNoise.independent = true;
    noiseCfg.positionSensorNoise.mu = 0;
    noiseCfg.positionSensorNoise.sigma = 2.0;
    noiseCfg.positionSensorNoise.lower = -2;
    noiseCfg.positionSensorNoise.upper = 2;
    
    noiseCfg.velocitySensorNoise.enabled = true;
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

    // printf("init_file:%s\n", initFileTemplate);
    // printf("out_file:%s\n", resultFileTemplate);

    if (!quietMode) {
        printf ("You have %d logical cores.\n", numCpuCores);
        printf("horizon = %d \t\t time steps = %d \t\t time step = %.4f \t\t control step = %.3f\n", h, steps, c.dt, c.ct);
        printf("amax = %.3f \t\t vmax = %.3f \t\t dmin = %f\n", c.amax, c.vmax, c.dmin);
        printf("rc = %.3f \t\t wc = %.3f \t\t mc = %.3f\n", param.rc, param.wc, param.mc);
        printf("rs = %.3f \t\t ws = %.3f \t\t ms = %.3f\n", param.rs, param.ws, param.ms);
        printf("rv = %.3f \t\t wv = %.3f \t\t wu = %.3f\n", param.rv, param.wv, param.wu);
        printf("predAt = %.3f \t\t predAv = %.3f \t\t predDist = %.3f\n", param.predAt, param.predAv, c.predAvD);
        printf("wt = %.3f \t\t wspc = %.3f \t\t wob = %.3f\n", param.wt, param.wspc, param.wob);
        printf("lR = %.6f \t\t delta = %.9f \t\t maxIters = %d \t\t precision = %.4f\n", 
               gdOptions.learningRate, gdOptions.delta, gdOptions.maxIters, gdOptions.precision);
    }
    
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
        simParams[i].resultFileTemplate = resultFileTemplate;
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
        //runMPC(simNum, initFileTemplate, resultFileTemplate, h, steps, c, param, gdOptions);
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
    double timeEnd = clock();
    double timeSeconds = (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC;
    if (!quietMode) {
        printf("Total running time = %f seconds\n", timeSeconds);
    }

    printf("Malloc: %ld\t Free: %ld\n", malloc_count, free_count);
}
