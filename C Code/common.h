#ifndef _COMMON_H_
#define _COMMON_H_

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
// #include "normal.h"

#ifdef _WIN32
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define PI 3.14159265358979323846

// BAIG WAS HERE
extern long free_count;
extern long malloc_count;

/**
 * Defines the set of obstacles(points).
 */
typedef struct pointObstacles{
    int n;
    double  *px;
    double  *py;
} pointObstacles;

/**
 * Defines the set of obstacles(rectangles).
 */
typedef struct Obstacles{
    int n;
    struct Rectangle *R;
} Obstacles;
/**
 * Defines A rectangle.
 */
typedef struct Rectangle{
    double cx;
    double cy;
    double theta;
    double l;
    double w;
} Rectangle;
/**
 * Defines system state.
 */
typedef struct State { 
    int numBirds;      // Number of birds.
    double *x;  // Array of x-coordinates of all birds.
    double *y;  // Array of y-coordinates of all birds.
    double *vx; // Array of vx components of velocity of all birds.
    double *vy; // Array of vy components of velocity of all birds.
    int *policy; //0 for AC(MPC) and 1 for BC(Reynbold's).
} State;

/**
 * Defines control action.
 */
typedef struct Action {
    int numBirds;      // Number of birds.
    double *ax; // Array of ax components of acceleration of all birds.
    double *ay; // Array of ay components of acceleration of all birds.
} Action;

typedef struct Direction{
    int numBirds;
    double *d_vx;
    double *d_vy;
    double *d_ux;
    double *d_uy;
}Direction;

/**
 * Defines system constraints.
 */
typedef struct Constraint {
    double vmax; // Maximum velocity
    double amax; // Maximum acceleration
    double dt;   // Time step. TODO: find a better place to house time step.
    double ct;  // Control step.
    double dmin;
    double predAvD; // Distance maintained by the agents with the predator. 
    double alpha;
    double beta;
    Obstacles *obs;
    pointObstacles *pointObs;
    int predator; // 1 If last agent is predator, else 0
    double pFactor; // maximum acc and max speed factor for predator  


} Constraint;

/**
 * Defines parameters used in fitness function.
 */
typedef struct FitnessParameter {
    double rc;  // Radius of local neighborhood for cohesion calculation.
    double msc;  
    double rs;  // Radius of local neighborhood for separation calculation.
    double ms;
    double mc;
    double rv;  // Radius of local neighborhood for velocity matching calculation.
    double wc;  // Weight of cohesion term in the fitness function.
    double ws;  // Weight of separation term in the fitness function.
    double wv;  // Weight of velocity matching term in the fitness function.
    double wu;  // Weight of control input term in the cost function.
    double wd;  // Weight of the average squared distance term in the cost function.
    double predAt; // Weight of the predator attack term. 
    double predAv; // Weight of the predator avoid term. 
    double wt; // Weight of th target term 
    double wspc; // Weight of collision-avoidance penalty term
    double wob; // Weight of obstacle avoidance penalty term
    int knn; // Number of nearest neighbors
} FitnessParameter;

/**
 * Defines the configurations for logging.
 */
typedef struct LogConfig {
    char *logFileX;         // Name of the log file to write x-coordinates of the birds.
    char *logFileY;         // Name of the log file to write y-coordinates of the birds.
    char *logFileVx;        // Name of the log file to write vx.
    char *logFileVy;        // Name of the log file to write vy.
    char *logFileAx;        // Name of the log file to write ax.
    char *logFileAy;        // Name of the log file to write ay.
    char *logFileFitness;   // Name of the log file to write fitness values.
    char *logFileLieDerivative; // Name of the log file to write lie-derivative.
    char *logFilePolicy;    //The control policy of the agent. AC/BC
} LogConfig;

/**
 * Defines options for gradient descent algorithm.
 */
typedef struct GradientDescentOptions {
    double learningRate;
    double delta;          // The delta for approximating gradient. fgrad(x) = (f(x + delta) - f(x - delta)) / (2 * delta);
    double precision;      // Stop gradient descent loop when ||current_x - previous_x|| <= precision.
    int maxIters;          // Maximum number of iterations.
} GradientDescentOptions;

/**
 * Defines parameters for noise.
 */
typedef struct NoiseParameters {
    bool enabled;       // A flag indicating whether noise is enabled.
    bool independent;   // A flag indicating whether noise is independent across agents.
    bool uniform;       // A flag indicating: if true, noise is uniformly distributed, if false, noise is normally distributed.
    double mu;          // If uniform is false, mean of Gaussian distribution.
    double sigma;       // If uniform is false, standard deviation of Gaussian distribution.
    double lower;       // If uniform is true, lower bound on the range of uniform distribution.
    double upper;       // If uniform is true, upper bound on the range of uniform distribution.
    
} NoiseParameters;

/**
 * Defines noise configurations.
 */
typedef struct NoiseConfig {
    NoiseParameters positionSensorNoise;
    NoiseParameters velocitySensorNoise;
    NoiseParameters actuatorNoise;
} NoiseConfig;

void appendToFile(char * name, double val);
void computeDistances(State *s, double distances[s->numBirds][s->numBirds]);
void computeSquaredDistances(State *s, double sqDistances[s->numBirds][s->numBirds]);
void computeSquaredDistancesInNeighborhood(State *s, double sqDistances[s->numBirds][s->numBirds], int birdIdx, int *neighborIds, int numNeighbors);
void computeNeighbors(State *s, int neighbors[s->numBirds][s->numBirds], int numNeighbors[s->numBirds], double sqDistances[s->numBirds][s->numBirds], double radius);
void getKNN(State *s, int neighbors[s->numBirds][s->numBirds], double sqDistances[s->numBirds][s->numBirds], FitnessParameter *param);
double pointToRectangle(Rectangle *R, double x, double y);
int compare ( const void *pa, const void *pb ); 
State dynamics(State *s, Action *a, Constraint *c);
State dynamics_max_braking(State *s, Action *a, Constraint *c);
double braking_distance(State *s, int i, int j, Constraint *c);
double braking_distance_worst(State *s, int i, int j, double pos[2], double vel[2], Constraint *c);
double braking_distance_propotional(State *s, int i, int j, double pos[2], double vel[2], Constraint *c);

void state_to_adjacency_matrix( int n, double distances[n][n], int M[n][n], double rs);
void print_matric(int n, int M[n][n]);
int is_connected(int n, int M[n][n]);


/**
 * Calculates the cost of a specified control action. The cost is defined as the 
 * sum of squares of control action magnitudes.
 * @param a The control action to calculate the cost.
 * @return Return the control action cost, which is the sum of squares of control action magnitudes.
 */
double controlCost(Action *a);
void logStates(State *states, int steps, LogConfig *logCfg);
void logActions(Action *actions, int steps, LogConfig *logCfg);
void logFitnesses(double *fitnesses, int steps, LogConfig *logCfg);
void logLD(double *LD, int steps, LogConfig *logCfg);


State cloneState(State *s) ;
void addNoiseToState(State *s, NoiseConfig *noise, uint32_t *seed, uint32_t kn[128], float fn[128], float wn[128]);

static void AddRect(Rectangle *r, double cx, double cy, double w, double l, double theta ){
    r->cx = cx;
    r->cy = cy;
    r->l = l;
    r->w = w;
    r->theta = theta;
}

static void initPointObstacles(pointObstacles *P, int n){
    P->n = n;
    P->px = (double *)malloc(n * sizeof(double)); 
    P->py = (double *)malloc(n * sizeof(double)); 
    //Baig was here
    malloc_count = malloc_count + 2;
}

static void initObstacles(Obstacles *O, int n){
    O->n = n;
    O->R = (Rectangle *)malloc(n * sizeof(Rectangle)); 
    //Baig was here
    malloc_count = malloc_count + 1;
}

/**
 * Initializes a state structure to hold data for a specified number of birds.
 * @param s The state to initialize.
 * @param numBirds The number of birds.
 * @remark Memory is allocated using malloc, so the caller is responsible for freeing the memory.
 *         This method makes no attempts to set initial values for the state variables.
 *         The caller should assume the state variables hold random values after initializing with this method.
 */
static void initState(State *s, int numBirds)
{
    s->numBirds = numBirds;
    s->x = (double *)malloc(numBirds * sizeof(double));
    s->y = (double *)malloc(numBirds * sizeof(double));
    s->vx = (double *)malloc(numBirds * sizeof(double));
    s->vy = (double *)malloc(numBirds * sizeof(double));
    s->policy = (int *)calloc(numBirds, sizeof(int)); //Zero initialized
    malloc_count += 5;
}

/**
 * Initializes a control action structure to hold data for a specified number of birds.
 * @param a The control action to initialize.
 * @param numBirds The number of birds.
 * @remark Memory is allocated using malloc, so the caller is responsible for freeing the memory.
 *         This method makes no attempts to set initial values for the fields of the struct.
 *         The caller should assume the fields hold random values after initializing with this method.
 */
static void initAction(Action *a, int numBirds)
{
    a->numBirds = numBirds;
    a->ax = (double *)malloc((numBirds) * sizeof(double));
    a->ay = (double *)malloc((numBirds) * sizeof(double));
    malloc_count += 2;
}

static void freeObstacles(Obstacles *O){
    free(O->R);
}

static void freePointObstacles(pointObstacles *P){
    free(P->px);
    free(P->py);
}
/**
 * Frees the memory allocated to the pointers in a specified state.
 * @param s The state to free memory.
 */
static void freeState(State *s)
{
    free(s->x);
    free(s->y);
    free(s->vx);
    free(s->vy);
    free(s->policy);

    // BAIG WAS HERE
    free_count += 5;

    s->numBirds = 0;
}

/**
 * Frees the memory allocated to the pointers in a control action.
 * @param a The control action to free memory.
 */
static void freeAction(Action *a)
{
    free(a->ax);
    free(a->ay);

    // BAIG WAS HERE
    free_count += 2;
}

/**
 * Frees the memory allocated to the pointers in a Direction.
 * @param d The Direction to free memory.
 */
static void freeDirection(Direction *d)
{
    free(d->d_vx);
    free(d->d_vy);    
    free(d->d_ux);
    free(d->d_uy);

    // BAIG WAS HERE
    free_count += 4;
}

/**
 * Adds two control actions and stores the sum in another control action.
 * @param a1 The first control action to add.
 * @param a2 The second control action to add.
 * @param sum The control action that will hold the sum of a1 and a2.
 * @remark The sum must be initialized, e.g., by calling initAction, before calling this method.
 *         It is assumed that a1 and a2 have the same numBirds.
 *         sum->numBirds is not changed.
 */
static void addActions(Action *a1, Action *a2, Action *sum)
{
    for (int i = 0; i < a1->numBirds; i++) {
        sum->ax[i] = a1->ax[i] + a2->ax[i];
        sum->ay[i] = a1->ay[i] + a2->ay[i];
    }
}


/**
 * Returns a pseudo-random number in the specified range.
 * @param min Lower bound of the range, inclusive.
 * @param max Upper bound of the range, inclusive.
 * @returns A pseudo-random number within the specified range.
 */
static inline double randd(double min, double max) 
{
    return ((double) rand() / (double) RAND_MAX) * (max - min) + min;
}

/**
 * Generates an array of pseudo-random numbers in the specified range.
 * @param min Lower bound of the range, inclusive.
 * @param max Upper bound of the range, inclusive.
 * @param len The number of numbers to store to the specified result array.
 * @param r The array to store generated pseudo-random numbers.
 */
static void randdVector(double min, double max, int len, double r[len])
{
    for (int i = 0; i < len; i++) {
        r[i] = randd(min, max);
    }
}   

// /**
 // * Returns a pseudonormal number with specified mean and standard deviation.
 // * @param mu Mean of the normal distribution.
 // * @param sigma Standard deviation of the normal distribution.
 // * @param seed Pointer to a random seed for the random generator to use and change.
 // * @returns A pseudonormal number with the specified mean and standard deviation.
 // */
// static inline double normrnd(double mu, double sigma, int *seed)
// {
    // return r8_normal_ab(mu, sigma, seed);
// }

// /**
 // * Generates an array of pseudonormal number with specified mean and standard deviation.
 // * @param mu Mean of the normal distribution.
 // * @param sigma Standard deviation of the normal distribution.
 // * @param seed Pointer to a random seed for the random generator to use and change.
 // * @param len The number of numbers to store to the specified result array.
 // * @param r The array to store generated pseudonormal numbers.
 // */
// static void normrndVector(double mu, double sigma, int *seed, int len, double r[len])
// {
    // for (int i = 0; i < len; i++) {
        // r[i] = normrnd(mu, sigma, seed);
    // }
// } 

static uint32_t initRandomSeed()
{
    return (uint32_t)time(NULL);
}

/**
 * Returns the norm-2 of a vector of length 2.
 * @param x The first element of the vector to compute norm-2.
 * @param y The second element of the vector to compute norm-2.
 * @returns norm-2 value of the specified vector.
 */
static inline double norm(double x, double y) 
{
    return sqrt(x * x + y * y);
}

/**
 * Returns the dot product of two vectors.
 * @param x1 The first element of the first vector.
 * @param y1 The second element of the first vector.
 * @param x2 The first element of the second vector.
 * @param y2 The second element of the second vector.
 * @return Returns the dot product of the specified vectors.
 */
static inline double dot(double x1, double y1, double x2, double y2) 
{
    return x1*x2 + y1*y2;
}

/**
 * Normalizes a vector of length 2.
 * @param x The first element of the vector.
 * @param y The second element of the vector.
 */
static void normalize(double *x, double *y)
{
    double mag = norm(*x, *y);
    if (mag != 0) {
        *x /= mag;
        *y /= mag;
    }
}

/**
 * Sets the magnitude a vector of length 2.
 * @param x The first element of the vector.
 * @param y The second element of the vector.
 * @param mag The magnitude of the resultant vector.
 */
static void set_magnitude(double *x, double *y, double mag)
{
    normalize(x, y);
    *x *= mag;
    *y *= mag;
}

/**
 * Returns multivariate PDF. Immitates Matlab's mvnpdf function.
 * @param x The first variable.
 * @param y The second variable.
 * @param a
 * @param b
 */
static inline double mvnpdf(double x, double y, double a, double b) 
{
    return exp(-0.5 * (a * x * x + b * y * y));
}

/**
 * Ensures that the norm-2 of the specified vector is less than or equal to the specified magnitude.
 * @param x Pointer to the first element of the vector to trim.
 * @param y Pointer to the second element of the vector to trim.
 * @param mag Upper limit on the norm-2 of the vector. If norm-2 of the specified vector exceeds
 *            mag, then x and y will be scaled down.
 */
static void trim(double *x, double *y, double mag)
{
    double normXY = norm(*x, *y);
    if (normXY <= mag)
        return;
    *x *= mag / normXY;
    *y *= mag / normXY;
}

/**
 * Ensures that the norm-2 of the specified vector is less than or equal to the specified magnitude.
 * @param x Pointer to the first element of the vector to trim.
 * @param y Pointer to the second element of the vector to trim.
 * @param mag Upper limit on the norm-2 of the vector. If norm-2 of the specified vector exceeds
 *            mag, then x and y will be scaled down.
 */
static void trimActionSeqNaive(Action *a, int h, double amax)
{
    int n = a[0].numBirds;
    for (int j = 0; j < h; j++) {
        for (int k = 0; k < n; k++) { 
            double normA = norm(a[j].ax[k], a[j].ay[k]);
            if (normA > amax) {
                a[j].ax[k] *= amax / normA;
                a[j].ay[k] *= amax / normA;
            }
        }
    }
}

/**
 * Ensures that the norm-2 of the specified vector is less than or equal to the specified magnitude.
 * @param x Pointer to the first element of the vector to trim.
 * @param y Pointer to the second element of the vector to trim.
 * @param mag Upper limit on the norm-2 of the vector. If norm-2 of the specified vector exceeds
 *            mag, then x and y will be scaled down.
 */
static void trimActionSeq(Action *a, int h, double amax)
{
    int n = a[0].numBirds;

    double * maxAgentAcc;   
    maxAgentAcc = (double *)malloc(h*sizeof(double));
    memset(maxAgentAcc, 0, h*sizeof(double));

    for (int j = 0; j < h; j++) {
        double temp = 0;
        for (int k = 0; k < n; k++) { 
            temp = norm(a[j].ax[k], a[j].ay[k]);
            if (temp > maxAgentAcc[j]){
                maxAgentAcc[j] = temp;
            }            
        }
    }

    for (int j = 0; j < h; j++) {
        if (maxAgentAcc[j] > amax){
            for (int k = 0; k < n; k++) { 
                a[j].ax[k] *= amax / maxAgentAcc[j];
                a[j].ay[k] *= amax / maxAgentAcc[j];   
            }
        }
    } 
    free(maxAgentAcc);
}

/**
 * Concaternates two strings.
 * @param s1 The first string.
 * @param s2 The second string.
 * @return A string resulted from appending s2 to the end of s1.
 * @remark Caller is responsible for freeing the memory after using the returned string.
 */
static char* concat(char *s1, char *s2){
    char *result = (char*)malloc(strlen(s1) + strlen(s2)+1);
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

/**
 * Converts a string to lower case.
 * @param s The string to convert to lower case.
 */
static inline void str2lower(char *s)
{
    char *p = s;
    for ( ; *p; ++p) *p = tolower(*p);
}

/**
 * Trims leading and trailing white spaces off a string.
 * @param str The string to trim.
 * @return Returns a pointer to a substring of the original string.
 * @remark This function returns a pointer to a substring of the original string.
 * If the given string was allocated dynamically, the caller must not overwrite
 * that pointer with the returned value, since the original pointer must be
 * deallocated using the same allocator with which it was allocated.  The return
 * value must NOT be deallocated using free() etc.
 * @reference https://stackoverflow.com/questions/122616/how-do-i-trim-leading-trailing-whitespace-in-a-standard-way
 */
static char *strtrim(char *str)
{
  char *end;

  // Trim leading space
  while (isspace((unsigned char)*str)) str++;

  if (*str == 0)  // All spaces?
    return str;

  // Trim trailing space
  end = str + strlen(str) - 1;
  while (end > str && isspace((unsigned char)*end)) end--;

  // Write new null terminator
  *(end + 1) = 0;

  return str;
}


int getNumberOfCores();

#endif
