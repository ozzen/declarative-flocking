#ifndef _CONF_H_
#define _CONF_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "common.h"
#include "reynolds.h"

// BAIG WAS HERE
long free_count;
long malloc_count;

/**
To be added later.
**/
int readPointObstacles(char *filename, pointObstacles *O);

/**
* To be added later.
**/
int readRectObstacles(char *filename, Obstacles *O);

/**
 * Reads parameters for global MPC controller from a file.
 * @param filename Name of the configuration file to read parameters from.
 * @param c The struct to read system constraints into.
 * @param param The struct to read the fitness parameters into.
 * @param options The struct to read gradient descent parameters into.
 * @param h The prediction horizon.
 * @param steps The number of simulation time steps.
 */
int readGlobalMpcConfig(char *filename, Constraint *c, FitnessParameter *param, GradientDescentOptions *options, int *h, int *steps);

/**
 * Reads parameters for Reynold's model.
 * @param filename Name of the configuration file to read parameters from.
 * @param c The struct to read system constraints into.
 * @param param The struct to read the Reynold's model parameters into.
 * @param steps The number of simulation time steps.
 */
int readReynoldsConfig(char *filename, Constraint *c, ReynoldsParameters *param, int *steps);

/**
 * Reads the number of birds from the first line in a specified configuration file.
 * @param filename Name of the configuration file to read the number of birds.
 * @return Returns -1 if there are errors reading the specified file. 
 * Otherwise returns the number of birds from the first line in the specified file.
 */
int readNumBirds(char *filename);

/**
 * Reads configuration (initial position and velocity) of the birds from a specified configuration file.
 * @param filename Name of the configuration file to read.
 * @param x Pointer to the array to write the x-coordinates of the initial positions to.
 *          The array must be large enough to hold the number of birds specified in param numBirds.
 * @param y Pointer to the array to write the y-coordinates of the initial positions to.
 *          The array must be large enough to hold the number of birds specified in param numBirds.
 * @param vx Pointer to the array to write the vx components of the initial velocities to.
 *          The array must be large enough to hold the number of birds specified in param numBirds.
 * @param vy Pointer to the array to write the vy components of the initial velocities to.
 *          The array must be large enough to hold the number of birds specified in param numBirds.
 * @return Returns 0 if the operation finished successfully. Returns -1 otherwise.
 * @remark Use readNumBirds function to read the number of birds stored in the file first. 
 *         Use this information to create the arrays x, y, vx, vy with the correct size.
 *
 */
int readConf(char *filename, double *x, double *y, double *vx, double *vy, int numBirds);

int readParams(char *filename, FitnessParameter *param);

int readConfTarget(char *filename, pointObstacles *T);

/**
 * Reads configuration (initial position and velocity) of the birds from a specified configuration file.
 * @param filename Name of the configuration file to read.
 * @param s Pointer to the State. 
 * @return Returns 0 if the operation finished successfully. Returns -1 otherwise.
 * @remark Use readNumBirds function to read the number of birds stored in the file first. 
 *         Use this information to create the arrays x, y, vx, vy with the correct size.
 *
 */
int readConfState(char *filename, State *s);

/**
 * Writes configuration (initial position and velocity) of the birds to a specified configuration file.
 * @param filename Name of the configuration file to write.
 * @param x Pointer to the array of x-coordinates of the initial positions.
 * @param y Pointer to the array of y-coordinates of the initial positions.
 * @param vx Pointer to the array of vx components of the initial velocities.
 * @param vy Pointer to the array of vy components of the initial velocities.
 * @return Returns 0 if the operation finished successfully. Returns -1 otherwise.
 */
int writeConf(char *filename, double *x, double *y, double *vx, double *vy, int numBirds);

int writeConfPoints(char *filename, pointObstacles *O, int numBirds);

int readConstraints(char *filename, Constraint *c);

#endif