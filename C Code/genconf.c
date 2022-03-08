/*
 * genconf.c
 *
 * Purpose: Generate random initial configurations of the birds and save them to text files.
 *
 * How to compile: 
 *     gcc genconf.c conf.c -O3 -o genconf
 * On Windows, use MinGW or cygwin to compile:
 *     gcc genconf.c conf.c -O3 -o genconf.exe
 * How to use: run 
 *     ./genconf
 * to print usage.
 *
 */

#include <stdbool.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <getopt.h> 
#include "common.h"
#include "conf.h"

/**
 * 
 * 
 * 
 */
double braking_distance_new(double P_ij[2], double v_ij[2], Constraint *c)
{  

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
 * Initializes a random initial configuration.
 * @param x Pointer to an array to hold the initial x-coordinates.
 * @param y Pointer to an array to hold the initial y-coordinates.
 * @param vx Pointer to an array to hold the initial velocities in the x-direction.
 * @param vy Pointer to an array to hold the initial velocities in the y-direction.
 * @param numBirds Number of birds.
 * @param dmin Minimum distance between any two birds.
 * @param xmin Minimum x-coordinate of the bounding box.
 * @param xmax Maximum x-coordinate of the bounding box.
 * @param ymin Minimum y-coordinate of the bounding box.
 * @param ymax Maximum y-coordinate of the bounding box.
 * @param vmin Minimum initial velocity.
 * @param vmax Maximum initial velocity.
 * @param targets Struct containig the target locations for each agent.
 */
void init_circle(double *x, double *y, double *vx, double *vy, 
          int numBirds, double dmin, double xmin, double xmax, 
          double ymin, double ymax, double vminx, double vmaxx,
          double vminy, double vmaxy, pointObstacles *targets) {
    double delta_theta = (2 * PI) / numBirds;
    double theta = 0;
    double radius = (2 * dmin * numBirds) / (2 * PI);

    for (int i = 0; i < numBirds; i++){
        x[i] = radius * cos(theta);
        y[i] = radius * sin(theta);
        vx[i] = 0;
        vy[i] = 0;

        targets->px[i] = radius * cos(theta + PI);
        targets->py[i] = radius * sin(theta + PI);

        theta += delta_theta;
    }
}

/** 
 * Initializes a random initial configuration.
 * @param x Pointer to an array to hold the initial x-coordinates.
 * @param y Pointer to an array to hold the initial y-coordinates.
 * @param vx Pointer to an array to hold the initial velocities in the x-direction.
 * @param vy Pointer to an array to hold the initial velocities in the y-direction.
 * @param numBirds Number of birds.
 * @param dmin Minimum distance between any two birds.
 * @param xmin Minimum x-coordinate of the bounding box.
 * @param xmax Maximum x-coordinate of the bounding box.
 * @param ymin Minimum y-coordinate of the bounding box.
 * @param ymax Maximum y-coordinate of the bounding box.
 * @param vmin Minimum initial velocity.
 * @param vmax Maximum initial velocity.
 */
void init(double *x, double *y, double *vx, double *vy, 
          int numBirds, double dmin, double xmin, double xmax, 
          double ymin, double ymax, double vminx, double vmaxx,
          double vminy, double vmaxy) {
    
    Constraint c;
    FitnessParameter param;
    char *confFile = "C Code/gmpc.conf";
    readConstraints(confFile, &c);
    readParams(confFile, &param);

    // x[numBirds - 1] = 0;
    // y[numBirds - 1] = 50;
    // vx[numBirds - 1] = 0;
    // vy[numBirds - 1] = 0;
    int isCollision;
    do {
        isCollision = 0;
        for (int j = 0; j < numBirds; j++) {
            x[j] = randd(xmin, xmax);
            y[j] = randd(ymin, ymax);

            vx[j] = randd(vminx, vmaxx);
            vy[j] = randd(vminy, vmaxy);

            // for (int ki = 0; ki < j-1; ki++){
            //     double P_ij[2] = {x[ki]-x[j]};
            //     double v_ij[2] = {y[]-y[bj]};            
            // }
        }
        for (int bi = 0; bi < numBirds; bi++) {
            for (int bj = bi + 1; bj < numBirds; bj++) {
                // double d_brake = braking_distance(x[bi], y[bi], x[bj], y[bj], vx[bi], vy[bi], vx[bj], vy[bj], &c);
                double P_ij[2] = {x[bi]-x[bj], y[bi]-y[bj]};
                double v_ij[2] = {vx[bi]-vx[bj], vy[bi]-vy[bj]};
                double d_brake = braking_distance_new(P_ij, v_ij, &c);
                if (norm(x[bi]-x[bj], y[bi]-y[bj]) - d_brake < c.dmin)
                {
                    isCollision = 1;
                    break;
                }
            }
            if (isCollision)
                break;
        }
        // Check implemented for connectivity of flock
        State s;
        initState(&s, numBirds);
        s.numBirds = numBirds;
        s.x = x;
        s.y = y;
        s.vx = vx;
        s.vy = vy;

        //Proximity-net is connected. Connected components = 1.
        double sqDistances[numBirds][numBirds];
        int M[numBirds][numBirds];
        // int M[5][5] = {
        //     {0,0,0,0,1},
        //     {0,0,1,0,0},
        //     {0,1,0,1,0},
        //     {0,0,1,0,1},
        //     {1,0,0,1,0},
        // };
        computeSquaredDistances(&s, sqDistances);
        state_to_adjacency_matrix(numBirds, sqDistances, M, param.rs);
        // print_matric(numBirds, M);

        int connected = is_connected(numBirds, M);
        // printf("connected: %d\n", connected);
        if (connected == 0){
            isCollision = 1;
        }
    } while (isCollision);
}

void printUsage()
{
    printf("genconf: Generate a random initial configuration (position and velocity) of the birds.\n\n");
    printf("Usage: genconf [-xXyYvVdnbqh] file\n");
    printf("    -x xmin   Minimum x-coordinate of the bounding box.\n");
    printf("    -X xmax   Maximum x-coordinate of the bounding box.\n");
    printf("    -y ymin   Minimum y-coordinate of the bounding box.\n");
    printf("    -Y ymax   Maximum y-coordinate of the bounding box.\n");
    printf("    -v vmin   Minimum initial velocity.\n");
    printf("    -V vmax   Maximum initial velocity.\n");
    printf("    -d dmin   Minimum distance between two birds.\n");
    printf("    -n n      Number of birds.\n");
    printf("    -b nconf  Number of configurations (batch mode).\n");
    printf("    -q        Quiet mode.\n");
    printf("    -h        Print help.\n");
    printf("    file      File name to write results to. In batch mode, use '%%d' to specify the configuration number.\n\n");
    printf("Example 1: Generate one random initial configuration\n");
    printf("    genconf -x -1 -y -1 -n 100 init_conf.txt\n\n");
    printf("Example 2: Generate 50 random initial configurations and save them to init_conf_1.txt, init_conf_2.txt, ...\n");
    printf("    genconf -x -1 -y -1 -n 100 -b 50 init_conf_%%d.txt\n");
}

/**
 * Checks input from command line arguments for errors.
 * @param xmin Minimum x-coordinate of the bounding box.
 * @param xmax Maximum x-coordinate of the bounding box.
 * @param ymin Minimum y-coordinate of the bounding box.
 * @param ymax Maximum y-coordinate of the bounding box.
 * @param vmin Minimum initial velocity.
 * @param vmax Maximum initial velocity.
 * @param numBirds Number of birds.
 * @param dmin Minimum distance between any two birds.
 * @param filename Name of the file to write configuration to.
 * @return Returns 0 if no errors. Otherwise returns the number of errors.
 */
int checkInput(double xmin, double xmax, double ymin, double ymax, double vminx, double vmaxx, 
                double vminy, double vmaxy, int numBirds, double dmin, char *filename)
{
    int ret = 0;
    if (xmin >= xmax)
    {
        ret += 1;
        printf("ERROR: xmin must be less than xmax\n");
    }
    if (ymin >= ymax)
    {
        ret += 1;
        printf("ERROR: ymin must be less than ymax\n");
    }
    if (vminx >= vmaxx)
    {
        ret += 1;
        printf("ERROR: vminx must be less than vmaxx\n");
    }
    if (vminy >= vmaxy)
    {
        ret += 1;
        printf("ERROR: vminy must be less than vmaxy\n");
    }
    if (numBirds < 2)
    {
        ret += 1;
        printf("ERROR: number of birds must be greater than 1\n");
    }
    if (dmin < 0)
    {
        ret += 1;
        printf("ERROR: dmin must be non-negative\n");
    }
    return ret;
}

int main(int argc, char *argv[]) 
{
    bool quietMode = false;
    int numConfs = 1; // Number of configurations to generate
    double xmin = 0;
    double xmax = 1;
    double ymin = 0;
    double ymax = 1;
    double vminx = 0;
    double vmaxx = 2;    
    double vminy = 0;
    double vmaxy = 2;
    double dmin = 2.0;
    int numBirds = 2;
    char *filename;
    char *confFile = "conf.conf";
    int opt;
    while ((opt = getopt(argc, argv, "x:X:y:Y:u:U:v:V:d:n:b:qh")) != -1) {
        switch (opt) {
        case 'x': xmin = atof(optarg); break;
        case 'X': xmax = atof(optarg); break;
        case 'y': ymin = atof(optarg); break;
        case 'Y': ymax = atof(optarg); break;
        case 'u': vminx = atof(optarg); break;
        case 'U': vmaxx = atof(optarg); break;
        case 'v': vminy = atof(optarg); break;
        case 'V': vmaxy = atof(optarg); break;
        case 'd': dmin = atof(optarg); break;
        case 'n': numBirds = atoi(optarg); break;
        case 'b': numConfs = atoi(optarg); break;
        case 'q': quietMode = true; break;
        case 'h': printUsage(); exit(EXIT_SUCCESS);
        default:
            printUsage();
            exit(EXIT_FAILURE);
        }
    }
    
    if (optind < argc) {
        filename = argv[optind];
    }
    else {
        printUsage();
        exit(EXIT_FAILURE);
    }
    
    if (!quietMode) {
        printf("Bounding box: xmin = %f, xmax = %f, ymin = %f, ymax = %f\n", 
               xmin, xmax, ymin, ymax);
        printf("Velocity bound: vminx = %f, vmaxx = %f\n", vminx, vmaxx);
        printf("Minimum distance between two birds: %f\n", dmin);
        printf("Number of birds: %d\n", numBirds);
        printf("Output file: %s\n", filename);
    }
    
    int errors = checkInput(xmin, xmax, ymin, ymax, vminx, vmaxx, vminy, vmaxy, numBirds, dmin, filename);
    
    if (errors != 0) {
        exit(EXIT_FAILURE);
    }
    
    double x[numBirds];
    double y[numBirds];
    double vx[numBirds];
    double vy[numBirds];

    pointObstacles targets;
    initPointObstacles(&targets,  numBirds);

    
    /* Intializes random number generator */
    time_t t;
    srand((unsigned) time(&t));
    
    if (numConfs < 1)
        numConfs = 1;
    
    char fname[500];
    char tname[500];
    char *targetfilename = "experiments_circle/src/target_conf_%d.txt";

    
    for (int i = 0; i < numConfs; i++) {
        init(x, y, vx, vy, numBirds, dmin, xmin, xmax, ymin, ymax, vminx, vmaxx, vminy, vmaxy);
        // init_circle(x, y, vx, vy, numBirds, dmin, xmin, xmax, ymin, ymax, vminx, vmaxx, vminy, vmaxy, &targets);
        sprintf(fname, filename, (i+1));
        sprintf(tname, targetfilename, (i+1));
        if (writeConf(fname, x, y, vx, vy, numBirds) != 0) {
            printf("Error writing configuration to file.\n");
            exit(EXIT_FAILURE);
        }        
        if (writeConfPoints(tname, &targets, numBirds) != 0) {
            printf("Error writing configuration to file.\n");
            exit(EXIT_FAILURE);
        }
    }
    freePointObstacles(&targets);
}
