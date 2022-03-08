#include <stdio.h>
#include "conf.h"

int readPointObstacles(char *filename, pointObstacles *P){
    FILE *fp = fopen(filename, "r");
    if (fp){
        char line[50];
        char strVal[50];
        double data[2];

        fgets(line, 50, fp);
        sscanf(line, "%s", strVal);
        P->n = atoi(strVal);

        initPointObstacles(P, P->n); //Allocates Memory

        for(int i = 0; i < P->n; i++){
            for (int j = 0; j<2; j++){
                fgets(line, 50, fp);
                sscanf(line, "%s", strVal);
                data[j] = atof(strVal);
            }
            P->px[i] = data[0];
            P->py[i] = data[1];
        }
        fclose(fp);
        return 0;
    }
    return -1;
}

int readRectObstacles(char *filename, Obstacles *O){
    FILE *fp = fopen(filename, "r");
    if (fp){
        char line[50];
        char strVal[50];
        double data[5];


        fgets(line, 50, fp);
        sscanf(line, "%s", strVal);
        O->n = atoi(strVal);
        initObstacles(O, O->n); //Allocates Memory

        for(int i = 0; i < O->n; i++){
            for (int j = 0; j<5; j++){
                fgets(line, 50, fp);
                sscanf(line, "%s", strVal);
                data[j] = atof(strVal);
            }
            Rectangle R;
            AddRect(&R, data[0], data[1], data[2], data[3], data[4]);
            O->R[i] = R;
        }
        fclose(fp);
        return 0;
    }
    return -1;
}

int readParams(char *filename, FitnessParameter *param)
{
    FILE *fp = fopen(filename, "r");
    if (fp)
    {
        char line[500];
        char *trimmedLine = NULL;
        char paramName[100];
        char strVal[200];
        
        while (fgets(line, 500, fp) != NULL) {
            trimmedLine = strtrim(line);
            if (trimmedLine[0] == '#') // Comment line
                continue;
            sscanf(trimmedLine, "%s %s ", paramName, strVal);
            str2lower(paramName);
            if (0 == strcmp(paramName, "rs")) {
                param->rs = atof(strVal); 
                continue;
            }
        }
        fclose(fp);
        return 0;
    }
    return -1;
}

int readConstraints(char *filename, Constraint *c)
{
    FILE *fp = fopen(filename, "r");
    if (fp)
    {
        char line[500];
        char *trimmedLine = NULL;
        char paramName[100];
        char strVal[200];
        
        while (fgets(line, 500, fp) != NULL) {
            trimmedLine = strtrim(line);
            if (trimmedLine[0] == '#') // Comment line
                continue;
            sscanf(trimmedLine, "%s %s ", paramName, strVal);
            str2lower(paramName);
            if (0 == strcmp(paramName, "amax")) {
                c->amax = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "vmax")) {
                c->vmax = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "dt")) {
                c->dt = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "ct")) {
                c->ct = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "dmin")) {
                c->dmin = atof(strVal); 
                continue;
            }
        }
        fclose(fp);
        return 0;
    }
    return -1;
}

int readGlobalMpcConfig(char *filename, Constraint *c, FitnessParameter *param, GradientDescentOptions *options, int *h, int *steps)
{
    FILE *fp = fopen(filename, "r");
    if (fp)
    {
        char line[500];
        char *trimmedLine = NULL;
        char paramName[100];
        char strVal[200];
        //double dblVal;
        //size_t len = 0;
        //ssize_t read;
        
        while (fgets(line, 500, fp) != NULL) {
            trimmedLine = strtrim(line);
            if (trimmedLine[0] == '#') // Comment line
                continue;
            sscanf(trimmedLine, "%s %s ", paramName, strVal);
            str2lower(paramName);
            if (0 == strcmp(paramName, "predator")) {
                c->predator = atoi(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "pfactor")) {
                c->pFactor = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "amax")) {
                c->amax = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "vmax")) {
                c->vmax = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "dt")) {
                c->dt = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "ct")) {
                c->ct = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "dmin")) {
                c->dmin = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "predavd")) {
                c->predAvD = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "alpha")) {
                c->alpha = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "beta")) {
                c->beta = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "knn")) {
                param->knn = atoi(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "rc")) {
                param->rc = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "rs")) {
                param->rs = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "rv")) {
                param->rv = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "wc")) {
                param->wc = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "ws")) {
                param->ws = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "wv")) {
                param->wv = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "wu")) {
                param->wu = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "wd")) {
                param->wd = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "mc")) {
                param->mc = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "ms")) {
                param->ms = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "predat")) {
                param->predAt = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "predav")) {
                param->predAv = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "wt")) {
                param->wt = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "wspc")) {
                param->wspc = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "wob")) {
                param->wob = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "learningrate")) {
                options->learningRate = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "delta")) {
                options->delta = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "precision")) {
                options->precision = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "maxiters")) {
                options->maxIters = atoi(strVal); 
                continue;
            }
            if ((0 == strcmp(paramName, "h")) || (0 == strcmp(paramName, "horizon"))) {
                *h = atoi(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "steps")) {
                *steps = atoi(strVal); 
                continue;
            }            
        }
        fclose(fp);
        return 0;
    }
    return -1;
}

int readReynoldsConfig(char *filename, Constraint *c, ReynoldsParameters *param, int *steps)
{
    FILE *fp = fopen(filename, "r");
    if (fp)
    {
        char line[500];
        char *trimmedLine = NULL;
        char paramName[100];
        char strVal[200];
        //double dblVal;
        //size_t len = 0;
        //ssize_t read;
        
        while (fgets(line, 500, fp) != NULL) {
            trimmedLine = strtrim(line);
            if (trimmedLine[0] == '#') // Comment line
                continue;
            sscanf(trimmedLine, "%s %s ", paramName, strVal);
            str2lower(paramName);
            if (0 == strcmp(paramName, "amax")) {
                c->amax = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "vmax")) {
                c->vmax = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "dt")) {
                c->dt = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "rc")) {
                param->rc = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "rs")) {
                param->rs = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "rv")) {
                param->rv = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "wc")) {
                param->wc = atof(strVal); 
                continue;
            }
            if (0 == strcmp(paramName, "ws")) {
                param->ws = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "wv")) {
                param->wv = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "wu")) {
                param->wu = atof(strVal);
                continue;
            }
            if (0 == strcmp(paramName, "steps")) {
                *steps = atoi(strVal); 
                continue;
            }
        }
        fclose(fp);
        return 0;
    }
    return -1;
}

int readNumBirds(char *filename)
{
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp) {
        int n;
        // First line is the number of birds
        fscanf(fp, "%d ", &n);
        fclose(fp);
        return n;
    }
    return -1;
}

int readConf(char *filename, double *x, double *y, double *vx, double *vy, int numBirds)
{
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp) {
        int n;
        // First line is the number of birds
        fscanf(fp, "%d ", &n);
        n = numBirds;
        
        float t1,t2,t3,t4;
        
        int i = 0;
        while (fscanf(fp, "%f %f %f %f ", &t1, &t2, &t3, &t4) != EOF) {
            //sscanf(s, "%f %f %f %f", &t1, &t2, &t3, &t4);
            //printf("%f %f %f %f\n", t1, t2, t3, t4);
            x[i] = t1;
            y[i] = t2;
            vx[i] = t3;
            vy[i] = t4;
            i++;
            if (i >= n) {
                break;
            }
        }
        fclose(fp);
        return 0;
    }
    return -1;
}

int readConfState(char *filename, State *s)
{
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp) {
        int n;
        // First line is the number of birds
        fscanf(fp, "%d ", &n);
        n = s->numBirds;
        
        float t1,t2,t3,t4;
        
        
        int i = 0;
        while (fscanf(fp, "%f %f %f %f ", &t1, &t2, &t3, &t4) != EOF) {
            //sscanf(s, "%f %f %f %f", &t1, &t2, &t3, &t4);
            //printf("%f %f %f %f\n", t1, t2, t3, t4);
            s->x[i] = t1;
            s->y[i] = t2;
            s->vx[i] = t3;
            s->vy[i] = t4;
            s->policy[i] = 0;
            i++;
            if (i >= n) {
                break;
            }
        }
        fclose(fp);
        return 0;
    }
    return -1;
}

int readConfTarget(char *filename, pointObstacles *T)
{
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp) {
        int n;
        // First line is the number of birds
        fscanf(fp, "%d ", &n);
        // n = T->n;
        float t1,t2;

        int i = 0;
        while (fscanf(fp, "%f %f ", &t1, &t2) != EOF) {
            //sscanf(s, "%f %f %f %f", &t1, &t2, &t3, &t4);
            //printf("%f %f %f %f\n", t1, t2, t3, t4);
            T->px[i] = t1;
            T->py[i] = t2;
            i++;
            if (i >= n) {
                break;
            }
        }
        fclose(fp);
        return 0;
    }
    return -1;
}

int writeConf(char *filename, double *x, double *y, double *vx, double *vy, int numBirds)
{
    FILE *fp;
    fp = fopen(filename, "w+");
    if (fp) {
        fprintf(fp, "%d\n", numBirds);
        for (int i = 0; i < numBirds; i++) {
            fprintf(fp, "%f %f %f %f\n", x[i], y[i], vx[i], vy[i]);
        }
        fclose(fp);
        return 0;
    }
    return -1;
}

int writeConfPoints(char *filename, pointObstacles *O, int numBirds)
{
    FILE *fp;
    fp = fopen(filename, "w+");
    if (fp) {
        fprintf(fp, "%d\n", numBirds);
        for (int i = 0; i < numBirds; i++) {
            fprintf(fp, "%f %f\n", O->px[i], O->py[i]);
        }
        fclose(fp);
        return 0;
    }
    return -1;
}