#ifndef _REYNOLDS_H_
#define _REYNOLDS_H_

typedef struct ReynoldsParameters {
    double rc;  // Radius of local neighborhood for cohesion calculation.
    double rs;  // Radius of local neighborhood for separation calculation.
    double rv;  // Radius of local neighborhood for velocity matching (alignment) calculation.
    double wc;  // Weight of cohesion term in the total steering force.
    double ws;  // Weight of separation term in the total steering force.
    double wv;  // Weight of velocity matching (alignment) term in the total steering force.
    double wu;  // Weight of the control input.
} ReynoldsParameters;

bool steer_total_distributed(State *s, Action *total, int i, FitnessParameter *param, Constraint *c);
bool steer_max_avoidance(State *s, Action *total, int i, ReynoldsParameters *param, Constraint *c);
bool steer_max_braking(State *s, Action *total, int i, Constraint *c);
bool steer_separation_reynolds(State *s, Action *total, int i, FitnessParameter *param, Constraint *c);
bool steer_cohesion_reynolds(State *s, Action *total, int i, FitnessParameter *param, Constraint *c);
#endif