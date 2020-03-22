#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#include "include/object.h"
#include <vector> // for 2D vectors
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void RK4(vector<object> &bodies,vector<object> &C,double dt);
void write_in_OBS(object m,object OBS,ofstream& file);
void get_cartesian(double GM,double t,double t0,object &m);
double orbital_period(double GM, double a);
void write_hamiltonian(vector<object> &bodies,double t,ofstream& file);
void update_COM(object &COM,vector<object> &bodies);
void get_keplerian(double GM,object &m,object COM);
double fRand(double fMin, double fMax);
double compute_GM(vector<object> &bodies);
void check_body_mass(object &b);
void get_relative_cartesian(double GM,object &obs,double t,double t0,object &m);

#endif // FUNCTIONS_H_INCLUDED
