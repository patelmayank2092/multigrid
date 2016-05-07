#ifndef SOLVER_H
#define SOLVER_H

#include<iostream>
#include<vector>
#include<cmath>
#include"Grid.h"

class Solver
{
private:
int level;
double ngp_;
int Vcycle;

public:

std::vector<Grid*> lev_Vec;
std::vector<double> u_initial;

void display_u();
void display_u_app(int );
void display_frc(int );
void display_res(int );

double get_n_Vcycle();


Solver(int l_level, int n_Vcycle, std::vector<double> &u);
Solver(int l_level, int n_Vcycle);

void RBGS(int l_level);
void pre_smoothing(int l_level);
void post_smoothing(int l_level);

void residual(int l_level);

void restriction(int l_level);

void prolongation(int l_level);

std::vector<double> get_u_app(int l_level);

void Simulation();

~Solver();
};

#endif


	

