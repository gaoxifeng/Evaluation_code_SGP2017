#pragma once
#include "global_types.h"
#include "Eigen/Dense"
#include <omp.h>
using namespace Eigen;


void initialization_parameters();

void initializeVectorT(vector<int> &Ts,int t,int n);
void initializeVectorT(vector<bool> &Ts,bool t,int n);
void initializeVectorT(vector<double> &Ts,double t,int n);

bool insideVectorT(vector<int> Ts,int t);
//basic functions
void set_exclusion(vector<int> large_set,vector<int> small_set,vector<int> &result_set);
int set_contain(vector<int> large_set,int elm);
bool set_contain(vector<int> large_set,vector<int> small_set);
bool set_contain(vector<int> large_set,vector<int> small_set,int num);//contain num elements of small_set
void set_cross(vector<int> set1,vector<int> set2,vector<int> &result_set);
void set_cross_ids(vector<int> set1,vector<int> set2,vector<int> &result_set);
void set_redundent_clearn(vector<int> &set);
//math
double matrix_determinant(vector<double> &C1,vector<double> &C2,vector<double> &C3,vector<double> &C4);
double sign_of_value(double x);
//vector
void append_vector(vector<int> &v1,vector<int> v2);
void interpolation_vector(vector<double> v1,vector<double> v2,vector<double> &v12,double w);
void reverse_vector(vector<int> &vs);
//geometry
void point_line_projection(vector<double> v1,vector<double> v2,vector<double> v,double &t);
void point_line_projection(double v1[3],double v2[3],double v[3],double &t);
void triangle_coordinates(vector<double> v1,vector<double> v2,vector<double> v3,vector<double> v,vector<double> &ws);
void triangle_coordinates(double v1[3],double v2[3],double v3[3],double v[3],vector<double> &ws);
void triangle_area(vector<double> v1,vector<double> v2,vector<double> v3,double &area);
void triangle_area(double v1[3],double v2[3],double v3[3],double &area);
void segment_sphere_intersection(double center[3],double v2[3],double v3[],double radius);
double cal_volume_Tet(double v0[3],double v1[3],double v2[3],double v3[3]);

