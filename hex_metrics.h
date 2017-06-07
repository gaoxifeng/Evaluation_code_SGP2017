#pragma once
#include "io.h"
#include "global_types.h"
#include <vector>
#include "global_functions.h"
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
class hex_metrics
{
public:
	vector<double> all_metrics;
	double minimum_diagonal;
	double average_diagonal;
	double maximum_dimension;
	double average_dimension;
	double minimum_distortion;
	double average_distortion;
	double maximum_edge_ratio;
	double average_edge_ratio;
	double minimum_jacobian;
	double average_jacobian;
	double maximum_maximum_edge_ratio;
	double average_maximum_edge_ratio;
	double maximum_condition;
	double average_condition;
	double maximum_mean_condition;
	double average_mean_condition;
	double maximum_oddy;
	double average_oddy;
	double minimum_relative_size_squared;
	double average_relative_size_squared;
	double minimum_scaled_jacobian;
	double average_scaled_jacobian;
	double minimum_shape;
	double average_shape;
	double minimum_shape_size;
	double average_shape_size;
	double minimum_shear;
	double average_shear;
	double minimum_shear_size;
	double average_shear_size;
	double maximum_skew;
	double average_skew;
	double minimum_stretch;
	double average_stretch;
	double maximum_taper;
	double average_taper;
	double minimum_volume;
	double average_volume;
	double maximum_volume;
	double maximum_gradient;
	double average_gradient;
	double minimum_jacobian_t1;
	double average_jacobian_t1;
	double minimum_jacobian_t2;
	double average_jacobian_t2;
	double average_volume_24;
public:
	vector<Vertex> Vs;
	vector<Edge> Es;
	vector<Cube_F> Fs;
	vector<Hex> Hs;
	vector<MatrixXd> Jaco_matrices;
	vector<double> Jaco;

	vector<MatrixXd> Jaco_matrices_D;
	vector<double> Jaco_D;
public:
	h_io hio;

public:
	hex_metrics(void);
	void measurements(vector<Vertex> &Vs,vector<Edge> &Es,vector<Cube_F> &Fs,vector<Hex> &Hs);
	void initialization();
	void printout_metrics();

	double metric_for_an_element(int metric_id, int hid,vector<Vertex> &Vs,vector<Hex> &Hs);
	double acceptable_metrics_for_an_element(vector<Vertex> &Vs,vector<Hex> &Hs,int hid);


	double hex_diagonal(vector<vector<double>> coordinates);
	double hex_dimension(vector<vector<double>> coordinates);
	double hex_distortion(vector<Vertex> &Vs,vector<Hex> &Hs,int hid,vector<vector<double>> coordinates);
	double hex_edge_ratio(vector<vector<double>> coordinates);
	double hex_jacobian(vector<vector<double>> coordinates);
	double hex_maximum_edge_ratio(vector<vector<double>> coordinates);
	double hex_condition(vector<vector<double>> coordinates);
	double hex_mean_condition(vector<vector<double>> coordinates);
	double hex_oddy(vector<vector<double>> coordinates);
	double hex_relative_size_squared(vector<vector<double>> coordinates);
	double hex_scaled_jacobian(vector<vector<double>> coordinates);
	double hex_shape(vector<vector<double>> coordinates);
	double hex_shape_and_size(vector<vector<double>> coordinates);
	double hex_shear(vector<vector<double>> coordinates);
	double hex_shear_and_size(vector<vector<double>> coordinates);
	double hex_skew(vector<vector<double>> coordinates);
	double hex_stretch(vector<vector<double>> coordinates);
	double hex_taper(vector<vector<double>> coordinates);
	double hex_volume(vector<vector<double>> coordinates);
	void hex_volume2();

	double hex_gradient_2(vector<Vertex> &Vs,vector<Cube_F> &Fs,vector<Hex> &Hs);
	double hex_gradient(vector<Vertex> &Vs,vector<Edge> &Es,vector<Cube_F> &Fs,vector<Hex> &Hs);
	double hex_jacobian_t1(vector<Vertex> &Vs,vector<Hex> &Hs,int hid, vector<vector<double>> coordinates);
	double hex_jacobian_t2(vector<Vertex> &Vs,vector<Hex> &Hs,int hid, vector<vector<double>> coordinates);

//utilities
	double hex_edge_length(int max_min, vector<vector<double>> coordinates);
	double diag_length(int max_min, vector<vector<double>> coordinates);
	void make_hex_edges(vector<vector<double>> coordinates, vector<vector<double>> &edgevectors);
	void calc_hex_efg( int efg_index, vector<double> &efg, vector<vector<double>> coordinates);
	double jacobian_fun(vector<double> v0,vector<double> v1,vector<double> v2,vector<double> v3);
	double condition_comp(vector<double> xxi, vector<double> xet, vector<double> xze);
	double oddy_comp( vector<double> xxi, vector<double> xet, vector<double> xze);
	int v_hex_get_weight(vector<double> &v1,vector<double> &v2,vector<double> &v3);
	double scaled_jacobian(vector<double> v0,vector<double> v1,vector<double> v2,vector<double> v3);
	double dot_cross_vectors(vector<double> v0,vector<double> v1,vector<double> v2);
	void jacobian_matrices(vector<Vertex> &Vs,vector<Hex> &Hs,int hid);
	void jacobian_matrices_derivative(vector<Vertex> &Vs,vector<Hex> &Hs,int hid);

//translate to tet-mesh (one hex to 24 tets)
	vector<Vertex> TVs;vector<Tet> TTs;
	void hex2tetmesh(vector<Vertex> &Vs,vector<Edge> &Es,vector<Cube_F> &Fs,vector<Hex> &Hs);
	~hex_metrics(void);
};

