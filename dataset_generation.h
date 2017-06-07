#pragma once
#include "global_types.h"
#include "hex_metrics.h"
#include "global_functions.h"
#include "io.h"
#include <time.h>       /* time */
class dataset_generation
{
public:
	vector<int> minimum_metrics;//diagonal,distortion,Jacobian,relative size squared,
					//scaled Jacobian,shape,shape size, shear, shear size, stretch, volume
	vector<int> average_metrics;//diagonal,distortion,jacobian,....

	vector<vector<double>> minimum_metrics_ranges;
	vector<vector<double>> average_metrics_ranges;

	vector<vector<int>> minimum_datasets;
	vector<vector<vector<int>>> minimum_metric_ids;
	vector<vector<int>> average_datasets;
	vector<vector<vector<int>>> average_metric_ids;

public:
	hex_metrics hm;
	h_io io;

public:
	dataset_generation(void);

	void minimum_metrics_generation(char *path_sub);
	double minimum_metric_datagen(int metric_id, vector<Vertex> Vs,double &noise_s, double range_min, double range_max,char *path_sub);
	double minimum_metric_datagen_single(int metric_id, vector<Vertex> Vs,double &noise_s, double range_min, double range_max,char *path_sub);
	bool minimum_acceptable_this_mesh(vector<double> mesh_metric_values);

	void average_metrics_generation(char *path_sub);
	double average_metric_datagen(int metric_id, vector<Vertex> Vs,double &noise_s, double range_max, double range_min,char *path_sub,int &start_id,vector<int> which_min_ids);
	bool average_acceptable_this_mesh(vector<double> mesh_metric_values);
	bool acceptable_min_element(vector<int> which_min_ids,vector<Vertex> &Vs,vector<Hex> &H_Hs,int hid);
	int percentage;double noise_ratio;
	vector<int> perturbe_same_average(vector<Vertex> Vs,double noise_s);

	void extract_edges();
	void extract_faces();
	//point in polyhedron
	void extract_bounding_polyhedron();
	Vertex bmin,bmax;double radius;
	char Inpolyhedron(Vertex v, vector<Vertex> &Vs);
	void RandomRay(Vertex &r);
	char SegTriInt(T_F T,Vertex q,Vertex r, Vertex p, vector<Vertex> &Vs);
	char SegTriCross(T_F T,Vertex q, Vertex r, vector<Vertex> &Vs);
	char InTri3D(T_F T, int m, Vertex p, vector<Vertex> &Vs);
	char InTri2D(Vertex Tp[3], Vertex pp);
	char SegPlaneInt(T_F T, Vertex q, Vertex r, Vertex p, int *m, vector<Vertex> &Vs);
	int PlaneCoeff(T_F T, Vertex &N, double *D, vector<Vertex> &Vs);
	void NormalVec(Vertex a, Vertex b, Vertex c, Vertex &N);
	int VolumeSign(Vertex v1,Vertex v2,Vertex v3,Vertex p);
	int AreaSign(Vertex v1,Vertex v2,Vertex v3);
	~dataset_generation(void);
};

