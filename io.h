#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "global_types.h"

using namespace std;
class h_io
{
public:
	
public:
	h_io(void);
	

	void read_hex_mesh_off(	vector<Vertex> &Vs, vector<Hex> &Hexs,char * fname);
	void write_hex_mesh_off(vector<Vertex> &Vs, vector<Hex> &Hexs,char * fname);
	void write_hex_mesh_metrics(vector<double> metrics,char * fname);

	void write_file_list(vector<string> path,char *fname);
	void write_hex_mesh_jacobianList(vector<double> minJs,vector<double> AveJs,char *fname);
	~h_io(void);
};

