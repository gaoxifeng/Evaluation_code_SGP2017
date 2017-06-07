#ifndef GLOBAL_TYPE_HEADER
#define GLOBAL_TYPE_HEADER
#include <vector>
using namespace std;

struct Vertex
{
	int index;
	int isboundary;
	double v[3];
	std::vector<int> neighborv;
	std::vector<int> neighbort;
	std::vector<int> neighborf;
	std::vector<int> neighborh;
	std::vector<int> neighbore;

	std::vector<int> boundingTFS;
	std::vector<int> one_ringv;

};
struct Edge
{
	int vid[2];
	int index;
	int where_location;
	std::vector<int> neighborf;
	std::vector<int> neighborh;
};
struct Cube_F
{
	int index;
	int is_boundary;
	Vertex cv[4];
	int eids[4];
	double center[3];

	std::vector<int> neighbor_ES;
	std::vector<int> neighbor_FS;
	std::vector<int> neighbor_CS;

	std::vector<int> neighbor_TFS;
};
struct Hex
{
	int index;
	int vid[8];
	std::vector<int> neighborE;
	std::vector<int> neighborF;
	std::vector<int> neighborH;

	double center[3];
	double size_gradient[3];
	double total_volume;
	double volumes[8];
};
struct T_F
{
	int index;
	int is_boundary;
	int v[3];

	std::vector<int> neighbor_ES;
	std::vector<int> neighbor_FS;
	std::vector<int> neighbor_CS;
};
struct Tet
{
	int index;
	int is_boundary;
	vector<int> vid;
};

const double EPS = 1.0e-7;
const double PI=3.1415926535898;
const double MIN_RANGE=-1e+20;
const double MAX_RANGE=+1e+20;

const int hex_tetra_table[8][4]=
{
	{0,3,4,1},
	{1,0,5,2},
	{2,1,6,3},
	{3,2,7,0},
	{4,7,5,0},
	{5,4,6,1},
	{6,5,7,2},
	{7,6,4,3},
};

#define DISTANCE(a,b,c)	{a=sqrt(pow((b)[0]-(c)[0],2)+pow((b)[1]-(c)[1],2)+pow((b)[2]-(c)[2],2));}
#define CROSSVECTOR3(a,b,c)       {(a)[0]=(b)[1]*(c)[2]-(b)[2]*(c)[1]; \
	(a)[1]=(b)[2]*(c)[0]-(b)[0]*(c)[2]; \
	(a)[2]=(b)[0]*(c)[1]-(b)[1]*(c)[0];}
#define DOTVECTOR3(b,c)       ((b)[0]*(c)[0]+(b)[1]*(c)[1]+(b)[2]*(c)[2]);//non-commutative
#define NORMALIZE(a,b)       {(a)[0]=(a)[0]/b; \
	(a)[1]=(a)[1]/b; \
	(a)[2]=(a)[2]/b;}

#define SQUARED_LENGTH(a)       ((a)[0]*(a)[0]+	(a)[1]*(a)[1]+(a)[2]*(a)[2]);

#define VECTOR_ADD(a,b)       {(a)[0] = (a)[0]+ (b)[0];\
	(a)[1] = (a)[1]+(b)[1];\
	(a)[2] = (a)[2]+(b)[2];}
#define VECTOR_ADD2(a,b,c)       {(a)[0] = (b)[0]+ (c)[0];\
	(a)[1] = (b)[1]+(c)[1];\
	(a)[2] = (b)[2]+(c)[2];}
#define VECTOR_MINUS(a,b)       {(a)[0] = (a)[0]- (b)[0];\
	(a)[1] = (a)[1]-(b)[1];\
	(a)[2] = (a)[2]-(b)[2];}
#define VECTOR_MINUS2(a,b,c)     {(a)[0] = (b)[0]- (c)[0];\
	(a)[1] = (b)[1]-(c)[1];\
	(a)[2] = (b)[2]-(c)[2];}
extern vector<Vertex> H_Vs; 
extern vector<Edge> H_Es; 
extern vector<Cube_F> H_Fs; 
extern vector<Hex> H_Hs; 

extern vector<T_F> H_TFs; 


extern vector<Vertex> H_Vs_gen; 
extern vector<Hex> H_Hs_gen; 

static int rst[8][3]={
	-1,-1,-1,
	1,-1,-1,
	1,1,-1,
	-1,1,-1,
	-1,-1,1,
	1,-1,1,
	1,1,1,
	-1,1,1,
};
static double RST[8][3]={
	-1,-1,-1,
	1,-1,-1,
	1,1,-1,
	-1,1,-1,
	-1,-1,1,
	1,-1,1,
	1,1,1,
	-1,1,1,
};
#endif