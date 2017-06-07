#include "global_functions.h"
#include "global_types.h"

void initialization_parameters()
{
};
void initializeVectorT(vector<int> &Ts,int t,int n)
{
	Ts.clear();
	for(int i=0;i<n;i++)
		Ts.push_back(t);
};
void initializeVectorT(vector<bool> &Ts,bool t,int n)
{
	Ts.clear();
	for(int i=0;i<n;i++)
		Ts.push_back(t);
};
void initializeVectorT(vector<double> &Ts,double t,int n)
{
	Ts.clear();
	for(int i=0;i<n;i++)
		Ts.push_back(t);
};


bool insideVectorT(vector<int> Ts,int t)
{
	for(int i=0;i<Ts.size();i++)
		if(Ts[i]==t)
			return true;
	return false;
};

//basic functions
void set_exclusion(vector<int> large_set,vector<int> small_set,vector<int> &result_set)
{
	for(int i=0;i<large_set.size();i++)
	{
		bool inside=false;
		for(int j=0;j<small_set.size();j++)
		{
			if(small_set[j]==large_set[i])
			{
				inside=true;
				break;
			}
		}
		if(!inside)
			result_set.push_back(large_set[i]);
	}
}
int set_contain(vector<int> large_set,int elm)
{
	for(int j=0;j<large_set.size();j++)
		if(elm==large_set[j])
			return j;

	return -1;
};
bool set_contain(vector<int> large_set,vector<int> small_set)
{
	for(int i=0;i<small_set.size();i++)
	{
		bool inside=false;
		for(int j=0;j<large_set.size();j++)
		{
			if(small_set[i]==large_set[j])
			{
				inside=true;
				break;
			}
		}
		if(!inside)
			return false;
	}
	return true;
};
bool set_contain(vector<int> large_set,vector<int> small_set,int num)
{
	for(int i=0;i<small_set.size();i++)
	{
		bool inside=false;
		for(int j=0;j<large_set.size();j++)
		{
			if(small_set[i]==large_set[j])
			{
				num--;
				break;
			}
		}
	}
	if(num==0)
		return true;
	return false;
}
void set_cross(vector<int> set1,vector<int> set2,vector<int> &result_set)
{
	result_set.clear();
	for(int i=0;i<set1.size();i++)
	{
		bool inside=false;
		for(int j=0;j<set2.size();j++)
		{
			if(set2[j]==set1[i])
			{
				inside=true;
				break;
			}
		}
		if(inside)
			result_set.push_back(set1[i]);
	}
}
void set_cross_ids(vector<int> set1,vector<int> set2,vector<int> &result_set)
{
	result_set.clear();
	for(int i=0;i<set1.size();i++)
	{
		bool inside=false;
		for(int j=0;j<set2.size();j++)
		{
			if(set2[j]==set1[i])
			{
				inside=true;
				break;
			}
		}
		if(inside)
			result_set.push_back(i);
	}
}
void set_redundent_clearn(vector<int> &set)
{
	vector<int> set_copy;
	for(int i=0;i<set.size();i++)
	{
		bool have=false;
		for(int j=i+1;j<set.size();j++)
			if(set[i]==set[j])
				have=true;
		if(!have)
			set_copy.push_back(set[i]);
	}
	set=set_copy;
}

//math
double matrix_determinant(vector<double> &C1,vector<double> &C2,vector<double> &C3,vector<double> &C4)
{
	Matrix4f M(4,4);
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			switch(i)
			{
			case 0:
				M(i,j)=C1[j];
				break;
			case 1:
				M(i,j)=C2[j];
				break;
			case 2:
				M(i,j)=C3[j];
				break;
			case 3:
				M(i,j)=C4[j];
				break;
			}
		}
	}
	return M.determinant();
}
double sign_of_value(double x)
{
	if(abs(x)<EPS)
		return 0.0;
	if(x<0)
		return -1.0;
	return 1.0;
}
//vector 
void append_vector(vector<int> &v1,vector<int> v2)
{
	for(int i=0;i<v2.size();i++)
		v1.push_back(v2[i]);
}
void interpolation_vector(vector<double> v1,vector<double> v2,vector<double> &v12,double w)
{
	for(int i=0;i<v1.size();i++)
	{
		v12.push_back(v1[i]+w*(v2[i]-v1[i]));
	}
}
void reverse_vector(vector<int> &vs)
{
	vector<int> tempvs=vs;
	vs.clear();
	for(int i=tempvs.size()-1;i>=0;i--)
		vs.push_back(tempvs[i]);
}
//geometry
void point_line_projection(vector<double> v1,vector<double> v2,vector<double> v,double &t)
{
	vector<double> vv1,v21;
	for(int i=0;i<v1.size();i++)
	{
		vv1.push_back(v[i]-v1[i]);
		v21.push_back(v2[i]-v1[i]);
	}
	double nv21_2=v21[0]*v21[0]+v21[1]*v21[1]+v21[2]*v21[2];
	if(abs(nv21_2)>=EPS)
		t=(vv1[0]*v21[0]+vv1[1]*v21[1]+vv1[2]*v21[2])/nv21_2;
	else
		t=-1;
}
void point_line_projection(double v1[3],double v2[3],double v[3],double &t)
{
	vector<double> vv1,v21;
	for(int i=0;i<3;i++)
	{
		vv1.push_back(v[i]-v1[i]);
		v21.push_back(v2[i]-v1[i]);
	}
	double nv21_2=v21[0]*v21[0]+v21[1]*v21[1]+v21[2]*v21[2];
	if(abs(nv21_2)>=EPS)
		t=(vv1[0]*v21[0]+vv1[1]*v21[1]+vv1[2]*v21[2])/nv21_2;
	else
		t=-1;
}
void triangle_coordinates(vector<double> v1,vector<double> v2,vector<double> v3,vector<double> v,vector<double> &ws)
{
	double area,area1,area2,area3;
	triangle_area(v2,v3,v,area1);
	ws.push_back(area1);
	triangle_area(v3,v1,v,area2);
	ws.push_back(area2);
	triangle_area(v1,v2,v,area3);
	ws.push_back(area3);
	area=ws[0]+ws[1]+ws[2];
	if(abs(area)>EPS)
	{
		ws[0]/=area;
		ws[1]/=area;
		ws[2]/=area;
	}else
	{
		ws[0]=ws[1]=ws[2]=1.0/3;
	}
}
void triangle_coordinates(double v1[3],double v2[3],double v3[3],double v[3],vector<double> &ws)
{
	for(int i=0;i<3;i++)
		printf("v1 %f, v2 %f, v3 %f, v %f\n",v1[i],v2[i],v3[i],v[i]);
	double area,area1,area2,area3;
	triangle_area(v2,v3,v,area1);
	ws.push_back(area1);
	triangle_area(v3,v1,v,area2);
	ws.push_back(area2);
	triangle_area(v1,v2,v,area3);
	ws.push_back(area3);

	area=ws[0]+ws[1]+ws[2];
	printf("a1 %f,a2 %f, a3 %f a %f\n",area1,area2,area3,area);
	ws[0]/=area;
	ws[1]/=area;
	ws[2]/=area;
}
void triangle_area(double v1[3],double v2[3],double v3[3],double &area)
{
	double dis1,dis2,dis3;
	DISTANCE(dis1,v1,v2);
	DISTANCE(dis2,v1,v3);
	DISTANCE(dis3,v3,v2);
	double len=(dis1+dis2+dis3)/2;
	area=sqrt(len*(len-dis1)*(len-dis2)*(len-dis3));	
}
void triangle_area(vector<double> v1,vector<double> v2,vector<double> v3,double &area)
{
	double dis1,dis2,dis3;
	DISTANCE(dis1,v1,v2);
	DISTANCE(dis2,v1,v3);
	DISTANCE(dis3,v3,v2);
	double len=(dis1+dis2+dis3)/2;
	area=sqrt(len*(len-dis1)*(len-dis2)*(len-dis3));
}
void segment_sphere_intersection(double center[3],double v2[3],double v3[],double radius)
{
	double x0=v2[0]-center[0];
	double x1=v2[1]-center[1];
	double x2=v2[2]-center[2];
	double a,b,c;//a for x^2, b for x, c is constant
	a=x0*x0+x1*x1+x2*x2;
	b=2*(center[0]*x0+center[1]*x1+center[2]*x2);
	c=center[0]*center[0]+center[0]*center[0]+center[0]*center[0]-radius*radius;
	double r1,r2;
	r1=(-b+sqrt(b*b-4*a*c))/(2*a);
	r2=(-b-sqrt(b*b-4*a*c))/(2*a);

	double v31[3],v32[3];
	v31[0]=center[0]+r1*(v2[0]-center[0]);
	v31[1]=center[1]+r1*(v2[1]-center[1]);
	v31[2]=center[2]+r1*(v2[2]-center[2]);
	v32[0]=center[0]+r2*(v2[0]-center[0]);
	v32[1]=center[1]+r2*(v2[1]-center[1]);
	v32[2]=center[2]+r2*(v2[2]-center[2]);

	double dis1,dis2;
	DISTANCE(dis1,v2,v31);DISTANCE(dis2,v2,v32);
	v3[0]=v31[0];
	v3[1]=v31[1];
	v3[2]=v31[2];
	if(dis1>dis2)
	{
		v3[0]=v32[0];
		v3[1]=v32[1];
		v3[2]=v32[2];
	}
}
double cal_volume_Tet(double v0[3],double v1[3],double v2[3],double v3[3])
{
	double v1v0[3],v2v0[3],v3v0[3];
	for(int i=0;i<3;i++)
	{
		v1v0[i]=v1[i]-v0[i];
		v2v0[i]=v2[i]-v0[i];
		v3v0[i]=v3[i]-v0[i];
	}

	double norm1=sqrt(v1v0[0]*v1v0[0]+v1v0[1]*v1v0[1]+v1v0[2]*v1v0[2]);
	double norm2=sqrt(v2v0[0]*v2v0[0]+v2v0[1]*v2v0[1]+v2v0[2]*v2v0[2]);
	double norm3=sqrt(v3v0[0]*v3v0[0]+v3v0[1]*v3v0[1]+v3v0[2]*v3v0[2]);

	double volume=v1v0[0]*(v2v0[1]*v3v0[2]-v2v0[2]*v3v0[1])-v1v0[1]*(v2v0[0]*v3v0[2]-v2v0[2]*v3v0[0])+v1v0[2]*(v2v0[0]*v3v0[1]-v2v0[1]*v3v0[0]); 
	return volume/6;
}