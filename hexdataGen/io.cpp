#include "io.h"

h_io::h_io(void)
{

}


void h_io::read_hex_mesh_off(vector<Vertex> &Vs, vector<Hex> &Hexs,char * fname)
{
	Vs.clear();Hexs.clear();
	char file[300];
	std::fstream f(fname,std::ios::in);
	char s[1024],sread[1024];
	int vnum,tnum;	double x,y,z;
	f.getline(s,1023);
	sscanf(s,"%s",&sread);
	f.getline(s,1023);
	sscanf(s,"%d%d%f",&vnum,&tnum,&x);
	for(int i=0;i<vnum;i++)
	{
		f.getline(s,1023);int temp=-1;
		sscanf(s,"%lf %lf %lf",&x,&y,&z);
		Vertex v;
		v.v[0]=x;v.v[1]=y;v.v[2]=z;
		v.index=i;
		Vs.push_back(v);
	}
	int o,p,q;
	for(int i=0;i<tnum;i++)
	{
		f.getline(s,1023);
		int num,a,b,c,d;
		sscanf(s,"%d %d %d %d %d %d %d %d %d",&num,&vnum,&o,&p,&q,&a,&b,&c,&d);
				
		Hex h;
		h.vid[0]=vnum;
		h.vid[1]=o;
		h.vid[2]=p;
		h.vid[3]=q;
		h.vid[4]=a;
		h.vid[5]=b;
		h.vid[6]=c;
		h.vid[7]=d;
		h.index=i;
		Hexs.push_back(h);

		Vs[vnum].neighborh.push_back(i);
		Vs[o].neighborh.push_back(i);
		Vs[p].neighborh.push_back(i);
		Vs[q].neighborh.push_back(i);
		Vs[a].neighborh.push_back(i);
		Vs[b].neighborh.push_back(i);
		Vs[c].neighborh.push_back(i);
		Vs[d].neighborh.push_back(i);
	}

	f.close();
}
void h_io::write_hex_mesh_off(vector<Vertex> &Vs, vector<Hex> &Hexs,char * fname)
{
	char file[300];
	fstream f(fname,ios::out);

	f<<"OFF"<<endl;
	f<<Vs.size()<<" "<<Hexs.size()<<" "<<0<<endl;
	for(int i=0;i<Vs.size();i++)
		f<<setprecision(10)<<Vs[i].v[0]<<" "<<setprecision(10)<<Vs[i].v[1]<<" "<<setprecision(10)<<Vs[i].v[2]<<endl;

	for(int i=0;i<Hexs.size();i++)
	{
		f<<"10 "<<Hexs[i].vid[0]<<" "<<Hexs[i].vid[1]<<" "<<Hexs[i].vid[2]<<" "<<Hexs[i].vid[3];
		f<<" "<<Hexs[i].vid[4]<<" "<<Hexs[i].vid[5]<<" "<<Hexs[i].vid[6]<<" "<<Hexs[i].vid[7]<<" 0 0"<<endl;
	}

	f.close();
}
void h_io::write_hex_mesh_metrics(vector<double> metrics,char * fname)
{
	char file[300];
	fstream f(fname,ios::out);
	for(int i=0;i<metrics.size();i++)
	{
		f<<metrics[i]<<endl;
	}

	f.close();
}
void h_io::write_file_list(vector<string> path,char *fname)
{
	fstream f(fname,ios::out);

	for(int i=0;i<path.size();i++)
	{
		f<<path[i]<<endl;
	}
	f.close();
}
void h_io::write_hex_mesh_jacobianList(vector<double> minJs,vector<double> AveJs,char *fname)
{
	fstream f(fname,ios::out);

	for(int i=0;i<minJs.size();i++)
	{
		f<<minJs[i]<<" ";
	}
	f<<endl;
	for(int i=0;i<AveJs.size();i++)
	{
		f<<AveJs[i]<<" ";
	}
	f.close();
}
h_io::~h_io(void)
{
}
