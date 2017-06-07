#include "io.h"
#include "dataset_generation.h"
h_io io;

dataset_generation *datagen;
void datasets(int which);
void out_metrics();

char path_IH[300];// "D:/data/original_hex_meshes/hanger.off";
char path_OH[300];// "D:/data/sig_data/hanger/"; Note that, there must be two folders under path_OH - "min" and "ave"
char Choices[300];// "MIN";

int main(int argc, char *argv[])
{
	if (argc < 4) { cout << "not enough parameters!" << endl; return 0; }
	sprintf(Choices, "%s", argv[1]);
	sprintf(path_IH, "%s", argv[2]);
	sprintf(path_OH, "%s", argv[3]);

	if(strcmp (Choices,"AVE")==0)
	{
		datasets(0);
	}else if(strcmp (Choices,"MIN")==0)
	{
		datasets(1);
	}else if(strcmp (Choices,"METRIC")==0)
	{
		out_metrics();
	}
	return 0;
}
void datasets(int which)
{
	io.read_hex_mesh_off(H_Vs,H_Hs, path_IH);
	
	datagen=new dataset_generation();
	datagen->extract_edges();
	datagen->extract_faces();
	datagen->extract_bounding_polyhedron();

	if(which==0)
		datagen->average_metrics_generation(path_OH);
	else if(which==1)
		datagen->minimum_metrics_generation(path_OH);
}
void out_metrics()
{
	char filename[300];
	datagen=new dataset_generation();
	hex_metrics h_metric;
	int num=63;
	for(int i=1;i<num;i++)
	{
		cout<<i<<endl;

		char filename[300];
		sprintf(filename,"%s%i%s",path_IH,i,".off");
		io.read_hex_mesh_off(H_Vs,H_Hs,filename);

		datagen->extract_edges();
		datagen->extract_faces();
		h_metric.hex2tetmesh(H_Vs,H_Es,H_Fs,H_Hs);
		h_metric.measurements(H_Vs,H_Es,H_Fs,H_Hs);

		sprintf(filename,"%s%d%s",path_OH,i,".metrics");
		io.write_hex_mesh_metrics(h_metric.all_metrics,filename);
	}
}
