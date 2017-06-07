#include "dataset_generation.h"


dataset_generation::dataset_generation(void)
{
	minimum_metrics.push_back(0);//diagonal
	minimum_metrics.push_back(4);//distortion
	minimum_metrics.push_back(8);//jacobian
	minimum_metrics.push_back(18);//relative size squared
	minimum_metrics.push_back(20);//scaled jacobian
	minimum_metrics.push_back(22);//shape
	minimum_metrics.push_back(24);//shape size
	minimum_metrics.push_back(26);//shear
	minimum_metrics.push_back(28);//shear size
	minimum_metrics.push_back(32);//stretch
	minimum_metrics.push_back(36);//volume

	average_metrics.push_back(1);//diagonal
	//average_metrics.push_back(3);//dimension
	average_metrics.push_back(5);//distortion
	average_metrics.push_back(9);//jacobian
	average_metrics.push_back(19);//relative size squared
	average_metrics.push_back(21);//scaled jacobian
	average_metrics.push_back(23);//shape
	average_metrics.push_back(25);//shape size
	average_metrics.push_back(27);//shear
	average_metrics.push_back(29);//shear size
	average_metrics.push_back(33);//stretch
	average_metrics.push_back(37);//volume

	percentage=1;
	noise_ratio=0;
}


void dataset_generation::minimum_metrics_generation(char *path_sub)
{
	hm.measurements(H_Vs,H_Es,H_Fs,H_Hs);
	hm.printout_metrics();


	bool Only_varying_Single_V=true;
	if(Only_varying_Single_V)
	{
		vector<vector<float>> scaled_jacobians;
		for(int i=0;i<H_Hs.size();i++)
		{
			vector<vector<double>> coordinates;
			for(int j = 0; j < 8; j++)
			{
				vector<double> coord;
				coord.push_back(H_Vs[H_Hs[i].vid[j]].v[0]);	
				coord.push_back(H_Vs[H_Hs[i].vid[j]].v[1]);	
				coord.push_back(H_Vs[H_Hs[i].vid[j]].v[2]);	
				coordinates.push_back(coord);
			}
			double jacobi, min_jacobi =1, lengths;
			double len1_sq, len2_sq, len3_sq; 
			vector<double> xxi, xet, xze;

			hm.calc_hex_efg(1, xxi, coordinates);
			hm.calc_hex_efg(2, xet, coordinates);
			hm.calc_hex_efg(3, xze, coordinates);

			jacobi =  hm.dot_cross_vectors(xxi,xet,xze);

			len1_sq = SQUARED_LENGTH(xxi);
			len2_sq = SQUARED_LENGTH(xet);
			len3_sq = SQUARED_LENGTH(xet);

			lengths = sqrt( len1_sq * len2_sq * len3_sq );
			min_jacobi = jacobi / lengths;
			vector<float> s_jacobians;
			for(int j = 0; j < 8; j++)
			{
				int v0,v1,v2,v3;
				v0=hex_tetra_table[j][0];v1=hex_tetra_table[j][1];
				v2=hex_tetra_table[j][2];v3=hex_tetra_table[j][3];

				jacobi = hm.scaled_jacobian(coordinates[v0],coordinates[v1],coordinates[v2],coordinates[v3]);
				s_jacobians.push_back(jacobi);
			}
			scaled_jacobians.push_back(s_jacobians);
		}
		float max_ave=0;int Start_vid=-1;
		for(int i=0;i<H_Vs.size();i++)
		{
			if(H_Vs[i].isboundary!=1)
			{
				float average_J=0;
				for(int j=0;j<H_Vs[i].neighborh.size();j++)
				{
					int nhid=H_Vs[i].neighborh[j];
					for(int k=0;k<scaled_jacobians[nhid].size();k++)
						average_J+=scaled_jacobians[nhid][k];
				}
				average_J/=8*H_Vs[i].neighborh.size();
				if(average_J>max_ave)
				{
					max_ave=average_J;
					Start_vid=i;
				}
			}
		}
		cout<<"Chosen disturbing Vertex: "<<Start_vid<<" "<<max_ave<<endl;
		for(int i=0;i<H_Vs.size();i++)
		{
			if(i!=Start_vid)
				H_Vs[i].isboundary=1;
		}
	}

	int Total_series=20;
	for(int i=0;i<minimum_metrics.size();i++)
	{
		vector<double> a_metric_range;
		for(int j=0;j<Total_series;j++)
			a_metric_range.push_back(hm.all_metrics[minimum_metrics[i]]*(1-(double)j/(Total_series-1)));
		minimum_metrics_ranges.push_back(a_metric_range);

		vector<vector<int>> a_metric_ids;
		for(int j=0;j<Total_series;j++)
		{
			vector<int> a_slot;
			a_metric_ids.push_back(a_slot);
		}
		minimum_metric_ids.push_back(a_metric_ids);
	}


	double len=0;
	for(int i=0;i<H_Es.size();i++)
	{
		double len_;
		DISTANCE(len_,H_Vs[H_Es[i].vid[0]].v,H_Vs[H_Es[i].vid[1]].v);
		len+=len_;
	}
	len/=H_Es.size();

	double noise_s=0.1*len,noise_s_cur;
	
	for(int i = 0 ;i < minimum_metrics.size(); i++)
	{
		hm.measurements(H_Vs,H_Es,H_Fs,H_Hs);

		noise_s_cur=noise_s;
		cout<<"Metric "<<i<<endl;
		for(int j=0;j<Total_series-1;j++)
		{
			cout<<"!!!!!!!!!!!!!!!!!step "<<j<<endl;
			int before_num = minimum_datasets.size();
			minimum_metric_datagen(minimum_metrics[i], H_Vs,noise_s_cur, minimum_metrics_ranges[i][j], minimum_metrics_ranges[i][j+1],path_sub);
			//minimum_metric_datagen_single(minimum_metrics[i], H_Vs,noise_s_cur, minimum_metrics_ranges[i][j], minimum_metrics_ranges[i][j+1],path_sub);
			cout<<"num samples "<<minimum_datasets.size()-before_num<<endl;
		}
	}
}
double dataset_generation::minimum_metric_datagen(int metric_id, vector<Vertex> Vs,double &noise_s, double range_max, double range_min,char *path_sub)
{
	char fname[300],metricname[300];
	sprintf(fname,"%s%s",path_sub,"/min/mesh_");

	vector<Vertex> Temp_Vs;

	int howmany=0,threshold_how = 1,threshold = 2,iterations=0;
	while(howmany<threshold_how)
	{
		Temp_Vs=Vs;

		iterations++;
		int lower=0,higher=0;

		int Start_vid = rand()%Vs.size();

		for(int i=0;i<Vs.size();i++)
		{
			Vertex vtx=Vs[(i+Start_vid)%Vs.size()],temp_vtx=Vs[(i+Start_vid)%Vs.size()];
			if(vtx.isboundary==1)
				continue;

			bool thisone=false;
			bool no_negative=true;
			double	k,t,d;
			double minvalue;

			vector<double> cur_mini_metrics_values;
			for(int j=0;j<minimum_metrics.size();j++)
			{
				cur_mini_metrics_values.push_back(hm.all_metrics[minimum_metrics[j]]);
			}
			int loop_=0;
			while(true)
			{
				loop_++;
				for(int j=0;j<minimum_metrics.size();j++)
				{
					cur_mini_metrics_values[j]=hm.all_metrics[minimum_metrics[j]];
				}
				minvalue=hm.all_metrics[metric_id];

				k=(1.0*rand()/(RAND_MAX+1.0)-0.5)/5*noise_s;
				t=(1.0*rand()/(RAND_MAX+1.0)-0.5)/5*noise_s;
				d=(1.0*rand()/(RAND_MAX+1.0)-0.5)/5*noise_s;
				vtx.v[0]+=k*10;
				vtx.v[1]+=t*10;
				vtx.v[2]+=d*10;

				Vs[(i+Start_vid)%Vs.size()]=vtx;

				if(Inpolyhedron(vtx,Vs) == 'o')
				{
					Vs[(i+Start_vid)%Vs.size()]=temp_vtx;
					break;
				}

				for(int m=0;m<vtx.neighborh.size();m++)
				{
					int hid=vtx.neighborh[m];

					double cur_value = hm.metric_for_an_element(metric_id,hid,Vs,H_Hs);
					if(cur_value<minvalue)
						minvalue=cur_value;

					for(int n=0;n<minimum_metrics.size();n++)
					{
						cur_value = hm.metric_for_an_element(minimum_metrics[n],hid,Vs,H_Hs);
						if(cur_value<cur_mini_metrics_values[n])
							cur_mini_metrics_values[n]=cur_value;
					}

					if(!hm.acceptable_metrics_for_an_element(Vs,H_Hs,hid))
						no_negative=false;
					
				}
				if(!no_negative)
					break;
				if(range_min<=minvalue&&minvalue<range_max)
				{
					if(minimum_acceptable_this_mesh(cur_mini_metrics_values))
					{
						//printf("metric %f\n",minvalue);
						thisone=true;
						Vs[(i+Start_vid)%Vs.size()]=vtx;
						
						hm.measurements(Vs,H_Es,H_Fs,H_Hs);
						sprintf(metricname,"%s%d%s",fname,minimum_datasets.size(),".metrics");
						io.write_hex_mesh_metrics(hm.all_metrics,metricname);
						
						char path_[300];
						sprintf(path_, "%s%d%s", path_sub, minimum_datasets.size(), ".off");
						io.write_hex_mesh_off(Vs, H_Hs, path_);


						howmany++;
						if(howmany>=threshold_how)
							break;
					}else
						;//cout<<"false"<<endl;
				}else if(minvalue<range_min)
				{
					lower++;
					break;
				}else if(minvalue>range_max)
				{
					higher++;
					if(loop_>threshold)
						break;
					else 
						continue;
				}
				if(howmany>=threshold_how)
					break;
			}
			Vs[(i+Start_vid)%Vs.size()]=temp_vtx;
			
			if(howmany>=threshold_how)
				break;
		}
		if(howmany<threshold_how&&lower<higher)
		{
			noise_s*=5;
			printf("higher noise: %f\n",noise_s);
		}
		else if(howmany<threshold_how&&lower>higher)
		{
			noise_s*=0.2;
			printf("lower noise: %f\n",noise_s);
		}
		if(iterations>=3)
			break;
	}

	return false;
}
double dataset_generation::minimum_metric_datagen_single(int metric_id, vector<Vertex> Vs,double &noise_s, double range_max, double range_min,char *path_sub)
{
	char fname[300],metricname[300];
	sprintf(fname,"%s%s",path_sub,"/min/mesh_");

	vector<Vertex> Temp_Vs;
	
	int Start_vid = rand()%Vs.size();

	int howmany=0,threshold_how = 1,threshold = 2,iterations=0;
	while(howmany<threshold_how)
	{
		Temp_Vs=Vs;

		iterations++;
		int lower=0,higher=0;

		bool found_already=false;
		for(int i=0;i<Vs.size();i++)
		{
			Vertex vtx=Vs[(i+Start_vid)%Vs.size()],temp_vtx=Vs[(i+Start_vid)%Vs.size()];
			if(vtx.isboundary==1)
				continue;
			if(found_already)
				break;
			else
				found_already=true;

			bool thisone=false;
			bool no_negative=true;
			double	k,t,d;
			double minvalue;

			vector<double> cur_mini_metrics_values;
			for(int j=0;j<minimum_metrics.size();j++)
			{
				cur_mini_metrics_values.push_back(hm.all_metrics[minimum_metrics[j]]);
			}
			int loop_=0;
			while(true)
			{
				loop_++;
				for(int j=0;j<minimum_metrics.size();j++)
				{
					cur_mini_metrics_values[j]=hm.all_metrics[minimum_metrics[j]];
				}
				minvalue=hm.all_metrics[metric_id];

				k=(1.0*rand()/(RAND_MAX+1.0)-0.5)/5*noise_s;
				t=(1.0*rand()/(RAND_MAX+1.0)-0.5)/5*noise_s;
				d=(1.0*rand()/(RAND_MAX+1.0)-0.5)/5*noise_s;
				vtx.v[0]+=k*10;
				vtx.v[1]+=t*10;
				vtx.v[2]+=d*10;

				Vs[(i+Start_vid)%Vs.size()]=vtx;

				if(Inpolyhedron(vtx,Vs) == 'o')
				{
					Vs[(i+Start_vid)%Vs.size()]=temp_vtx;
					break;
				}

				for(int m=0;m<vtx.neighborh.size();m++)
				{
					int hid=vtx.neighborh[m];

					double cur_value = hm.metric_for_an_element(metric_id,hid,Vs,H_Hs);
					if(cur_value<minvalue)
						minvalue=cur_value;

					for(int n=0;n<minimum_metrics.size();n++)
					{
						cur_value = hm.metric_for_an_element(minimum_metrics[n],hid,Vs,H_Hs);
						if(cur_value<cur_mini_metrics_values[n])
							cur_mini_metrics_values[n]=cur_value;
					}

					if(!hm.acceptable_metrics_for_an_element(Vs,H_Hs,hid))
						no_negative=false;

				}
				if(!no_negative)
					break;
				if(range_min<=minvalue&&minvalue<range_max)
				{
					if(minimum_acceptable_this_mesh(cur_mini_metrics_values))
					{
						//printf("metric %f\n",minvalue);
						thisone=true;
						Vs[(i+Start_vid)%Vs.size()]=vtx;

						hm.measurements(Vs,H_Es,H_Fs,H_Hs);
						sprintf(metricname,"%s%d%s",fname,minimum_datasets.size(),".metrics");
						io.write_hex_mesh_metrics(hm.all_metrics,metricname);
						char path_[300];
						sprintf(path_, "%s%d%s", path_sub, minimum_datasets.size(), ".off");
						io.write_hex_mesh_off(Vs, H_Hs, path_);


						howmany++;
						if(howmany>=threshold_how)
							break;
					}else
						;//cout<<"false"<<endl;
				}else if(minvalue<range_min)
				{
					lower++;
					break;
				}else if(minvalue>range_max)
				{
					higher++;
					if(loop_>threshold)
						break;
					else 
						continue;
				}
				if(howmany>=threshold_how)
					break;
			}
			Vs[(i+Start_vid)%Vs.size()]=temp_vtx;

			if(howmany>=threshold_how)
				break;
		}
		if(howmany<threshold_how&&lower<higher)
		{
			noise_s*=5;
			printf("higher noise: %f\n",noise_s);
		}
		else if(howmany<threshold_how&&lower>higher)
		{
			noise_s*=0.2;
			printf("lower noise: %f\n",noise_s);
		}
		if(iterations>=3)
			break;
	}

	return false;
}
bool dataset_generation::minimum_acceptable_this_mesh(vector<double> mesh_metric_values)
{
	vector<int> metric_ids;
	for(int i = 0;i < minimum_metrics_ranges.size();i++)
	{
		for(int j=0;j<minimum_metrics_ranges[i].size()-1;j++)
		{
			if(mesh_metric_values[i]<=minimum_metrics_ranges[i][j]&&mesh_metric_values[i]>=minimum_metrics_ranges[i][j+1])
			{	
				metric_ids.push_back(j);
				break;
			}
		}
	}
	for(int i=0;i<metric_ids.size();i++)
	{
		if(!minimum_metric_ids[i][metric_ids[i]].size())
		{
			for(int j=0;j<metric_ids.size();j++)
				minimum_metric_ids[j][metric_ids[j]].push_back(minimum_datasets.size());
			minimum_datasets.push_back(metric_ids);

// 			for(int j=0;j<metric_ids.size();j++)
// 				cout<<" "<<metric_ids[j];
// 			cout<<endl;
			return true;
		}
	}

	bool allthesame=false;
	for(int i =0; i<minimum_metric_ids[1][metric_ids[1]].size();i++)
	{
		int data_id=minimum_metric_ids[1][metric_ids[1]][i];
		bool thesame=true;
		for(int j=0;j<metric_ids.size();j++)
		{
			if(minimum_datasets[data_id][j] != metric_ids[j])
				thesame=false;	
		}
		if(thesame)
			allthesame=true;
	}
	if(!allthesame)
	{
		minimum_datasets.push_back(metric_ids);

// 		for(int j=0;j<metric_ids.size();j++)
// 			cout<<" "<<metric_ids[j];
// 		cout<<endl;
		return true;
	}

	return false;
}

void dataset_generation::average_metrics_generation(char *path_sub)
{
	double len=0;
	for(int i=0;i<H_Es.size();i++)
	{
		double len_;
		DISTANCE(len_,H_Vs[H_Es[i].vid[0]].v,H_Vs[H_Es[i].vid[1]].v);
		len+=len_;
	}
	len/=H_Es.size();
	double noise_s=0.1*len,noise_s_cur;

	int Total_series=20;

	for(int i=0;i<average_metrics.size();i++)
	{
		vector<vector<int>> a_metric_ids;
		for(int j=0;j<Total_series;j++)
		{
			vector<int> a_slot;
			a_metric_ids.push_back(a_slot);
		}
		average_metric_ids.push_back(a_metric_ids);
	}

	hm.measurements(H_Vs,H_Es,H_Fs,H_Hs);
	for(int i=0;i<minimum_metrics.size();i++)
	{
		vector<double> a_metric_range;
		for(int j=0;j<Total_series/2+1;j++)
			a_metric_range.push_back(hm.all_metrics[minimum_metrics[i]]*(1-(double)j/(Total_series/2)));
		minimum_metrics_ranges.push_back(a_metric_range);
	}

	vector<Vertex> Temp_Vs;
	vector<Hex> Temp_Hs;
	char fname[300];

	int start_id=0;
	for(int S=0; S<58;S++)
	{
		cout<<"Mesh "<<S<<endl;

		sprintf(fname,"%s%s%d%s",path_sub,"/min/mesh_",S+1,".off");
		io.read_hex_mesh_off(Temp_Vs,Temp_Hs,fname);
		for(int i=0;i<Temp_Vs.size();i++)
		{
			H_Vs[i].v[0]=Temp_Vs[i].v[0];
			H_Vs[i].v[1]=Temp_Vs[i].v[1];
			H_Vs[i].v[2]=Temp_Vs[i].v[2];
		}

		hm.measurements(H_Vs,H_Es,H_Fs,H_Hs);

		vector<int> which_min_ids;
		for(int i=0;i<minimum_metrics_ranges.size();i++)
		{
			bool found=false;
			for(int j=0;j<minimum_metrics_ranges[i].size()-1;j++)
			{
				if(hm.all_metrics[minimum_metrics[i]]<=minimum_metrics_ranges[i][j]&&
					hm.all_metrics[minimum_metrics[i]]>minimum_metrics_ranges[i][j+1])
				{
					which_min_ids.push_back(j);
					found=true;
					break;
				}
			}
			if(!found)
				which_min_ids.push_back(0);
		}
		//hm.printout_metrics();

		average_metrics_ranges.clear();
		for(int i=0;i<average_metrics.size();i++)
		{
			vector<double> a_metric_range;
			for(int j=0;j<Total_series;j++)
				a_metric_range.push_back(hm.all_metrics[average_metrics[i]]*(1-(double)j/(Total_series-1)));
			average_metrics_ranges.push_back(a_metric_range);
		}

		for(int i = 0 ;i < average_metrics.size(); i++)
		{
			noise_s_cur=noise_s;
			cout<<"Metric "<<i<<endl;
			for(int j=0;j<Total_series-1;j++)
			{
				cout<<"!!!!!!!!!!!!!!!!!step "<<j<<endl;
				int before_num = average_datasets.size();
				average_metric_datagen(average_metrics[i], H_Vs,noise_s_cur, average_metrics_ranges[i][j], average_metrics_ranges[i][j+1],path_sub,start_id,which_min_ids);
				cout<<"num samples "<<average_datasets.size()-before_num<<endl;
			}
		}
	}
}
double dataset_generation::average_metric_datagen(int metric_id, vector<Vertex> Vs,double &noise_s, double range_max, double range_min,char *path_sub,int &start_id,vector<int> which_min_ids)
{
	char fname[300],metricname[300];
	sprintf(fname,"%s%s",path_sub,"/ave/mesh_");

	double cur_mean_value=0;
	int iteration=0,threshold_iter=1,threshold=2;int fail=0,failin=0,failinsame=0,samefail=0;
	int stillvalid=0;

	int Num_b=0;
	for(int j=0;j<Vs.size();j++)
	{
		if(Vs[j].isboundary==1)
			continue;
		Num_b++;
	}
	int Num_cur=0, Start_vid=0,Num_Div=5;
	vector<Vertex> Temp_Vs;
	Temp_Vs=Vs;
	srand (time(NULL));
	Num_cur = Num_b/(rand()%Num_Div+1);
	while(true)
	{
		Start_vid = rand()%Vs.size();
		stillvalid=0;


		double	k,t,d;
		long l=0;
		int counter=0,V_Ind=0;
		while(counter<=Num_cur)
		{
			int cur_id=((V_Ind++)+Start_vid)%Vs.size();
			Vertex vtx=Vs[cur_id],temp_vtx=Vs[cur_id];
			if(vtx.isboundary==1)
				continue;
			counter++;

			srand (time(NULL));

			k=(1.0*rand()/(RAND_MAX+1.0)-0.5)/5*noise_s;
			t=(1.0*rand()/(RAND_MAX+1.0)-0.5)/5*noise_s;
			d=(1.0*rand()/(RAND_MAX+1.0)-0.5)/5*noise_s;
			vtx.v[0]+=k;
			vtx.v[1]+=t;
			vtx.v[2]+=d;

			Vs[cur_id]=vtx;

			if(Inpolyhedron(vtx,Vs) == 'o')
			{
				Vs[cur_id]=temp_vtx;
				continue;
			}

			for(int m=0;m<vtx.neighborh.size();m++)
			{
				int hid=vtx.neighborh[m];

				if(!hm.acceptable_metrics_for_an_element(Vs,H_Hs,hid)||!acceptable_min_element(which_min_ids,Vs,H_Hs,hid))
				{
					Vs[cur_id] = temp_vtx;
					break;
				}
			}
			stillvalid=true;
		}
		if(!stillvalid)
			return false;
		hm.measurements(Vs,H_Es,H_Fs,H_Hs);
		vector<double> cur_ave_metrics_values;
		for(int j=0;j<average_metrics.size();j++)
		{
			cur_ave_metrics_values.push_back(hm.all_metrics[average_metrics[j]]);
		}
		cur_mean_value = hm.all_metrics[metric_id];

		if(range_min<=cur_mean_value&&cur_mean_value<=range_max)
		{
			bool havenot=false;
			#pragma omp critical(dataupdate)
			{
				havenot=average_acceptable_this_mesh(cur_ave_metrics_values);
			}
			if(havenot)
			{
				hm.measurements(Vs,H_Es,H_Fs,H_Hs);
				sprintf(metricname,"%s%d%s",fname,start_id,".metrics");
				io.write_hex_mesh_metrics(hm.all_metrics,metricname);
				
				char path_[300];
				sprintf(path_, "%s%d%s", path_sub, start_id++, ".off");
				io.write_hex_mesh_off(Vs, H_Hs, path_);

				iteration++;
			}else if(samefail>=2)
			{	
				failinsame++;
			}else
				failinsame++;
		}else
			failin++;
		if(failin>threshold*Num_Div||failinsame>threshold*Num_Div)
		{
			failin=0;failinsame=0;
			if(range_min > cur_mean_value)
			{
				noise_s*=0.5;
				fail++;

				Vs=Temp_Vs;
			}
			else if(range_max < cur_mean_value)
			{
				noise_s*=2;
				fail++;
				Temp_Vs=Vs;
			}else 
			{
				fail++;
				noise_s*=2;
				Temp_Vs=Vs;
			}
		}
		if(iteration>=threshold_iter)
			break;
		if(fail>2*threshold)
			break;
	}
	return true;
}
bool dataset_generation::average_acceptable_this_mesh(vector<double> mesh_metric_values)
{
	vector<int> metric_ids;
	for(int i = 0;i < average_metrics_ranges.size();i++)
	{
		for(int j=0;j<average_metrics_ranges[i].size()-1;j++)
		{
			if(mesh_metric_values[i]<=average_metrics_ranges[i][j]&&mesh_metric_values[i]>=average_metrics_ranges[i][j+1])
			{	
				metric_ids.push_back(j);
				break;
			}
		}
	}
	for(int i=0;i<metric_ids.size();i++)
	{
		if(!average_metric_ids[i][metric_ids[i]].size())
		{
			for(int j=0;j<metric_ids.size();j++)
				average_metric_ids[j][metric_ids[j]].push_back(average_datasets.size());
			average_datasets.push_back(metric_ids);

			return true;
		}
	}

	bool allthesame=false;
	for(int i =0; i<average_metric_ids[1][metric_ids[0]].size();i++)
	{
		int data_id=average_metric_ids[1][metric_ids[0]][i];
		bool thesame=true;
		for(int j=0;j<metric_ids.size();j++)
		{
			if(average_datasets[data_id][j] != metric_ids[j])
				thesame=false;	
		}
		if(thesame)
			allthesame=true;
	}
	if(!allthesame)
	{
		average_datasets.push_back(metric_ids);
		return true;
	}

	return false;
}
bool dataset_generation::acceptable_min_element(vector<int> which_min_ids,vector<Vertex> &Vs,vector<Hex> &H_Hs,int hid)
{
	vector<vector<double>> coordinates;
	for(int j = 0; j < 8; j++)
	{
		vector<double> coord;
		coord.push_back(Vs[H_Hs[hid].vid[j]].v[0]);	
		coord.push_back(Vs[H_Hs[hid].vid[j]].v[1]);	
		coord.push_back(Vs[H_Hs[hid].vid[j]].v[2]);	
		coordinates.push_back(coord);
	}
	for(int i=0;i<minimum_metrics.size();i++)
	{
		double mvalue=hm.metric_for_an_element(minimum_metrics[i], hid, Vs, H_Hs);
		if(mvalue<minimum_metrics_ranges[i][which_min_ids[i]+1])
			return false;
	}
	return true;
}
vector<int> dataset_generation::perturbe_same_average(vector<Vertex> Vs,double noise_s)
{
	double	k,t,d;

	for(int j=0;j<Vs.size();j++)
	{
		Vertex vtx=Vs[j],temp_vtx=Vs[j];
		if(vtx.isboundary==1)
			continue;

		if(j%percentage!=0)
			continue;

		k=(1.0*rand()/(RAND_MAX+1.0)-0.5)/5*noise_s;
		t=(1.0*rand()/(RAND_MAX+1.0)-0.5)/5*noise_s;
		d=(1.0*rand()/(RAND_MAX+1.0)-0.5)/5*noise_s;
		vtx.v[0]+=k;
		vtx.v[1]+=t;
		vtx.v[2]+=d;

		Vs[j]=vtx;

		for(int m=0;m<vtx.neighborh.size();m++)
		{
			int hid=vtx.neighborh[m];

			if(!hm.acceptable_metrics_for_an_element(Vs,H_Hs,hid))
			{
				Vs[j] = temp_vtx;
				break;
			}
		}
	}
	hm.measurements(Vs,H_Es,H_Fs,H_Hs);
	
	vector<double> cur_ave_metrics_values;
	for(int j=0;j<average_metrics.size();j++)
	{
		cur_ave_metrics_values.push_back(hm.all_metrics[average_metrics[j]]);
	}

	vector<int> metric_ids;
	for(int m = 0;m < average_metrics_ranges.size();m++)
	{
		for(int n=0;n<average_metrics_ranges[m].size()-1;n++)
		{
			if(cur_ave_metrics_values[m]<=average_metrics_ranges[m][n]&&cur_ave_metrics_values[m]>=average_metrics_ranges[m][n+1])
			{	
				metric_ids.push_back(n);
				break;
			}
		}
	}
	return metric_ids;
}

void dataset_generation::extract_edges()
{
	H_Es.clear();
	vector<Edge> Temp_HEs;
	for(int i=0;i<H_Hs.size();i++)
	{
		Edge e;
		//up
		e.vid[0]=H_Hs[i].vid[0];
		e.vid[1]=H_Hs[i].vid[1];
		Temp_HEs.push_back(e);

		e.vid[0]=H_Hs[i].vid[1];
		e.vid[1]=H_Hs[i].vid[2];
		Temp_HEs.push_back(e);

		e.vid[0]=H_Hs[i].vid[2];
		e.vid[1]=H_Hs[i].vid[3];
		Temp_HEs.push_back(e);

		e.vid[0]=H_Hs[i].vid[3];
		e.vid[1]=H_Hs[i].vid[0];
		Temp_HEs.push_back(e);
		//down
		e.vid[0]=H_Hs[i].vid[4];
		e.vid[1]=H_Hs[i].vid[5];
		Temp_HEs.push_back(e);

		e.vid[0]=H_Hs[i].vid[5];
		e.vid[1]=H_Hs[i].vid[6];
		Temp_HEs.push_back(e);

		e.vid[0]=H_Hs[i].vid[6];
		e.vid[1]=H_Hs[i].vid[7];
		Temp_HEs.push_back(e);

		e.vid[0]=H_Hs[i].vid[7];
		e.vid[1]=H_Hs[i].vid[4];
		Temp_HEs.push_back(e);
		//connect
		e.vid[0]=H_Hs[i].vid[0];
		e.vid[1]=H_Hs[i].vid[4];
		Temp_HEs.push_back(e);

		e.vid[0]=H_Hs[i].vid[1];
		e.vid[1]=H_Hs[i].vid[5];
		Temp_HEs.push_back(e);

		e.vid[0]=H_Hs[i].vid[2];
		e.vid[1]=H_Hs[i].vid[6];
		Temp_HEs.push_back(e);

		e.vid[0]=H_Hs[i].vid[3];
		e.vid[1]=H_Hs[i].vid[7];
		Temp_HEs.push_back(e);
	}

	for(int i=0;i<H_Vs.size();i++)
		H_Vs[i].neighborv.clear();

	int E_N=0;
	for(int i=0;i<Temp_HEs.size();i++)
	{
		bool havesame=false;

		int id1=Temp_HEs[i].vid[0];int id2=Temp_HEs[i].vid[1];

		for(int j=0;j<H_Vs[id1].neighborv.size();j++)
		{
			if(H_Vs[id1].neighborv[j]==id2)
			{
				havesame=true;
				break;
			}
		}
		if(!havesame)
		{
			Temp_HEs[i].index=E_N++;

			for(int j=0;j<H_Vs[Temp_HEs[i].vid[0]].neighborh.size();j++)
			{
				int h1=H_Vs[Temp_HEs[i].vid[0]].neighborh[j];
				for(int k=0;k<H_Vs[Temp_HEs[i].vid[1]].neighborh.size();k++)
				{
					int h2=H_Vs[Temp_HEs[i].vid[1]].neighborh[k];
					if(h1==h2)
						Temp_HEs[i].neighborh.push_back(h1);
				}
			}

			Temp_HEs[i].where_location=-1;
			H_Es.push_back(Temp_HEs[i]);
			H_Vs[id1].neighbore.push_back(Temp_HEs[i].index);
			H_Vs[id1].neighborv.push_back(id2);

			H_Vs[id2].neighbore.push_back(Temp_HEs[i].index);
			H_Vs[id2].neighborv.push_back(id1);
		}
	}

	for(int i=0;i<H_Es.size();i++)
	{
		for(int j=0;j<H_Es[i].neighborh.size();j++)
			H_Hs[H_Es[i].neighborh[j]].neighborE.push_back(i);
	}

}
void dataset_generation::extract_faces()
{
	H_Fs.clear();
	vector<Cube_F> Temp_HFs;
	int F_N=0;
	for(int i=0;i<H_Hs.size();i++)
	{
		Cube_F hf[6];

		hf[0].cv[0]=H_Vs[H_Hs[i].vid[0]];
		hf[0].cv[1]=H_Vs[H_Hs[i].vid[1]];
		hf[0].cv[2]=H_Vs[H_Hs[i].vid[2]];
		hf[0].cv[3]=H_Vs[H_Hs[i].vid[3]];

		hf[1].cv[0]=H_Vs[H_Hs[i].vid[4]];
		hf[1].cv[1]=H_Vs[H_Hs[i].vid[5]];
		hf[1].cv[2]=H_Vs[H_Hs[i].vid[6]];
		hf[1].cv[3]=H_Vs[H_Hs[i].vid[7]];

		hf[2].cv[0]=H_Vs[H_Hs[i].vid[0]];
		hf[2].cv[1]=H_Vs[H_Hs[i].vid[1]];
		hf[2].cv[2]=H_Vs[H_Hs[i].vid[5]];
		hf[2].cv[3]=H_Vs[H_Hs[i].vid[4]];

		hf[3].cv[0]=H_Vs[H_Hs[i].vid[0]];
		hf[3].cv[1]=H_Vs[H_Hs[i].vid[4]];
		hf[3].cv[2]=H_Vs[H_Hs[i].vid[7]];
		hf[3].cv[3]=H_Vs[H_Hs[i].vid[3]];

		hf[4].cv[0]=H_Vs[H_Hs[i].vid[3]];
		hf[4].cv[1]=H_Vs[H_Hs[i].vid[2]];
		hf[4].cv[2]=H_Vs[H_Hs[i].vid[6]];
		hf[4].cv[3]=H_Vs[H_Hs[i].vid[7]];

		hf[5].cv[0]=H_Vs[H_Hs[i].vid[1]];
		hf[5].cv[1]=H_Vs[H_Hs[i].vid[2]];
		hf[5].cv[2]=H_Vs[H_Hs[i].vid[6]];
		hf[5].cv[3]=H_Vs[H_Hs[i].vid[5]];

		for(int j=0;j<6;j++)
		{
			bool have=false;
			for(int m=0;m<4;m++)
			{
				for(int n=0;n<H_Vs[hf[j].cv[m].index].neighborf.size();n++)
				{
					int F_id=H_Vs[hf[j].cv[m].index].neighborf[n];

					bool all_have=true;
					for(int p=0;p<4;p++)
					{
						bool exist_v=false;

						for(int q=0;q<4;q++)
							if(hf[j].cv[p].index==H_Fs[F_id].cv[q].index)
								exist_v=true;
						if(!exist_v)
							all_have=false;
					}
					if(all_have)
					{
						have=true;
						H_Fs[F_id].neighbor_CS.push_back(i);
					}
				}
			}
			if(!have)
			{
				hf[j].index=F_N++;

				for(int k=0;k<4;k++)
				{
					int id1=hf[j].cv[k].index;int id2=hf[j].cv[(k+1)%4].index;
					bool found=false;
					for(int m=0;m<H_Vs[id1].neighbore.size();m++)
					{
						int edge1=H_Vs[id1].neighbore[m];
						for(int n=0;n<H_Vs[id2].neighbore.size();n++)
						{
							int edge2=H_Vs[id2].neighbore[n];
							if(edge1==edge2)
							{
								hf[j].eids[k]=edge1;found=true;
							}
							if(found)
								break;
						}
						if(found)
							break;
					}
				}
				hf[j].is_boundary=-1;

				H_Fs.push_back(hf[j]);
				H_Vs[hf[j].cv[0].index].neighborf.push_back(hf[j].index);
				H_Vs[hf[j].cv[1].index].neighborf.push_back(hf[j].index);
				H_Vs[hf[j].cv[2].index].neighborf.push_back(hf[j].index);
				H_Vs[hf[j].cv[3].index].neighborf.push_back(hf[j].index);

				H_Fs[H_Fs.size()-1].neighbor_CS.push_back(i);
			}
		}
	}


	for(int i=0;i<H_Fs.size();i++)
	{
		for(int j=0;j<H_Fs[i].neighbor_CS.size();j++)
			H_Hs[H_Fs[i].neighbor_CS[j]].neighborF.push_back(i);
	}

	for(int i=0;i<H_Hs.size();i++)
	{
		vector<int> f_ids=H_Hs[i].neighborF;
		H_Hs[i].neighborF.clear();
		for(int j=0;j<f_ids.size();j++)
		{
			bool havesame=false;
			for(int k=j+1;k<f_ids.size();k++)
				if(f_ids[j]==f_ids[k])
					havesame=true;
			if(!havesame)
			{
				H_Hs[i].neighborF.push_back(f_ids[j]);
			}
		}
	}

	for(int i=0;i<H_Vs.size();i++)
		H_Vs[i].isboundary=-1;

	for(int i=0;i<H_Fs.size();i++)
	{
		if(H_Fs[i].neighbor_CS.size()==1)
		{
			for(int j=0;j<4;j++)
			{
				int eid=H_Fs[i].eids[j];
				H_Es[eid].where_location=1;
				int vid=H_Fs[i].cv[j].index;
				H_Vs[vid].isboundary=1;
			}
			H_Fs[i].is_boundary=1;

		}
	}
	for(int i=0;i<H_Es.size();i++)
	{
		H_Es[i].where_location=-1;
		int v1,v2;
		v1=H_Es[i].vid[0];v2=H_Es[i].vid[1];
		if(H_Vs[v1].isboundary==1&&H_Vs[v2].isboundary==1)
			H_Es[i].where_location=1;

		for(int j=0;j<H_Vs[v1].neighborf.size();j++)
		{
			int f1=H_Vs[v1].neighborf[j];
			for(int k=0; k<H_Vs[v2].neighborf.size();k++)
			{
				int f2=H_Vs[v2].neighborf[k];
				if(f1==f2)
				{
					bool already=false;
					for(int m=0;m<H_Es[i].neighborf.size();m++)
						if(H_Es[i].neighborf[m]==f1)
							already=true;
					if(!already)
						H_Es[i].neighborf.push_back(f1);
				}
			}
		}
	}
}
void dataset_generation::extract_bounding_polyhedron()
{
	H_TFs.clear();

	for(int i = 0; i < H_Fs.size();i++)
	{
		T_F tf1,tf2;

		tf1.index = 2*i;
		tf1.v[0] =H_Fs[i].cv[0].index;
		tf1.v[1] =H_Fs[i].cv[1].index;
		tf1.v[2] =H_Fs[i].cv[2].index;

		tf2.index = 2*i+1;
		tf2.v[0] =H_Fs[i].cv[2].index;
		tf2.v[1] =H_Fs[i].cv[3].index;
		tf2.v[2] =H_Fs[i].cv[0].index;

		H_Fs[i].neighbor_TFS.push_back(tf1.index);
		H_Fs[i].neighbor_TFS.push_back(tf2.index);

		H_TFs.push_back(tf1);
		H_TFs.push_back(tf2);
	}

	for(int i = 0 ; i < H_Vs.size(); i++ )
	{
		vector<int> allfs,allvs;
		for(int j = 0; j< H_Vs[i].neighborh.size();j++)
		{
			int nhid= H_Vs[i].neighborh[j];
			for(int k = 0 ; k< H_Hs[nhid].neighborF.size(); k++)
				allfs.push_back(H_Hs[nhid].neighborF[k]);
			for(int k = 0 ; k< 8; k++)
				allvs.push_back(H_Hs[nhid].vid[k]);
		}
		set_redundent_clearn(allfs);
		set_redundent_clearn(allvs);
		vector<int> v_vector, boundingvs;
		v_vector.push_back(i);
		set_exclusion(allvs,v_vector,boundingvs);
		H_Vs[i].one_ringv = boundingvs;

		vector<int> boundingfs;
		set_exclusion(allfs,H_Vs[i].neighborf,boundingfs);

		for(int j = 0; j< boundingfs.size();j++)
		{
			H_Vs[i].boundingTFS.push_back(H_Fs[boundingfs[j]].neighbor_TFS[0]);
			H_Vs[i].boundingTFS.push_back(H_Fs[boundingfs[j]].neighbor_TFS[1]);
		}
	}
}



//point in polyhedron
char dataset_generation::Inpolyhedron(Vertex v, vector<Vertex> &Vs)
{
	bmin.v[0] = bmin.v[1] = bmin.v[2] = MAX_RANGE;
	bmax.v[0] = bmax.v[1] = bmax.v[2] = MIN_RANGE ;
	for(int i=0;i < v.one_ringv.size();i++)
	{
		for(int j =0;j<3;j++)
		{
			if(bmin.v[j]>Vs[v.one_ringv[i]].v[j])
				bmin.v[j]=Vs[v.one_ringv[i]].v[j];
			if(bmax.v[j]<Vs[v.one_ringv[i]].v[j])
				bmax.v[j]=Vs[v.one_ringv[i]].v[j];
		}
	}
	DISTANCE(radius,bmin.v,bmax.v);
	radius *= 2;

	Vertex r,p;
	int f, k = 0, crossings = 0;
	char code ='?';
	for(int i =0;i<3;i++)
	{
		if(v.v[i]>bmax.v[i]||v.v[i]<bmin.v[i])
			return 'o';
	}

	bool degenerate_ray = true;
	while(degenerate_ray)
	{
		degenerate_ray = false;

		crossings = 0;
		RandomRay(r);
		VECTOR_ADD2(r.v,r.v,v.v);

		for(f = 0; f <v.boundingTFS.size();f++)
		{
			code = '?';

			int tfid= v.boundingTFS[f];
			int v1 =H_TFs[tfid].v[0], v2 = H_TFs[tfid].v[1], v3 = H_TFs[tfid].v[2];
			for(int k = 0; k< 3; k++)
			{
				if(v.v[k] < Vs[v1].v[k]&& v.v[k] < Vs[v2].v[k] &&v.v[k] < Vs[v3].v[k]
				&&r.v[k] < Vs[v1].v[k]&& r.v[k] < Vs[v2].v[k] &&r.v[k] < Vs[v3].v[k])
				{	code = '0'; break;}
				else if(v.v[k] > Vs[v1].v[k]&& v.v[k] > Vs[v2].v[k] &&v.v[k] > Vs[v3].v[k]
				&&r.v[k] > Vs[v1].v[k]&& r.v[k] > Vs[v2].v[k] &&r.v[k] > Vs[v3].v[k])
				{	code = '0'; break;}
			}

			if(code =='0')
				code = '0';
			else
				code = SegTriInt(H_TFs[v.boundingTFS[f]], v, r, p,Vs);

			if(code == 'p' || code == 'v' || code == 'e')
			{
				degenerate_ray = true;
				break;
			}else if(code == 'f')
			{
				crossings++;
			}else if(code == 'V' || code == 'E' || code == 'F')
				return code;
			else if (code == '0')
				;
		}
	}

	if(crossings%2 ==1)
		return 'i';
	else
		return 'o';
}
void dataset_generation::RandomRay(Vertex &r)
{
	double theta=1.0*rand()/(RAND_MAX+1.0)*PI;
	double phi=1.0*rand()/(RAND_MAX+1.0)*2*PI;
	r.v[0]= radius*sin(theta)*cos(phi);
	r.v[1]= radius*sin(theta)*sin(phi);
	r.v[2]= radius*cos(theta);
}
char dataset_generation::SegTriInt(T_F T,Vertex q,Vertex r, Vertex p, vector<Vertex> &Vs)
{
	int code,m;
	code = SegPlaneInt (T, q, r, p, &m,Vs);
	if(code == 'q')
		return InTri3D (T,m,q,Vs);
	else if(code =='r')
		return InTri3D(T,m,r,Vs);
	else if(code =='p')
		return 'p';//InPlane(T,m,q,r,p);
	else 
		return SegTriCross(T,q,r,Vs);
}
char dataset_generation::SegTriCross(T_F T,Vertex q, Vertex r, vector<Vertex> &Vs)
{
	double vol0,vol1,vol2;
	vol0 = VolumeSign(q,Vs[T.v[0]],Vs[T.v[1]],r);
	vol1 = VolumeSign(q,Vs[T.v[1]],Vs[T.v[2]],r);
	vol2 = VolumeSign(q,Vs[T.v[2]],Vs[T.v[0]],r);
	if((vol0>0 && vol1>0 && vol2 > 0 )||(vol0<0 && vol1<0 && vol2 <0 ))
		return 'f';
	if((vol0>0 || vol1>0 || vol2 > 0 )&&(vol0<0 || vol1<0 || vol2 <0 ))
		return '0';
	else if(abs(vol0)==0 || abs(vol1)==0 || abs(vol2)==0 )
		return '0';//error
	else if((abs(vol0)==0 && abs(vol1)==0)||
		(abs(vol0)==0 && abs(vol2)==0)||
		(abs(vol1)==0 && abs(vol2)==0))
		return 'v';
	else if(abs(vol0)==0 || abs(vol1)==0|| abs(vol2)==0)
		return 'e';
	else
		return '0';//error
}
char dataset_generation::InTri3D(T_F T, int m, Vertex p, vector<Vertex> &Vs)
{
	int i,j,k;Vertex pp, Tp[3];
	j=0;
	for(i=0;i<3;i++)
	{
		if(i!=m)
		{
			pp.v[j] = p.v[i];
			for( k =0; k<3;k++)
				Tp[k].v[j]=Vs[T.v[k]].v[i];
			j++;
		}
	}
	return InTri2D(Tp,pp);
}
char dataset_generation::InTri2D(Vertex Tp[3], Vertex pp)
{
	double area0,area1,area2;
	area0=AreaSign(pp,Tp[0],Tp[1]);
	area1=AreaSign(pp,Tp[1],Tp[2]);
	area2=AreaSign(pp,Tp[2],Tp[0]);

	if((area0==0&&area1>0&&area2>0)||
		(area1==0&&area0>0&&area2>0)||
		(area2==0&&area0>0&&area1>0))
		return 'E';
	if((area0==0&&area1<0&&area2>0)||
		(area1==0&&area0<0&&area2>0)||
		(area2==0&&area0<0&&area1>0))
		return 'E';
	if((area0>0&&area1>0&&area2>0)||
		(area1<0&&area0<0&&area2<0))
		return 'F';
	if((area0==0&&area1==0&&area2==0))
		return '0';//error
	if((area0==0&&area1==0)||
		(area0==0&&area2==0)||
		(area1==0&&area2==0))
		return 'V';
	else
		return '0';
}
char dataset_generation::SegPlaneInt(T_F T, Vertex q, Vertex r, Vertex p, int *m, vector<Vertex> &Vs)
{
	Vertex N, rq; double D, num, denom, t;
	*m = PlaneCoeff(T, N, &D,Vs);
	num = D - DOTVECTOR3(q.v,N.v);
	VECTOR_MINUS2(rq.v,r.v,q.v);//SubVec, need to check
	denom = DOTVECTOR3(rq.v,N.v);
	if(denom == 0.0)
	{
		if(num == 0.0)
			return 'p';
		else
			return '0';
	}
	else 
		t = num/denom;

	int i; 
	for (i =0; i<3; i++)
		p.v[i] = q.v[i] +t*(r.v[i] - q.v[i]);

	if(t>0.0 && t<1.0)
		return '1';
	else if (num == 0)
		return 'q';
	else if( num == denom)
		return 'r';
	else 
		return '0';
}
int dataset_generation::PlaneCoeff(T_F T, Vertex &N, double *D, vector<Vertex> &Vs)
{
	int i, m=0; double t, biggest =0.0;
	NormalVec(Vs[T.v[0]],Vs[T.v[1]],Vs[T.v[2]],N);
	*D=DOTVECTOR3(Vs[T.v[0]].v,N.v);
	for (i=0;i<3;i++)
	{
		t = fabs(N.v[i]);
		if(t>biggest)
		{
			biggest=t;
			m=i;
		}
	}
	return m;
}
void dataset_generation::NormalVec(Vertex a, Vertex b, Vertex c, Vertex &N)
{
	N.v[0]=(c.v[2]-a.v[2])*(b.v[1]-a.v[1])-(b.v[2]-a.v[2])*(c.v[1]-a.v[1]);
	N.v[1]=(b.v[2]-a.v[2])*(c.v[0]-a.v[0])-(b.v[0]-a.v[0])*(c.v[2]-a.v[2]);
	N.v[2]=(b.v[0]-a.v[0])*(c.v[1]-a.v[1])-(b.v[1]-a.v[1])*(c.v[0]-a.v[0]);
}
int dataset_generation::VolumeSign(Vertex v1,Vertex v2,Vertex v3,Vertex p)
{
	double vol;
	double ax,ay,az,bx,by,bz,cx,cy,cz;
	ax=v1.v[0]-p.v[0];
	ay=v1.v[1]-p.v[1];
	az=v1.v[2]-p.v[2];

	bx=v2.v[0]-p.v[0];
	by=v2.v[1]-p.v[1];
	bz=v2.v[2]-p.v[2];

	cx=v3.v[0]-p.v[0];
	cy=v3.v[1]-p.v[1];
	cz=v3.v[2]-p.v[2];

	vol = ax*(by*cz-bz*cy)+
		ay*(bz*cx-bx*cz)+
		az*(bx*cy-by*cx);
	if(vol>0) return 1;
	else if(vol<0) return -1;
	else
		return 0;
}
int dataset_generation::AreaSign(Vertex v1,Vertex v2,Vertex v3)
{
	double area2;
	area2 = (v2.v[0]-v1.v[0])*(v3.v[1]-v1.v[1]) -
		(v3.v[0]-v1.v[0])*(v2.v[1]-v1.v[1]) ;
	if(area2>0) return 1;
	else if(area2<0) return -1;
	else return 0;
}
dataset_generation::~dataset_generation(void)
{
}
