#include "hex_metrics.h"


hex_metrics::hex_metrics(void)
{
}
void hex_metrics::measurements(vector<Vertex> &Vs,vector<Edge> &Es,vector<Cube_F> &Fs,vector<Hex> &Hs)
{
	initialization();
	all_metrics.clear();

	for(int i=0;i<Hs.size();i++)
	{
		vector<vector<double>> coordinates;
		for(int j = 0; j < 8; j++)
		{
			vector<double> coord;
			coord.push_back(Vs[Hs[i].vid[j]].v[0]);	
			coord.push_back(Vs[Hs[i].vid[j]].v[1]);	
			coord.push_back(Vs[Hs[i].vid[j]].v[2]);	
			coordinates.push_back(coord);
		}

		double cur_volume = hex_volume(coordinates);
		if(minimum_volume>cur_volume)
			minimum_volume=cur_volume;
		if(maximum_volume<cur_volume)
			maximum_volume=cur_volume;
		average_volume+=cur_volume;
	}
	average_volume/=Hs.size();

	for(int i=0;i<Hs.size();i++)
	{
		vector<vector<double>> coordinates;
		for(int j = 0; j < 8; j++)
		{
			vector<double> coord;
			coord.push_back(Vs[Hs[i].vid[j]].v[0]);	
			coord.push_back(Vs[Hs[i].vid[j]].v[1]);	
			coord.push_back(Vs[Hs[i].vid[j]].v[2]);	
			coordinates.push_back(coord);
		}
		double cur_diagonal = hex_diagonal(coordinates);
		if(minimum_diagonal>cur_diagonal)
			minimum_diagonal=cur_diagonal;
		average_diagonal+=cur_diagonal;

		double cur_dimension = hex_dimension(coordinates);
		if(maximum_dimension<cur_dimension)
			maximum_dimension=cur_dimension;
		average_dimension+=cur_dimension;

		double cur_distortion = hex_distortion(Vs,Hs,i,coordinates);
		if(minimum_distortion>cur_distortion)
			minimum_distortion=cur_distortion;
		average_distortion+=cur_distortion;

		double cur_edge_ratio = hex_edge_ratio(coordinates);
		if(maximum_edge_ratio<cur_edge_ratio)
			maximum_edge_ratio=cur_edge_ratio;
		average_edge_ratio+=cur_edge_ratio;

		double cur_jacobian = hex_jacobian(coordinates);
		if(minimum_jacobian>cur_jacobian)
			minimum_jacobian=cur_jacobian;
		average_jacobian+=cur_jacobian;

		double cur_max_edge_ratio = hex_maximum_edge_ratio(coordinates);
		if(maximum_maximum_edge_ratio<cur_max_edge_ratio)
			maximum_maximum_edge_ratio=cur_max_edge_ratio;
		average_maximum_edge_ratio+=cur_max_edge_ratio;

		double cur_condition = hex_condition(coordinates);
		if(maximum_condition<cur_condition)
			maximum_condition=cur_condition;
		average_condition+=cur_condition;

		double cur_mean_condition = hex_mean_condition(coordinates);
		if(maximum_mean_condition<cur_mean_condition)
			maximum_mean_condition=cur_mean_condition;
		average_mean_condition+=cur_mean_condition;

		double cur_oddy = hex_oddy(coordinates);
		if(maximum_oddy<cur_oddy)
			maximum_oddy=cur_oddy;
		average_oddy+=cur_oddy;

		double cur_relative_size_squared = hex_relative_size_squared(coordinates);
		if(minimum_relative_size_squared >cur_relative_size_squared )
			minimum_relative_size_squared =cur_relative_size_squared ;
		average_relative_size_squared +=cur_relative_size_squared ;

		double cur_scaled_jacobian = hex_scaled_jacobian(coordinates);
		if(minimum_scaled_jacobian>cur_scaled_jacobian)
			minimum_scaled_jacobian=cur_scaled_jacobian;
		average_scaled_jacobian+=cur_scaled_jacobian;

		double cur_shape = hex_shape(coordinates);
		if(minimum_shape>cur_shape)
			minimum_shape=cur_shape;
		average_shape+=cur_shape;

		double cur_shape_and_size = hex_shape_and_size(coordinates);
		if(minimum_shape_size>cur_shape_and_size )
			minimum_shape_size=cur_shape_and_size ;
		average_shape_size+=cur_shape_and_size ;

		double cur_shear = hex_shear(coordinates);
		if(minimum_shear>cur_shear)
			minimum_shear=cur_shear;
		average_shear+=cur_shear;

		double cur_shear_and_size = hex_shear_and_size(coordinates);
		if(minimum_shear_size>cur_shear_and_size )
			minimum_shear_size=cur_shear_and_size ;
		average_shear_size+=cur_shear_and_size ;

		double cur_skew = hex_skew(coordinates);
		if(maximum_skew<cur_skew)
			maximum_skew=cur_skew;
		average_skew+=cur_skew;

		double cur_stretch = hex_stretch(coordinates);
		if(minimum_stretch>cur_stretch)
			minimum_stretch=cur_stretch;
		average_stretch+=cur_stretch;

		double cur_taper = hex_taper(coordinates);
		if(maximum_taper<cur_taper)
			maximum_taper=cur_taper;
		average_taper+=cur_taper;

// 		double cur_gradient = hex_gradient();
// 		if(maximum_gradient>cur_gradient)
// 			maximum_gradient=cur_gradient;
// 		average_gradient+=cur_gradient;

		double cur_jacobian_t1 = hex_jacobian_t1(Vs,Hs,i, coordinates);
		if(minimum_jacobian_t1>cur_jacobian_t1)
			minimum_jacobian_t1=cur_jacobian_t1;
		average_jacobian_t1+=cur_jacobian_t1;

		double cur_jacobian_t2 = hex_jacobian_t2(Vs,Hs,i, coordinates);
		if(minimum_jacobian_t2>cur_jacobian_t2)
			minimum_jacobian_t2=cur_jacobian_t2;
		average_jacobian_t2+=cur_jacobian_t2;
	}
	average_diagonal/=Hs.size(); average_dimension/=Hs.size();
	average_distortion/=Hs.size();average_edge_ratio/=Hs.size();
	average_jacobian/=Hs.size();average_maximum_edge_ratio/=Hs.size();
	average_condition/=Hs.size();average_mean_condition/=Hs.size();
	average_oddy/=Hs.size();average_relative_size_squared/=Hs.size();
	average_scaled_jacobian/=Hs.size();average_shape/=Hs.size();
	average_shape_size/=Hs.size();average_shear/=Hs.size();
	average_shear_size/=Hs.size();average_skew/=Hs.size();
	average_stretch/=Hs.size();average_taper/=Hs.size();
	//average_gradient/=Hs.size();
	average_jacobian_t1/=Hs.size();average_jacobian_t2/=Hs.size();

	//hex_gradient(Vs,Es,Fs,Hs);
	hex_gradient_2(Vs,Fs,Hs);
	hex2tetmesh(H_Vs,H_Es,H_Fs,H_Hs);
	hex_volume2();

	all_metrics.push_back(minimum_diagonal);
	all_metrics.push_back(average_diagonal);
	all_metrics.push_back(maximum_dimension);
	all_metrics.push_back(average_dimension);
	all_metrics.push_back(minimum_distortion);
	all_metrics.push_back(average_distortion);
	all_metrics.push_back(maximum_edge_ratio);
	all_metrics.push_back(average_edge_ratio);
	all_metrics.push_back(minimum_jacobian);
	all_metrics.push_back(average_jacobian);
	all_metrics.push_back(maximum_maximum_edge_ratio);
	all_metrics.push_back(average_maximum_edge_ratio);
	all_metrics.push_back(maximum_condition);
	all_metrics.push_back(average_condition);
	all_metrics.push_back(maximum_mean_condition);
	all_metrics.push_back(average_mean_condition);
	all_metrics.push_back(maximum_oddy);
	all_metrics.push_back(average_oddy);
	all_metrics.push_back(minimum_relative_size_squared);
	all_metrics.push_back(average_relative_size_squared);
	all_metrics.push_back(minimum_scaled_jacobian);
	all_metrics.push_back(average_scaled_jacobian);
	all_metrics.push_back(minimum_shape);
	all_metrics.push_back(average_shape);
	all_metrics.push_back(minimum_shape_size);
	all_metrics.push_back(average_shape_size);
	all_metrics.push_back(minimum_shear);
	all_metrics.push_back(average_shear);
	all_metrics.push_back(minimum_shear_size);
	all_metrics.push_back(average_shear_size);
	all_metrics.push_back(maximum_skew);
	all_metrics.push_back(average_skew);
	all_metrics.push_back(minimum_stretch);
	all_metrics.push_back(average_stretch);
	all_metrics.push_back(maximum_taper);
	all_metrics.push_back(average_taper);
	all_metrics.push_back(minimum_volume);
	all_metrics.push_back(average_volume);
	all_metrics.push_back(maximum_volume);
	all_metrics.push_back(maximum_gradient);
	all_metrics.push_back(average_gradient);
	all_metrics.push_back(minimum_jacobian_t1);
	all_metrics.push_back(average_jacobian_t1);
	all_metrics.push_back(minimum_jacobian_t2);
	all_metrics.push_back(average_jacobian_t2);
	all_metrics.push_back(average_volume_24);
}
void hex_metrics::initialization()
{
	double Min_Value=10000000,Max_value=-10000000;
	minimum_diagonal = Min_Value;
	average_diagonal = 0;
	maximum_dimension = Max_value;
	average_dimension = 0;
	minimum_distortion = Min_Value;
	average_distortion = 0;
	maximum_edge_ratio = Max_value;
	average_edge_ratio = 0;
	minimum_jacobian = Min_Value;
	average_jacobian = 0;
	maximum_maximum_edge_ratio= Max_value;
	average_maximum_edge_ratio= 0;
	maximum_condition = Max_value;
	average_condition = 0;
	maximum_mean_condition= Max_value;
	average_mean_condition= 0;
	maximum_oddy = Max_value;
	average_oddy = 0;
	minimum_relative_size_squared = Min_Value;
	average_relative_size_squared = 0;
	minimum_scaled_jacobian = Min_Value;
	average_scaled_jacobian = 0;
	minimum_shape = Min_Value;
	average_shape = 0;
	minimum_shape_size = Min_Value;
	average_shape_size = 0;
	minimum_shear = Min_Value;
	average_shear = 0;
	minimum_shear_size = Min_Value;
	average_shear_size = 0;
	maximum_skew = Max_value;
	average_skew = 0;
	minimum_stretch = Min_Value;
	average_stretch = 0;
	maximum_taper = Max_value;
	average_taper = 0;
	minimum_volume = MAX_RANGE;
	average_volume = 0;
	maximum_volume = Max_value;
	maximum_gradient = Max_value;
	average_gradient = 0;
	minimum_jacobian_t1 = Min_Value;
	average_jacobian_t1 = 0;
	minimum_jacobian_t2 = Min_Value;
	average_jacobian_t2 = 0;
	average_volume_24 = 0;
}
void hex_metrics::printout_metrics()
{
	cout<<"minimum_diagonal "<< minimum_diagonal<<endl;
	cout<<"average_diagonal "<< average_diagonal<<endl;
	cout<<"maximum_dimension "<<maximum_dimension<<endl;
	cout<<"average_dimension "<<average_dimension<<endl;
	cout<<"minimum_distortion "<<minimum_distortion<<endl;
	cout<<"average_distortion "<<average_distortion<<endl;
	cout<<"maximum_edge_ratio "<<maximum_edge_ratio<<endl;
	cout<<"average_edge_ratio "<<average_edge_ratio<<endl;
	cout<<"minimum_jacobian "<<minimum_jacobian<<endl;
	cout<<"average_jacobian "<<average_jacobian<<endl;
	cout<<"maximum_maximum_edge_ratio "<<maximum_maximum_edge_ratio<<endl;
	cout<<"average_maximum_edge_ratio "<<average_maximum_edge_ratio<<endl;
	cout<<"maximum_condition "<<maximum_condition<<endl;
	cout<<"average_condition "<<average_condition<<endl;
	cout<<"maximum_mean_condition "<<maximum_mean_condition<<endl;
	cout<<"average_mean_condition "<<average_mean_condition<<endl;
	cout<<"maximum_oddy "<<maximum_oddy<<endl;
	cout<<"average_oddy "<<average_oddy<<endl;
	cout<<"minimum_relative_size_squared "<<minimum_relative_size_squared<<endl;
	cout<<"average_relative_size_squared "<<average_relative_size_squared<<endl;
	cout<<"minimum_scaled_jacobian "<<minimum_scaled_jacobian<<endl;
	cout<<"average_scaled_jacobian "<<average_scaled_jacobian<<endl;
	cout<<"minimum_shape "<<minimum_shape<<endl;
	cout<<"average_shape "<<average_shape<<endl;
	cout<<"minimum_shape_size "<<minimum_shape_size<<endl;
	cout<<"average_shape_size "<<average_shape_size<<endl;
	cout<<"minimum_shear "<<minimum_shear<<endl;
	cout<<"average_shear "<<average_shear<<endl;
	cout<<"minimum_shear_size "<<minimum_shear_size<<endl;
	cout<<"average_shear_size "<<average_shear_size<<endl;
	cout<<"maximum_skew "<<maximum_skew<<endl;
	cout<<"average_skew "<<average_skew<<endl;
	cout<<"minimum_stretch "<<minimum_stretch<<endl;
	cout<<"average_stretch "<<average_stretch<<endl;
	cout<<"maximum_taper "<<maximum_taper<<endl;
	cout<<"average_taper "<<average_taper<<endl;
	cout<<"minimum_volume "<<minimum_volume<<endl;
	cout<<"average_volume "<<average_volume<<endl;
	cout<<"maximum_volume "<<maximum_volume<<endl;
	cout<<"maximum_gradient "<<maximum_gradient<<endl;
	cout<<"average_gradient "<<average_gradient<<endl;
	cout<<"minimum_jacobian_t1 "<<minimum_jacobian_t1<<endl;
	cout<<"average_jacobian_t1 "<<average_jacobian_t1<<endl;
	cout<<"minimum_jacobian_t2 "<<minimum_jacobian_t2<<endl;
	cout<<"average_jacobian_t2 "<<average_jacobian_t2<<endl;
	cout<<"average_volume_24 "<<average_volume_24<<endl;
}
double hex_metrics::metric_for_an_element(int metric_id, int hid,vector<Vertex> &Vs,vector<Hex> &Hs)
{
	vector<vector<double>> coordinates;
	for(int j = 0; j < 8; j++)
	{
		vector<double> coord;
		coord.push_back(Vs[Hs[hid].vid[j]].v[0]);	
		coord.push_back(Vs[Hs[hid].vid[j]].v[1]);	
		coord.push_back(Vs[Hs[hid].vid[j]].v[2]);	
		coordinates.push_back(coord);
	}

	switch (metric_id)
	{
	case 0:case 1:
		return hex_diagonal(coordinates);
		break;
	case 2:case 3:
		return hex_dimension(coordinates);
		break;
	case 4:case 5:
		return hex_distortion(Vs,Hs,hid,coordinates);
		break;
	case 6:case 7:
		return hex_edge_ratio(coordinates);
		break;
	case 8:case 9:
		return hex_jacobian(coordinates);
		break;
	case 10:case 11:
		return hex_maximum_edge_ratio(coordinates);
		break;
	case 12:case 13:
		return hex_condition(coordinates);
		break;
	case 14:case 15:
		return hex_mean_condition(coordinates);
		break;
	case 16:case 17:
		return hex_oddy(coordinates);
		break;
	case 18:case 19:
		return hex_relative_size_squared(coordinates);
		break;
	case 20:case 21:
		return hex_scaled_jacobian(coordinates);
		break;
	case 22:case 23:
		return hex_shape(coordinates);
		break;
	case 24:case 25:
		return hex_shape_and_size(coordinates);
		break;
	case 26:case 27:
		return hex_shear(coordinates);
		break;
	case 28:case 29:
		return hex_shear_and_size(coordinates);
		break;
	case 30:case 31:
		return hex_skew(coordinates);
		break;
	case 32:case 33:
		return hex_stretch(coordinates);
		break;
	case 34:case 35:
		return hex_taper(coordinates);
		break;
	case 36:case 37:case 38:
		return hex_volume(coordinates);
		break;
	case 39:case 40:
		return 0;//hex_gradient();
		break;
	case 41:case 42:
		return hex_jacobian_t1(Vs,Hs,hid,coordinates);
		break;
	case 43:case 44:
		return hex_jacobian_t2(Vs,Hs,hid,coordinates);
		break;
	}
}
double hex_metrics::acceptable_metrics_for_an_element(vector<Vertex> &Vs,vector<Hex> &Hs,int hid)
{
	vector<vector<double>> coordinates;
	for(int j = 0; j < 8; j++)
	{
		vector<double> coord;
		coord.push_back(Vs[Hs[hid].vid[j]].v[0]);	
		coord.push_back(Vs[Hs[hid].vid[j]].v[1]);	
		coord.push_back(Vs[Hs[hid].vid[j]].v[2]);	
		coordinates.push_back(coord);
	}
	if(hex_jacobian(coordinates)>=0&&hex_scaled_jacobian(coordinates)>=0)
		return true;
	return false;
}
double hex_metrics::hex_diagonal(vector<vector<double>> coordinates)
{
	double min_diag = diag_length( 0, coordinates ); 
	double max_diag = diag_length( 1, coordinates );

	double diagonal = min_diag;
	if(max_diag > 0)
		return diagonal = diagonal / max_diag;
	else
		cout<<"ERROR at hex_diagonal"<<endl;
}
#define SQR(x) ((x) * (x))
double hex_metrics::hex_dimension(vector<vector<double>> coordinates)
{
	double gradop[9][4];

	double x1 = coordinates[0][0];
	double x2 = coordinates[1][0];
	double x3 = coordinates[2][0];
	double x4 = coordinates[3][0];
	double x5 = coordinates[4][0];
	double x6 = coordinates[5][0];
	double x7 = coordinates[6][0];
	double x8 = coordinates[7][0];

	double y1 = coordinates[0][1];
	double y2 = coordinates[1][1];
	double y3 = coordinates[2][1];
	double y4 = coordinates[3][1];
	double y5 = coordinates[4][1];
	double y6 = coordinates[5][1];
	double y7 = coordinates[6][1];
	double y8 = coordinates[7][1];

	double z1 = coordinates[0][2];
	double z2 = coordinates[1][2];
	double z3 = coordinates[2][2];
	double z4 = coordinates[3][2];
	double z5 = coordinates[4][2];
	double z6 = coordinates[5][2];
	double z7 = coordinates[6][2];
	double z8 = coordinates[7][2];

	double z24 = z2 - z4;
	double z52 = z5 - z2;
	double z45 = z4 - z5;
	gradop[1][1] = ( y2*(z6-z3-z45) + y3*z24 + y4*(z3-z8-z52)
		+ y5*(z8-z6-z24) + y6*z52 + y8*z45 ) / 12.0;

	double z31 = z3 - z1;
	double z63 = z6 - z3;
	double z16 = z1 - z6;
	gradop[2][1] = ( y3*(z7-z4-z16) + y4*z31 + y1*(z4-z5-z63)
		+ y6*(z5-z7-z31) + y7*z63 + y5*z16 ) / 12.0;

	double z42 = z4 - z2;
	double z74 = z7 - z4;
	double z27 = z2 - z7;
	gradop[3][1] = ( y4*(z8-z1-z27) + y1*z42 + y2*(z1-z6-z74)
		+ y7*(z6-z8-z42) + y8*z74 + y6*z27 ) / 12.0;

	double z13 = z1 - z3;
	double z81 = z8 - z1;
	double z38 = z3 - z8;
	gradop[4][1] = ( y1*(z5-z2-z38) + y2*z13 + y3*(z2-z7-z81)
		+ y8*(z7-z5-z13) + y5*z81 + y7*z38 ) / 12.0;

	double z86 = z8 - z6;
	double z18 = z1 - z8;
	double z61 = z6 - z1;
	gradop[5][1] = ( y8*(z4-z7-z61) + y7*z86 + y6*(z7-z2-z18)
		+ y1*(z2-z4-z86) + y4*z18 + y2*z61 ) / 12.0;

	double z57 = z5 - z7;
	double z25 = z2 - z5;
	double z72 = z7 - z2;
	gradop[6][1] = ( y5*(z1-z8-z72) + y8*z57 + y7*(z8-z3-z25)
		+ y2*(z3-z1-z57) + y1*z25 + y3*z72 ) / 12.0;

	double z68 = z6 - z8;
	double z36 = z3 - z6;
	double z83 = z8 - z3;
	gradop[7][1] = ( y6*(z2-z5-z83) + y5*z68 + y8*(z5-z4-z36)
		+ y3*(z4-z2-z68) + y2*z36 + y4*z83 ) / 12.0;

	double z75 = z7 - z5;
	double z47 = z4 - z7;
	double z54 = z5 - z4;
	gradop[8][1] = ( y7*(z3-z6-z54) + y6*z75 + y5*(z6-z1-z47)
		+ y4*(z1-z3-z75) + y3*z47 + y1*z54 ) / 12.0;

	double x24 = x2 - x4;
	double x52 = x5 - x2;
	double x45 = x4 - x5;
	gradop[1][2] = ( z2*(x6-x3-x45) + z3*x24 + z4*(x3-x8-x52)
		+ z5*(x8-x6-x24) + z6*x52 + z8*x45 ) / 12.0;

	double x31 = x3 - x1;
	double x63 = x6 - x3;
	double x16 = x1 - x6;
	gradop[2][2] = ( z3*(x7-x4-x16) + z4*x31 + z1*(x4-x5-x63)
		+ z6*(x5-x7-x31) + z7*x63 + z5*x16 ) / 12.0;

	double x42 = x4 - x2;
	double x74 = x7 - x4;
	double x27 = x2 - x7;
	gradop[3][2] = ( z4*(x8-x1-x27) + z1*x42 + z2*(x1-x6-x74)
		+ z7*(x6-x8-x42) + z8*x74 + z6*x27 ) / 12.0;

	double x13 = x1 - x3;
	double x81 = x8 - x1;
	double x38 = x3 - x8;
	gradop[4][2] = ( z1*(x5-x2-x38) + z2*x13 + z3*(x2-x7-x81)
		+ z8*(x7-x5-x13) + z5*x81 + z7*x38 ) / 12.0;

	double x86 = x8 - x6;
	double x18 = x1 - x8;
	double x61 = x6 - x1;
	gradop[5][2] = ( z8*(x4-x7-x61) + z7*x86 + z6*(x7-x2-x18)
		+ z1*(x2-x4-x86) + z4*x18 + z2*x61 ) / 12.0;

	double x57 = x5 - x7;
	double x25 = x2 - x5;
	double x72 = x7 - x2;
	gradop[6][2] = ( z5*(x1-x8-x72) + z8*x57 + z7*(x8-x3-x25)
		+ z2*(x3-x1-x57) + z1*x25 + z3*x72 ) / 12.0;

	double x68 = x6 - x8;
	double x36 = x3 - x6;
	double x83 = x8 - x3;
	gradop[7][2] = ( z6*(x2-x5-x83) + z5*x68 + z8*(x5-x4-x36)
		+ z3*(x4-x2-x68) + z2*x36 + z4*x83 ) / 12.0;

	double x75 = x7 - x5;
	double x47 = x4 - x7;
	double x54 = x5 - x4;
	gradop[8][2] = ( z7*(x3-x6-x54) + z6*x75 + z5*(x6-x1-x47)
		+ z4*(x1-x3-x75) + z3*x47 + z1*x54 ) / 12.0;

	double y24 = y2 - y4;
	double y52 = y5 - y2;
	double y45 = y4 - y5;
	gradop[1][3] = ( x2*(y6-y3-y45) + x3*y24 + x4*(y3-y8-y52)
		+ x5*(y8-y6-y24) + x6*y52 + x8*y45 ) / 12.0;

	double y31 = y3 - y1;
	double y63 = y6 - y3;
	double y16 = y1 - y6;
	gradop[2][3] = ( x3*(y7-y4-y16) + x4*y31 + x1*(y4-y5-y63)
		+ x6*(y5-y7-y31) + x7*y63 + x5*y16 ) / 12.0;

	double y42 = y4 - y2;
	double y74 = y7 - y4;
	double y27 = y2 - y7;
	gradop[3][3] = ( x4*(y8-y1-y27) + x1*y42 + x2*(y1-y6-y74)
		+ x7*(y6-y8-y42) + x8*y74 + x6*y27 ) / 12.0;

	double y13 = y1 - y3;
	double y81 = y8 - y1;
	double y38 = y3 - y8;
	gradop[4][3] = ( x1*(y5-y2-y38) + x2*y13 + x3*(y2-y7-y81)
		+ x8*(y7-y5-y13) + x5*y81 + x7*y38 ) / 12.0;

	double y86 = y8 - y6;
	double y18 = y1 - y8;
	double y61 = y6 - y1;
	gradop[5][3] = ( x8*(y4-y7-y61) + x7*y86 + x6*(y7-y2-y18)
		+ x1*(y2-y4-y86) + x4*y18 + x2*y61 ) / 12.0;

	double y57 = y5 - y7;
	double y25 = y2 - y5;
	double y72 = y7 - y2;
	gradop[6][3] = ( x5*(y1-y8-y72) + x8*y57 + x7*(y8-y3-y25)
		+ x2*(y3-y1-y57) + x1*y25 + x3*y72 ) / 12.0;

	double y68 = y6 - y8;
	double y36 = y3 - y6;
	double y83 = y8 - y3;
	gradop[7][3] = ( x6*(y2-y5-y83) + x5*y68 + x8*(y5-y4-y36)
		+ x3*(y4-y2-y68) + x2*y36 + x4*y83 ) / 12.0;

	double y75 = y7 - y5;
	double y47 = y4 - y7;
	double y54 = y5 - y4;
	gradop[8][3] = ( x7*(y3-y6-y54) + x6*y75 + x5*(y6-y1-y47)
		+ x4*(y1-y3-y75) + x3*y47 + x1*y54 ) / 12.0;

	//     calculate element volume and characteristic element aspect ratio
	//     (used in time step and hourglass control) - 

	double volume =  coordinates[0][0] * gradop[1][1]
	+ coordinates[1][0] * gradop[2][1]
	+ coordinates[2][0] * gradop[3][1]
	+ coordinates[3][0] * gradop[4][1]
	+ coordinates[4][0] * gradop[5][1]
	+ coordinates[5][0] * gradop[6][1]
	+ coordinates[6][0] * gradop[7][1]
	+ coordinates[7][0] * gradop[8][1];
	double aspect = .5*SQR(volume) /
		( SQR(gradop[1][1]) + SQR(gradop[2][1])
		+ SQR(gradop[3][1]) + SQR(gradop[4][1])
		+ SQR(gradop[5][1]) + SQR(gradop[6][1])
		+ SQR(gradop[7][1]) + SQR(gradop[8][1])
		+ SQR(gradop[1][2]) + SQR(gradop[2][2])
		+ SQR(gradop[3][2]) + SQR(gradop[4][2])
		+ SQR(gradop[5][2]) + SQR(gradop[6][2])
		+ SQR(gradop[7][2]) + SQR(gradop[8][2])
		+ SQR(gradop[1][3]) + SQR(gradop[2][3])
		+ SQR(gradop[3][3]) + SQR(gradop[4][3])
		+ SQR(gradop[5][3]) + SQR(gradop[6][3])
		+ SQR(gradop[7][3]) + SQR(gradop[8][3]) );

	return sqrt(aspect);
}
double hex_metrics::hex_distortion(vector<Vertex> &Vs,vector<Hex> &Hs,int hid,vector<vector<double>> coordinates)
{
	double weight=.577350269189626;
	for(int j=0;j<8;j++)
	{
		RST[j][0]=rst[j][0]/abs(rst[j][0])*weight;
		RST[j][1]=rst[j][1]/abs(rst[j][1])*weight;
		RST[j][2]=rst[j][2]/abs(rst[j][2])*weight;
	}

	jacobian_matrices(Vs,Hs,hid);
	jacobian_matrices_derivative(Vs,Hs,hid);
	double elev=0, min_J=10000;
	for(int j=0;j<Jaco.size();j++)
	{
		double absj=abs(Jaco[j]);
// 		if(min_J>absj)
// 			min_J=absj;
		elev+=absj;
	}
	for(int j=0;j<Jaco_D.size();j++)
	{
		double absj=abs(Jaco_D[j]);
		if(min_J>absj)
			min_J=absj;
	}
	//elev = hex_volume(coordinates);
	double distortion=min_J*8/elev;

	return distortion;
}
double hex_metrics::hex_edge_ratio(vector<vector<double>> coordinates)
{
	vector<vector<double>> edgevectors;
	make_hex_edges(coordinates, edgevectors);

	double a2 = SQUARED_LENGTH(edgevectors[0]);
	double b2 = SQUARED_LENGTH(edgevectors[1]);
	double c2 = SQUARED_LENGTH(edgevectors[2]);
	double d2 = SQUARED_LENGTH(edgevectors[3]);
	double e2 = SQUARED_LENGTH(edgevectors[4]);
	double f2 = SQUARED_LENGTH(edgevectors[5]);
	double g2 = SQUARED_LENGTH(edgevectors[6]);
	double h2 = SQUARED_LENGTH(edgevectors[7]);
	double i2 = SQUARED_LENGTH(edgevectors[8]);
	double j2 = SQUARED_LENGTH(edgevectors[9]);
	double k2 = SQUARED_LENGTH(edgevectors[10]);
	double l2 = SQUARED_LENGTH(edgevectors[11]);

	double mab,mcd,mef,Mab,Mcd,Mef;
	double mgh,mij,mkl,Mgh,Mij,Mkl;

	if ( a2 < b2 )
	{
		mab = a2;
		Mab = b2;
	}
	else // b2 <= a2
	{
		mab = b2;
		Mab = a2;
	}
	if ( c2 < d2 )
	{
		mcd = c2;
		Mcd = d2;
	}
	else // d2 <= c2
	{
		mcd = d2;
		Mcd = c2;
	}
	if ( e2 < f2 )
	{
		mef = e2;
		Mef = f2;
	}
	else // f2 <= e2
	{
		mef = f2;
		Mef = e2;
	}
	if ( g2 < h2 )
	{
		mgh = g2;
		Mgh = h2;
	}
	else // h2 <= g2
	{
		mgh = h2;
		Mgh = g2;
	}
	if ( i2 < j2 )
	{
		mij = i2;
		Mij = j2;
	}
	else // j2 <= i2
	{
		mij = j2;
		Mij = i2;
	}
	if ( k2 < l2 )
	{
		mkl = k2;
		Mkl = l2;
	}
	else // l2 <= k2
	{
		mkl = l2;
		Mkl = k2;
	}

	double m2;
	m2 = mab < mcd ? mab : mcd;
	m2 = m2  < mef ? m2  : mef;
	m2 = m2  < mgh ? m2  : mgh;
	m2 = m2  < mij ? m2  : mij;
	m2 = m2  < mkl ? m2  : mkl;

	double M2;
	M2 = Mab > Mcd ? Mab : Mcd;
	M2 = M2  > Mef ? M2  : Mef;
	M2 = M2  > Mgh ? M2  : Mgh;
	M2 = M2  > Mij ? M2  : Mij;
	M2 = M2  > Mkl ? M2  : Mkl;
	m2 = m2  < mef ? m2  : mef;

	M2 = Mab > Mcd ? Mab : Mcd;
	M2 = M2  > Mef ? M2  : Mef;

	double edge_ratio;
	if(m2 <= 0)
		cout<<"ERROR hex_edge_ratio"<<endl;
	else
		return sqrt( M2 / m2 );
}
double hex_metrics::hex_jacobian(vector<vector<double>> coordinates)
{
	double jacobian = MAX_RANGE;
	double current_jacobian; 
	vector<double> xxi, xet, xze;

	calc_hex_efg(1, xxi, coordinates );
	calc_hex_efg(2, xet, coordinates );
	calc_hex_efg(3, xze, coordinates );

	double crosse[3];
	CROSSVECTOR3(crosse,xet,xze);
	current_jacobian=DOTVECTOR3(xxi,crosse);
	current_jacobian = current_jacobian / 64.0;
	if ( current_jacobian < jacobian ) { jacobian = current_jacobian; }

	for(int j=0;j<8;j++)
	{

		int v0,v1,v2,v3;
		v0=hex_tetra_table[j][0];v1=hex_tetra_table[j][1];
		v2=hex_tetra_table[j][2];v3=hex_tetra_table[j][3];
		
		current_jacobian =jacobian_fun(coordinates[v0],coordinates[v1],coordinates[v2],coordinates[v3]); 
		if ( current_jacobian < jacobian ) 
		{ jacobian = current_jacobian; }
	}
	return jacobian;
}
double hex_metrics::hex_maximum_edge_ratio(vector<vector<double>> coordinates)
{
	double aspect;
	double aspect_12, aspect_13, aspect_23;
	vector<double> efg1,efg2,efg3;
	calc_hex_efg( 1, efg1, coordinates);
	calc_hex_efg( 2, efg2, coordinates);
	calc_hex_efg( 3, efg3, coordinates);

	double mag_efg1 = SQUARED_LENGTH(efg1);
	double mag_efg2 = SQUARED_LENGTH(efg2);
	double mag_efg3 = SQUARED_LENGTH(efg3);
	mag_efg1 =sqrt(mag_efg1);
	mag_efg2 =sqrt(mag_efg2);
	mag_efg3 =sqrt(mag_efg3);

	double max1 = mag_efg1,min1 = mag_efg1,
		max2 = mag_efg1,min2 = mag_efg1,
		max3 = mag_efg2,min3 = mag_efg2;
	if(max1<mag_efg2)
		max1=mag_efg2;
	if(min1>mag_efg2)
		min1=mag_efg2;
	aspect_12 = max1/min1;

	if(max2<mag_efg3)
		max2=mag_efg3;
	if(min2>mag_efg3)
		min2=mag_efg3;
	aspect_13 = max2/min2;

	if(max3<mag_efg3)
		max3=mag_efg3;
	if(min3>mag_efg3)
		min3=mag_efg3;
	aspect_23 = max3/min3;
	aspect =aspect_12;
	if(aspect<aspect_13)
		aspect=aspect_13;
	if(aspect<aspect_23)
		aspect=aspect_23;
	return aspect;
}
double hex_metrics::hex_condition(vector<vector<double>> coordinates)
{
	vector<double> xxi, xet, xze;
	xxi.push_back(0);xxi.push_back(0);xxi.push_back(0);
	xet.push_back(0);xet.push_back(0);xet.push_back(0);
	xze.push_back(0);xze.push_back(0);xze.push_back(0);

	double condition =0;
	for(int i = 0 ;i < 8;i++)
	{
		int v0,v1,v2,v3;
		v0=hex_tetra_table[i][0];v1=hex_tetra_table[i][1];
		v2=hex_tetra_table[i][2];v3=hex_tetra_table[i][3];

		VECTOR_MINUS2(xxi,coordinates[v1],coordinates[v0]);
		VECTOR_MINUS2(xet,coordinates[v2],coordinates[v0]);
		VECTOR_MINUS2(xze,coordinates[v3],coordinates[v0]);

		double current_condition = condition_comp( xxi, xet, xze );
		if(i == 0)
			condition = current_condition;
		else if ( current_condition > condition )
			{ condition = current_condition; }
	}
	if(condition>100)
		;//cout<<"here"<<endl;
	condition /= 3.;

	return condition;
}
double hex_metrics::hex_mean_condition(vector<vector<double>> coordinates)
{
	vector<double> xxi, xet, xze;
	xxi.push_back(0);xxi.push_back(0);xxi.push_back(0);
	xet.push_back(0);xet.push_back(0);xet.push_back(0);
	xze.push_back(0);xze.push_back(0);xze.push_back(0);

	double condition =0;
	for(int i = 0 ;i < 8;i++)
	{
		int v0,v1,v2,v3;
		v0=hex_tetra_table[i][0];v1=hex_tetra_table[i][1];
		v2=hex_tetra_table[i][2];v3=hex_tetra_table[i][3];

		VECTOR_MINUS2(xxi,coordinates[v1],coordinates[v0]);
		VECTOR_MINUS2(xet,coordinates[v2],coordinates[v0]);
		VECTOR_MINUS2(xze,coordinates[v3],coordinates[v0]);

		condition += condition_comp( xxi, xet, xze );	
	}
	if(condition>100)
		;//cout<<"here"<<endl;
	condition /= 24.0;

	return condition;
}
double hex_metrics::hex_oddy(vector<vector<double>> coordinates)
{
	double oddy = 0.0, current_oddy;
	vector<double> xxi, xet, xze;

	calc_hex_efg(1, xxi, coordinates );
	calc_hex_efg(2, xet, coordinates );
	calc_hex_efg(3, xze, coordinates );

	current_oddy = oddy_comp( xxi, xet, xze);
	if ( current_oddy > oddy ) { oddy = current_oddy; }

	for(int i = 0 ;i < 8;i++)
	{
		int v0,v1,v2,v3;
		v0=hex_tetra_table[i][0];v1=hex_tetra_table[i][1];
		v2=hex_tetra_table[i][2];v3=hex_tetra_table[i][3];

		VECTOR_MINUS2(xxi,coordinates[v1],coordinates[v0]);
		VECTOR_MINUS2(xet,coordinates[v2],coordinates[v0]);
		VECTOR_MINUS2(xze,coordinates[v3],coordinates[v0]);

		current_oddy = oddy_comp( xxi, xet, xze );
		 if ( current_oddy > oddy )
			{ oddy = current_oddy; }
	}
	return oddy;
}
double hex_metrics::hex_relative_size_squared(vector<vector<double>> coordinates)
{
	double size = 0;
	double tau; 

	vector<double> xxi, xet, xze;
	xxi.push_back(0);xxi.push_back(0);xxi.push_back(0);
	xet.push_back(0);xet.push_back(0);xet.push_back(0);
	xze.push_back(0);xze.push_back(0);xze.push_back(0);

	double det, det_sum = 0;

	v_hex_get_weight( xxi, xet, xze );

	//This is the average relative size 
	double detw = dot_cross_vectors(xxi,xet,xze);

	for(int i = 0; i < 8; i++)
	{
		int v0,v1,v2,v3;
		v0=hex_tetra_table[i][0];v1=hex_tetra_table[i][1];
		v2=hex_tetra_table[i][2];v3=hex_tetra_table[i][3];

		VECTOR_MINUS2(xxi,coordinates[v1],coordinates[v0]);
		VECTOR_MINUS2(xet,coordinates[v2],coordinates[v0]);
		VECTOR_MINUS2(xze,coordinates[v3],coordinates[v0]);

		det = dot_cross_vectors(xxi,xet,xze);
		det_sum += det;  
	}

	{
		tau = det_sum / ( 8*detw);
		if(tau>1.0/tau)
			tau = 1.0/tau;
		size = tau*tau; 
	}
	return size;
}
double hex_metrics::hex_scaled_jacobian(vector<vector<double>> coordinates)
{
	double jacobi, min_jacobi =1, lengths;
	double len1_sq, len2_sq, len3_sq; 
	vector<double> xxi, xet, xze;

	calc_hex_efg(1, xxi, coordinates);
	calc_hex_efg(2, xet, coordinates);
	calc_hex_efg(3, xze, coordinates);

	jacobi =  dot_cross_vectors(xxi,xet,xze);

	len1_sq = SQUARED_LENGTH(xxi);
	len2_sq = SQUARED_LENGTH(xet);
	len3_sq = SQUARED_LENGTH(xet);

	lengths = sqrt( len1_sq * len2_sq * len3_sq );
	min_jacobi = jacobi / lengths;

	for(int i = 0; i < 8; i++)
	{
		int v0,v1,v2,v3;
		v0=hex_tetra_table[i][0];v1=hex_tetra_table[i][1];
		v2=hex_tetra_table[i][2];v3=hex_tetra_table[i][3];

		jacobi = scaled_jacobian(coordinates[v0],coordinates[v1],coordinates[v2],coordinates[v3]);
		if(min_jacobi>jacobi)
			min_jacobi = jacobi;
	}

	return min_jacobi;
}
double hex_metrics::hex_shape(vector<vector<double>> coordinates)
{
	double det, shape;
	double min_shape = 1.0; 
	static const double two_thirds = 2.0/3.0;

	vector<double> xxi, xet, xze;
	xxi.push_back(0);xxi.push_back(0);xxi.push_back(0);
	xet.push_back(0);xet.push_back(0);xet.push_back(0);
	xze.push_back(0);xze.push_back(0);xze.push_back(0);

	for(int i = 0; i < 8; i++)
	{
		int v0,v1,v2,v3;
		v0=hex_tetra_table[i][0];v1=hex_tetra_table[i][1];
		v2=hex_tetra_table[i][2];v3=hex_tetra_table[i][3];

		VECTOR_MINUS2(xxi,coordinates[v1],coordinates[v0]);
		VECTOR_MINUS2(xet,coordinates[v2],coordinates[v0]);
		VECTOR_MINUS2(xze,coordinates[v3],coordinates[v0]);

		det = dot_cross_vectors(xxi,xet,xze);
		double dot1,dot2,dot3;
		dot1 = DOTVECTOR3(xxi,xxi);
		dot2 = DOTVECTOR3(xet,xet);
		dot3 = DOTVECTOR3(xze,xze);

		shape = 3 * pow( det, two_thirds) / (dot1+dot2+dot3);
		if ( shape < min_shape ) 
			{ min_shape = shape; }
	}
	return min_shape;
}
double hex_metrics::hex_shape_and_size(vector<vector<double>> coordinates)
{
	double size = hex_relative_size_squared(coordinates);
	double shape = hex_shape(coordinates );

	return size * shape;
}
double hex_metrics::hex_shear(vector<vector<double>> coordinates)
{
	double shear;
	double min_shear = 1.0; 
	vector<double> xxi, xet, xze;
	xxi.push_back(0);xxi.push_back(0);xxi.push_back(0);
	xet.push_back(0);xet.push_back(0);xet.push_back(0);
	xze.push_back(0);xze.push_back(0);xze.push_back(0);

	double det, len1_sq, len2_sq, len3_sq, lengths;

	for(int i = 0; i < 8; i++)
	{
		int v0,v1,v2,v3;
		v0=hex_tetra_table[i][0];v1=hex_tetra_table[i][1];
		v2=hex_tetra_table[i][2];v3=hex_tetra_table[i][3];

		VECTOR_MINUS2(xxi,coordinates[v1],coordinates[v0]);
		VECTOR_MINUS2(xet,coordinates[v2],coordinates[v0]);
		VECTOR_MINUS2(xze,coordinates[v3],coordinates[v0]);

		len1_sq = SQUARED_LENGTH(xxi);
		len2_sq = SQUARED_LENGTH(xet);
		len3_sq = SQUARED_LENGTH(xze);

		lengths = sqrt( len1_sq * len2_sq * len3_sq );
		det =  dot_cross_vectors(xxi,xet,xze);

		shear = det / lengths;
		if (shear < min_shear ) 
		{ min_shear =shear; }
	}
	return min_shear;
}
double hex_metrics::hex_shear_and_size(vector<vector<double>> coordinates)
{
	double size = hex_relative_size_squared(coordinates );
	double shear = hex_shear(coordinates );

	return shear * size; 
}
double hex_metrics::hex_skew(vector<vector<double>> coordinates)
{
	double skew_1, skew_2, skew_3;

	vector<double> xxi, xet, xze;

	calc_hex_efg(1, xxi, coordinates );
	calc_hex_efg(2, xet, coordinates );
	calc_hex_efg(3, xze, coordinates );

	double slengh=SQUARED_LENGTH(xxi);
	NORMALIZE(xxi, sqrt(slengh));
	slengh=SQUARED_LENGTH(xet);
	NORMALIZE(xet, sqrt(slengh));
	slengh=SQUARED_LENGTH(xze);
	NORMALIZE(xze, sqrt(slengh));

	skew_1=DOTVECTOR3(xxi,xet);
	skew_2=DOTVECTOR3(xxi,xze);
	skew_3=DOTVECTOR3(xet,xze);

	skew_1=abs(skew_1);
	skew_2=abs(skew_2);
	skew_3=abs(skew_3);

	double max_skew = skew_1;
	if(max_skew<skew_2)
		max_skew=skew_2;
	if(max_skew<skew_3)
		max_skew=skew_3;
	return max_skew;
}
double hex_metrics::hex_stretch(vector<vector<double>> coordinates)
{
	static const double HEX_STRETCH_SCALE_FACTOR = sqrt(3.0);

	double min_edge = hex_edge_length( 0, coordinates );
	double max_diag = diag_length( 1, coordinates );  

	return HEX_STRETCH_SCALE_FACTOR * (min_edge/max_diag);
}
double hex_metrics::hex_taper(vector<vector<double>> coordinates)
{
	vector<double> efg1,efg2,efg3,efg12,efg13,efg23;	
	calc_hex_efg( 1, efg1, coordinates);
	calc_hex_efg( 2, efg2, coordinates);
	calc_hex_efg( 3, efg3, coordinates);
	calc_hex_efg( 12, efg12, coordinates);
	calc_hex_efg( 13, efg13, coordinates);
	calc_hex_efg( 23, efg23, coordinates);

	double sl1,sl2,sl3,sl12,sl13,sl23;
	sl1=SQUARED_LENGTH(efg1);
	sl2=SQUARED_LENGTH(efg2);
	sl3=SQUARED_LENGTH(efg3);
	sl12=SQUARED_LENGTH(efg12);
	sl13=SQUARED_LENGTH(efg13);
	sl23=SQUARED_LENGTH(efg23);

	double taper_1,taper_2,taper_3;
	if(sl1>sl2)
		taper_1 = sqrt(sl12/sl2);
	else
		taper_1 = sqrt(sl12/sl1);
	if(sl1>sl3)
		taper_2 = sqrt(sl13/sl3);
	else
		taper_2 = sqrt(sl13/sl1);
	if(sl2>sl3)
		taper_3 = sqrt(sl23/sl3);
	else
		taper_3 = sqrt(sl23/sl2);

	double taper = taper_1;
	if(taper<taper_2)
		taper=taper_2;
	if(taper<taper_3)
		taper=taper_3;
	return taper;
}
double hex_metrics::hex_volume(vector<vector<double>> coordinates)
{
	vector<double> xxi, xet, xze;

	calc_hex_efg(1, xxi, coordinates );
	calc_hex_efg(2, xet, coordinates );
	calc_hex_efg(3, xze, coordinates );
	
	return dot_cross_vectors(xxi,xet,xze)/64.0;

// 	double volumes[4];
// 	volumes[0]=cal_volume_Tet(coordinates[1],coordinates[0],coordinates[5],coordinates[2]);
// 	volumes[1]=cal_volume_Tet(coordinates[3],coordinates[2],coordinates[0],coordinates[6]);
// 	volumes[2]=cal_volume_Tet(coordinates[7],coordinates[0],coordinates[3],coordinates[6]);
// 	volumes[3]=cal_volume_Tet(coordinates[4],coordinates[7],coordinates[5],coordinates[0]);
// 
// 	if(volumes[0]+volumes[1]+volumes[2]+volumes[3]<0)
// 	{
// 		return -(volumes[0]+volumes[1]+volumes[2]+volumes[3])/6;
// 	}else
// 	{
// 		return  (volumes[0]+volumes[1]+volumes[2]+volumes[3])/6;
// 
}
void hex_metrics::hex_volume2()
{
	for(int i=0;i<TTs.size();i++)
	{
// 		vector<double[3]> coordinates;
// 		for(int j = 0; j < 4; j++)
// 		{
// 			double coord[3];
// 			coord[0]=(Vs[TTs[i].vid[j]].v[0]);	
// 			coord[1]=(Vs[TTs[i].vid[j]].v[1]);	
// 			coord[2]=(Vs[TTs[i].vid[j]].v[2]);	
// 			coordinates.push_back(coord);
// 		}
		average_volume_24+=abs(cal_volume_Tet(TVs[TTs[i].vid[0]].v,TVs[TTs[i].vid[1]].v,TVs[TTs[i].vid[2]].v,TVs[TTs[i].vid[3]].v));
	}
	average_volume_24*=24;
	average_volume_24/=TTs.size();

}
double hex_metrics::hex_gradient_2(vector<Vertex> &Vs,vector<Cube_F> &Fs,vector<Hex> &Hs)
{
	vector<double> volumes;
	for(int i=0;i<Hs.size();i++)
	{
		vector<vector<double>> coordinates;
		for(int j = 0; j < 8; j++)
		{
			vector<double> coord;
			coord.push_back(Vs[Hs[i].vid[j]].v[0]);	
			coord.push_back(Vs[Hs[i].vid[j]].v[1]);	
			coord.push_back(Vs[Hs[i].vid[j]].v[2]);	
			coordinates.push_back(coord);
		}
		volumes.push_back(hex_volume(coordinates));
	}
	vector<double> gradients;
	for(int i=0;i<Hs.size();i++)
	{
		vector<int> f6s=Hs[i].neighborF,fs_left;
		vector<int> x2,y2,z2,v41,v42;
		x2.push_back(f6s[0]);
		for(int j=0;j<4;j++)
			v41.push_back(Fs[x2[0]].cv[j].index);
		for(int j=1;j<f6s.size();j++)
		{
			v42.clear();
			for(int k=0;k<4;k++)
				v42.push_back(Fs[f6s[j]].cv[k].index);
			vector<int> vscross;
			set_cross(v41,v42,vscross);
			if(!vscross.size())
			{
				x2.push_back(f6s[j]);
				break;
			}
		}

		set_exclusion(f6s,x2,fs_left);
		f6s=fs_left;
		y2.push_back(f6s[0]);
		v41.clear();
		for(int j=0;j<4;j++)
			v41.push_back(Fs[y2[0]].cv[j].index);
		for(int j=1;j<f6s.size();j++)
		{
			v42.clear();
			for(int k=0;k<4;k++)
				v42.push_back(Fs[f6s[j]].cv[k].index);
			vector<int> vscross;
			set_cross(v41,v42,vscross);
			if(!vscross.size())
			{
				y2.push_back(f6s[j]);
				break;
			}
		}

		set_exclusion(f6s,y2,fs_left);
		f6s=fs_left;
		z2.push_back(f6s[0]);
		v41.clear();
		for(int j=0;j<4;j++)
			v41.push_back(Fs[z2[0]].cv[j].index);
		for(int j=1;j<f6s.size();j++)
		{
			v42.clear();
			for(int k=0;k<4;k++)
				v42.push_back(Fs[f6s[j]].cv[k].index);
			vector<int> vscross;
			set_cross(v41,v42,vscross);
			if(!vscross.size())
			{
				z2.push_back(f6s[j]);
				break;
			}
		}

		vector<int> xhs,yhs,zhs;
		for(int j=0;j<x2.size();j++)
			for(int k=0;k<Fs[x2[j]].neighbor_CS.size();k++)
				if(Fs[x2[j]].neighbor_CS[k]!=i)
					xhs.push_back(Fs[x2[j]].neighbor_CS[k]);
		if(xhs.size()==1)
			Hs[i].size_gradient[0]=(volumes[xhs[0]]-volumes[i])/2;
		else if(xhs.size()==2)
			Hs[i].size_gradient[0]=(volumes[xhs[0]]-volumes[xhs[1]])/2;
		else
			Hs[i].size_gradient[0]=0;

		for(int j=0;j<y2.size();j++)
			for(int k=0;k<Fs[y2[j]].neighbor_CS.size();k++)
				if(Fs[y2[j]].neighbor_CS[k]!=i)
					yhs.push_back(Fs[y2[j]].neighbor_CS[k]);
		if(yhs.size()==1)
			Hs[i].size_gradient[1]=(volumes[yhs[0]]-volumes[i])/2;
		else if(yhs.size()==2)
			Hs[i].size_gradient[1]=(volumes[yhs[0]]-volumes[yhs[1]])/2;
		else
			Hs[i].size_gradient[1]=0;

		for(int j=0;j<z2.size();j++)
			for(int k=0;k<Fs[z2[j]].neighbor_CS.size();k++)
				if(Fs[z2[j]].neighbor_CS[k]!=i)
					zhs.push_back(Fs[z2[j]].neighbor_CS[k]);
		if(zhs.size()==1)
			Hs[i].size_gradient[2]=(volumes[zhs[0]]-volumes[i])/2;
		else if(zhs.size()==2)
			Hs[i].size_gradient[2]=(volumes[zhs[0]]-volumes[zhs[1]])/2;
		else
			Hs[i].size_gradient[2]=0;


		double scalar=sqrt(Hs[i].size_gradient[0]*Hs[i].size_gradient[0]+
			Hs[i].size_gradient[1]*Hs[i].size_gradient[1]+
			Hs[i].size_gradient[2]*Hs[i].size_gradient[2]);
		gradients.push_back(scalar);
		average_gradient+=scalar;
		if(scalar>maximum_gradient)
		{	
			maximum_gradient=scalar;
			//cout<<i<<endl;
		}
	}
	average_gradient=average_gradient/Hs.size();
	return 1.0;
}
double hex_metrics::hex_gradient(vector<Vertex> &Vs,vector<Edge> &Es,vector<Cube_F> &Fs,vector<Hex> &Hs)
{
	for(int i=0;i<Hs.size();i++)
	{
		Hs[i].center[0]=Hs[i].center[1]=Hs[i].center[2]=0;
		for(int j=0;j<8;j++)
		{
			Hs[i].center[0]+=Vs[Hs[i].vid[j]].v[0];
			Hs[i].center[1]+=Vs[Hs[i].vid[j]].v[1];
			Hs[i].center[2]+=Vs[Hs[i].vid[j]].v[2];
		}
		Hs[i].center[0]/=8;
		Hs[i].center[1]/=8;
		Hs[i].center[2]/=8;
	}
	for(int i=0;i<Hs.size();i++)
	{
		std::vector<int> hs;
		for(int j=0;j<6;j++)
		{
			for(int k=0;k<Fs[Hs[i].neighborF[j]].neighbor_CS.size();k++)
			{
				int h=Fs[Hs[i].neighborF[j]].neighbor_CS[k];
				if(h!=i)
				{
					bool alreadyhave=false;
					for(int m=0;m<hs.size();m++)
					{
						if(hs[m]==h)
							alreadyhave=true;
					}
					if(!alreadyhave)
						hs.push_back(h);
				}
			}
		}

		vector<vector<double>> coordinates;
		for(int k = 0; k < 8; k++)
		{
			vector<double> coord;
			coord.push_back(Vs[Hs[i].vid[k]].v[0]);	
			coord.push_back(Vs[Hs[i].vid[k]].v[1]);	
			coord.push_back(Vs[Hs[i].vid[k]].v[2]);	
			coordinates.push_back(coord);
		}
		double curvolume=hex_volume(coordinates);
		std::vector<double> hvectors[3],scalar_vs;
		for(int j=0;j<hs.size();j++)
		{
			for(int k=0;k<3;k++)
			{
				hvectors[k].push_back(Hs[hs[j]].center[k]-Hs[i].center[k]);
			}

			coordinates.clear();
			for(int k = 0; k < 8; k++)
			{
				vector<double> coord;
				coord.push_back(Vs[Hs[hs[j]].vid[k]].v[0]);	
				coord.push_back(Vs[Hs[hs[j]].vid[k]].v[1]);	
				coord.push_back(Vs[Hs[hs[j]].vid[k]].v[2]);	
				coordinates.push_back(coord);
			}
			double hidvolume = hex_volume(coordinates);
			double minused_volume=pow(hidvolume,double(1.0/3))-pow(curvolume,double(1.0/3));
			scalar_vs.push_back(minused_volume);
		}

		MatrixXf A(hs.size(),3);
		VectorXf b(hs.size()),x(3);
		for(int j=0;j<hs.size();j++)
		{
			A(j,0)=hvectors[0][j];
			A(j,1)=hvectors[1][j];
			A(j,2)=hvectors[2][j];
			b(j)=scalar_vs[j];
		}
		x=A.colPivHouseholderQr().solve(b);
		for(int j=0;j<3;j++)
			Hs[i].size_gradient[j]=x(j);
	}

	vector<double> gradients;
	for(int i=0;i<Hs.size();i++)
	{
		double scalar=sqrt(Hs[i].size_gradient[0]*Hs[i].size_gradient[0]+
			Hs[i].size_gradient[1]*Hs[i].size_gradient[1]+
			Hs[i].size_gradient[2]*Hs[i].size_gradient[2]);
		gradients.push_back(scalar);
		average_gradient+=scalar;
		if(scalar>maximum_gradient)
		{	maximum_gradient=scalar;
		cout<<i<<endl;}
	}
	vector<int> hs_ids;
	vector<Hex> tempHs;
	for(int i=0;i<6;i++)
	{
		for(int j=0;j<Fs[Hs[301].neighborF[i]].neighbor_CS.size();j++)
			hs_ids.push_back(Fs[Hs[301].neighborF[i]].neighbor_CS[j]);
	}
	hs_ids.push_back(301);
	set_redundent_clearn(hs_ids);
	for(int i=0;i<hs_ids.size();i++)
		tempHs.push_back(Hs[hs_ids[i]]);

	average_gradient=average_gradient/Hs.size();
	return 1.0;
}
double hex_metrics::hex_jacobian_t1(vector<Vertex> &Vs,vector<Hex> &Hs,int hid, vector<vector<double>> coordinates)
{
	double weight=.57735;
	for(int j=0;j<8;j++)
	{
		RST[j][0]=rst[j][0]/abs(rst[j][0])*weight;
		RST[j][1]=rst[j][1]/abs(rst[j][1])*weight;
		RST[j][2]=rst[j][2]/abs(rst[j][2])*weight;
	}

	jacobian_matrices(Vs,Hs,hid);

	double min_J_divided=10000;
	for(int j=0;j<Jaco_matrices.size();j++)
	{
		MatrixXd Jacoi=Jaco_matrices[j];

		double len1=Jacoi.row(0).norm();
		double len2=Jacoi.row(1).norm();
		double len3=Jacoi.row(2).norm();

		double absj1=abs(Jacoi.determinant()*(27/pow((len1+len2+len3),3)));//old tao
		if(min_J_divided>absj1)
			min_J_divided=absj1;
	}
	return min_J_divided;
}
double hex_metrics::hex_jacobian_t2(vector<Vertex> &Vs,vector<Hex> &Hs,int hid, vector<vector<double>> coordinates)
{
	double weight=.57735;
	for(int j=0;j<8;j++)
	{
		RST[j][0]=rst[j][0]/abs(rst[j][0])*weight;
		RST[j][1]=rst[j][1]/abs(rst[j][1])*weight;
		RST[j][2]=rst[j][2]/abs(rst[j][2])*weight;
	}

	jacobian_matrices(Vs,Hs,hid);
	double minimum_J_divided2=10000;
	for(int j=0;j<Jaco_matrices.size();j++)
	{
		MatrixXd Jacoi=Jaco_matrices[j];

		double len1=Jacoi.row(0).norm();
		double len2=Jacoi.row(1).norm();
		double len3=Jacoi.row(2).norm();

		if(len1<len2)
			len1 = len2;
		if(len1<len3)
			len1 = len3;

		double absj2=abs(Jacoi.determinant()*(1.0/pow((double)len1,(double)3)));//old tao	
		if(minimum_J_divided2>absj2)
			minimum_J_divided2=absj2;
	}
	return minimum_J_divided2;
}

double hex_metrics::hex_edge_length(int max_min, vector<vector<double>> coordinates)
{
	double temp[3], edge[12];
	int i;

	//lengths^2 of edges
	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[1][i] - coordinates[0][i];
		temp[i] = temp[i] * temp[i];
	}
	edge[0] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[2][i] - coordinates[1][i];
		temp[i] = temp[i] * temp[i];
	}
	edge[1] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[3][i] - coordinates[2][i];
		temp[i] = temp[i] * temp[i];
	}
	edge[2] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[0][i] - coordinates[3][i];
		temp[i] = temp[i] * temp[i];
	}
	edge[3] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[5][i] - coordinates[4][i];
		temp[i] = temp[i] * temp[i];
	}
	edge[4] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[6][i] - coordinates[5][i];
		temp[i] = temp[i] * temp[i];
	}
	edge[5] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[7][i] - coordinates[6][i];
		temp[i] = temp[i] * temp[i];
	}
	edge[6] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[4][i] - coordinates[7][i];
		temp[i] = temp[i] * temp[i];
	}
	edge[7] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[4][i] - coordinates[0][i];
		temp[i] = temp[i] * temp[i];
	}
	edge[8] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[5][i] - coordinates[1][i];
		temp[i] = temp[i] * temp[i];
	}
	edge[9] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[6][i] - coordinates[2][i];
		temp[i] = temp[i] * temp[i];
	}
	edge[10] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[7][i] - coordinates[3][i];
		temp[i] = temp[i] * temp[i];
	}
	edge[11] = sqrt( temp[0] + temp[1] + temp[2] );

	double _edge = edge[0];

	if ( max_min == 0)
	{
		for( i = 1; i<12; i++) 
			if(_edge>edge[i])
				_edge=edge[i];
		return (double)_edge;
	}  
	else
	{
		for( i = 1; i<12; i++) 
			if(_edge<edge[i])
				_edge=edge[i];
		return (double)_edge;
	}
}
double hex_metrics::diag_length(int max_min, vector<vector<double>> coordinates)
{
	double temp[3], diag[4];
	int i;
	//lengths^2  f diag nals
	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[6][i] - coordinates[0][i];
		temp[i] = temp[i] * temp[i];
	}
	diag[0] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[4][i] - coordinates[2][i];
		temp[i] = temp[i] * temp[i];
	}
	diag[1] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[7][i] - coordinates[1][i];
		temp[i] = temp[i] * temp[i];
	}
	diag[2] = sqrt( temp[0] + temp[1] + temp[2] );

	for (i = 0; i < 3; i++ )
	{
		temp[i] = coordinates[5][i] - coordinates[3][i];
		temp[i] = temp[i] * temp[i];
	}
	diag[3] = sqrt( temp[0] + temp[1] + temp[2] );

	double diagonal = diag[0];
	if ( max_min == 0 )  //Return min diagonal
	{ 
		for( i = 1; i<4; i++)
			if(diagonal > diag[i] )
				diagonal = diag[i];
		return diagonal;
	}
	else          //Return max diagonal
	{
		for( i = 1; i<4; i++)
			if(diagonal < diag[i] )
				diagonal = diag[i];
		return diagonal;  
	}
}
void hex_metrics::make_hex_edges(vector<vector<double>> coordinates, vector<vector<double>> &edgevectors)
{
	vector<double> edge;
	edge.push_back(coordinates[1][0] - coordinates[0][0]);
	edge.push_back(coordinates[1][1] - coordinates[0][1]);
	edge.push_back(coordinates[1][2] - coordinates[0][2]);
	edgevectors.push_back(edge);
	edge.clear();

	edge.push_back(coordinates[2][0] - coordinates[1][0]);
	edge.push_back(coordinates[2][1] - coordinates[1][1]);
	edge.push_back(coordinates[2][2] - coordinates[1][2]);
	edgevectors.push_back(edge);
	edge.clear();

	edge.push_back(coordinates[3][0] - coordinates[2][0]);
	edge.push_back(coordinates[3][1] - coordinates[2][1]);
	edge.push_back(coordinates[3][2] - coordinates[2][2]);
	edgevectors.push_back(edge);
	edge.clear();

	edge.push_back(coordinates[0][0] - coordinates[3][0]);
	edge.push_back(coordinates[0][1] - coordinates[3][1]);
	edge.push_back(coordinates[0][2] - coordinates[3][2]);
	edgevectors.push_back(edge);				   
	edge.clear();
											   
	edge.push_back(coordinates[5][0] - coordinates[4][0]);
	edge.push_back(coordinates[5][1] - coordinates[4][1]);
	edge.push_back(coordinates[5][2] - coordinates[4][2]);
	edgevectors.push_back(edge);				   
	edge.clear();
											   
	edge.push_back(coordinates[6][0] - coordinates[5][0]);
	edge.push_back(coordinates[6][1] - coordinates[5][1]);
	edge.push_back(coordinates[6][2] - coordinates[5][2]);
	edgevectors.push_back(edge);				   
	edge.clear();
										   
	edge.push_back(coordinates[7][0] - coordinates[6][0]);
	edge.push_back(coordinates[7][1] - coordinates[6][1]);
	edge.push_back(coordinates[7][2] - coordinates[6][2]);
	edgevectors.push_back(edge);				   
	edge.clear();
											   
	edge.push_back(coordinates[4][0] - coordinates[7][0]);
	edge.push_back(coordinates[4][1] - coordinates[7][1]);
	edge.push_back(coordinates[4][2] - coordinates[7][2]);
	edgevectors.push_back(edge);				   
	edge.clear();
											   
	edge.push_back(coordinates[4][0] - coordinates[0][0]);
	edge.push_back(coordinates[4][1] - coordinates[0][1]);
	edge.push_back(coordinates[4][2] - coordinates[0][2]);
	edgevectors.push_back(edge);				   
	edge.clear();
											   
	edge.push_back(coordinates[5][0] - coordinates[1][0]);
	edge.push_back(coordinates[5][1] - coordinates[1][1]);
	edge.push_back(coordinates[5][2] - coordinates[1][2]);
	edgevectors.push_back(edge);				   
	edge.clear();
											   
	edge.push_back(coordinates[6][0] - coordinates[2][0]);
	edge.push_back(coordinates[6][1] - coordinates[2][1]);
	edge.push_back(coordinates[6][2] - coordinates[2][2]);
	edgevectors.push_back(edge);				   
	edge.clear();
										   
	edge.push_back(coordinates[7][0] - coordinates[3][0]);
	edge.push_back(coordinates[7][1] - coordinates[3][1]);
	edge.push_back(coordinates[7][2] - coordinates[3][2]);
	edgevectors.push_back(edge);
	edge.clear();

}
void hex_metrics::calc_hex_efg( int efg_index, vector<double> &efg,vector<vector<double>> coordinates)
{

	//double efg[3];
	efg.push_back(0);
	efg.push_back(0);
	efg.push_back(0);

	switch(efg_index) {
	case 1:
		VECTOR_ADD(efg, coordinates[1]);
		VECTOR_ADD(efg, coordinates[2]);
		VECTOR_ADD(efg, coordinates[5]);
		VECTOR_ADD(efg, coordinates[6]);
		VECTOR_MINUS(efg, coordinates[0]);
		VECTOR_MINUS(efg, coordinates[3]);
		VECTOR_MINUS(efg, coordinates[4]);
		VECTOR_MINUS(efg, coordinates[7]);
		break;
	case 2:
		VECTOR_ADD(efg, coordinates[2]);
		VECTOR_ADD(efg, coordinates[3]);
		VECTOR_ADD(efg, coordinates[6]);
		VECTOR_ADD(efg, coordinates[7]);
		VECTOR_MINUS(efg, coordinates[0]);
        VECTOR_MINUS(efg, coordinates[1]);
		VECTOR_MINUS(efg, coordinates[4]);
		VECTOR_MINUS(efg, coordinates[5]);
		break;
	case 3:
		VECTOR_ADD(efg, coordinates[4]);
		VECTOR_ADD(efg, coordinates[5]);
		VECTOR_ADD(efg, coordinates[6]);
		VECTOR_ADD(efg, coordinates[7]);
		VECTOR_MINUS(efg, coordinates[0]);
		VECTOR_MINUS(efg, coordinates[1]);
		VECTOR_MINUS(efg, coordinates[2]);
		VECTOR_MINUS(efg, coordinates[3]);
		break;
	case 12:
		VECTOR_ADD(efg, coordinates[0]);
		VECTOR_ADD(efg, coordinates[2]);
		VECTOR_ADD(efg, coordinates[4]);
		VECTOR_ADD(efg, coordinates[6]);
		VECTOR_MINUS(efg, coordinates[1]);
		VECTOR_MINUS(efg, coordinates[3]);
		VECTOR_MINUS(efg, coordinates[5]);
		VECTOR_MINUS(efg, coordinates[7]);
		break;
	case 13:
		VECTOR_ADD(efg, coordinates[0]);
		VECTOR_ADD(efg, coordinates[3]);
		VECTOR_ADD(efg, coordinates[5]);
		VECTOR_ADD(efg, coordinates[6]);
		VECTOR_MINUS(efg, coordinates[1]);
		VECTOR_MINUS(efg, coordinates[2]);
		VECTOR_MINUS(efg, coordinates[4]);
		VECTOR_MINUS(efg, coordinates[7]);
		break;

	case 23:
		VECTOR_ADD(efg, coordinates[0]);
		VECTOR_ADD(efg, coordinates[1]);
		VECTOR_ADD(efg, coordinates[6]);
		VECTOR_ADD(efg, coordinates[7]);
		VECTOR_MINUS(efg, coordinates[2]);
		VECTOR_MINUS(efg, coordinates[3]);
		VECTOR_MINUS(efg, coordinates[4]);
		VECTOR_MINUS(efg, coordinates[5]);
		break;

	case 123://probably error
		VECTOR_ADD(efg, coordinates[0]);
		VECTOR_ADD(efg, coordinates[2]);
		VECTOR_ADD(efg, coordinates[5]);
		VECTOR_ADD(efg, coordinates[7]);
		VECTOR_MINUS(efg, coordinates[1]);
		VECTOR_MINUS(efg, coordinates[5]);
		VECTOR_MINUS(efg, coordinates[6]);
		VECTOR_MINUS(efg, coordinates[2]);
		break;
	}

	//return efg;
}
double hex_metrics::jacobian_fun(vector<double> v0,vector<double> v1,vector<double> v2,vector<double> v3)
{//need to verify further whether correct or not
	double v1v0[3],v2v0[3],v3v0[3];
	for(int i=0;i<3;i++)
	{
		v1v0[i]=v1[i]-v0[i];
		v2v0[i]=v2[i]-v0[i];
		v3v0[i]=v3[i]-v0[i];
	}

	return v1v0[0]*(v2v0[1]*v3v0[2]-v2v0[2]*v3v0[1])-v1v0[1]*(v2v0[0]*v3v0[2]-v2v0[2]*v3v0[0])+v1v0[2]*(v2v0[0]*v3v0[1]-v2v0[1]*v3v0[0]); 	
}
double hex_metrics::condition_comp( vector<double> xxi, vector<double> xet, vector<double> xze)
{
	double det = dot_cross_vectors(xxi,xet,xze);
	double dot1=DOTVECTOR3(xxi,xxi);
	double dot2=DOTVECTOR3(xet,xet);
	double dot3=DOTVECTOR3(xze,xze);
	double term1 = dot1+ dot2 + dot3;
	double cross1[3],cross2[3],cross3[3];
	CROSSVECTOR3(cross1,xxi,xet);
	CROSSVECTOR3(cross2,xet,xze);
	CROSSVECTOR3(cross3,xze,xxi);

	dot1 = DOTVECTOR3(cross1,cross1);dot2 = DOTVECTOR3(cross2,cross2);
	dot3= DOTVECTOR3(cross3,cross3);
	double term2 =  dot1+ dot2+ dot3;
	return sqrt( term1 * term2 ) / det;
}
double hex_metrics::oddy_comp( vector<double> xxi, vector<double> xet, vector<double> xze)
{
	static const double third=1.0/3.0;

	double g11, g12, g13, g22, g23, g33, rt_g;

	g11 = DOTVECTOR3(xxi,xxi);
	g12 = DOTVECTOR3(xxi,xet); 
	g13 = DOTVECTOR3(xxi,xze); 
	g22 = DOTVECTOR3(xet,xet); 
	g23 = DOTVECTOR3(xet,xze); 
	g33 = DOTVECTOR3(xze,xze);
	rt_g = dot_cross_vectors(xxi,xet,xze);

	double oddy_metric;
	{
		double norm_G_squared = g11*g11 + 2.0*g12*g12 + 2.0*g13*g13 + g22*g22 + 2.0*g23*g23 +g33*g33;

		double norm_J_squared = g11 + g22 + g33;

		oddy_metric = ( norm_G_squared - third*norm_J_squared*norm_J_squared ) / pow( rt_g, 4.*third );
	}
	return oddy_metric;
}
int hex_metrics::v_hex_get_weight(vector<double> &v1,vector<double> &v2,vector<double> &v3)
{
	if ( average_volume == 0)
		return 0;

	v1[0] =1;v1[1]=v1[2]=0;
	v2[0] =0;v2[1]=1; v2[2]=0;
	v3[0] =0;v3[1]=0; v3[2]=1;

	double scale = pow((double)average_volume/ dot_cross_vectors(v1,v2,v3), (double)0.33333333333); 
	for(int i = 0; i<3;i++)
	{
		v1[i]=v1[i]*scale;
		v2[i]=v2[i]*scale;
		v3[i]=v3[i]*scale;
	}

	return 1;
}
double hex_metrics::scaled_jacobian(vector<double> v0,vector<double> v1,vector<double> v2,vector<double> v3)
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

	for(int i=0;i<3;i++)
	{
		if(abs(norm1)>=EPS)
			v1v0[i]=v1v0[i]/norm1;
		else
			v1v0[i]=0;
		if(abs(norm2)>=EPS)
			v2v0[i]=v2v0[i]/norm2;
		else
			v2v0[i]=0;
		if(abs(norm2)>=EPS)
			v3v0[i]=v3v0[i]/norm3;
		else
			v3v0[i]=0;
	}

	double scaled_jacobian=v1v0[0]*(v2v0[1]*v3v0[2]-v2v0[2]*v3v0[1])-v1v0[1]*(v2v0[0]*v3v0[2]-v2v0[2]*v3v0[0])+v1v0[2]*(v2v0[0]*v3v0[1]-v2v0[1]*v3v0[0]); 
	return scaled_jacobian;
}
double hex_metrics::dot_cross_vectors(vector<double> v0,vector<double> v1,vector<double> v2)
{
	double v12[3];
	CROSSVECTOR3(v12,v1,v2);
	return DOTVECTOR3(v0,v12);
}
void hex_metrics::jacobian_matrices(vector<Vertex> &Vs,vector<Hex> &Hs,int hid)
{
	Jaco.clear();
	Jaco_matrices.clear();

	vector<double> coords;
	for(int i=0;i<8;i++)
	{
		coords.push_back(Vs[Hs[hid].vid[i]].v[0]);
		coords.push_back(Vs[Hs[hid].vid[i]].v[1]);
		coords.push_back(Vs[Hs[hid].vid[i]].v[2]);
	}
	//jacobian
	for(int i=0;i<8;i++)
	{
		MatrixXd Jacoi(3,3);
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<3;k++)
			{
				Jacoi(j,k)=0;
				for(int m=0;m<8;m++)
				{
					if(j==0)
					{
						Jacoi(j,k)+=1.0/8*rst[m][0]*(1+RST[i][1]*rst[m][1])*(1+RST[i][2]*rst[m][2])*Vs[Hs[hid].vid[m]].v[k];
					}
					else if(j==1)
					{
						Jacoi(j,k)+=1.0/8*rst[m][1]*(1+RST[i][0]*rst[m][0])*(1+RST[i][2]*rst[m][2])*Vs[Hs[hid].vid[m]].v[k];						
					}
					else if(j==2)
					{
						Jacoi(j,k)+=1.0/8*rst[m][2]*(1+RST[i][0]*rst[m][0])*(1+RST[i][1]*rst[m][1])*Vs[Hs[hid].vid[m]].v[k];
					}
				}
			}
		}
		Jaco_matrices.push_back(Jacoi);
		Jaco.push_back(Jacoi.determinant());
	}
}
void hex_metrics::jacobian_matrices_derivative(vector<Vertex> &Vs,vector<Hex> &Hs,int hid)
{
	Jaco_D.clear();
	Jaco_matrices_D.clear();

	vector<double> coords;
	for(int i=0;i<8;i++)
	{
		coords.push_back(Vs[Hs[hid].vid[i]].v[0]);
		coords.push_back(Vs[Hs[hid].vid[i]].v[1]);
		coords.push_back(Vs[Hs[hid].vid[i]].v[2]);
	}
	//jacobian
	for(int i=0;i<8;i++)
	{
		MatrixXd Jacoi_D(3,3);
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<3;k++)
			{
				Jacoi_D(j,k)=0;
				for(int m=0;m<8;m++)
				{
					if(j==0)
					{
						Jacoi_D(j,k)+=1.0/8*rst[m][0]*(1+rst[i][1]*rst[m][1])*(1+rst[i][2]*rst[m][2])*Vs[Hs[hid].vid[m]].v[k];
					}
					else if(j==1)
					{
						Jacoi_D(j,k)+=1.0/8*rst[m][1]*(1+rst[i][0]*rst[m][0])*(1+rst[i][2]*rst[m][2])*Vs[Hs[hid].vid[m]].v[k];						
					}
					else if(j==2)
					{
						Jacoi_D(j,k)+=1.0/8*rst[m][2]*(1+rst[i][0]*rst[m][0])*(1+rst[i][1]*rst[m][1])*Vs[Hs[hid].vid[m]].v[k];
					}
				}
			}
		}
		Jaco_matrices_D.push_back(Jacoi_D);
		Jaco_D.push_back(Jacoi_D.determinant());
	}
}
void hex_metrics::hex2tetmesh(vector<Vertex> &Vs,vector<Edge> &Es,vector<Cube_F> &Fs,vector<Hex> &Hs)
{
	TVs.clear();TTs.clear();
	for(int i=0;i<Fs.size();i++)
	{
		Fs[i].center[2]=Fs[i].center[1]=Fs[i].center[0]=0;
		for(int j=0;j<4;j++)
		{
			for(int k=0;k<3;k++)
			{
				Fs[i].center[k]+=Vs[Fs[i].cv[j].index].v[k];
			}
		}
		Fs[i].center[2]/=4;
		Fs[i].center[1]/=4;
		Fs[i].center[0]/=4;
	}
	for(int i=0;i<Hs.size();i++)
	{
		Hs[i].center[2]=Hs[i].center[1]=Hs[i].center[0]=0;
		for(int j=0;j<8;j++)
		{
			for(int k=0;k<3;k++)
			{
				Hs[i].center[k]+=Vs[Hs[i].vid[j]].v[k];
			}
		}
		Hs[i].center[2]/=8;
		Hs[i].center[1]/=8;
		Hs[i].center[0]/=8;
	}

	TVs=Vs;
	vector<Vertex> FVs,HVs;
	//face v
	for(int i=0;i<Fs.size();i++)
	{
		Vertex v;
		v.index=TVs.size();
		v.v[0]=Fs[i].center[0];
		v.v[1]=Fs[i].center[1];
		v.v[2]=Fs[i].center[2];
		FVs.push_back(v);
		TVs.push_back(v);
	}
	//hex v
	for(int i=0;i<Hs.size();i++)
	{
		Vertex v;
		v.index=TVs.size();
		v.v[0]=Hs[i].center[0];
		v.v[1]=Hs[i].center[1];
		v.v[2]=Hs[i].center[2];
		HVs.push_back(v);
		TVs.push_back(v);
	}
	//tets for a hex
	for(int i=0;i<Hs.size();i++)
	{
		for(int j=0;j<6;j++)
		{
			int fid=Hs[i].neighborF[j];
			for(int k=0;k<4;k++)
			{
				int v1=Es[Fs[fid].eids[k]].vid[0],v2=Es[Fs[fid].eids[k]].vid[1];
				Tet t;
				t.vid.push_back(HVs[i].index);
				t.vid.push_back(v1);
				t.vid.push_back(v2);
				t.vid.push_back(FVs[fid].index);
				TTs.push_back(t);
			}
		}
	}
}
hex_metrics::~hex_metrics(void)
{
}
