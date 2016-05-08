#include "cvl.h"
//#define debug

/*Class COE*/
COE::COE(int number_layer, src s)
{
	num_layer = number_layer;
	d = number_layer*s.dz;
	sigma_max = (m + 1) / (sqrt(eps_r) * 150 * PI * s.dz);
	fstream myfile;
	myfile.open(filename, ios::out);
	myfile.close();
}

//the distance of every field points of H and E in CPML layer.
float COE::set_distance(src s, float xcoor)
{
	float bench_x;
	//for bottom left & right corner
	if (xcoor >= 0.f && xcoor < num_layer){
		bench_x = num_layer;
		distance = sqrtf(
			powf(fabsf(bench_x - xcoor), 2)
			);
		return distance*s.dz;
	}
	if (xcoor > (num_layer + s.size_x) && xcoor <= (2 * num_layer + s.size_x)){
		bench_x = num_layer + s.size_x;
		distance = sqrtf(
			powf(fabsf(bench_x - xcoor), 2)
			);
		return distance*s.dz;
	}
	//for FDTD area
	if (xcoor >= num_layer && xcoor <= num_layer + s.size_x){
		return 0;
	}
}

float COE::set_sigma(src s, float xcoor)
{
	sigma = sigma_max*powf(set_distance(s, xcoor) / d, m);
	return sigma;
}

float COE::set_kappa(src s, float xcoor)
{
	//kappa=1;
	kappa = 1 + ((kappa_max-1)*powf(set_distance(s, xcoor) / d, m));
	return kappa;
}

float COE::set_alpha(src s, float xcoor)
{
	alpha = set_sigma(s, xcoor) / (eps0*set_kappa(s, xcoor));
	return alpha;
}

float COE::set_c(src s, float xcoor)
{
	kappa = set_kappa(s, xcoor);
	alpha = set_alpha(s, xcoor);
	c = (1 / kappa)
		*(exp(-alpha*s.dt) - 1);
	return c;
}

float COE::set_coe(src s, float xcoor)
{
	cvl_coe = exp(-set_alpha(s, xcoor)*s.dt);
	return cvl_coe;
}

void COE::checkout(src s)
{
	//for coefficient
	cout << "convolution class checkout:" << endl;
	cout << "\tm = " << m << endl;
	cout << "\teps_r = " << eps_r << endl;
	cout << "\teps0 = " << eps0 << endl;
	cout << "\tsigma_max = " << sigma_max << endl;
	cout << "\tkappa_max = " << kappa_max << endl;
	cout << "\td = " << d << endl;
	cout << endl;

	int i, j;
	fstream myfile;
	myfile.open(filename, ios::app);
	for (j = 0; j < 17; j++){
		myfile << set_distance(s, j) << "\t";
	}
	myfile << endl;
	for (j = 0; j < 17; j++){
		myfile << set_kappa(s, j) << "\t";
	}
	myfile << endl;
	for (j = 0; j < 17; j++){
		myfile << set_sigma(s, j) << "\t";
	}
	myfile << endl;
	for (j = 0; j < 17; j++){
		myfile << set_alpha(s, j) << "\t";
	}
	myfile << endl;
	for (j = 0; j < 17; j++){
		myfile << set_c(s, j) << "\t";
	}
	myfile << endl;
	myfile.close();
}

COE::~COE()
{
}
