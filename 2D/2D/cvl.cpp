#include "cvl.h"
//#define debug

/*Class COE*/
COE::COE(int number_layer, src s)
{
	num_layer = number_layer;
	d = number_layer*s.dz;
	sigma_max = (m + 1) / (sqrt(eps_r) * 150 * PI * s.dz);
}

//the distance of every field points of H and E in CPML layer.
float COE::set_distance(src source, float ycoor, float xcoor)
{
	float bench_x, bench_y;
	//for bottom left & right corner
	if (ycoor >= 0 && ycoor < 3.0)
	{
		if (xcoor >= 0.0 && xcoor < 3.0){
			bench_x = 3.0f;
			bench_y = 3.0f;
		}
		if (xcoor > 13.0f && xcoor <= 16.0f){
			bench_x = 13.0f;
			bench_y = 3.0f;
		}
		distance = sqrtf(
			powf(fabsf(bench_x - xcoor), 2) + sqrtf(powf(fabsf(bench_y - ycoor), 2))
			);
		return distance;
	}
	//for bottom and top center
	if (xcoor >= 3 && xcoor <= 13)
	{
		//for bottom center area
		if (ycoor >= 0.f && ycoor < 3.0f){
			bench_y = 3.0f;
		}
		//for top center area
		if (ycoor > 13.f && ycoor <= 16.f){
			bench_y = 13.0f;
		}
		distance = fabsf(bench_y - ycoor);
		return distance;
	}
	//for left and right center
	if (ycoor >= 3.f && ycoor <= 13.f){
		//for left center area
		if (xcoor >= 0.f && xcoor < 3.f){
			bench_x = 3.f;
		}
		//for right center area
		if (xcoor > 13.f && xcoor <= 16.f){
			bench_x = 13.f;
		}
		distance = fabsf(bench_x - xcoor);
		return distance;
	}
	//for top's left and right corner
	if (ycoor>13.f && ycoor <=16.f)
	{
		//left corner
		if (xcoor>=0.f && xcoor <3.f){
			bench_x = 3.f;
			bench_y = 13.f;
		}
		if (xcoor > 13.f && xcoor<=16.f){
			bench_y = 13.f;
			bench_x = 13.f;
		}
		distance = sqrtf(
			powf(fabsf(bench_x - xcoor), 2) + sqrtf(powf(fabsf(bench_y - ycoor), 2))
			);
		return distance;
	}
	//for FDTD area
	return 0;
}

float COE::set_sigma(src source, float ycoor, float xcoor)
{
	//		sigma_full[i] = sigma_max * powf((distance_full[i] / d), m);
	return 0;
}

float COE::set_kappa(src source, float ycoor, float xcoor)
{
	return 1;
//	kappa = 1 + (kappa_max - 1)*powf(set_distance(s, ycoor, xcoor) / d, m);
//	return kappa;
}

float COE::set_alpha(src source, float ycoor, float xcoor)
{
	//		alpha_full[i] = sigma_full[i] / (eps0 * kappa_full[i]);
	return 0;
}

float COE::set_c(src source, float ycoor, float xcoor)
{
	//		c_full[i] = (1 / kappa_full[i])*(exp(-alpha_full[i] * source.dt) - 1);
	return 0;
}

float COE::set_coe(src source, float ycoor, float xcoor)
{
	//		cvl_full_coe[i] = exp(-alpha_full[i] * source.dt);
	return 0;
}

void COE::checkout()
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

	cout << "\tside grid number x(side_width):\t" << side_width << endl;
	cout << "\tside grid number y(side_height):\t" << side_height << endl;
	cout << "\tup & down grid number x(ud_width):\t" << ud_width << endl;
	cout << "\tup & down grid number y(ud_height):\t" << ud_height << endl;
	cout << endl;
	//for variables
	//int i;
	//cout << "c_full" << endl;
	//for (i = 0; i < num_layer; i++)
	//{
	//	cout << c_full[i] << endl;
	//}
}

COE::~COE()
{
}


/******************* Class COV ************************/
COV::COV()
{
}

COV::~COV()
{
	delete left;
	left = NULL;
	delete right;
	right = NULL;
	delete top;
	top = NULL;
	delete btm;
	btm = NULL;
}

/**************************** Convolutional Terms *******************************/
HXZ::HXZ(COE c, src s)
{
	full = new area(s.size_x + 2 * c.num_layer + 1, s.size_y + 2 * c.num_layer, "hy.txt");
	//	left = new area(c.num_layer, s.size_y, "hxzl.txt");
	//	right = new area(c.num_layer, s.size_y, "hxzr.txt");
	//	top = new area(s.size_x + 1, c.num_layer, "hxzt.txt");
	//	btm = new area(s.size_x + 1, c.num_layer, "hxzb.txt");
}

HYZ::HYZ(COE c, src s)
{
	full = new area(s.size_x + 2 * c.num_layer, s.size_y + 2 * c.num_layer + 1, "hx.txt");
	//	left = new area(c.num_layer, s.size_y + 2 * c.num_layer + 1, "hyzl.txt");
	//	right = new area(c.num_layer, s.size_y + 2 * c.num_layer + 1, "hyzr.txt");
	//	top = new area(s.size_x, c.num_layer, "hyzt.txt");
	//	btm = new area(s.size_x, c.num_layer, "hyzb.txt");
}

EXY::EXY(COE c, src s)
{
	full = new area(s.size_x + 2 * c.num_layer + 1, s.size_y + 2 * c.num_layer + 1, "ez.txt");
	//	left = new area(c.num_layer, s.size_y + 2 * c.num_layer + 1, "exyl.txt");
	//	right = new area(c.num_layer, s.size_y + 2 * c.num_layer + 1, "exyr.txt");
	//	top = new area(s.size_x + 1, c.num_layer, "exyt.txt");
	//	btm = new area(s.size_x + 1, c.num_layer, "exyb.txt");
}

EYX::EYX(COE c, src s)
{
	full = new area(s.size_x + 2 * c.num_layer + 1, s.size_y + 2 * c.num_layer + 1, "ez.txt");
	//	left = new area(c.num_layer, s.size_y + 2 * c.num_layer + 1, "eyxl.txt");
	//	right = new area(c.num_layer, s.size_y + 2 * c.num_layer + 1, "eyxr.txt");
	//	top = new area(s.size_x + 1, c.num_layer, "eyxt.txt");
	//	btm = new area(s.size_x + 1, c.num_layer, "eyxb.txt");
}