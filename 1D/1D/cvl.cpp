#include "cvl.h"


cvl::cvl(src source)
{
	Hzx = (float *)malloc((source.size_space)*sizeof(float));
	memset(Hzx, 0, (source.size_space)*sizeof(float));
	Hxz = (float *)malloc((source.size_space)*sizeof(float));
	memset(Hxz, 0, (source.size_space)*sizeof(float));
	Ezy = (float*)malloc((source.size_space + 1)*sizeof(float));
	memset(Ezy, 0, (source.size_space + 1)*sizeof(float));
	Eyz = (float*)malloc((source.size_space + 1)*sizeof(float));
	memset(Eyz, 0, (source.size_space + 1)*sizeof(float));

	dz = source.dz;
	dt = source.dt;
	size_space = source.size_space;
	d = 10.f * source.dz;
}


cvl::~cvl()
{
}

void cvl::set_sigma(float z, float z0)
{
	sigma = sigma_max*powf(fabsf(z - z0), m) / powf(d, m);
}

void cvl::set_kappa(float z, float z0)
{
	kappa = 1 + (kappa_max - 1) * (powf(fabsf(z - z0), m)) / powf(d, m);
}

void cvl::set_alpha(float z, float z0)
{
	set_sigma(z, z0);
	set_kappa(z, z0);
	alpha = sigma / (eps0 * kappa);
}

void cvl::set_c()
{
	c = (1 / kappa) * (exp(-alpha*dt) - 1);
}

void cvl::cmp_H(int current_timestep, int z0, float* cvl_H, float *Ex)
{
	int i;
	for (i = 0; i < size_space; i++)
	{
		set_alpha(i, z0);
		set_c();
		
		cvl_H[i] = exp(-alpha*dt)*cvl_H[i] + c*(Ex[i + 1] - Ex[i]) / dz;
	}
}

void cvl::cmp_H(int current_timestep, int z0, float* cvl_H)
{
	int i;
	for (i = 0; i < size_space; i++)
	{
		set_alpha(i, z0);
		set_c();

		cvl_H[i] = exp(-alpha * dt)*cvl_H[i];
	}
}

void cvl::cmp_E(int current_timestep, int z0, float *cvl_E)
{
	int i;
	for (i = 1; i < (size_space + 1) - 1; i++)
	{
		set_alpha(i, z0);
		set_c();

		cvl_E[i] = exp(-alpha * dt) * cvl_E[i];
	}
}

void cvl::cmp_E(int current_timestep, int z0, float *cvl_E, float *Hy)
{
	int i;
	for (i = 1; i < (size_space + 1) - 1; i++)
	{
		set_alpha(i, z0);
		set_c();

		cvl_E[i] = exp(-alpha * dt) * cvl_E[i] + c * (Hy[i] - Hy[i - 1]) / dz;
	}
}