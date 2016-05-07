#include "computation.h"

void cmp_ez(
	E e, H h,
	EXY exy, EYX eyx, COE c, src s, int time)
{
	int i, j;
	float xcoor, ycoor, kappa;
	area *local = e.ez;

	for (i = 0; i < e.ez->height; i++){
		for (j = 0; j < e.ez->width; j++){
			//transfer the indexes of vector of the coordinates in field
			xcoor = j;
			ycoor = i;

			//FDTD area
			//if the point is on the boundaries
			if (i == 0 || i == (e.ez->height - 1)
				|| j == 0 || j == (e.ez->width - 1)){
				e.ez->p.at(i*e.ez->width + j) = 0.f;
				continue;
			}
			kappa = c.set_kappa(s, ycoor, xcoor);
			e.ez->p.at(i*e.ez->width + j) +=
				e.coe_E*(
				(h.hy->p.at(i*h.hy->width + j) - h.hy->p.at((i - 1)*h.hy->width + j)) / (kappa*s.dz)
				- (h.hx->p.at(i*h.hx->width + j) - h.hx->p.at(i*h.hx->width + j - 1)) / (kappa*s.dz)
				)
				+ e.coe_E*(exy.full->p.at(i*exy.full->width + j) - eyx.full->p.at(i*eyx.full->width + j));
		}
	}
}

void cmp_hx(H h, E e, HYZ hyz, COE c, src s)
{
	int i, j;
	float xcoor, ycoor;

	for (i = 0; i < h.hx->height; i++)
	{
		for (j = 0; j < h.hx->width; j++)
		{
			//transfer the indexes of vector of the coordinates in field
			xcoor = j + 0.5f;
			ycoor = i;
			h.hx->p.at(i*h.hx->width + j) +=
				(-h.coe_h)*(
				e.ez->p.at(i*e.ez->width + j + 1) - e.ez->p.at(i*e.ez->width + j)
				) / (c.set_kappa(s, ycoor, xcoor)*s.dz)
				+ h.coe_h*(-hyz.full->p.at(i*hyz.full->width + j));
		}
	}
}

void cmp_hy(H h, E e, HXZ hxz, COE c, src s)
{
	int i, j;
	float xcoor, ycoor;
	area *local;
	local = h.hy;
	for (i = 0; i < h.hy->height; i++)
	{
		for (j = 0; j < h.hy->width; j++)
		{
			//transfer the indexes of vector to the coordinates in field
			xcoor = j;
			ycoor = i + 0.5f;
			local->p.at(i*local->width + j) +=
				h.coe_h*(
				e.ez->p.at((i + 1)*e.ez->width + j) - e.ez->p.at(i*e.ez->width + j)
				) / (c.set_kappa(s, ycoor, xcoor)*s.dz)
				+ h.coe_h*(hxz.full->p.at(i*hxz.full->width + j));
		}
	}
}

void cmp_hxz(E e, HXZ hxz, COE c, src s)
{
	int i, j;
	float xcoor, ycoor, kappa;
	area *local;
	local = hxz.full;
	for (i = 0; i < local->height; i++)
	{
		for (j = 0; j < local->height; j++)
		{
			xcoor = j;
			ycoor = i + 0.5f;
			local->p.at(i*local->width + j) = c.set_c(s, ycoor, xcoor)*(
				e.ez->p.at((i + 1)*e.ez->width + j) - e.ez->p.at(i*e.ez->width + j)
				)
				+ c.set_coe(s, ycoor, xcoor)*local->p.at(i*local->width + j);
		}
	}
}

void cmp_hyz(E e, HYZ hyz, COE c, src s)
{
	int i, j;
	float xcoor, ycoor;
	area *local = hyz.full;
	for (i = 0; i < local->height; i++){
		for (j = 0; j < local->width; j++){
			xcoor = j + 0.5f;
			ycoor = i;
			local->p.at(i*local->width + j) = c.set_c(s, ycoor, xcoor)*(
				e.ez->p.at(i*e.ez->width + j + 1) - e.ez->p.at(i*e.ez->width + j)
				) +
				c.set_coe(s, ycoor, xcoor)*local->p.at(i*local->width + j);
		}
	}
}

void cmp_exy(H h, EXY exy, COE c, src s)
{
	int i, j;
	float xcoor, ycoor;
	area *local = exy.full;
	for (i = 0; i < local->height; i++){
		for (j = 0; j < local->width; j++){
			if (i == 0 || i == (local->height - 1)
				|| j == 0 || j == (local->width - 1)){
				local->p.at(i*local->width + j) = 0.f;
				continue;
			}
			xcoor = i;
			ycoor = j;
			local->p.at(i*local->width + j) = c.set_c(s, ycoor, xcoor)*(
				h.hy->p.at(i*h.hy->width + j) - h.hy->p.at((i - 1)*h.hy->width + j)
				) +
				c.set_c(s, ycoor, xcoor)*local->p.at(i*local->width + j);
		}
	}
}

void cmp_eyx(H h, EYX eyx, COE c, src s)
{
	int i, j;
	float xcoor, ycoor;
	area *local = eyx.full;
	for (i = 0; i < local->height; i++){
		for (j = 0; j < local->width; j++){
			if (i == 0 || i == (local->height - 1)
				|| j == 0 || j == (local->width - 1)){
				local->p.at(i*local->width + j) = 0.f;
				continue;
			}
			xcoor = j;
			ycoor = i;
			local->p.at(i*local->width + j) = c.set_c(s, ycoor, xcoor)*(
				h.hx->p.at(i*h.hx->width + j) - h.hx->p.at(i*h.hx->p.at(i*h.hx->width + j - 1))
				) +
				c.set_coe(s, ycoor, xcoor)*local->p.at(i*local->width + j);
		}
	}
}