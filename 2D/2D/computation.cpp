#include "computation.h"

void cmp_ez(
	E e, H h,
	EXY exy, EYX eyx, COE c, src s, int time)
{
	int i, j;
	//
	float kappa;
	float xcoor, ycoor;

	for (i = 0; i < e.ez->height; i++)
	{
		for (j = 0; j < e.ez->width; j++)
		{
			//cout <<"time: "<<time<<", ez("<<i<<", "<<j<<") judge: ";
			//if the point is on the boundaries
			if (i == 0 || i == (e.ez->height - 1)
				|| j == 0 || j == (e.ez->width - 1))
			{
				e.ez->p.at(i*e.ez->width + j) = 0.f;
				//cout << "ez boundary" << endl;
				continue;
			}

			//transfer the indexes of vector of the coordinates in field
			xcoor = j;
			ycoor = i;
			kappa = c.set_kappa(s, ycoor, xcoor);
			//computation
			e.ez->p.at(i*e.ez->width + j) +=
				e.coe_E*(
				(h.hy->p.at(i*h.hy->width + j) - h.hy->p.at((i - 1)*h.hy->width + j)) / (kappa*s.dz)
				- (h.hx->p.at(i*h.hx->width + j) - h.hx->p.at(i*h.hx->width + j - 1)) / (kappa*s.dz)
				)
				+ e.coe_E*(exy.full->p.at(i*exy.full->width + j) - eyx.full->p.at(i*eyx.full->width + j));
			//cout << "get in computation" << endl;
			if (0&&time == 0){
				cout << "ez(" << xcoor << ", " << ycoor << "): " << e.ez->p.at(i*e.ez->width + j) << endl;
			}
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
			//computation
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
	for (i = 0; i < h.hy->height; i++)
	{
		for (j = 0; j < h.hy->width; j++)
		{
			//transfer the indexes of vector to the coordinates in field
			xcoor = j;
			ycoor = i + 0.5f;
			h.hy->p.at(i*h.hy->width + j) +=
				h.coe_h*(
				e.ez->p.at((i + 1)*e.ez->width + j) - e.ez->p.at(i*e.ez->width + j)
				) / (c.set_kappa(s, ycoor, xcoor)*s.dz)
				+ h.coe_h*(hxz.full->p.at(i*hxz.full->width + j));
		}
	}
}