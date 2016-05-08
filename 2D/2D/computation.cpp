#include "computation.h"

void cmp_ez(
	area *e, area *hy, area *hx,
	area *exy, area *eyx, COE c, src s, int time)
{
	int i, j;
	float xcoor, ycoor, kappa;
	float eps = e->epsilon;

	for (i = 0; i < e->height; i++){
		for (j = 0; j < e->width; j++){
			//transfer the indexes of vector of the coordinates in field
			xcoor = j;
			ycoor = i;

			//if the point is on the boundaries
			if (i == 0 || i == (e->height - 1)
				|| j == 0 || j == (e->width - 1)){
				e->p.at(i*e->width + j) = 0.f;
				continue;
			}
			kappa = c.set_kappa(s, ycoor, xcoor);
			e->p.at(i*e->width + j) +=
				(s.dt / (eps*kappa*s.dz))*(
				(hy->p.at(i*hy->width + j) - hy->p.at((i - 1)*hy->width + j))
				- (hx->p.at(i*hx->width + j) - hx->p.at(i*hx->width + j - 1))
				) +
				(s.dt / eps)*(exy->p.at(i*exy->width + j) - eyx->p.at(i*eyx->width + j));
			if (time == 6){
				if (i == 2 && j == 8){
					cout << e->p.at(i*e->width + j) << endl;
					cout << "part 1: " << (s.dt / (eps*kappa*s.dz))*(
						(hy->p.at(i*hy->width + j) - hy->p.at((i - 1)*hy->width + j))
						- (hx->p.at(i*hx->width + j) - hx->p.at(i*hx->width + j - 1))
						) << "\t"
						<< (s.dt / eps)<<"\t"
						<<(exy->p.at(i*exy->width + j) - eyx->p.at(i*eyx->width + j))
						<< endl;
				}
			}
		}
	}
}

void cmp_hx(area *hx, area* e, area *hyz, COE c, src s)
{
	int i, j;
	float xcoor, ycoor, kappa, mu = hx->mu;

	for (i = 0; i < hx->height; i++)
	{
		for (j = 0; j < hx->width; j++)
		{
			//transfer the indexes of vector of the coordinates in field
			xcoor = j + 0.5f;
			ycoor = i;
			kappa = c.set_kappa(s, ycoor, xcoor);
			hx->p.at(i*hx->width + j) = hx->p.at(i*hx->width + j) -
				(s.dt / (mu*kappa*s.dz))*(
				e->p.at(i*e->width + j + 1) - e->p.at(i*e->width + j)
				)
				+ (s.dt / mu)*(-hyz->p.at(i*hyz->width + j));
		}
	}
}

void cmp_hy(area *hy, area *e, area *hxz, COE c, src s)
{
	int i, j;
	float xcoor, ycoor, kappa, mu = hy->mu;

	for (i = 0; i < hy->height; i++)
	{
		for (j = 0; j < hy->width; j++)
		{
			//transfer the indexes of vector to the coordinates in field
			xcoor = j;
			ycoor = i + 0.5f;
			kappa = c.set_kappa(s, ycoor, xcoor);
			hy->p.at(i*hy->width + j) +=
				(s.dt / (mu*kappa*s.dz))*(
				e->p.at((i + 1)*e->width + j) - e->p.at(i*e->width + j)
				)
				+ (s.dt / mu)*(hxz->p.at(i*hxz->width + j));
		}
	}
}

void cmp_hxz(area *hxz, area* e, COE c, src s)
{
	int i, j;
	float xcoor, ycoor, kappa;

	for (i = 0; i < hxz->height; i++)
	{
		for (j = 0; j < hxz->width; j++)
		{
			xcoor = j;
			ycoor = i + 0.5f;
			kappa = c.set_kappa(s, ycoor, xcoor);
			hxz->p.at(i*hxz->width + j) =
				c.set_c(s, ycoor, xcoor)*(
				e->p.at((i + 1)*e->width + j) - e->p.at(i*e->width + j)
				) / s.dz
				+ c.set_coe(s, ycoor, xcoor)*hxz->p.at(i*hxz->width + j);
		}
	}
}

void cmp_hyz(area *hyz, area *e, COE c, src s, int time)
{
	int i, j;
	float xcoor, ycoor, kappa;

	float past;
	for (i = 0; i < hyz->height; i++){
		for (j = 0; j < hyz->width; j++){
			xcoor = j + 0.5f;
			ycoor = i;
			kappa = c.set_kappa(s, ycoor, xcoor);
			past = hyz->p.at(i*hyz->width + j);
			hyz->p.at(i*hyz->width + j) = c.set_c(s, ycoor, xcoor)*(
				e->p.at(i*e->width + j + 1) - e->p.at(i*e->width + j)
				) / s.dz +
				c.set_coe(s, ycoor, xcoor)*hyz->p.at(i*hyz->width + j);
		}
	}
}

void cmp_exy(area *exy, area *hy, COE c, src s)
{
	int i, j;
	float xcoor, ycoor;

	for (i = 0; i < exy->height; i++){
		for (j = 0; j < exy->width; j++){
			if (i == 0 || i == (exy->height - 1)
				|| j == 0 || j == (exy->width - 1)){
				exy->p.at(i*exy->width + j) = 0.f;
				continue;
			}
			xcoor = j;
			ycoor = i;
			exy->p.at(i*exy->width + j) =
				c.set_c(s, ycoor, xcoor)*(
				hy->p.at(i*hy->width + j) - hy->p.at((i - 1)*hy->width + j)
				) / s.dz +
				c.set_coe(s, ycoor, xcoor)*exy->p.at(i*exy->width + j);
		}
	}
}

void cmp_eyx(area *eyx, area *hx, COE c, src s, int time)
{
	int i, j;
	float xcoor, ycoor;

	for (i = 0; i < eyx->height; i++){
		for (j = 0; j < eyx->width; j++){
			if (i == 0 || i == (eyx->height - 1)
				|| j == 0 || j == (eyx->width - 1)){
				eyx->p.at(i*eyx->width + j) = 0.f;
				continue;
			}
			xcoor = j;
			ycoor = i;
			if (i*eyx->width + j + 1 >= eyx->grid_num){
				cout << i << "\t" << j << endl;
			}
			eyx->p.at(i*eyx->width + j) =
				c.set_c(s, ycoor, xcoor)*
				(hx->p.at(i*hx->width + j) - hx->p.at(i*hx->width + j - 1)) / s.dz
				+ c.set_coe(s, ycoor, xcoor)*(eyx->p.at(i*eyx->width + j));
		}
	}
}