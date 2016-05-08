#include "computation.h"

void cmp_ex(
	E e, H h,
	area *ezy, COE c, src s, int time)
{
	int j;
	float xcoor, kappa;
	area *local = e.ex;

	for (j = 0; j < local->length; j++){
		//transfer the indexes of vector of the coordinates in field
		xcoor = j;

		//FDTD area
		//if the point is on the boundaries
		if (j == 0 || j == (local->length - 1)){
			local->p.at(j) = 0.f;
			continue;
		}
		kappa = c.set_kappa(s, xcoor);
		local->p.at(j) +=
			(s.dt / (e.epsilon*kappa*s.dz))*(
			-(h.hy->p.at(j) - h.hy->p.at(j - 1))
			)
			+ (s.dt / e.epsilon)*(-ezy->p.at(j));
	}
}

void cmp_hy(H h, E e, area *hzx, COE c, src s, int time)
{
	int j;
	float xcoor, kappa;
	area *local;
	local = h.hy;
	for (j = 0; j < local->length; j++)
	{
		//transfer the indexes of vector to the coordinates in field
		xcoor = j + 0.5;
		kappa = c.set_kappa(s, xcoor);
		local->p.at(j) +=
			-(s.dt / (h.mu*kappa*s.dz))*(
			e.ex->p.at(j + 1) - e.ex->p.at(j)
			)
			- (s.dt / h.mu)*(hzx->p.at(j));
	}
}

void cmp_hzx(E e, area *hzx, COE c, src s, int time)
{
	int j;
	float xcoor, kappa;
	area *local;
	local = hzx;
	for (j = 0; j < local->length; j++)
	{
		xcoor = j + 0.5;
		float past = local->p.at(2);
		local->p.at(j) = c.set_c(s, xcoor)*(
			e.ex->p.at(j + 1) - e.ex->p.at(j)
			) / s.dz
		+c.set_coe(s, xcoor)*local->p.at(j);
		if (time == 7){
			if (j == 2){
				cout << "hzx(" << j << "): " << local->p.at(j) << endl;
				cout << "1 part: "
					<< c.set_coe(s, xcoor)*past << "\t"
					<< c.set_coe(s, xcoor) << "\t"
					<< past << "\t"
					<< "part 2: "
					<< endl;
			}
		}
	}
}

void cmp_ezy(H h, area *exy, COE c, src s)
{
	int i, j;
	float xcoor;
	area *local = exy;
	for (j = 0; j < local->length; j++){
		if (j == 0 || j == (local->length - 1)){
			if (j == 0){
				local->p.at(j) = local->p.at(j + 1);
			}
			if (j == (local->length - 1)){
				local->p.at(j) = local->p.at(j - 1);
			}
			continue;
		}
		xcoor = j;
		local->p.at(j) = c.set_c(s, xcoor)*(
			h.hy->p.at(j) - h.hy->p.at(j - 1)
			) / s.dz +
			c.set_coe(s, xcoor)*local->p.at(j);
	}
}
