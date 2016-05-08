#include <iostream>
#include "H.h" 
#include "E.h"
#include "cvl.h"
#include "source.h"

using namespace std;

void cmp_ez(area *e, area *hy, area *hx, area *exy, area *eyx, COE, src, int);

void cmp_hx(area *hx, area *e, area *hyz, COE, src);

void cmp_hy(area *hy, area *e, area *hxz, COE, src);

void cmp_exy(area *exy, area *hy, COE, src);

void cmp_eyx(area *eyx, area *hx, COE, src, int);

void cmp_hyz(area *hyz, area *e, COE, src, int);

void cmp_hxz(area *hxz, area *e, COE, src);