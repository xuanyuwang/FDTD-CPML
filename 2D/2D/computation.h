#include <iostream>
#include "H.h" 
#include "E.h"
#include "cvl.h"
#include "source.h"

using namespace std;

void cmp_ez(E, H, EXY, EYX, COE, src,int);

void cmp_hx(H, E, HYZ, COE, src);

void cmp_hy(H, E, HXZ, COE, src);

void cmp_exy(EXY, H, EXY, COE, src);

void cmp_eyx(EYX, H, EYX, COE, src);

void cmp_hyz(HYZ, E, HYZ, COE, src);

void cmp_hxz(HXZ, E, HXZ, COE, src);