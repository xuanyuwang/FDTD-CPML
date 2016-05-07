#include <iostream>
#include "H.h" 
#include "E.h"
#include "cvl.h"
#include "source.h"

using namespace std;

void cmp_ez(E, H, EXY, EYX, COE, src,int);

void cmp_hx(H, E, HYZ, COE, src);

void cmp_hy(H, E, HXZ, COE, src);

void cmp_exy(H, EXY, COE, src);

void cmp_eyx( H, EYX, COE, src);

void cmp_hyz( E, HYZ, COE, src);

void cmp_hxz( E, HXZ, COE, src);