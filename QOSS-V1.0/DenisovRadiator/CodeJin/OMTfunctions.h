#pragma once

#ifndef OMTFUNCTIONS_H

#define OMTFUNCTIONS_H

#include<complex>
#include<cmath>
#include"../util/Constant_Var.h"


using namespace std;

double sig(double _mag, double _z0, double _lc, double _ls, double _z);
double dsig(double _mag, double _z0, double _lc, double _ls, double _z);

complex<double> GAMMA0(double _beta1, double _beta2, int _l1, int _l2, int _m1, int _m2, double _z, double _dis1, double _dis2);

complex<double> GAMMAZ(double _beta1, double _beta2, int _l1, int _l2, int _m1, int _m2, double _z, double _a, double _dis1, double _ddis1, double _dis2, double _ddis2);

complex<double> CQR(double freq, double _a, double _beta1, double _beta2, int _l1, int _l2, double _z, double _dis1, double _ddis1, double _dis2, double _ddis2, double _kz1, double _kz2, double _Chi1, double _Chi2, int _m1, int _m2);


#endif

