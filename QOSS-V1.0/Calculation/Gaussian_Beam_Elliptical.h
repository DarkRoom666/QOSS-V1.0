#pragma once
//����һ��������Բ��˹�����ֲ���ר�ú���
//�ø�˹��������Z�������򴫲�������λ����Z=0��
//by XD 20180305

#ifndef GAUSSIAN_BEAM_ELLIPTICAL_H
#define GAUSSIAN_BEAM_ELLIPTICAL_H

#include <iostream>
#include <cmath>
#include <complex>//�����࣬C++�Դ�
#include <string>
#include "../util/Constant_Var.h"
#include "../util/Vector3.h"

using namespace std;

//���ز�������Z����������ֵW(Wx��Wy��ʽһ��)
double Gauss_Omega_Elliptical(double frequency0, double w0, double z0);

//���ز�������Z���������ʰ뾶R(Wx��Wy��ʽһ��)
double Gauss_R_Elliptical(double frequency0, double w0, double z0);

//������λֵ
double Phi_Elliptical(double frequency0, double w0x, double w0y, double z0);

//����Բ��˹�����������P��x,y,z���ĸ���������ֵ
complex<double> Gauss_E_Elliptical(double frequency0, double w0x, double w0y, Vector3 P);

#endif