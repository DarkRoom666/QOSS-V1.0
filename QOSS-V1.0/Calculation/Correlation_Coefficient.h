#pragma once

//����һ�����㳡�ֲ���ʸ���������˹���ϵ����ר�ú���
//����Ϊ���������򳡷���EX��EY��������ʽ��
//ÿ�μ���һ��Ҫ��˳������EX��EY
//by XD 20180321

#ifndef CORRELATION_COEFFICIENT_H
#define CORRELATION_COEFFICIENT_H

#include <iostream>
#include <cmath>
#include <complex>//�����࣬C++�Դ�
#include <string>
#include <vector>
#include "../util/Vector3.h"
#include "../util/Constant_Var.h"
#include "Gaussian_Beam_Circular.h"

using namespace std;

//���Լ������ƥ���˹���ϵ��ʱ�Ķ�Ӧ����±��Լ���˹��������
void Correlation_Coefficient_Omega(double frequency0, double ds,int N0,
	vector <vector <complex<double>>> EX,
	vector <vector <complex<double>>> EY,//N0*N0
	int &index, int &MaxE_i, int &MaxE_j, double &W0);
//indexΪ0,������ΪEX��indexΪ1,������ΪEY
//MaxE_i��MaxE_j�����������ǿ�����±�
//WO�Ǿ���ƥ��õ��ĸ�˹��������

//���Լ������ƥ���˹�������ϵ��
double Scalar_Correlation_Coefficient(double frequency0, double ds, int N0,
	vector <vector <complex<double>>> EX,
	vector <vector <complex<double>>> EY);

//���Լ������ƥ���˹ʸ�����ϵ��
double Vector_Correlation_Coefficient(double frequency0, double ds, int N0,
	vector <vector <complex<double>>> EX,
	vector <vector <complex<double>>> EY);

//���������Ŀ������W0�ĸ�˹�������ϵ��
double Scalar_Correlation_Coefficient_TargetW0(double frequency0, double ds, int N0,
	vector <vector <complex<double>>> EX,
	vector <vector <complex<double>>> EY,
	double Target_W0);

//���������Ŀ������W0�ĸ�˹ʸ�����ϵ��
double Vector_Correlation_Coefficient_TargetW0(double frequency0, double ds, int N0,
	vector <vector <complex<double>>> EX,
	vector <vector <complex<double>>> EY,
	double Target_W0);

#endif
