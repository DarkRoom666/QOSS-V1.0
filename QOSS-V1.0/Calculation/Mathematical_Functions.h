#pragma once
//����һ������׳��Լ���ר�ú���
#ifndef MATHEMATICAL FUNCTIONS_H
#define MATHEMATICAL FUNCTIONS_H

#include <iostream>
#include <cmath>
#include <complex>//�����࣬C++�Դ�
#include <vector>//vector���ڴ���������������  
#include <string>

#include "../util/Constant_Var.h"

using namespace std;

double fact(int n);//�Ѳ��ԣ��׳˺���û����


//����M��Bessel������һ�׵���
double D1_J(int m, double z);


//���س��ֲ�ֵH_phi(������ʽ�ģ�Ĭ��Ϊ����)
complex<double> H_phi(double Kz0, double Kr0, int m0, double r0, double Phi0);

	
complex<double> H_r(double Kz0, double Kr0, int m0, double r0, double Phi0);


//���س��ֲ�ֵE_phi(������ʽ)
complex<double> E_phi(double Kz0, double Kr0, int m0, double r0,
	                  double Phi0, double Omega0, double Beta0);

	
//���س��ֲ�ֵE_r(������ʽ)
complex<double> E_r(double Kz0, double Kr0, int m0, double r0,
	double Phi0, double Omega0, double Beta0);


//���س��ֲ�ֵH_z(������ʽ)
complex<double> H_z(double Kr0, int m0, double r0, double Phi0);
	

//������������TE���ֲ�

//���س��ֲ�ֵH_phi_LeftRo(������ʽ)
complex<double> H_phi_LeftRo(double Kz0, double Kr0, int m0, double r0, double Phi0);


//���س��ֲ�ֵH_r_LeftRo(������ʽ)
complex<double> H_r_LeftRo(double Kz0, double Kr0, int m0, double r0, double Phi0);

	
//���س��ֲ�ֵE_phi_LeftRo(������ʽ)
complex<double> E_phi_LeftRo(double Kz0, double Kr0, int m0, double r0,
	                         double Phi0, double Omega0, double Beta0); 
	
//���س��ֲ�ֵE_r_LeftRo(������ʽ)
complex<double> E_r_LeftRo(double Kz0, double Kr0, int m0, double r0,
	                       double Phi0, double Omega0, double Beta0);

//���س��ֲ�ֵH_z_LeftRo(������ʽ)
complex<double> H_z_LeftRo(double Kr0, int m0, double r0, double Phi0);


#endif // 
//20171026�Ѿ���֤