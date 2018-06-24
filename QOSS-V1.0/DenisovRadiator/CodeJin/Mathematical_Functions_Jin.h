#pragma once
//����һ������׳��Լ���ר�ú���

#ifndef MATHEMATICAL_FUNCTIONS_JIN_H
#define MATHEMATICAL_FUNCTIONS_JIN_H

#include <iostream>
#include <cmath>
#include <complex>//�����࣬C++�Դ�
#include <vector>//vector���ڴ���������������  
#include <string>
#include "../util/Constant_Var.h"

using namespace std;

//Add by Jin Ming
double rootdbessel(int m, int n);	//Find the nth root of mth order derivate bessel function of the first kind

double rootbessel(int m, int n);	//Find the nth root of mth order bessel function of the first kind

double kzmnTE(int m, int n, double lambda, double radius);

#endif // 
//20171026�Ѿ���֤