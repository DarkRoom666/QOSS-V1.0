#pragma once
//����һ������λ�ֲ����������λ�����ĳ��򣬲���Matlab����Version5.m�еĲ���
//������λ��λ�ǻ���
//by XD 20180308
//20180310ͨ������

#ifndef PHASEUNWRAP_H
#define PHASEUNWRAP_H

#include <iostream>
#include <cmath>
#include <math.h>
#include <string>
#include "../util/Constant_Var.h"
#include <vector>//vector���ڴ���������������
#include "../util/Vector3.h"


using namespace std;

//һά��λ���۵�����ά��λ���۵�������Ҫ�õ�����
//�����Phase�������۵����ģ���arg����atan2�����õ�����Χ�ǣ�-Pi��Pi]
vector <double> Unwrap_1D(int N0,vector <double> Phase);

//��ά��λ���۵�����,��N0*N0��λ�����������
vector <vector<double>> Unwrap_2D(int N0, vector <vector <double>> Phase);

#endif
