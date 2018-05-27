#pragma once
//����һ��������ο����ų��ֲ��Ĵ������򣬲��������������ɢ��˹�����ֽ�����㷨
//����Ϊ���γ��ֲ������Ըþ��в����Ǹ��ֲ�ƽ����Ϊ�ֲ�����ϵ
//Ĭ��������Ey����������ʽ����������
//���������󵥵����������ֲ��Ͷ����ĳ��ֲ�
//by XD 20180303
//20180307������ȷ��by XD
//20180314�ٴβ���б��˹��������ȷ,by XD
//test
#ifndef BEAMPROPAGATION_H
#define BEAMPROPAGATION_H

#include <iostream>
#include <cmath>
#include <math.h>
#include <complex>//�����࣬C++�Դ�
#include <string>
#include "../util/Constant_Var.h"
#include <vector>//vector���ڴ���������������
#include "../util/Vector3.h"
#include "Gaussian_Beam_Elliptical.h"


using namespace std;

//�ֽ���泡�����طֽ����Actual_SplitTimes�ͷֽ���ĸ�˹������Ϣ����Split_Info
void Field_Split(double frequency0, double ds,int N0,
	                        vector <vector <complex <double>>> &E0,
	                        int &Actual_SplitTimes, vector <vector <double>> &Split_Info);

//��������Field_Split�����õ�Actual_SplitTimes��Split_Info�Ļ����ϣ�
//������泡������������P�ĳ�ǿ
complex<double> Calculate_SinglePoint(double frequency0, double ds,int N0,
	                        int Actual_SplitTimes, vector <vector <double>> &Split_Info,
	                        Vector3 P);

#endif
