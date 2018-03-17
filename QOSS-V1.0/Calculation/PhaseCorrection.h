#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>//��������������
#include "Gaussian_Beam_Circular.h"//Բ��˹��������
#include "BeamPropagation.h"//���泡��ż������ʹ�ø�˹�����ֽ�������㷨
//vector���ڴ��������������ݣ�����BeamPropagation.h������
#include "PhaseUnwrap.h"//��λ������۵�����
#include "../util/Vector3.h"
//Vector3.h���Ѿ�������Vector3.h��
//����������ļ����Ѿ����ù�"Constant_Val.h"�������治��Ҫ�����include������
//#include "Constant_Val.h"//��ѧ����H�ļ�

using namespace std;
/*
*	created by XD 2018/02/05
*   ��λĬ��Ϊ��
*   function ������λ�������ļ�
*/

#ifndef PHASECORRECTION_H // ȫ��д
#define PHASECORRECTION_H

namespace calculation//������Calculation�����ռ���

{ // ����ռ� ���������ظ�

	class PhaseCorrection
	{
	public:

		PhaseCorrection
			(double _Frequency = 44.0e9,//Ƶ��
				double _Lamda = C_Speed / 44.0e9,//����
				double _Target_W0 = 0.01,//Ŀ���˹�������
				int _N_Mirror = 101,//����N*N�����������
				int _N_InField = 101,//������λ������������泡����ɢ��N_Input_Field*N_Input_Field
				double _InField_ds=0.001
				);

		~PhaseCorrection();

		//�������ã������У�
		void Set(double _Frequency, double _Target_W0, int _N_Mirror, int _N_InField,
			vector <vector <Vector3>> &_Initial_Mirror,
			vector <vector <complex<double>>> &_InField_EX,
			vector <vector <complex<double>>> &_InField_EY,
			double _InField_ds);

		//��������ϵ�������ã��ڴ�֮ǰҪ������Set������
		void SetCoordinate(Vector3 _InField_Position, Vector3 _InField_X, Vector3 _InField_Y, Vector3 _InField_Z,
			Vector3 _Mirror_Position, Vector3 _Mirror_X, Vector3 _Mirror_Y, Vector3 _Mirror_Z,
			Vector3 _Reflect_Point, Vector3 _Reflect_Point_Normal);

		//������Ҫ�����ľ���
		void SetMirror(vector <vector <Vector3>> &_Mirror );

		//��������
		vector <vector <Vector3>> Single_Correction();


	private:
		double Frequency;//Ƶ��
		double Lamda;//����
		double Target_W0;//Ŀ���˹����������
		double D_in;//������泡��������ľ��루���ĵ㵽���ĵ㣩
		double D_out;//������ڵ�������ľ��루���ĵ㵽���ĵ㣩
		//�����ֶ�������Ϊ����ʱ��Ҫ����ȷ��
		int N_Mirror;//���������������N_Mirror*N_Mirror��
	    int N_InField;//������泡������������N_Input_Field*N_Input_Field��
		vector <vector <Vector3>> Initial_Mirror;//�����ʼ����ĵ����
		//�����ĵ����궼�Ǿֲ�����ϵMirror_X, Mirror_Y, Mirror_Z�µ�
		vector <vector <complex<double>>> InField_EX,InField_EY;
		//�������볡�ĳ��ֲ�(InField_X��InField_Y�������)����-x,-y��x,y��
		vector <vector <complex<double>>> InField_EX1,InField_EY1;
		//����������볡�ĳ��ֲ�(InField_X1��InField_Y1�������)
		double InField_ds;//�������볡�ĵ���
		Vector3 InField_Position, InField_X, InField_Y, InField_Z;//���볡����ϵ���ĵ�λ���Լ��ֲ�����ϵ
		Vector3 InField_X1, InField_Y1, InField_Z1;//���볡���ĵ�λ���Լ��ֲ�����ϵ��������ΪInField_Y1��
		Vector3 Mirror_Position, Mirror_X, Mirror_Y, Mirror_Z;//��������ϵ���ĵ�λ���Լ��ֲ�����ϵ
		Vector3 OutField_Position, OutField_X, OutField_Y, OutField_Z;//���䳡������ϵ���ĵ�λ���Լ��ֲ�����ϵ
		Vector3 Reflect_Point, Reflect_Point_Normal;//���巴��㼰�õ�ķ���
		double Theta;//�������ڷ����淴����ϵķ���Ƕ�

		//��������ϵ�ı�ʾ����ȫ������ϵ����

		

	};
}

#endif // TRACELIGHT_H