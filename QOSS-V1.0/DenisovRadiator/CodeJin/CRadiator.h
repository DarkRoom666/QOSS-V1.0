#pragma once
//define Radiator Structure
//�����������FDTD������ģ�ͣ�
//length Structure
#include<omp.h>
#include<cmath>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cstdio>
#include<cstdlib>
#include<complex>
#include"../DenisovRadiator/CodeJin/GEMS_Memory.h"
#include"../DenisovRadiator/CodeJin/Global_Vars.h"

using namespace::std;

#ifndef CRADIATOR_H  // ȫ��д
#define CRADIATOR_H

class CRadiator {

	public:
		CRadiator(void);
		CRadiator(
		int _Kind = 1,	//1�ǵͽ�TEģʽ;2�Ǹ߽�TEģʽ��3�Ǿ��β���ģʽ
		int _Type = 1,	//ģʽ TE 1 ��TM 2��
		int _Rotation = 0,	//�������� ����Ϊ1������Ϊ2��
		int _m = 0,
		int _n = 1,
		int Nx = 100,
		int Ny = 100,
		float _lambda = 3e-3,	//����
		float _Radius = 20e-3,	//Բ�����뾶
		float _Lc = 20e-3,		//Բ��ģ�пڸ߶�
		float _Lp = 10e-3,		//Բ��ģ�����θ߶�
		float _wa = 8e-3,		//���β�������
		float _wb = 4e-3,		//���β����̱�
		float _Hp = 8e-3,		//�����θ߶�
		float _Ha = 20e-3,		//����θ߶�
		float _aa = 20e-3,		//������泤��
		float _bb = 15e-3		//�������̱�
	);//Ĭ�Ϲ��캯��

	~CRadiator(void) {
		// ���з������ڴ�ǵ��ͷ�
		// if(p) {
		//   delete p;
		//   p = NULL;
		// }	
		//position
		Free_1D(px); px = NULL;
		Free_1D(py); py = NULL;
		Free_1D(pz); pz = NULL;
		Free_3D(EpsMap, Nz_model);	EpsMap = NULL;
		Free_3D(Esig, Nz_model);	Esig = NULL;
	};
	//����ͽ׷���������
	void ResetRadiator_RoundL(float _lambda, float _Radius, int _Nx, int _Ny, int _N_spa,float _Lc, float _Lp);
	//���õͽ׷���������
	void SetFirstMirror_RoundL(float _F);
	//����߽׷���������
	void ResetRadiator_RoundH(int _Rotation,float _lambda, float _Radius, int _Nx, int _Ny, int _N_spa, float _Lc, float _Lp);

	//ȷ���ͽ�Բ�����������
	void SetRadiator_RoundL( float _Lc, float _Lp);
	//ȷ���߽�Բ�����������
	void SetRadiator_RoundH(int _Rotation, float lambda, float Radius, float Lc, float Lp);
	
	
	//���������������(�ǹ��Ͳ��֣��ͽ�Բ�������)
	void GenerateCellArray_RoundL(void);
	//���������������(�ǹ��Ͳ��֣��߽�Բ�������)
	void GenerateCellArray_RoundH(float _phic);
	//��������������ݵ��м亯��
	float RFxy(float x, float y, float _phic);
	//����������������м亯��
	bool RPECStructure(float x, float y, float z, float _phic);

	//���������������(�ǹ��Ͳ��֣����β�������)
	void GenerateCellArray_Rec(void);
	//���������������(�ǹ��Ͳ��֣����η��������в���)
	void GenerateCellArray_RecArray(int Nunit);

	//ȷ���ز�������������
	void SetRadiator_Rectan(int Type, int m, int n, float lambda, float wa, float wb, float Hp, float Ha, float aa, float bb );
	//����Logfile
	void SetLogfile(fstream* logAdd);
private:	//data

	int Kind;		//��������ʽ
	int Type;		//����������ģʽ	
	int Rotation;	//Բ����������
	int m;
	int n;
	int Nx;			//����ģʽ��X�����������
	int Ny;			//����ģʽ��Y�����������
	float lambda;	//�粨��
	float Radius;	//Բ�����뾶
	float Lc;		//Բ��ģ�пڸ߶�
	float Lp;		//Բ��ģ�����θ߶�
	float phic;
	float wa;		//���β�������
	float wb;		//���β����̱�
	float Hp;		//�����θ߶�
	float Ha;		//����θ߶�
	float aa;		//������泤��
	float bb;		//�������̱�
	float ds;		//�ʷ־���

	fstream* Logfile;

public:
	int Nx_model;	//Model Size
	int Ny_model;	//Model Size	
	int Nz_model;	//Model Size

	int Px_model;	//Start Position of Excitation Port
	int Py_model;	//Start Position of Excitation Port
	int Pz_model;	//Start Position of Excitation Port
	//position
	float* px;
	float* py;
	float* pz;
	//uniform mesh	//����������ܲ�һ��Ŷ��Ҫ���ֿ�
	float dx;		//
	float dy;
	float dz;
	//non uniform mesh �ݲ���
	float* dx_nu;
	float* dy_nu;
	float* dz_nu;

	int*** EpsMap;	//Mark	
	float*** Esig;	//Passed into non-conformal FDTD

	int N_spa;

	float MirrorHeight; //First Mirror Height

	//define geometry functions

};
#endif // TRACELIGHT_H

