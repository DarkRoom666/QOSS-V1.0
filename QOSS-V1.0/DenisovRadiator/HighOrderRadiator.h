#pragma once
//����ļ����������߽�Բ��ģ��������һ�����ĵ�������
//����ļ����ɸ������㷴���ĺ���
//����ļ����ɶ���߽�Բ��ģʽ�Ķ˿ڹ�����
//����-�Ķ�-��д��ά��
//��λm
//2018��6��
//jinmingaps@163.com
//���ȵ�λ m
//�Ƕȵ�λ rad
//ע�� Ĭ���п�λ�õĸ߶���0�߶�λ�ã�ˮƽλ��������ԭ��
//ע�� Ĭ��������

#ifndef HIGHORDERRADIATOR_H
#define HIGHORDERRADIATOR_H
#include <cmath>
#include <vector>
#include "../util/Vector3.h"
#include "../DenisovRadiator/CodeJin/OMTfunctions.h"
#include "../DenisovRadiator/CodeJin/Mathematical_Functions_Jin.h"

using namespace std;
class HighOrderRadiator {

public:

	HighOrderRadiator(void);
	~HighOrderRadiator(void);

	//ע�⣺
	//�������- ���ط������������ڽ�ģ��
	// ���ã� SetRadiatorBasicParas(3) + SetRadiatorTurbulenceParas(10); ���ܣ�GetRadiatorLattice1() + GetRadiatorLattice2();

	//�������- ����һ�����������ڽ�ģ��
	// ���ã� SetFirstMirrorParas(5); ���ܣ�GetFirstMirrorLattice();

	//�������- ����������ڲ�����׷��
	// ���ã� SetRadiatorBasicParas(4) + SetRadiatorTurbulenceParass(10); ���ܣ�CalculateRefiectionPoint();

	//�������- ���ع���׷���ö˿���ʼ����
	// ���ã� SetModeParas(5)�����ܣ�GetModeRays();

	//�趨����
	void SetRadiatorBasicParas(double _Ra, double _zcut, double _phic);
	void SetRadiatorBasicParas(double _Ra, double _lcut, double _zcut, double _phic);//����
	void SetRadiatorTurbulenceParas(double _delBeta1, double _delBeta2, double _z0_1, double _z0_2, double _ls_1, double _ls_2, double _lc_1, double _lc_2, double _mag1, double _mag2);
	//����һ��������
	void SetFirstMirrorParas(double _F1, double _F2, double _Rc, double _Phi, double _Lcut);
	//���ù���׷�ٵ�Դ
	void SetModeParas(double _freq, double _m, double _n, double _Ra, double _zcut);

	//����һ�����ĵ��� ��ά����
	int GetFirstMirrorLattice(vector<vector<Vector3>> &_lattice, double _ds);
	//���ط������ĵ��� һά���� ������ϳ�һ��ԲͰ
	int GetRadiatorLattice1(vector<Vector3> &_lattice, double _ds);
	int GetRadiatorLattice2(vector<Vector3> &_lattice, double _ds);
	//���ط���������˿ڴ��Ĺ���
	int GetModeRays(vector<Vector3> &_pos, vector<Vector3> &_dir, int _Nr, int _Nphi);
	
	//����׷���ú��� 
	inline double R(double _phi, double _t);
	inline Vector3 Position(double _phi, double _t);
	inline Vector3 DelPhiVec(double _phi, double _t);
	inline Vector3 DelTVec(double _phi, double _t);
	inline Vector3 Normal(double _phi, double _t);
	//���뷶Χ 0~2pi
	bool OnRadiatorSurface(double _phi, double _t);
	
	//XD��Ӽ��㷴��㺯��
	bool CalculateReflectionPoint(Vector3 StartPoint, Vector3 IncidentVector, Vector3 &ReflectionPoint);

private:
	//������һ��������ļ��β���
	double F1;
	double F2;
	double Rc;
	double Phi;
	double Lcut;
	//�����Ƿ���������ļ��β���
	double Ra;
	double delBeta1;
	double delBeta2;
	int l1; 
	int l2;
	double z0_1;
	double z0_2;
	double ls_1;
	double ls_2;
	double lc_1;
	double lc_2;
	double mag1;
	double mag2;
	double Zcut;
	double phic;

	double frequency;
	int m;
	int n;
};

#endif


