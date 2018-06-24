#pragma once
//����ļ����������߽�Բ��ģ��������һ�����ĵ�������
//����-��д��ά��
//��λm
//2018��6��
//jinmingaps@163.com
//���ȵ�λ m
//�Ƕȵ�λ rad
//ע�� Ĭ���п�λ�õĸ߶���0�߶�λ�ã�ˮƽλ��������ԭ��
#ifndef HIGHORDERRADIATOR_H
#define HIGHORDERRADIATOR_H
#include <cmath>
#include <vector>
#include "../util/Vector3.h"
#include "../DenisovRadiator/CodeJin/OMTfunctions.h"

using namespace std;
class HighOrderRadiator {

public:

	HighOrderRadiator(void);
	~HighOrderRadiator(void);
	//�趨����
	void SetRadiatorBasicParas(double _Ra, double _zcut, double _phic);
	void SetRadiatorTurbulenceParas(double _delBeta1, double _delBeta2, double _z0_1, double _z0_2, double _ls_1, double _ls_2, double _lc_1, double _lc_2, double _mag1, double _mag2);
	void SetFirstMirrorParas(double _F1, double _F2, double _Rc, double _Phi, double _Lcut);

	//����һ�����ĵ��� ��ά����
	int GetFirstMirrorLattice(vector<vector<Vector3>> &_lattice, double _ds);
	//���ط������ĵ��� һά���� ������ϳ�һ��ԲͰ
	int GetRadiatorLattice1(vector<Vector3> &_lattice, double _ds);
	int GetRadiatorLattice2(vector<Vector3> &_lattice, double _ds);
	
	//����׷���ú��� 
	inline double R(double _phi, double _t);
	inline Vector3 Position(double _phi, double _t);
	inline Vector3 DelPhiVec(double _phi, double _t);
	inline Vector3 DelTVec(double _phi, double _t);
	inline Vector3 Normal(double _phi, double _t);
	//���뷶Χ 0~2pi
	bool OnRadiatorSurface(double _phi, double _t);
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

};

#endif


