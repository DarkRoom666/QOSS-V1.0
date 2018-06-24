#pragma once
//这个文件给出构建高阶圆电模辐射器和一级镜的点云数据
//金铭-编写和维护
//单位m
//2018年6月
//jinmingaps@163.com
//长度单位 m
//角度单位 rad
//注意 默认切口位置的高度在0高度位置，水平位置在坐标原点
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
	//设定参数
	void SetRadiatorBasicParas(double _Ra, double _zcut, double _phic);
	void SetRadiatorTurbulenceParas(double _delBeta1, double _delBeta2, double _z0_1, double _z0_2, double _ls_1, double _ls_2, double _lc_1, double _lc_2, double _mag1, double _mag2);
	void SetFirstMirrorParas(double _F1, double _F2, double _Rc, double _Phi, double _Lcut);

	//返回一级镜的点阵 二维数组
	int GetFirstMirrorLattice(vector<vector<Vector3>> &_lattice, double _ds);
	//返回辐射器的点阵 一维数列 两个面合成一个圆桶
	int GetRadiatorLattice1(vector<Vector3> &_lattice, double _ds);
	int GetRadiatorLattice2(vector<Vector3> &_lattice, double _ds);
	
	//光线追踪用函数 
	inline double R(double _phi, double _t);
	inline Vector3 Position(double _phi, double _t);
	inline Vector3 DelPhiVec(double _phi, double _t);
	inline Vector3 DelTVec(double _phi, double _t);
	inline Vector3 Normal(double _phi, double _t);
	//输入范围 0~2pi
	bool OnRadiatorSurface(double _phi, double _t);
private:
	//这里是一级镜所需的几何参数
	double F1;
	double F2;
	double Rc;
	double Phi;
	double Lcut;
	//这里是辐射器所需的几何参数
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


