#pragma once
//这个文件给出构建高阶圆电模辐射器和一级镜的点云数据
//这个文件还可给出计算反射点的函数
//这个文件还可定义高阶圆电模式的端口光线束
//金铭-夏冬-编写和维护
//单位m
//2018年6月
//jinmingaps@163.com
//长度单位 m
//角度单位 rad
//注意 默认切口位置的高度在0高度位置，水平位置在坐标原点
//注意 默认右旋！

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

	//注意：
	//功能配对- 返回辐射器点阵用于建模：
	// 设置： SetRadiatorBasicParas(3) + SetRadiatorTurbulenceParas(10); 功能：GetRadiatorLattice1() + GetRadiatorLattice2();

	//功能配对- 返回一级镜点阵用于建模：
	// 设置： SetFirstMirrorParas(5); 功能：GetFirstMirrorLattice();

	//功能配对- 计算辐射器内部光线追踪
	// 设置： SetRadiatorBasicParas(4) + SetRadiatorTurbulenceParass(10); 功能：CalculateRefiectionPoint();

	//功能配对- 返回光线追踪用端口起始波束
	// 设置： SetModeParas(5)；功能：GetModeRays();

	//设定参数
	void SetRadiatorBasicParas(double _Ra, double _zcut, double _phic);
	void SetRadiatorBasicParas(double _Ra, double _lcut, double _zcut, double _phic);//重载
	void SetRadiatorTurbulenceParas(double _delBeta1, double _delBeta2, double _z0_1, double _z0_2, double _ls_1, double _ls_2, double _lc_1, double _lc_2, double _mag1, double _mag2);
	//设置一级镜参数
	void SetFirstMirrorParas(double _F1, double _F2, double _Rc, double _Phi, double _Lcut);
	//设置光线追踪的源
	void SetModeParas(double _freq, double _m, double _n, double _Ra, double _zcut);

	//返回一级镜的点阵 二维数组
	int GetFirstMirrorLattice(vector<vector<Vector3>> &_lattice, double _ds);
	//返回辐射器的点阵 一维数列 两个面合成一个圆桶
	int GetRadiatorLattice1(vector<Vector3> &_lattice, double _ds);
	int GetRadiatorLattice2(vector<Vector3> &_lattice, double _ds);
	//返回辐射器入射端口处的光线
	int GetModeRays(vector<Vector3> &_pos, vector<Vector3> &_dir, int _Nr, int _Nphi);
	
	//光线追踪用函数 
	inline double R(double _phi, double _t);
	inline Vector3 Position(double _phi, double _t);
	inline Vector3 DelPhiVec(double _phi, double _t);
	inline Vector3 DelTVec(double _phi, double _t);
	inline Vector3 Normal(double _phi, double _t);
	//输入范围 0~2pi
	bool OnRadiatorSurface(double _phi, double _t);
	
	//XD添加计算反射点函数
	bool CalculateReflectionPoint(Vector3 StartPoint, Vector3 IncidentVector, Vector3 &ReflectionPoint);

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

	double frequency;
	int m;
	int n;
};

#endif


