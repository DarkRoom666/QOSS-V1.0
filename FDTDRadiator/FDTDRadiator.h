#pragma once
#include <stdexcept>
#include <cmath>
#include <ctime>
#include <fstream>
#include <vector>
#include <complex>
#include "Position3D.h"	//公用
#include "Matrix4D.h"	//公用
#include "Vector3D.h"	//公用
#include <string>
#include "HuygensBoxData.h"
#include "../CUDAPhysicalOptics/CUDAPhysicalOptics.h"
//注意这里不些引用

using namespace std;

//注意 这不是一个Thread哦，可以写不同的运行函数！

//This Part Set the FDTD Computation for a Radiator
class _declspec(dllexport) FDTDRadiator //整个类都写成了接口
{
public:
	FDTDRadiator();
	~FDTDRadiator();

	//FDTDcomputation
	void SetUpModelRadiator(string _filename);
	void SetUpExcRadiator(string _filename);
	void SetFreqList(void);
	void SetUpCommonFDTD(double _freq, double _ompNum, int _N_spa, int _timemode, int _huygensmode);
	void runCommonFDTD(void);
	//ApertureComputation


	//This Function Set Up the LowOrder Vlasov Radiator Computation
	void SetUpLowOrderVlasovRadiator(int _WG_m, int _WG_n, double _Frequency, double _Radius, double _F, int _Ns, int _OmpNum);
	void SetUpAperturePlane(Position3D _AperturePosition, Vector3D _ApertureDirection, Vector3D _UDirection, Vector3D _VDirection, int _Nu, int _Nv, double _Lu, double _Lv);
	void SetInt(int input);
	void run();
	void SetpFun(void(*pFun)(int));
	void SetReturnInt(void(*returnint)(int, void*), void *_user);  // 注册回调函数
	void SetReturnFloat(void(*returnFloat)(float, void*), void*_user);// 注册回调函数
	void WriteApertureDataToFile(const char* _filename);
	void GetProFieldE(vector<vector<complex<double>>>& _Eu, vector<vector<complex<double>>>& _Ev, int _Nu, int _Nv);
	void GetProFieldH(vector<vector<complex<double>>>& _Hu, vector<vector<complex<double>>>& _Hv, int _Nu, int _Nv);
	void GetProPowerRatio(double& _PowerRatio);
	void LoadProFieldE(const char* _filename, vector<vector<complex<double>>>& _Eu, vector<vector<complex<double>>>& _Ev, int _Nu, int _Nv);
	void LoadProFieldH(const char* _filename, vector<vector<complex<double>>>& _Hu, vector<vector<complex<double>>>& _Hv, int _Nu, int _Nv);

private:
	// 回调函数指针
	void(*pFun)(int);
	//回调函数指针
	void(*returnInt)(int, void*);
	void(*returnFloat)(float, void*);

	void *user; // 回调函数的类指针

	void SelectDT(float _dx, float _dy, float _dz); 
public:
	//For test
	int number;

	//Waveguidemode
	int SourceKind;	//1：低阶TE模式；2：高阶TE模式；3：矩形波导模式
	int SourceType;	//1：TE模式；2：TM模式；
	int Rotation;	//0：不存在；1：左旋；2：右旋；
	int WG_m;		//TE mn
	int WG_n;		//TE mn
	double Frequency;	//频率 Hz
	double Radius;	//圆波导半径
	double F;		//Focus lengh for the First Mirror Low order TE case
	double WG_A;	//矩波导长边
	double WG_B;	//矩波导短边
	double WG_H;	//波导段长度
	double AP_A;	//开口段长边
	double AP_B;	//开口段短边
	double AP_H;	//开口段长度
	int Ns;			//采样点数
	int Ns1;

	//出射口面参数
	Position3D AperturePosition;
	Vector3D ApertureDirection;
	Vector3D UDirection;
	Vector3D VDirection;
	int Nu;
	int Nv;
	float Lu;
	float Lv;
	vector<vector<complex<float>>> Eu;
	vector<vector<complex<float>>> Ev;
	vector<vector<complex<float>>> Hu;//添加磁场
	vector<vector<complex<float>>> Hv;//添加电场
	
	//fstream logfile
	char* FILENAME;

	//圆电模式辐射器结构参量
	double PHId;
	double THETAd;
	double Rcd;
	double Lcd;
	float PHI;
	float THETA;
	float Rc;
	float Lc;
	float Lp;
	float ds;

	//建立画模型用  全局变量
	int Domian_Nx;
	int Domian_Ny;
	int Domain_Nz;
	int Domian_num_pml;

	float cdt;
	float Requireddt;
	int NN;

	double PowerRatio;

	//入射端口
	complex<float>** Ex_Port;	
	complex<float>** Ey_Port;	
	complex<float>** Hx_Port;	
	complex<float>** Hy_Port;

	//模型点阵
	float*** Eps;
	
	//端口
	int Nx_exc;			
	int Ny_exc;
	int Nz_exc;
	int Shiftx;
	int Shifty;
	//模型
	int Nx_model;
	int Ny_model;
	int Nz_model;
	float dx;
	float dy;
	float dz;
	float cx;//模型的中心位置	这就是huygensBox的中心位置！
	float cy;//模型的中心位置
	float cz;//模型的中心位置

	int OmpNum;
	int N_spa;	//10 by default
	int timemode; //0 single frequency computation; 1 multi-frequency computation
	int huygensmode; //0 five faces huygens Box; 1. cylinder faces to be established

	HuygensBox Huygens;
};

