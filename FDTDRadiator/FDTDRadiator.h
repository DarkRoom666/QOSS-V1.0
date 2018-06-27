#pragma once
#include <stdexcept>
#include <cmath>
#include <ctime>
#include <fstream>
#include <vector>
#include <complex>
#include "Position3D.h"	//����
#include "Matrix4D.h"	//����
#include "Vector3D.h"	//����
#include <string>
#include "HuygensBoxData.h"
#include "../CUDAPhysicalOptics/CUDAPhysicalOptics.h"
//ע�����ﲻЩ����

using namespace std;

//ע�� �ⲻ��һ��ThreadŶ������д��ͬ�����к�����

//This Part Set the FDTD Computation for a Radiator
class _declspec(dllexport) FDTDRadiator //�����඼д���˽ӿ�
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
	void SetReturnInt(void(*returnint)(int, void*), void *_user);  // ע��ص�����
	void SetReturnFloat(void(*returnFloat)(float, void*), void*_user);// ע��ص�����
	void WriteApertureDataToFile(const char* _filename);
	void GetProFieldE(vector<vector<complex<double>>>& _Eu, vector<vector<complex<double>>>& _Ev, int _Nu, int _Nv);
	void GetProFieldH(vector<vector<complex<double>>>& _Hu, vector<vector<complex<double>>>& _Hv, int _Nu, int _Nv);
	void GetProPowerRatio(double& _PowerRatio);
	void LoadProFieldE(const char* _filename, vector<vector<complex<double>>>& _Eu, vector<vector<complex<double>>>& _Ev, int _Nu, int _Nv);
	void LoadProFieldH(const char* _filename, vector<vector<complex<double>>>& _Hu, vector<vector<complex<double>>>& _Hv, int _Nu, int _Nv);

private:
	// �ص�����ָ��
	void(*pFun)(int);
	//�ص�����ָ��
	void(*returnInt)(int, void*);
	void(*returnFloat)(float, void*);

	void *user; // �ص���������ָ��

	void SelectDT(float _dx, float _dy, float _dz); 
public:
	//For test
	int number;

	//Waveguidemode
	int SourceKind;	//1���ͽ�TEģʽ��2���߽�TEģʽ��3�����β���ģʽ
	int SourceType;	//1��TEģʽ��2��TMģʽ��
	int Rotation;	//0�������ڣ�1��������2��������
	int WG_m;		//TE mn
	int WG_n;		//TE mn
	double Frequency;	//Ƶ�� Hz
	double Radius;	//Բ�����뾶
	double F;		//Focus lengh for the First Mirror Low order TE case
	double WG_A;	//�ز�������
	double WG_B;	//�ز����̱�
	double WG_H;	//�����γ���
	double AP_A;	//���ڶγ���
	double AP_B;	//���ڶζ̱�
	double AP_H;	//���ڶγ���
	int Ns;			//��������
	int Ns1;

	//����������
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
	vector<vector<complex<float>>> Hu;//��Ӵų�
	vector<vector<complex<float>>> Hv;//��ӵ糡
	
	//fstream logfile
	char* FILENAME;

	//Բ��ģʽ�������ṹ����
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

	//������ģ����  ȫ�ֱ���
	int Domian_Nx;
	int Domian_Ny;
	int Domain_Nz;
	int Domian_num_pml;

	float cdt;
	float Requireddt;
	int NN;

	double PowerRatio;

	//����˿�
	complex<float>** Ex_Port;	
	complex<float>** Ey_Port;	
	complex<float>** Hx_Port;	
	complex<float>** Hy_Port;

	//ģ�͵���
	float*** Eps;
	
	//�˿�
	int Nx_exc;			
	int Ny_exc;
	int Nz_exc;
	int Shiftx;
	int Shifty;
	//ģ��
	int Nx_model;
	int Ny_model;
	int Nz_model;
	float dx;
	float dy;
	float dz;
	float cx;//ģ�͵�����λ��	�����huygensBox������λ�ã�
	float cy;//ģ�͵�����λ��
	float cz;//ģ�͵�����λ��

	int OmpNum;
	int N_spa;	//10 by default
	int timemode; //0 single frequency computation; 1 multi-frequency computation
	int huygensmode; //0 five faces huygens Box; 1. cylinder faces to be established

	HuygensBox Huygens;
};

