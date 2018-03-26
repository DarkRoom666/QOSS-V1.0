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
//ɾȥ����Ҫ�� -20180321

using namespace std;

//This Part Set the FDTD Computation for a Radiator
class _declspec(dllexport) FDTDRadiator //�����඼д���˽ӿ�
{
public:
	FDTDRadiator();
	~FDTDRadiator();
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
	void LoadProFieldE(const char* _filename, vector<vector<complex<double>>>& _Eu, vector<vector<complex<double>>>& _Ev, int _Nu, int _Nv);

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

	//����ģ�͵Ĳ�������
	///int*** EFlag;
	///float *** EpsR;
	///float *** EpsI;

	float cdt;
	float Requireddt;
	int NN;
	int N_spa = 10;
	int OmpNum = 4;

	//fstream logfile;


};

