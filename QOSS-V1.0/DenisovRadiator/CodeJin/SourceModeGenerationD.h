#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <math.h>//��������������
#include "Mathematical_Functions.h"//�����׳ˡ�������������һ�׵����Ȼ�����ѧС����
#include "Vector3.h"
//Vector3.h���Ѿ�������Vector3.h��
//����������ļ����Ѿ����ù�"Constant_Val.h"�������治��Ҫ�����include������
//#include "Constant_Val.h"//��ѧ����H�ļ�
using namespace std;
/*
*	created by XD 2017/10/20
	modified by JM 2017 12
	ΪDenisov�ṩ����ģʽ���� �����о���
*   ��λĬ��Ϊ��
*   function ����Բ����TEmnģʽ�Լ����β����ĳ��ֲ�
*/

#ifndef SOURCEMODEGENERATIOND_H  // ȫ��д
#define SOURCEMODEGENERATIOND_H

namespace calculation
{ // ����ռ� ���������ظ�

	class SourceModeGenerationD //�������ֺ�ͷ�ļ�������ͬ
	{
	public:
		SourceModeGenerationD
		   (int _SourceKind=1,//�������࣬1�ǵͽ�TEģʽ��2�Ǹ߽�TEģʽ��3�Ǿ��β���ģʽ
			int _SourceType=1,//����TEģʽ����TMģʽ��TEΪ1��TMΪ2
            int _Rotation=0,//�������������ڵ�Ϊ0������Ϊ1������Ϊ2
			int _m=0 ,//TEmn
			int _n=1 ,//TEmn
			double _Frequency=44.0e9,//Ƶ��
			double _Radius=64.0e-3,//�����뾶
			double _Height=0,//���β����ĳ���
			double _Width=0,//���β����Ŀ��,����Ҫ���ڿ��
			int _N=250//������
		);// Ĭ�Ϲ��캯��

		// ���м̳� һ��Ҫ��virtual
		~SourceModeGenerationD() {
			// ���з������ڴ�ǵ��ͷ�
			// if(p) {
			//   delete p;
			//   p = NULL;
			// }
		};//��������
	
		void SetSource_Circular(int sourcekind, int sourcetype, int rotation,
			                    int m0, int n0, double frequency, double radius);
		//����Բ����ģʽ����

		void SetSource_Rectangle(int sourcekind, int sourcetype, int m0, int n0,
			                     double frequency, double height, double width);
		//���þ��β���ģʽ����

		void SetOutputProperty(int n1);
		//��������ľ��γ��ֲ�������ʼ��
		//n1��ָ�趨��������Ŀ��ֻ����ż������Ϊ�����㲻�ܰ���ԭ�㣬����ʽ�޷�����

		bool SetWaveguideMode(int m0, int n0);
		//�������M��N�ĺ���
		
		bool FieldCalculation_Circular();
		//����Բ����ģʽ�����,nֻ����ż������Ϊ�����㲻�ܰ���ԭ�㣬����ʽ�޷�����

		bool FieldCalculation_CircularT();
		//����FDTD������Ҫ�����򳡷ֲ�

		bool FieldCalculation_Rectangle();
		//���þ��β���ģʽ�����,nֻ����ż������Ϊ�����㲻�ܰ���ԭ�㣬����ʽ�޷�����
		//�����20171106

		bool SetFieldAmplitude(double K);//���ó�ǿ��һ������ֵ��EX EY EZ�е����ֵ��

		bool GetCircularWaveguideProperty(double &phi1,double &theta1,double &Rc,double &Lc);
		//���Բ�������粼��Ԩ���뽹ɢԲ�뾶�Ȳ��������
		//phi1�ǲ���Ԩ�ǣ�theta1�Ƿ���ʸ���ں����ͶӰ�벨�����ߵļн�
		//Rc�ǽ�ɢԲ�İ뾶
		//20171109��ɲ���֤��ȷ

		bool GetCircularConverterLauncherTriangle_TE0n(Vector3 &P1,Vector3 &P2, Vector3 &P3);
		//���S����άʾ��ͼ�е����ڻ������������ε�������

		bool GetCircularSystemOpticalPath_TE0n_2Mirrors(Vector3 &P4,
			Vector3 &V1, double &Length1, double &theta1,
			Vector3 &V2, double &Length2, double &theta2,
			Vector3 &V3);
		//���S����άʾ��ͼ�е����ڻ������������εĹ�·���ú���(˫����)

		bool GetCircularSystemOpticalPath_TE0n_3Mirrors(Vector3 &P4,
			Vector3 &V1, double &Length1, double &theta1,
			Vector3 &V2, double &Length2, double &theta2,
			Vector3 &V3, double &Length3, double &theta3,
			Vector3 &V4);
		//���S����άʾ��ͼ�е����ڻ������������εĹ�·���ú���(������)

		bool GetCircularSystemOpticalPath_TE0n_4Mirrors(Vector3 &P4,
			Vector3 &V1, double &Length1, double &theta1,
			Vector3 &V2, double &Length2, double &theta2,
			Vector3 &V3, double &Length3, double &theta3,
			Vector3 &V4, double &Length4, double &theta4,
			Vector3 &V5);
		//���S����άʾ��ͼ�е����ڻ������������εĹ�·���ú���(�ľ���)

		void GetEX(vector <vector <complex <double>>> &EX0);//���EX
		void GetHX(vector <vector <complex <double>>> &HX0);//���HX
		void GetEY(vector <vector <complex <double>>> &EY0);//���EY
		void GetHY(vector <vector <complex <double>>> &HY0);//���HY
		void GetEZ(vector <vector <complex <double>>> &EZ0);//���EZ
		void GetHZ(vector <vector <complex <double>>> &HZ0);//���HZ
		bool GetJ_R(vector<complex<double>> &_HPhi, vector<complex<double>> &_Hz, int _Nphi);//���HPhi




		// ���з����ڴ� �����д��ֵ��������� = 
		//TraceLight(const TraceLight& traceLight) {
//
	//	}
	//	TraceLight operator = (const TraceLight& traceLight) {
//
	//	}

		/*
		* �򵥺������ܽ��� + ����˵�� �Լ����ز���˵��
		*/
	//	void fun() {};

	private:

		int SourceKind;//�������࣬1�ǵͽ�TEģʽ��2�Ǹ߽�TEģʽ��3�Ǿ��β���ģʽ
		int SourceType;//����TEģʽ����TMģʽ��TEΪ1��TMΪ2
		int Rotation;//�������������ڵ�Ϊ0������Ϊ1������Ϊ2
		int m;//TEmn
		int n;//TEmn
		double Frequency;//Ƶ��
		double Radius;//�����뾶
		double Height;//���β����ĳ��ȣ�һ�㳤�ȴ��ڿ��
		double Width;//���β����Ŀ��
		int N;//��ɢ����������
		double ds;//��ɢ��������
		double Lc;//Cut �߶�
		vector <vector <complex <double>>>  EX,HX,EY,HY,EZ,HZ;
		
		vector <complex<double>> HPhi_r;
		vector <complex<double>> Hz_r;
		//���������������ʽ��6����������õ����ֲ���

	};

}

#endif // TRACELIGHT_H