/*
*	created by liyun 2018/1/23
*   function ���㷴���� ���뾵�ӵ�ָ�� ���������ߺͽ����Լ��Ƿ��ཻ
*   ���� ԴEx_In,Ey_In, Դλ�á��������򣬷�����
*   ����RayTracing������׷�ٲ���
*   version 1.0
*/
#ifndef PVVA_H
#define PVVA_H

#include <iostream>
#include <vector>
#include <cassert>
#include <vector>
#include <list>
#include <complex>

#include "../util/Vector3.h"
#include "../util/Vector2.h"
#include "../util/GEMS_Memory.h"
#include "../util/Constant_Var.h"
#include "FFTDI.h"

#include "VTK/include/Field.h"
#include "VTK/include/Mirror.h"
#include "RayTracing.h"



using namespace std;

namespace calculation
{
	class PVVA
	{
	public:
		PVVA();
		~PVVA();

		void Initialization(); // ���ݳ�ʼ��

		void setUnit(double factor); // ���õ�λ

		void setSource(const Field* _field); // ����Դ

		void setFre(double _fre);

		void setMirror(Mirror * mirror);

		void CalZ0Theta(); // ����z0 �� theta

						   // �õ�Դ���е������Լ���������
		void getSourcePoint(Vector3 & interPoint, Vector3 & n_Point) const;

		// Դ���е������Լ���������
		void setSourcePoint(const Vector3 & interPoint, const Vector3 & n_Point);

		static GraphTrans getSourceGraphTrans(const Vector3 & n_Point);

		bool CheckModle(Vector3 & interPoint, Vector3 & n_Point); // check�Ƿ��ܼ���

		void getField(Field * _field); // ���ؽ��

		void CalPlane(double dis); //����ƽ�洫�� ������ƽ��λ��

		void Reflect(); // ���㷴�� �ı���Exyz_R����λ�ͼ��� �����˷�������

		void CalAmplitude(); // ���㷴��ǰ���ͶӰ��� ����Exyz_R���ȱ任

		void InterVal(); // �Է�������Ĳ����в�ֵ��ʹ�������׼������

		void Result(double dis = 1.0); // �õ������

		void getPlane(Vector3 ** &_org, Vector3 ** &_n) const; // �õ�ƽ���͵�ķ���

	private:

		void Poynting(); //���������Ӧ͢ʸ��

		void AllocateMemory();  //�����ڴ�

		void updateSource_n(Vector3 new_n); // ����Դ�ķ������Լ���ת��Ϣ


		//�����Ĺ������ A * B(����)
		complex<double> ConjugateMul(const complex<double> &A, const complex<double> &B) const;
		//��������� A * B
		complex<double> ComMul(const complex<double> &A, const complex<double> &B) const;

		double CalDistance(const Vector3 &a, const Vector3 &b) const;
	
		/************************************************
		* ֱ����ƽ���ཻ��� origֱ����� dir ���߷���
		* ����t t>0 �ཻ��ֱ��ǰ�� t<0 �ཻ��ֱ�ߺ�
		* ƽ��Ϊ�㷨ʽ n*(P-P0) = 0
		* P0Ϊƽ���ϵ�һ�� nΪƽ�淨��
		***************************************************/
		double IntersectPlane(const Vector3 &orig, const Vector3 &dir, const Vector3 &Plane_org,
			const Vector3 &Plane_n, Vector3 &intersection);

		double IntersectPoint(const Vector3 &orig, const Vector3 &dir, const Vector3 &v0,
			const Vector3 &E1, const Vector3 &E2, Vector3 &intersection);

		// ���㷴��糡 ���뷨���� ����糡ExEyEz 
		void CalReflectExyz(const Vector3 &n, const complex<double> &Ex, const complex<double> &Ey,
			const complex<double> &Ez, complex<double> &Ex_out, complex<double> &Ey_out, complex<double> &Ez_out);

		// �����ĸ���Χ�ǵ����  ��Ч�ڼ���������ABC �� ������ADC�������
		double CalSquare(const Vector3 &A, const Vector3 &B, const Vector3 &C, const Vector3 &D) const;
		double CalSquare(const Vector3 &A, const Vector3 &B, const Vector3 &C) const;
		
		// ��ֵ�õġ����жϵ��Ƿ����������ڲ�
		bool InterVal_PointinTriangle(const Vector2 &A, const Vector2 &B,
			const Vector2 &C, const Vector2 &P);

	private:
		int N, M; // ������εĳ��� N*M
		double f; // Ƶ�� Ĭ�� 10.65e9
		double lamda; // ����
		//double k;
		double ds;
		double z0; //Ŀ�����

		double unit_factor; // ����ģ�͵ĵ�λת��Ϊ�׵�ϵ��

		Mirror * mirror;

		GraphTrans SourceGraphTrans;// Դ����ת����

		complex<double> ** Px, ** Py, ** Pz;	// ��Ӧ͢ʸ�� 
		complex<double> ** Ex1, ** Ey1, ** Ez1;	// ����calculation�����ĵ��ĵ�ų�
		complex<double> ** Hx1, ** Hy1, ** Hz1;	// 

		complex<double> ** Ex_In, ** Ey_In;	// Դ 

		complex<double> ** Ex_R, ** Ey_R, ** Ez_R;	// ������������ĵ��ĵ糡

		//Org_Sourceƽ���е� n_Sourceƽ�淨���� R_Source����������
		Vector3 Org_Source, n_Source, R_Source; 
		double theta_Source;   // Դ��-Z��ļн�

		Vector3 **Plane;    // ƽ���������꣨�������꣩
		Vector3 **n_Plane;  // ƽ�洫�������ķ�����
		Vector3 **R_Plane;  // �����ĸ��������
		Vector3 **Rn_Plane; // �����ĸ���ķ�����
	};
}



#endif