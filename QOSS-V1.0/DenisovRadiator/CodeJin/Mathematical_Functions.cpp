#include "Mathematical_Functions.h"

double fact(int n)
	{
		double s = 1.0;
		if (n > 0) {
			for (int i = 1; i <= n; i++)
			{
				s = s*i;
			}
		}
		else s = 1.0;
		return s;
	}

double D1_J(int m, double z)
{
	double a = 0.5*((_jn((m - 1), z)) - ((_jn((m + 1), z))));
	//����ѧ�����ߵļ�㹫ʽ
	return a;
}

complex<double> H_phi(double Kz0, double Kr0, int m0, double r0, double Phi0)
{
	complex <double> a1(0, -m0*Phi0);//-j*m0*Phi0
	complex <double> a2(0, 1);//j
	complex <double> a;
	a = -Kz0*m0 / Kr0 / Kr0 / r0*_jn(m0, Kr0*r0)*exp(a1) / a2 / 2.0;

	return a;

	//���س��ֲ�ֵH_phi(������ʽ)
}

complex<double> H_r(double Kz0, double Kr0, int m0, double r0, double Phi0)
{
	complex <double> a1(0, -m0*Phi0);//-j*m0*Phi0
	complex <double> a;
	a = -Kz0 / Kr0 *D1_J(m0, Kr0*r0)*exp(a1) / 2.0;

	return a;
}

complex<double> E_phi(double Kz0, double Kr0, int m0, double r0, double Phi0, double Omega0, double Beta0)
{
	complex <double> a;
	a = -Omega0*Mu0 / Beta0*H_r(Kz0, Kr0, m0, r0, Phi0);
	return a;
}

complex<double> E_r(double Kz0, double Kr0, int m0, double r0, double Phi0, double Omega0, double Beta0)
{
	complex <double> a;
	a = Omega0*Mu0 / Beta0*H_phi(Kz0, Kr0, m0, r0, Phi0);
	return a;
}

complex<double> H_z(double Kr0, int m0, double r0, double Phi0)
{
	complex <double> a1(0, -m0*Phi0);//-j*m0*Phi0
	complex <double> a;
	a = _jn(m0, Kr0*r0)*exp(a1) / 2.0;
	return a;

	//E_z=0
}

complex<double> H_phi_LeftRo(double Kz0, double Kr0, int m0, double r0, double Phi0)
{
	complex <double> a1(0, m0*Phi0);//j*m0*Phi0
	complex <double> a2(0, 1);//j
	complex <double> a;
	a = Kz0*m0 / Kr0 / Kr0 / r0*_jn(m0, Kr0*r0)*exp(a1) / 2.0 / a2;

	return a;
}

complex<double> H_r_LeftRo(double Kz0, double Kr0, int m0, double r0, double Phi0)
{
	complex <double> a1(0, m0*Phi0);//j*m0*Phi0
	complex <double> a;
	a = -Kz0 / Kr0 * D1_J(m0, Kr0*r0)*exp(a1) / 2.0;

	return a;
}

complex<double> E_phi_LeftRo(double Kz0, double Kr0, int m0, double r0, double Phi0,
	                         double Omega0, double Beta0)
{
	complex <double> a;
	a = -Omega0*Mu0 / Beta0*H_r_LeftRo(Kz0, Kr0, m0, r0, Phi0);
	return a;
}

complex<double> E_r_LeftRo(double Kz0, double Kr0, int m0, double r0, double Phi0, 
	                       double Omega0, double Beta0)
{
	complex <double> a;
	a = Omega0*Mu0 / Beta0*H_phi_LeftRo(Kz0, Kr0, m0, r0, Phi0);
	return a;
}

complex<double> H_z_LeftRo(double Kr0, int m0, double r0, double Phi0)
{
	complex <double> a1(0, m0*Phi0);//j*m0*Phi0
	complex <double> a;
	a = _jn(m0, Kr0*r0)*exp(a1) / 2.0;
	return a;

	//E_z_LeftRo=0
}


double rootdbessel(int m, int n)	//Find the nth root of mth order derivate bessel function of the first kind
{
	//���Ķ���д��Բ����ģʽ��������-SourceModeGeneration.cpp��
	//��ȡ��������Ϊ�����ĺ��� ����
	const int MaxLoopNumber = 10000;//���ѭ������
	const double step0 = 0.1; //��ɨ�Ĳ���
	const double step1 = step0 / MaxLoopNumber;//��һ��ϸɨ�Ĳ���

	int n0 = 0;//��n0���������ĳ�ʼ��
	int i0 = 0;//�ҵ��˵�n�����������������
	for (int i = 1; i <= MaxLoopNumber; i++)//���ⲻ��0��ʼ
	{
		if ((D1_J(m, i*step0) >= 0.0 && D1_J(m, (i + 1)*step0) <= 0.0)
			|| (D1_J(m, i*step0) <= 0.0 && D1_J(m, (i + 1)*step0) >= 0.0))
		{
			n0 = n0 + 1;//��ʼ����
			i0 = i;//���¸�����
		}
		if (n0 == n) break;//�ҵ���n����������ֱ������ѭ������ʱ����������[i0*ds��(i0+1)*ds]������

	}

	//��ʼ��һ��ϸɨ[i0*step0��(i0+1)*step0]����,�����С
	int j0 = 0;//�ҵ��˵�n�����������������
	for (int j = 0; j <= MaxLoopNumber; j++)//��0��ʼ
	{

		if ((D1_J(m, j*step1 + i0*step0) >= 0.0 && D1_J(m, (j + 1)*step1 + i0*step0) <= 0.0)
			|| (D1_J(m, j*step1 + i0*step0) <= 0.0 && D1_J(m, (j + 1)*step1 + i0*step0) >= 0.0))
		{
			j0 = j;//�ҵ���������䣬�����¸������±�
			break;//�ҵ���n����������ֱ������ѭ������ʱ����������[i*ds��(i+1)*ds]������
		}
	}

	double Xmn = 0.5*((j0*step1 + i0*step0) + ((j0 + 1)*step1 + i0*step0));
	//20171026�Ѿ���֤����ȫ����MathCad������ �Ķ�Good��
	return Xmn;
}

double rootbessel(int m, int n)	//Find the nth root of mth order bessel function of the first kind
{
	//��rootdbessel���������ϸ�д	����
	const int MaxLoopNumber = 10000;//���ѭ������
	const double step0 = 0.1; //��ɨ�Ĳ���
	const double step1 = step0 / MaxLoopNumber;//��һ��ϸɨ�Ĳ���

	int n0 = 0;//��n0���������ĳ�ʼ��
	int i0 = 0;//�ҵ��˵�n�����������������
	for (int i = 1; i <= MaxLoopNumber; i++)//���ⲻ��0��ʼ
	{
		if ((jn(m, i*step0) >= 0.0 && jn(m, (i + 1)*step0) <= 0.0)
			|| (jn(m, i*step0) <= 0.0 && jn(m, (i + 1)*step0) >= 0.0))
		{
			n0 = n0 + 1;//��ʼ����
			i0 = i;//���¸�����
		}
		if (n0 == n) break;//�ҵ���n����������ֱ������ѭ������ʱ����������[i0*ds��(i0+1)*ds]������

	}

	//��ʼ��һ��ϸɨ[i0*step0��(i0+1)*step0]����,�����С
	int j0 = 0;//�ҵ��˵�n�����������������
	for (int j = 0; j <= MaxLoopNumber; j++)//��0��ʼ
	{

		if ((jn(m, j*step1 + i0*step0) >= 0.0 && jn(m, (j + 1)*step1 + i0*step0) <= 0.0)
			|| (jn(m, j*step1 + i0*step0) <= 0.0 && jn(m, (j + 1)*step1 + i0*step0) >= 0.0))
		{
			j0 = j;//�ҵ���������䣬�����¸������±�
			break;//�ҵ���n����������ֱ������ѭ������ʱ����������[i*ds��(i+1)*ds]������
		}
	}

	double Xmn = 0.5*((j0*step1 + i0*step0) + ((j0 + 1)*step1 + i0*step0));
	//20171026�Ѿ���֤����ȫ����MathCad������ Good��
	return Xmn;
}

double kzmnTE(int m, int n, double lambda, double radius) {
	double kz;
	double k0 = 2 * Pi / lambda;
	double X0 = rootdbessel(m,n);
	double kr = X0 / radius;
	kz = sqrt(k0*k0 - kr*kr);
	return kz;
}