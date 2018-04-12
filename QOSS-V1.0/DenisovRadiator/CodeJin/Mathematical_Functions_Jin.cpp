#include "../DenisovRadiator/CodeJin/Mathematical_Functions_Jin.h"


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