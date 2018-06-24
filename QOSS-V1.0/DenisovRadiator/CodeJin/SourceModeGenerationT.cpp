#include "SourceModeGenerationT.h"

calculation::SourceModeGenerationT::SourceModeGenerationT(
	int _SourceKind, int _SourceType, int _Rotation, int _m, int _n,
	double _Frequency, double _Radius, double _Height, double _Width,
	int _N) :
	SourceKind(_SourceKind), SourceType(_SourceType), Rotation(_Rotation),
	m(_m), n(_n), Frequency(_Frequency), Radius(_Radius),
	Height(_Height), Width(_Width),
	N(_N)
{
	ds= _Radius * 2.0 / _N;//�������,Ĭ��ΪԲ��������
	
	EX.resize(_N); HX.resize(_N);
	EY.resize(_N); HY.resize(_N);
	EZ.resize(_N); HZ.resize(_N);
	for (int i = 0; i < _N; i++)
	{
		EX[i].resize(_N); HX[i].resize(_N);
		EY[i].resize(_N); HY[i].resize(_N);
		EZ[i].resize(_N); HZ[i].resize(_N);
	}//����6�������ľ���ߴ磺N*N
}//���캯������

//������������Ҫ������

//����Բ�����Ĳ���
void calculation::SourceModeGenerationT::SetSource_Circular(int _sourcekind, int _sourcetype,
	int _rotation, int _m0, int _n0, double _frequency, double _radius)
{
	SourceKind = _sourcekind;//ֻ����1����2
	SourceType = _sourcetype;
	Rotation = _rotation;
	m = _m0;
	n = _n0;
	Frequency = _frequency;
	Radius = _radius;
	ds = 2.0*Radius / N ;
	Height = 0;
	Width = 0;
}

void calculation::SourceModeGenerationT::SetSource_Rectangle(int sourcekind, int sourcetype, 
	int m0, int n0, double frequency, double height, double width)

{
	SourceKind = sourcekind;//ֻ����3
	SourceType = sourcetype;
	m = m0;
	n = n0;
	Frequency = frequency;
	Radius = 0;
	Height = height;
	Width = width;
	ds = Height / (N - 1);
}

void calculation::SourceModeGenerationT::SetOutputProperty(int n1)
{
	N = n1;//������
	if (SourceKind==1 || SourceKind == 2)//�����Բ�����Ļ�
	ds = Radius * 2.0 / N;//Բ�����Ĳ������
	else  ds = Height / N ;//���β����Ĳ�����࣬�Գ���Ϊ׼

	EX.resize(N); HX.resize(N);
	EY.resize(N); HY.resize(N);
	EZ.resize(N); HZ.resize(N);
	for (int i = 0; i < N; i++)
	{
		EX[i].resize(N); HX[i].resize(N); 
		EY[i].resize(N); HY[i].resize(N);
		EZ[i].resize(N); HZ[i].resize(N);
	}//����6�������ľ���ߴ磺N*N

	//20171030����Ҫ���S���һ��д�ĶԲ���
}


bool calculation::SourceModeGenerationT::SetWaveguideMode(int m0, int n0)
{
	if ((m0 < 0) || (n0 < 0)) return false;
	m = m0;
	n = n0;
	return true;
}

//���㲢���Բ����ģʽ���ֲ�-�����ų�-������д 201712 Ϊ FDTD�ṩ���룬�迼�ǲ���λ�õĿռ���λ
bool calculation::SourceModeGenerationT::FieldCalculation_CircularT()
{
	//��ʼ���㳡�ֲ�
	double Lamda = C_Speed / Frequency;
	double K0 = 2 * Pi / Lamda;
	double Omega = 2 * Pi * Frequency;

	//֮ǰ����ȡ��Сdelta���ݹ�ʽ�󵼵ķ�ʽ���Ȳ���꣬
	//�����ȡ����ʦ���������bessel�����ļ���չ��ʽ�����󵼣������üбƷ��õ��������ֵ
	//20171026����ѧ���������ټ��㵼������ʽ������

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
	//20171026�Ѿ���֤����ȫ����MathCad������ Good��

	double Kr = Xmn / Radius;
	if (Kr > K0) return false;//��ֹ

	double Kz = pow((K0*K0 - Kr*Kr), 0.5);
	double Lamda_c = 2.0 * Pi / Kr;
	double Beta = 2.0 * Pi*Frequency / C_Speed*pow((1.0 - Lamda*Lamda / Lamda_c / Lamda_c), 0.5);

	//���濪ʼ��д������ֲ�
	vector <double> X(N);//X���꣬����ΪN
	vector <double> Y(N);//Y���꣬����ΪN
	vector <vector <double>> R(N), JUDGE(N), PHI(N);
	vector <vector <complex <double>>> ER(N), HR(N), EPhi(N), HPhi(N);//��������
	for (int i = 0; i < N; i++)
	{
		R[i].resize(N); JUDGE[i].resize(N); PHI[i].resize(N);
		ER[i].resize(N); HR[i].resize(N); EPhi[i].resize(N); HPhi[i].resize(N);
	}//����R�Լ���ؾ���ĳߴ磺N*N

	ds = Radius * 2.0 / N;

	//First Generate Ex Hy
	for (int i = 0; i < N; i++)
	{
		X[i] = ds*(i - (N - 1.0) / 2.0);
		for (int j = 0; j < N; j++)
		{
			Y[j] = ds*(j - (N) / 2.0);
			R[i][j] = pow((X[i] * X[i] + Y[j] * Y[j]), 0.5);

			JUDGE[i][j] = 0.0;
			if (R[i][j]<Radius) JUDGE[i][j] = 1.0;//���ڰ뾶��Χ��JUDGE��Ϊ1������Ϊ0

			if ((X[i] == 0) && (Y[j] == 0)) PHI[i][j] = 0.0;
			else PHI[i][j] = atan2(Y[j], X[i]);//�����Ӧ�ļн�phi
			//atan2 ��������� 2Pi�ռ�

			if (JUDGE[i][j] != 0.0)//����ڲ����ڲ�
			{
				if (Rotation != 1)//������
				{
					ER[i][j] = E_r(Kz, Kr, m, R[i][j], PHI[i][j], Omega, Beta);
					HR[i][j] = H_r(Kz, Kr, m, R[i][j], PHI[i][j]);
					EPhi[i][j] = E_phi(Kz, Kr, m, R[i][j], PHI[i][j], Omega, Beta);
					HPhi[i][j] = H_phi(Kz, Kr, m, R[i][j], PHI[i][j]);
				}
				else //���������
				{
					ER[i][j] = E_r_LeftRo(Kz, Kr, m, R[i][j], PHI[i][j], Omega, Beta);
					HR[i][j] = H_r_LeftRo(Kz, Kr, m, R[i][j], PHI[i][j]);
					EPhi[i][j] = E_phi_LeftRo(Kz, Kr, m, R[i][j], PHI[i][j], Omega, Beta);
					HPhi[i][j] = H_phi_LeftRo(Kz, Kr, m, R[i][j], PHI[i][j]);
				}

				EX[i][j] = cos(PHI[i][j])*ER[i][j] - sin(PHI[i][j])*EPhi[i][j];
				HY[i][j] = sin(PHI[i][j])*HR[i][j] + cos(PHI[i][j])*HPhi[i][j];
			}
			else//����ڲ����ⲿ
			{
				ER[i][j] = 0; HR[i][j] = 0; EPhi[i][j] = 0; HPhi[i][j] = 0;
				EX[i][j] = 0;								HY[i][j] = 0; 
			}//��Բ������Χ�ⳡֵȫΪ0
		}//j
	}//i

	//Then Generate Ey Hx
	for (int i = 0; i < N; i++)
	{
		X[i] = ds*(i - N / 2.0);
		for (int j = 0; j < N; j++)
		{
			Y[j] = ds*(j - (N - 1.0) / 2.0);
			R[i][j] = pow((X[i] * X[i] + Y[j] * Y[j]), 0.5);

			JUDGE[i][j] = 0.0;
			if (R[i][j]<Radius) JUDGE[i][j] = 1.0;//���ڰ뾶��Χ��JUDGE��Ϊ1������Ϊ0

			if ((X[i] == 0) && (Y[j] == 0)) PHI[i][j] = 0.0;
			else PHI[i][j] = atan2(Y[j], X[i]);//�����Ӧ�ļн�phi
											   //atan2 ��������� 2Pi�ռ�

			if (JUDGE[i][j] != 0.0)//����ڲ����ڲ�
			{
				if (Rotation != 1)//������
				{
					ER[i][j] = E_r(Kz, Kr, m, R[i][j], PHI[i][j], Omega, Beta);
					HR[i][j] = H_r(Kz, Kr, m, R[i][j], PHI[i][j]);
					EPhi[i][j] = E_phi(Kz, Kr, m, R[i][j], PHI[i][j], Omega, Beta);
					HPhi[i][j] = H_phi(Kz, Kr, m, R[i][j], PHI[i][j]);
				}
				else //���������
				{
					ER[i][j] = E_r_LeftRo(Kz, Kr, m, R[i][j], PHI[i][j], Omega, Beta);
					HR[i][j] = H_r_LeftRo(Kz, Kr, m, R[i][j], PHI[i][j]);
					EPhi[i][j] = E_phi_LeftRo(Kz, Kr, m, R[i][j], PHI[i][j], Omega, Beta);
					HPhi[i][j] = H_phi_LeftRo(Kz, Kr, m, R[i][j], PHI[i][j]);
				}

				EY[i][j] = sin(PHI[i][j])*ER[i][j] + cos(PHI[i][j])*EPhi[i][j];
				HX[i][j] = cos(PHI[i][j])*HR[i][j] - sin(PHI[i][j])*HPhi[i][j];

			}
			else//����ڲ����ⲿ
			{
				ER[i][j] = 0; HR[i][j] = 0; EPhi[i][j] = 0; HPhi[i][j] = 0;
							  EY[i][j] = 0; HX[i][j] = 0;
			}//��Բ������Χ�ⳡֵȫΪ0
		}//j
	}//i

	return true;//������
}


//���㲢���Բ����ģʽ���ֲ� - ������
bool calculation::SourceModeGenerationT::FieldCalculation_Circular()
{
	//��ʼ���㳡�ֲ�
	double Lamda = C_Speed / Frequency;
	double K0 = 2 * Pi / Lamda;
	double Omega= 2 * Pi * Frequency;

	//֮ǰ����ȡ��Сdelta���ݹ�ʽ�󵼵ķ�ʽ���Ȳ���꣬
	//�����ȡ����ʦ���������bessel�����ļ���չ��ʽ�����󵼣������üбƷ��õ��������ֵ
	//20171026����ѧ���������ټ��㵼������ʽ������
	
	const int MaxLoopNumber = 10000;//���ѭ������
	const double step0 = 0.1; //��ɨ�Ĳ���
	const double step1 = step0 / MaxLoopNumber;//��һ��ϸɨ�Ĳ���

	int n0 = 0;//��n0���������ĳ�ʼ��
	int i0 = 0;//�ҵ��˵�n�����������������
	for (int i = 1; i <= MaxLoopNumber; i++)//���ⲻ��0��ʼ
	{
		if ((D1_J(m,i*step0)>=0.0 && D1_J(m, (i+1)*step0)<=0.0)
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
		
		if ((D1_J(m, j*step1 +i0*step0) >= 0.0 && D1_J(m, (j+1)*step1 + i0*step0) <= 0.0)
			|| (D1_J(m, j*step1 + i0*step0) <= 0.0 && D1_J(m, (j + 1)*step1 + i0*step0) >= 0.0))
		{
			j0 = j;//�ҵ���������䣬�����¸������±�
			break;//�ҵ���n����������ֱ������ѭ������ʱ����������[i*ds��(i+1)*ds]������
		}
	}

	double Xmn = 0.5*((j0*step1 + i0*step0) + ((j0 + 1)*step1 + i0*step0));
	//20171026�Ѿ���֤����ȫ����MathCad������

	double Kr = Xmn / Radius;
	if (Kr > K0) return false;//��ֹ

	double Kz = pow((K0*K0 - Kr*Kr), 0.5);
	double Lamda_c = 2.0 * Pi / Kr;
	double Beta= 2.0 * Pi*Frequency/C_Speed*pow((1.0 - Lamda*Lamda/ Lamda_c/Lamda_c), 0.5);
	
	//���濪ʼ��д������ֲ�
	vector <double> X(N);//X���꣬����ΪN
	vector <double> Y(N);//Y���꣬����ΪN
	vector <vector <double>> R(N), JUDGE(N), PHI(N);
	vector <vector <complex <double>>> ER(N), HR(N), EPhi(N), HPhi(N);//��������
	for (int i = 0; i < N; i++)
	{
		R[i].resize(N); JUDGE[i].resize(N); PHI[i].resize(N);
		ER[i].resize(N); HR[i].resize(N); EPhi[i].resize(N); HPhi[i].resize(N);
	}//����R�Լ���ؾ���ĳߴ磺N*N
	
	for (int i = 0; i < N; i++)
	{
		X[i] = ds*(i - (N - 1.0) / 2.0);
		for (int j = 0; j < N; j++)
		{
			Y[j] = ds*(j - (N - 1.0) / 2.0);
			R[i][j] = pow((X[i] * X[i] + Y[j] * Y[j]), 0.5);
			
			JUDGE[i][j] = 0.0;
			if (R[i][j]<Radius) JUDGE[i][j] = 1.0;//���ڰ뾶��Χ��JUDGE��Ϊ1������Ϊ0
			
			if ((X[i] == 0) && (Y[j] == 0)) PHI[i][j] = 0.0;
			else PHI[i][j] = atan2(Y[j] , X[i]);//�����Ӧ�ļн�phi
			
			//atan�����Χ-pi/2��pi/2 ��Ƿ�Χ��-PI~ PI��

			EZ[i][j] = 0;//TEģʽ��EzΪ0

			if (JUDGE[i][j] != 0.0)//����ڲ����ڲ�
			{
				if (Rotation != 1)//������
				{
					ER[i][j] = E_r(Kz, Kr, m, R[i][j], PHI[i][j], Omega, Beta);
					HR[i][j] = H_r(Kz, Kr, m, R[i][j], PHI[i][j]);
					EPhi[i][j] = E_phi(Kz, Kr, m, R[i][j], PHI[i][j], Omega, Beta);
					HPhi[i][j] = H_phi(Kz, Kr, m, R[i][j], PHI[i][j]);
					HZ[i][j] = H_z(Kr, m, R[i][j], PHI[i][j]);
				}
				else //���������
				{
					ER[i][j] = E_r_LeftRo(Kz, Kr, m, R[i][j], PHI[i][j], Omega, Beta);
					HR[i][j] = H_r_LeftRo(Kz, Kr, m, R[i][j], PHI[i][j]);
					EPhi[i][j] = E_phi_LeftRo(Kz, Kr, m, R[i][j], PHI[i][j], Omega, Beta);
					HPhi[i][j] = H_phi_LeftRo(Kz, Kr, m, R[i][j], PHI[i][j]);
					HZ[i][j] = H_z_LeftRo(Kr, m, R[i][j], PHI[i][j]);
				}

				EX[i][j] = cos(PHI[i][j])*ER[i][j] - sin(PHI[i][j])*EPhi[i][j];
				HX[i][j] = cos(PHI[i][j])*HR[i][j] - sin(PHI[i][j])*HPhi[i][j];
				EY[i][j] = sin(PHI[i][j])*ER[i][j] + cos(PHI[i][j])*EPhi[i][j];
				HY[i][j] = sin(PHI[i][j])*HR[i][j] + cos(PHI[i][j])*HPhi[i][j];
			}
			else//����ڲ����ⲿ
			{
				ER[i][j]=0; HR[i][j]=0; EPhi[i][j] = 0; HPhi[i][j] = 0; HZ[i][j] = 0;
				EX[i][j] = 0; HX[i][j] = 0; EY[i][j] = 0; HY[i][j] = 0;
			}//��Բ������Χ�ⳡֵȫΪ0
		}
	}
	return true;//������
}

//������β������ֲ�
bool calculation::SourceModeGenerationT::FieldCalculation_Rectangle()
{
	double Lamda = C_Speed / Frequency;
	double K0 = 2 * Pi / Lamda;
	double Omega = 2 * Pi * Frequency;
	double Lamda_c = 2.0 / pow((m*m*Pi*Pi/Height/ Height + n*n*Pi*Pi/Width/Width), 0.5);
	double Kc = 2.0 * Pi / Lamda_c;

	if (Kc > K0) return false;//Kc����K0ʱ��������ֹ������������ȥ��
	
	double Beta = 2.0 * Pi*Frequency / C_Speed*pow((1.0 - Lamda*Lamda / Lamda_c / Lamda_c), 0.5);
	complex <double> a1(0, 1);//��������j

	//���濪ʼ��д���β���������ֲ�
	vector <double> X(N);//X���꣬����ΪN
	vector <double> Y(N);//Y���꣬����ΪN
	vector <vector <double>> JUDGE(N);
	for (int i = 0; i < N; i++)
	{
		JUDGE[i].resize(N);
	}//����JUDGE�Լ���ؾ���ĳߴ磺N*N

	for (int i = 0; i < N; i++)
	{
		X[i] = ds*i;//�����X��Y����ʼ��Ϊ(0,0)���յ�Ϊ(a,b)
		for (int j = 0; j < N; j++)
		{
			Y[j] = ds*j;
			JUDGE[i][j] = 0;
			if ((abs((X[i])-ds*(N - 1.0) / 2.0) < 0.5*Height) &&
				(abs((Y[j]- ds*(N - 1.0) / 2.0)) < 0.5*Width))
			JUDGE[i][j] = 1;//�ڲ����ڲ���ʱ��JUDGEֵΪ1

			if (JUDGE[i][j] != 0.0)//����ڲ����ڲ�
			{
				if (SourceType == 1)//TEģʽ
				{
					HX[i][j] = a1*Beta / Kc / Kc*(m*Pi / Height)*sin(m*Pi / Height*X[i])*cos(n*Pi / Width*Y[i]);
					HY[i][j] = a1*Beta / Kc / Kc*(n*Pi / Width)*cos(m*Pi / Height*X[i])*sin(n*Pi / Width*Y[i]);
					HZ[i][j] = cos(m*Pi / Height*X[i])*cos(n*Pi / Width*Y[i]);
					EX[i][j] = Omega*Mu0 / Beta*HY[i][j];
					EY[i][j] = -Omega*Mu0 / Beta*HX[i][j];
					EZ[i][j] = 0;
				}
				if (SourceType == 2)//TMģʽ
				{
					EX[i][j] = -a1*Beta/Kc/Kc*(m*Pi/Height)*cos(m*Pi / Height*X[i])*sin(n*Pi / Width*Y[i]);
					EY[i][j] = -a1*Beta/Kc/Kc*(n*Pi/Width)*sin(m*Pi / Height*X[i])*cos(n*Pi / Width*Y[i]);
					EZ[i][j] = sin(m*Pi / Height*X[i])*sin(n*Pi / Width*Y[i]);
					HX[i][j] = -Omega*Eps0 / Beta*EY[i][j];
					HY[i][j] = Omega*Eps0 / Beta*EX[i][j];
					HZ[i][j] = 0;
				}
			}
			else//����ھ��β����ⲿ
			{
				EX[i][j] = 0; EY[i][j] = 0; EZ[i][j] = 0;
				HX[i][j] = 0; HY[i][j] = 0; HZ[i][j] = 0;
			}//�ھ��β�����Χ�ⳡֵȫΪ0

		}
	}
	
  return true;//������
}

//���ó�ǿ��һ������ֵ
bool calculation::SourceModeGenerationT::SetFieldAmplitude(double K)
{
	if (K<=0.0) return false;//K����Ϊ��ֵ
	double MAX_EX=0.0, MAX_EY = 0.0,MAX_EZ = 0.0,MAX = 0.0;//��ʼ������
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			{
			if (abs(EX[i][j]) >= MAX_EX) MAX_EX = abs(EX[i][j]);
			if (abs(EY[i][j]) >= MAX_EY) MAX_EY = abs(EY[i][j]);
			if (abs(EZ[i][j]) >= MAX_EZ) MAX_EZ = abs(EZ[i][j]);
			}
    }

	if (MAX_EX >= MAX_EY) MAX = MAX_EX; else MAX = MAX_EY;
	if (MAX_EZ >= MAX) MAX = MAX_EZ;//�õ��糡��ǿ���ֵMAX
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			EX[i][j] = EX[i][j] / MAX*K; HX[i][j] = HX[i][j] / MAX*K;
			EY[i][j] = EY[i][j] / MAX*K; HY[i][j] = HY[i][j] / MAX*K;
			EZ[i][j] = EZ[i][j] / MAX*K; HZ[i][j] = HZ[i][j] / MAX*K;
		}
	}//���˵糡�ֲ����ֵΪָ��ֵK
	return true;
}

bool calculation::SourceModeGenerationT::GetCircularWaveguideProperty(double &phi1, double &theta1, double &Rc, double &Lc)
{
	if (SourceKind == 3) return false;//�ú����������ھ��β���

	//�ٴ�ʹ�������Ѿ�����֤��������������һ�׵������ķ���
	double Lamda = C_Speed / Frequency;
	double K0 = 2 * Pi / Lamda;
	double Omega = 2 * Pi * Frequency;

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
	//20171026�Ѿ���֤����ȫ����MathCad������

	double Kr = Xmn / Radius;
	if (Kr > K0) return false;//��ֹ

	double Kz = pow((K0*K0 - Kr*Kr), 0.5);
	double Lamda_c = 2.0 * Pi / Kr;
	double Beta = 2.0 * Pi*Frequency / C_Speed*pow((1.0 - Lamda*Lamda / Lamda_c / Lamda_c), 0.5);

	//���濪ʼ���㲼��Ԩ�ǵ���ز���
	phi1 = asin(Kr / K0);//phi1�ǲ���Ԩ��
	theta1 = acos(m / Xmn); //theta1�Ƿ���ʸ���ں����ͶӰ�벨�����ߵļн�
	Rc = Radius*m / Xmn;//Rc�ǽ�ɢԲ�İ뾶

	//Modified by Ming Jin 20171213 
	if (SourceKind == 1) {//TE mode
		Lc = 2 * Radius*tan(Pi * 0.5 - phi1);	//2*a*cotphi
	}
	else if (SourceKind == 2) {//TM mode
		Lc = 2 * Pi * Radius * Radius * Kz / Xmn * sqrt(1 - pow(m / Xmn, 2)) / acos(m / Xmn);
	}

	return true;
}

bool calculation::SourceModeGenerationT::GetCircularConverterLauncherTriangle_TE0n(
	Vector3 &P1, Vector3 &P2, Vector3 &P3)
{
	if ((SourceKind !=1) && (SourceType != 1) && (m!=0)) return false;

	P1= Vector3(-Radius, 0, 0);
	P2= Vector3(Radius, 0, 0);
	double phi; double theta; double Rc;
	if (GetCircularWaveguideProperty(phi, theta, Rc,Lc))
		P3 = Vector3(Radius, 0, 2.0*Radius / (tan(phi)));
	return true;
}

bool calculation::SourceModeGenerationT::GetCircularSystemOpticalPath_TE0n_2Mirrors(
	Vector3 & P4, Vector3 & V1, double & Length1, double & theta1,
	                 Vector3 & V2, double & Length2, double & theta2, Vector3 & V3)
{
	if ((SourceKind != 1) && (SourceType != 1) && (m != 0)) return false;
	
	P4 = Vector3(Radius, 0, 0);
	double phi; double theta; double Rc;
	if (GetCircularWaveguideProperty(phi, theta, Rc,Lc))//�����S���Ҽ��һ��
	V1 = Vector3(-1, 0, 1 / (tan(phi)));
	Length1 = (1 + 1 * 1.25)*Radius / (sin(phi));
	theta1 = Pi - 2.0*phi;
	V2 = Vector3(1, 0, 1 / (tan(phi)));
	Length2 = 4 * 1.25*Radius / (cos(0.5*theta1));
	theta2 = 0.5*theta1;
	V3 = Vector3(-1, 0, 0);
	return true;
}
//˫�����·��ֵ�������

bool calculation::SourceModeGenerationT::GetCircularSystemOpticalPath_TE0n_3Mirrors(
	Vector3 & P4, Vector3 & V1, double & Length1, double & theta1,
	                 Vector3 & V2, double & Length2, double & theta2,
	                 Vector3 & V3, double & Length3, double & theta3, Vector3 & V4)
{
	if ((SourceKind != 1) && (SourceType != 1) && (m != 0)) return false;

	P4 = Vector3(Radius, 0, 0);
	double phi; double theta; double Rc;
	if (GetCircularWaveguideProperty(phi, theta, Rc,Lc))//�����S���Ҽ��һ��
		V1 = Vector3(-1, 0, 1 / (tan(phi)));
	Length1 = (1 + 1 * 1.25)*Radius / (sin(phi));
	theta1 = Pi - 2.0*phi;
	V2 = Vector3(1, 0, 1 / (tan(phi)));
	Length2 = 4 * 1.25*Radius / (cos(0.5*theta1));
	theta2 = 0.75*theta1;
	V3 = Vector3(-1, 0, 1*tan(0.25*theta1));
	Length3= 7 * 1.25*Radius / (cos(0.25*theta1));
	theta3 = 0.25*theta1;
	V4 = Vector3(1, 0, 0);
	return true;
}
//�������·��ֵ�������

bool calculation::SourceModeGenerationT::GetCircularSystemOpticalPath_TE0n_4Mirrors(
	Vector3 & P4, Vector3 & V1, double & Length1, double & theta1,
	                 Vector3 & V2, double & Length2, double & theta2,
	                 Vector3 & V3, double & Length3, double & theta3,
	                 Vector3 & V4, double & Length4, double & theta4, Vector3 & V5)
{
	if ((SourceKind != 1) && (SourceType != 1) && (m != 0)) return false;

	P4 = Vector3(Radius, 0, 0);
	double phi; double theta; double Rc;
	if (GetCircularWaveguideProperty(phi, theta, Rc,Lc))//�����S���Ҽ��һ��
		V1 = Vector3(-1, 0, 1 / (tan(phi)));
	Length1 = (1 + 1 * 1.25)*Radius / (sin(phi));
	theta1 = Pi - 2.0*phi;
	V2 = Vector3(1, 0, 1 / (tan(phi)));
	Length2 = 4 * 1.25*Radius / (cos(0.5*theta1));
	theta2 = 0.75*theta1;
	V3 = Vector3(-1, 0, 1 * tan(0.25*theta1));
	Length3 = 7 * 1.25*Radius / (cos(0.25*theta1));
	theta3 = 0.5*theta1;
	V4 = Vector3(1, 0, 1 * tan(0.25*theta1));
	Length4 = 7 * 1.25*Radius / (cos(0.25*theta1));
	theta4 = 0.25*theta1;
	V4 = Vector3(-1, 0, 0);
	return true;
}
//�ľ����·��ֵ�������20171116


void calculation::SourceModeGenerationT::GetEX(vector<vector<complex<double>>>& EX0)
{
	EX0.resize(N);
	for (int i = 0; i < N; i++)
	{
		EX0[i].resize(N);
	}//���ɾ���ߴ磺N*N
	//���濪ʼ��ֵ
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			EX0[i][j] = EX[i][j];
		}
	}
}

void calculation::SourceModeGenerationT::GetHX(vector<vector<complex<double>>>& HX0)
{
	HX0.resize(N);
	for (int i = 0; i < N; i++)
	{
		HX0[i].resize(N);
	}//���ɾ���ߴ磺N*N
	 //���濪ʼ��ֵ
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			HX0[i][j] = HX[i][j];
		}
	}
}

void calculation::SourceModeGenerationT::GetEY(vector<vector<complex<double>>>& EY0)
{
	EY0.resize(N);
	for (int i = 0; i < N; i++)
	{
		EY0[i].resize(N);
	}//���ɾ���ߴ磺N*N
	 //���濪ʼ��ֵ
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			EY0[i][j] = EY[i][j];
		}
	}
}

void calculation::SourceModeGenerationT::GetHY(vector<vector<complex<double>>>& HY0)
{
	HY0.resize(N);
	for (int i = 0; i < N; i++)
	{
		HY0[i].resize(N);
	}//���ɾ���ߴ磺N*N
	 //���濪ʼ��ֵ
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			HY0[i][j] = HY[i][j];
		}
	}
}

void calculation::SourceModeGenerationT::GetEZ(vector<vector<complex<double>>>& EZ0)
{
	EZ0.resize(N);
	for (int i = 0; i < N; i++)
	{
		EZ0[i].resize(N);
	}//���ɾ���ߴ磺N*N
	 //���濪ʼ��ֵ
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			EZ0[i][j] = EZ[i][j];
		}
	}
}

void calculation::SourceModeGenerationT::GetHZ(vector<vector<complex<double>>>& HZ0)
{
	HZ0.resize(N);
	for (int i = 0; i < N; i++)
	{
		HZ0[i].resize(N);
	}//���ɾ���ߴ磺N*N
	 //���濪ʼ��ֵ
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			HZ0[i][j] = HZ[i][j];
		}
	}
}

