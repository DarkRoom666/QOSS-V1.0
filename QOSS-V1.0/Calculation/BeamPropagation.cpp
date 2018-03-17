#include "BeamPropagation.h"

//���ش���z���Ⱥ�ĸ�˹����Wֵ
void Field_Split(double frequency0, double ds, int N0,
	vector <vector <complex <double>>> &E0,
	int &Actual_SplitTimes, vector <vector <double>> &Split_Info)
{
	double lamda = C_Speed / frequency0;//���㲨��

	complex <double> j0(0, 1);//����������λj0

	//���������N0*N0����,�����S��������˵
	vector <vector <complex<double>>> E1(N0, vector<complex<double>>(N0, 0));
	vector <vector <double>> AbsE1(N0, vector<double>(N0, 0));
	vector <vector <double>> PhaseE1(N0, vector<double>(N0, 0));

	double Whole_Power0 = 0;//������泡��������,��ʼ��Ϊ0

	for (int i = 0; i < N0; i++)
		for (int j = 0; j < N0; j++)
		{
			{
				//complex<double> mycomplex(2.000, 2);
				//AbsE0[i][j] = abs(mycomplex);

				E1[i][j] = E0[i][j];

				AbsE1[i][j] = abs(E0[i][j]);

				PhaseE1[i][j] = arg(E0[i][j]);

				Whole_Power0 = Whole_Power0 + AbsE1[i][j] * AbsE1[i][j];
			}
		}//�õ���λ�ͷ�ֵ�ֲ�(E1��������)
	//20180306������ˣ�û����
	//20180307ȡ��E0��ֵ����������ֻ��E1

	//���濪ʼ�����ֲ��ֽ�Ϊһϵ�еĸ�˹����

	//����������ֽ���Ϣ
	int Max_SplitTimes = 300;//�������ֽ����

	Actual_SplitTimes = 0;//��ʵ�ʷֽ������ʼ��Ϊ0

	Split_Info.resize(Max_SplitTimes);
	for (int i = 0; i < Max_SplitTimes; i++)
	{
		Split_Info[i].resize(15);
	}//Split_Info�Ƿֽ�����Ϣ��Max_SplitTimes*15�������е��ηֽ����15��������Ϣ

	double Remained_Power_Ratio = 0.002;//����ʣ��������ֵ��ʣ����������0.2%ʱֹͣ�ֽ�

	//Ϊ������ԣ����һ��С����
	vector <double> Power_Remained(Max_SplitTimes);

	//20180307��һЩ�������ǰ���壬��ʡʱ��
	//P1����MaxP��Ϊ���ĵĳ��ֲ�����ϵ
	vector <vector <Vector3>> P1(N0, vector<Vector3>(N0, 0));

	for (int k = 0; k < Max_SplitTimes; k++)
	{
		//Ѱ�ҳ�ǿ���ֵ��λ��
		double MaxE = 0;//��ʼ��
		int MaxP_i = 0;//maxposition
		int MaxP_j = 0;

		//�����ҵ���ǿ���ֵ,Ȼ���ҵ���ǿ���ֵ���ڵ�ii,jj
		for (int i = 0; i < N0; i++)
			for (int j = 0; j < N0; j++)
			{
				{
					if (AbsE1[i][j] > MaxE)
					{
						MaxE = AbsE1[i][j];
						MaxP_i = i;
						MaxP_j = j;
					}

				}
			}

		//�������ݷ������ֵ��ĵ�����
		Vector3 Ni, Nj;

		//�ȼ���Ni
		double a = 0, b = 0;
		if (MaxP_i > 0 & MaxP_i < N0 - 1)
		{
			a = PhaseE1[MaxP_i - 1][MaxP_j];
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i - 1][MaxP_j])>0.75 * 2.0 * Pi)
				a = PhaseE1[MaxP_i - 1][MaxP_j] + 2.0 * Pi;
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i - 1][MaxP_j]) < -0.75 * 2.0 * Pi)
				a = PhaseE1[MaxP_i - 1][MaxP_j] - 2.0 * Pi;//������λ�۵�

			b = PhaseE1[MaxP_i + 1][MaxP_j];
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i + 1][MaxP_j]) > 0.75 * 2.0 * Pi)
				b = PhaseE1[MaxP_i + 1][MaxP_j] + 2.0 * Pi;
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i + 1][MaxP_j]) < -0.75 * 2.0 * Pi)
				b = PhaseE1[MaxP_i + 1][MaxP_j] - 2.0 * Pi;//������λ�۵�
		}

		if (MaxP_i == 0)
		{
			a = PhaseE1[MaxP_i][MaxP_j];
			b = PhaseE1[MaxP_i + 2][MaxP_j];
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i + 2][MaxP_j]) > 0.75 * 2.0 * Pi)
				b = PhaseE1[MaxP_i + 2][MaxP_j] + 2.0 * Pi;
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i + 2][MaxP_j]) < -0.75 * 2.0 * Pi)
				b = PhaseE1[MaxP_i + 2][MaxP_j] - 2.0 * Pi;//������λ�۵�
		}

		if (MaxP_i == N0 - 1)
		{
			a = PhaseE1[MaxP_i - 2][MaxP_j];
			b = PhaseE1[MaxP_i][MaxP_j];
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i - 2][MaxP_j]) > 0.75 * 2.0 * Pi)
				b = PhaseE1[MaxP_i][MaxP_j] - 2.0 * Pi;
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i - 2][MaxP_j]) < -0.75 * 2.0 * Pi)
				b = PhaseE1[MaxP_i][MaxP_j] + 2.0 * Pi;//������λ�۵�
		}

		if (((b - a) / 2.0 / Pi*lamda) < (2.0*ds))
			Ni = Vector3(cos(Pi / 2.0 - acos((b - a) / 2.0 / Pi*lamda / 2.0 / ds)),
				0,
				sin(Pi / 2.0 - acos((b - a) / 2.0 / Pi*lamda / 2.0 / ds)));

		if (((b - a) / 2.0 / Pi*lamda) >= (2.0*ds))
			Ni = Vector3(1, 0, 0);
		//Ni�������

		//�������Nj����
		double c = 0, d = 0;
		if (MaxP_j > 0 & MaxP_j < N0 - 1)
		{
			c = PhaseE1[MaxP_i][MaxP_j - 1];
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i][MaxP_j - 1])>0.75 * 2.0 * Pi)
				c = PhaseE1[MaxP_i][MaxP_j - 1] + 2.0 * Pi;
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i][MaxP_j - 1]) < -0.75 * 2.0 * Pi)
				c = PhaseE1[MaxP_i][MaxP_j - 1] - 2.0 * Pi;//������λ�۵�

			d = PhaseE1[MaxP_i][MaxP_j + 1];
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i][MaxP_j + 1]) > 0.75 * 2.0 * Pi)
				d = PhaseE1[MaxP_i][MaxP_j + 1] + 2.0 * Pi;
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i][MaxP_j + 1]) < -0.75 * 2.0 * Pi)
				d = PhaseE1[MaxP_i][MaxP_j + 1] - 2.0 * Pi;//������λ�۵�
		}

		if (MaxP_j == 0)
		{
			c = PhaseE1[MaxP_i][MaxP_j];
			d = PhaseE1[MaxP_i][MaxP_j + 2];
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i][MaxP_j + 2]) > 0.75 * 2.0 * Pi)
				d = PhaseE1[MaxP_i][MaxP_j + 2] + 2.0 * Pi;
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i][MaxP_j + 2]) < -0.75 * 2.0 * Pi)
				d = PhaseE1[MaxP_i][MaxP_j + 2] - 2.0 * Pi;//������λ�۵�
		}

		if (MaxP_j == N0 - 1)
		{
			c = PhaseE1[MaxP_i][MaxP_j - 2];
			d = PhaseE1[MaxP_i][MaxP_j];
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i][MaxP_j - 2]) > 0.75 * 2.0 * Pi)
				d = PhaseE1[MaxP_i][MaxP_j] - 2.0 * Pi;
			if ((PhaseE1[MaxP_i][MaxP_j] - PhaseE1[MaxP_i][MaxP_j - 2]) < -0.75 * 2.0 * Pi)
				d = PhaseE1[MaxP_i][MaxP_j] + 2.0 * Pi;//������λ�۵�
		}

		if (((d - c) / 2.0 / Pi*lamda) < (2.0*ds))
			Nj = Vector3(0,
				cos(Pi / 2.0 - acos((d - c) / 2.0 / Pi*lamda / 2.0 / ds)),
				sin(Pi / 2.0 - acos((d - c) / 2.0 / Pi*lamda / 2.0 / ds)));

		if (((d - c) / 2.0 / Pi*lamda) >= (2.0*ds))
			Nj = Vector3(0, 1, 0);
		//Nj�������

		//���濪ʼ�������ǿ��Ĵ���������������һ��
		Vector3 MaxDirection;
		MaxDirection = Vector3(-Ni.z / Ni.x, -Nj.z / Nj.y, 1);
		MaxDirection.Normalization();

		//�µľֲ�����ϵ������������
		Vector3 Ix, Iy, Iz;
		Iz = MaxDirection;
		Iy = Vector3(0,
			MaxDirection.z / pow((MaxDirection.y*MaxDirection.y + MaxDirection.z*MaxDirection.z),0.5),
			-MaxDirection.y / pow((MaxDirection.y*MaxDirection.y + MaxDirection.z*MaxDirection.z), 0.5));
		Ix = Iy.Cross(Iz);//���

		//P1����MaxP��Ϊ���ĵĳ��ֲ�����ϵ��ͬʱ��б
		for (int i = 0; i < N0; i++)
			for (int j = 0; j < N0; j++)
			{
				{
					P1[i][j] = Vector3(ds*(i - (N0 - 1) / 2), ds*(j - (N0 - 1) / 2), 0) -
						Vector3(ds*(MaxP_i - (N0 - 1) / 2), ds*(MaxP_j - (N0 - 1) / 2), 0);
					P1[i][j] = Vector3(P1[i][j].Dot(Ix),
						P1[i][j].Dot(Iy),
						P1[i][j].Dot(Iz));
				}
			}

		//���濪ʼ��������
		double Dx_limit, Dy_limit;
		if (MaxP_i > 0 & MaxP_i < N0 - 1)
		{
			Dx_limit = ds*(N0 - 1 - MaxP_i);
			if (Dx_limit > ds*(MaxP_i - 0))
				Dx_limit = ds*(MaxP_i - 0);//�������ȡС��
		}
		else  Dx_limit = 3 * ds;

		if (MaxP_j > 0 & MaxP_j < N0 - 1)
		{
			Dy_limit = ds*(N0 - 1 - MaxP_j);
			if (Dy_limit > ds*(MaxP_j - 0))
				Dy_limit = ds*(MaxP_j - 0);//�������ȡС��
		}
		else  Dy_limit = 3 * ds;

		double W0x_max = Vector3(Dx_limit / 2.5, 0, 0).Dot(Ix);
		W0x_max = fabs(W0x_max);//fabs��������С���������ֵ

		double W0y_max = Vector3(0, Dy_limit / 2.5, 0).Dot(Iy);
		W0y_max = fabs(W0y_max);//fabs��������С���������ֵ

		//���濪ʼ����X�����Y���������
		//�ӷ�ֵ�������߿�ȥ��Ѱ�ҵ�һ����ֵ
		int ii0 = 0;
		if (MaxP_i == 0) ii0 = 0;
		else
		{
			for (int k = 0; k < MaxP_i; k++)
			{
				ii0 = k;
				if (AbsE1[MaxP_i - k - 1][MaxP_j]>AbsE1[MaxP_i - k][MaxP_j])
					break;//һ�ҵ�������ѭ��
			}
		}

		int ii1 = MaxP_i;
		if (MaxP_i == N0 - 1) ii0 = N0 - 1;
		else
		{
			for (int k = MaxP_i; k < N0 - 1; k++)
			{
				ii1 = k;
				if (AbsE1[k + 1][MaxP_j]>AbsE1[k][MaxP_j])
					break;//һ�ҵ�������ѭ��
			}
		}

		int jj0 = 0;
		if (MaxP_j == 0) jj0 = 0;
		else
		{
			for (int k = 0; k < MaxP_j; k++)
			{
				jj0 = k;
				if (AbsE1[MaxP_i][MaxP_j - k - 1]>AbsE1[MaxP_i][MaxP_j - k])
					break;//һ�ҵ�������ѭ��
			}
		}

		int jj1 = MaxP_j;
		if (MaxP_j == N0 - 1) jj0 = N0 - 1;
		else
		{
			for (int k = MaxP_j; k < N0 - 1; k++)
			{
				jj1 = k;
				if (AbsE1[MaxP_i][k + 1]>AbsE1[MaxP_i][k])
					break;//һ�ҵ�������ѭ��
			}
		}
		//����Wx
		double Dx1 = fabs(Vector3(ii0*ds, 0, 0).Dot(Ix));
		double Dx2 = fabs(Vector3((ii1 - MaxP_i)*ds, 0, 0).Dot(Ix));
		double Wx1 = sqrt(Dx1*Dx1 / (1.0e-10 - log(AbsE1[MaxP_i - ii0][MaxP_j] / AbsE1[MaxP_i][MaxP_j])));
		double Wx2 = sqrt(Dx2*Dx2 / (1.0e-10 - log(AbsE1[ii1][MaxP_j] / AbsE1[MaxP_i][MaxP_j])));
		double Wx;
		if (fabs(Wx1) > 1.0e-5 & fabs(Wx2) > 1.0e-5)
		{
			Wx = Wx1;
			if (Wx1 > Wx2)
				Wx = Wx2;//ȡ��С����ֵ
		}
		else Wx = Wx1 + Wx2;

		if (Wx > W0x_max) Wx = W0x_max;

		//����Wy
		double Dy1 = fabs(Vector3(0, jj0*ds, 0).Dot(Iy));
		double Dy2 = fabs(Vector3(0, (jj1 - MaxP_j)*ds, 0).Dot(Iy));
		double Wy1 = sqrt(Dy1*Dy1 / (1.0e-10 - log(AbsE1[MaxP_i][MaxP_j - jj0] / AbsE1[MaxP_i][MaxP_j])));
		double Wy2 = sqrt(Dy2*Dy2 / (1.0e-10 - log(AbsE1[MaxP_i][jj1] / AbsE1[MaxP_i][MaxP_j])));
		double Wy;
		if (fabs(Wy1) > 1.0e-5 & fabs(Wy2) > 1.0e-5)
		{
			Wy = Wy1;
			if (Wy1 > Wy2)
				Wy = Wy2;//ȡ��С����ֵ
		}
		else Wy = Wy1 + Wy2;

		if (Wy > W0y_max) Wy = W0y_max;

		//��ȷ����λ���ֵ�ı���ϵ��
		double Delta_Phase = PhaseE1[MaxP_i][MaxP_j]
			- arg(Gauss_E_Elliptical(frequency0, Wx, Wy, P1[MaxP_i][MaxP_j]));
		double Amplitude_Ratio = MaxE 
			/ abs(Gauss_E_Elliptical(frequency0, Wx, Wy, P1[MaxP_i][MaxP_j]))
			/ (Iy.Dot(Vector3(0, 1, 0)));

		//���濪ʼ���˴ηֽ����Ϣ���뵽Split_Info��
		Split_Info[k][0] = MaxP_i; Split_Info[k][1] = MaxP_j;
		Split_Info[k][2] = Ix.x; Split_Info[k][3] = Ix.y; Split_Info[k][4] = Ix.z;
		Split_Info[k][5] = Iy.x; Split_Info[k][6] = Iy.y; Split_Info[k][7] = Iy.z;
		Split_Info[k][8] = Iz.x; Split_Info[k][9] = Iz.y; Split_Info[k][10] = Iz.z;
		Split_Info[k][11] = Wx; Split_Info[k][12] = Wy;
		Split_Info[k][13] = Amplitude_Ratio; Split_Info[k][14] = Delta_Phase;
		//Split_Info.push_back(Split_Info_Single);//��Split_Info_Single����Split_Info��

		Actual_SplitTimes = k + 1;//ʵ�ʷֽ����+1

		double Whole_Power1 = 0;//������º�Ŀ��泡������,��ʼ��Ϊ0

		//���濪ʼ��ԭʼ���ֲ�����ȥ�˴ηֽ���ĸ�˹����
		//������ʣ������������ֽ���ֵ�ж�
		for (int i = 0; i < N0; i++)
			for (int j = 0; j < N0; j++)
			{
				{
					//��ȥ���ɵĵ�һ����Բ��˹�ֲ�����˹�ֲ�����ԭ�г��ֲ�ƥ����λ�ͷ��ȣ�
					E1[i][j] = E1[i][j]
						- Amplitude_Ratio*exp(j0*Delta_Phase)
						*Gauss_E_Elliptical(frequency0, Wx, Wy, P1[i][j])
						*(Iy.Dot(Vector3(0,1,0)));//��˹������ǿ��������Y����ķ���

					AbsE1[i][j] = abs(E1[i][j]);
					PhaseE1[i][j] = arg(E1[i][j]);
					//�ڼ�ȥ����Բ��˹�����ֲ��Ļ����ϸ���E1(��ȥһ��)

					Whole_Power1 = Whole_Power1 + abs(E1[i][j])*abs(E1[i][j]);
				}
			}

		Power_Remained[k] = Whole_Power1 / Whole_Power0;

		if (Power_Remained[k] <= Remained_Power_Ratio)
			break;//ʣ������С��Remained_Power_Ratioʱ��ֹͣ�ֽ�

		if (isnan(Power_Remained[k]))
		{
			Actual_SplitTimes = k;
			break;//�ֽ����NaNʱ��ֱ������
		}
			

	}
}

	//�������P���ĳ�ǿֵ
complex<double> Calculate_SinglePoint(double frequency0, double ds, int N0,
	                                  int Actual_SplitTimes, vector <vector <double>> &Split_Info,
	                                  Vector3 P)
{   
	complex <double> j0(0, 1);//����������λj0

	complex<double> E_P=0;
	for (int k = 0; k < Actual_SplitTimes; k++)
	{
		int ii = Split_Info[k][0];
		int jj = Split_Info[k][1];

		Vector3 Gauss_Position,Axis_X, Axis_Y, Axis_Z;
		Gauss_Position=Vector3(ds*(ii - (N0 - 1) / 2), ds*(jj - (N0 - 1) / 2), 0);
		Axis_X = Vector3(Split_Info[k][2], Split_Info[k][3], Split_Info[k][4]);
		Axis_Y = Vector3(Split_Info[k][5], Split_Info[k][6], Split_Info[k][7]);
		Axis_Z = Vector3(Split_Info[k][8], Split_Info[k][9], Split_Info[k][10]);

		Vector3 P1;//���P����б��˹������������ϵ�µ�����
		P1 = Vector3((P - Gauss_Position).Dot(Axis_X),
			         (P - Gauss_Position).Dot(Axis_Y),
			         (P - Gauss_Position).Dot(Axis_Z));
		
		double Wx = Split_Info[k][11]; double Wy = Split_Info[k][12];
		double Amplitude_Ratio = Split_Info[k][13]; double Delta_Phase = Split_Info[k][14];

		E_P = E_P + Gauss_E_Elliptical(frequency0, Wx, Wy, P1) * Amplitude_Ratio * exp(j0 * Delta_Phase);
	}

	return E_P;
}
