#include "PhaseCorrection.h"

calculation::PhaseCorrection::PhaseCorrection(double _Frequency, double _Lamda, double _Target_W0,
	int _N_Mirror, int _N_InField, double _InField_ds) :
	Frequency(_Frequency), Lamda(_Lamda), Target_W0(_Target_W0),
	N_Mirror(_N_Mirror), N_InField(_N_InField), InField_ds(_InField_ds)
{
	D_in = 0;
	D_out = 0;
	Theta = 0;
	//��ʼ��
	/*
	Initial_Mirror.resize(N_Mirror); In_Field.resize(N_Input_Field);
	for (int i = 0; i < N_Mirror; i++)
	{
		Initial_Mirror[i].resize(N_Mirror);
	}//�����resize

	for (int i = 0; i < N_Input_Field; i++)
	{
		In_Field[i].resize(N_Input_Field);
	}//�����resize
	*/
}

calculation::PhaseCorrection::~PhaseCorrection()
{
}

//������������
void calculation::PhaseCorrection::Set(double _Frequency, 
	double _Target_W0, int _N_Mirror, int _N_InField,
	const vector <vector <Vector3>> &_Initial_Mirror,
	const vector <vector <complex<double>>> &_InField_EX,
	const vector <vector <complex<double>>> &_InField_EY,
	double _InField_ds)
{
	Frequency = _Frequency;
	Lamda = C_Speed / Frequency;
	Target_W0 = _Target_W0;

	N_Mirror = _N_Mirror;
	N_InField = _N_InField;
	InField_ds = _InField_ds;

	//���S���ݸ��ҵľ������N*N�ģ���˲���Ҫresize,ֱ��ʹ�þ�����
	Initial_Mirror = _Initial_Mirror;
	InField_EX = _InField_EX;
	InField_EY = _InField_EY;

	/*ofstream outfilex1("PhsCorExEy.txt");
	for (int i = 0; i < 201; i++)
	{
		for (int j = 0; j < 201; j++)
		{
			outfilex1 //<< abs(InField_EX[i][j]) << " " 
				//<< arg(InField_EX[i][j]) << " "
				<< abs(InField_EY[i][j]) << " "
				<< arg(InField_EY[i][j]) << "\n";
		}
	}
	outfilex1.close();*/

	InField_EX1.resize(N_InField); InField_EY1.resize(N_InField);
	for (int i = 0; i < N_InField; i++)
	{
		InField_EX1[i].resize(N_InField);
		InField_EY1[i].resize(N_InField);
	}//�����resize
}

//��������ϵ�������ã��ڴ�֮ǰҪ������Set������
void calculation::PhaseCorrection::SetCoordinate(Vector3 _InField_Position, Vector3 _InField_X, Vector3 _InField_Y, Vector3 _InField_Z,
	Vector3 _Mirror_Position, Vector3 _Mirror_X, Vector3 _Mirror_Y, Vector3 _Mirror_Z,
	Vector3 _Reflect_Point, Vector3 _Reflect_Point_Normal)
{
	InField_Position = _InField_Position;
	InField_X = _InField_X;InField_Y = _InField_Y;InField_Z = _InField_Z;

	Mirror_Position = _Mirror_Position;
	Mirror_X = _Mirror_X; Mirror_Y = _Mirror_Y; Mirror_Z = _Mirror_Z;

	Reflect_Point = _Reflect_Point;
	Reflect_Point_Normal = _Reflect_Point_Normal;

	D_in = (_Reflect_Point - InField_Position).Length();//������ȥ���䳡���ĵ㣬�õ�D_in
	//D_out = _D_out;

	//����������ڷ����淴����ϵķ���Ƕ�Thetaֵ
	Theta = acos(InField_Z.Dot(Reflect_Point_Normal));
	if (Theta > (Pi / 2.0)) Theta = Pi - Theta;//���뱣֤ThetaС��Pi/2

	OutField_Z = InField_Z - Reflect_Point_Normal.operator*(2.0*(InField_Z.Dot(Reflect_Point_Normal)));
	OutField_Y = OutField_Z.Cross(InField_Z);
	OutField_Y.Normalization();
	OutField_X = OutField_Y.Cross(OutField_Z);
	OutField_X.Normalization();//�õ���������ľֲ�����ϵ

	//������Ҫ�漰�������������任���⣬Ҫ�����S����
	//����Ҫ���������ά����ת��Ϊ���볡�ľֲ�����ϵ�º������˹�����ľֲ�����ϵ��
	//���Ҫ��֤��������Ey,�����ζ����Ҫ�����볡����������������任������InField_X1��InField_Y1����InField_Z����
	//�ҵ�InField_EX��InField_EY���ֵ�ĵ㣬ȷ��������
	double MaxEX = 0, MaxEY = 0;//��ʼ��
	//int MaxEX_i = 0, MaxEY_i = 0;
	//int MaxEX_j = 0, MaxEY_j = 0;
	int ii = 0, jj = 0;

	//�����ҵ���ǿ���ֵ,Ȼ���ҵ���ǿ���ֵ���ڵ�ii,jj
	for (int i = 0; i < N_InField; i++)
		for (int j = 0; j < N_InField; j++)
		{
			{
				if (abs(InField_EX[i][j]) > MaxEX)
				{
					MaxEX = abs(InField_EX[i][j]);
					//MaxEX_i = i;
					//MaxEX_j = j;
				}

				if (abs(InField_EY[i][j]) > MaxEY)
				{
					MaxEY = abs(InField_EY[i][j]);
					//MaxEY_i = i;
					//MaxEY_j = j;
				}
			}
		}
	//���濪ʼ�Ա�
	//ii = MaxEX_i; jj = MaxEX_j;



	//20180314Ȩ��֮�ƣ�ֻȡExEy�еĽϴ�ֵ��������ϵΪEY1�����ں�������
	if (MaxEY > MaxEX)
	{
		//ii = MaxEY_i; jj = MaxEY_j;
		InField_Z1 = InField_Z;//Z���겻��
		InField_Y1 = InField_Y;
		InField_X1 = InField_X;
		for (int i = 0; i < N_InField; i++)
			for (int j = 0; j < N_InField; j++)
			{
				{
					InField_EX1[i][j] = InField_EX[i][j];
					InField_EY1[i][j] = InField_EY[i][j];
				}
			}
	}

	if (MaxEX >= MaxEY)
	{
		InField_Z1 = InField_Z;//Z���겻��
		InField_Y1 = InField_X;
		InField_X1 = InField_Y1.Cross(InField_Z1);
		for (int i = 0; i < N_InField; i++)
			for (int j = 0; j < N_InField; j++)
			{
				{
					InField_EY1[i][j] = InField_EX[i][j];
					InField_EX1[i][j] = -InField_EY[i][j];//������
				}
			}
	}

	/*
	//����Ӧ��EX��EY�ֿ��㣬������Ȩ��֮��
	//��������һ�����InField_X1,InField_Y1��InField_Z1������
	InField_Z1 = InField_Z;//Z���겻��
	//��InField_Y1������Ϊ����������
	InField_Y1 = InField_X.operator*(abs(InField_EX[ii][jj])) 
		       + InField_Y.operator*(abs(InField_EY[ii][jj]));
	InField_Y1.Normalization();//��һ��

	InField_X1 = InField_Y1.Cross(InField_Z1);//InField_Y1���InField_Z1�õ�InField_X1

	//�������ɶ�Ӧ��InField_EX1��InField_EY1
	double RatioX_X1 = InField_X.Dot(InField_X1); double RatioX_Y1 = InField_X.Dot(InField_Y1);
	double RatioY_X1 = InField_Y.Dot(InField_X1); double RatioY_Y1 = InField_Y.Dot(InField_Y1);
	for (int i = 0; i < N_InField; i++)
		for (int j = 0; j < N_InField; j++)
		{
			{
				InField_EX1[i][j] = RatioX_X1*InField_EX[i][j] + RatioY_X1*InField_EY[i][j];
				InField_EY1[i][j] = RatioX_Y1*InField_EX[i][j] + RatioY_Y1*InField_EY[i][j];
			}
		}
	*/
}

//������Ҫ�����ľ���
void calculation::PhaseCorrection::SetMirror(vector<vector<Vector3>>& _Mirror)
{
	Initial_Mirror = _Mirror;
}

//��������,ֻ����һ�Σ��������������µľ���ΪInitial_Mirror������ظ���
vector<vector<Vector3>> calculation::PhaseCorrection::Single_Correction()
{
	//���ȣ����ݳ����˹������С�����䲨��Ϣ�õ�D_out
	//�ȷֽ����䳡��������䴫���������洦������ֵ
	int Infield_splittimes = 0;
	vector <vector <double>> Infield_splitinfo;
	Field_Split(Frequency, InField_ds, N_InField, InField_EY1, Infield_splittimes, Infield_splitinfo);
	//ֻ����������InField_EY1

	//�����������
	vector<Vector3> Target_LineX(N_InField), Target_LineY(N_InField);
	vector<double> Target_LineX_E(N_InField), Target_LineY_E(N_InField);
	for (int i = 0; i < N_InField; i++)
	{
		Target_LineX[i] = Vector3((i - (N_InField - 1) / 2.0)*InField_ds, 0, D_in);
		Target_LineX_E[i] = abs(Calculate_SinglePoint(Frequency, InField_ds, 
			N_InField, Infield_splittimes, Infield_splitinfo, Target_LineX[i]));

		Target_LineY[i] = Vector3(0, (i - (N_InField - 1) / 2.0)*InField_ds, D_in);
		Target_LineY_E[i] = abs(Calculate_SinglePoint(Frequency, InField_ds,
			N_InField, Infield_splittimes, Infield_splitinfo, Target_LineY[i]));

	}



	//�ҵ�Target_LineX��Target_LineY�����ֵ��
	double Target_LineX_MaxE = 0, Target_LineY_MaxE = 0;
	int Target_LineX_MaxPoint, Target_LineY_MaxPoint;
	for (int i = 0; i < N_InField; i++)
	{
		if (Target_LineX_E[i] > Target_LineX_MaxE)
		{
			Target_LineX_MaxE = Target_LineX_E[i];
			Target_LineX_MaxPoint = i;
		}
	}

	for (int j = 0; j < N_InField; j++)
	{
		if (Target_LineY_E[j] > Target_LineY_MaxE)
		{
			Target_LineY_MaxE = Target_LineY_E[j];
			Target_LineY_MaxPoint = j;
		}
	}

	int iix0 = 0;
	for (int k = 0; k < Target_LineX_MaxPoint; k++)
	{
		iix0 = k;
		if ((Target_LineX_E[Target_LineX_MaxPoint - k - 1]<= (1.0 / Universal_Constant_e*Target_LineX_MaxE))
	        & (Target_LineX_E[Target_LineX_MaxPoint - k]>(1.0 / Universal_Constant_e*Target_LineX_MaxE)))
			break;//һ�ҵ�������ѭ��
	}

	int iix1 = Target_LineX_MaxPoint;
	for (int k = Target_LineX_MaxPoint; k < N_InField - 1; k++)
	{
		iix1 = k;
		if ((Target_LineX_E[k + 1] <= (1.0 / Universal_Constant_e*Target_LineX_MaxE))
			& (Target_LineX_E[k]>(1.0 / Universal_Constant_e*Target_LineX_MaxE)))
			break;
	}

	int jjy0 = 0;
	for (int k = 0; k < Target_LineY_MaxPoint; k++)
	{
		jjy0 = k;
		if ((Target_LineY_E[Target_LineY_MaxPoint - k - 1] <= (1.0 / Universal_Constant_e*Target_LineY_MaxE))
			& (Target_LineY_E[Target_LineY_MaxPoint - k]>(1.0 / Universal_Constant_e*Target_LineY_MaxE)))
			break;//һ�ҵ�������ѭ��
	}

	int jjy1 = Target_LineY_MaxPoint;
	for (int k = Target_LineY_MaxPoint; k < N_InField - 1; k++)
	{
		jjy1 = k;
		if ((Target_LineY_E[k + 1] <= (1.0 / Universal_Constant_e*Target_LineY_MaxE))
			& (Target_LineY_E[k]>(1.0 / Universal_Constant_e*Target_LineY_MaxE)))
			break;//һ�ҵ�������ѭ��
	}
	//�õ��ڷ����洦�����䳡�����Ĺ���ֵ
	double W_TargetPoint = InField_ds*(iix0 + (iix1 - Target_LineX_MaxPoint)
		                    + jjy0 + (jjy1 - Target_LineY_MaxPoint))/4.0;

	//�ٶ��ڷ���㴦�����䲨�������ڷ��䲨�����������򴫲��ĸ�˹�����ڴ˴�����Ҳ��W_TargetPoint
	//�ɴ˿��Ը��ݸ�˹����������ʽ����D_out
	D_out = Pi*Target_W0 / Lamda*pow(W_TargetPoint*W_TargetPoint - Target_W0*Target_W0, 0.5);

	//���濪ʼ��˹�������򴫲��������ϵ���λ
	Vector3 Gauss_Position = Reflect_Point + OutField_Z.operator*(D_out);
	vector <vector<double>> Unwrapped_Gauss_Phase(N_Mirror, vector<double>(N_Mirror, 0));
	for (int i = 0; i < N_Mirror; i++)
		for (int j = 0; j < N_Mirror; j++)
		{
			{
				Vector3 Position1 = Mirror_Position
					+ Mirror_X.operator*(Initial_Mirror[i][j].x)
					+ Mirror_Y.operator*(Initial_Mirror[i][j].y)
					+ Mirror_Z.operator*(Initial_Mirror[i][j].z);

				Vector3 Position2((Position1.operator-(Gauss_Position)).Dot(OutField_X),
					(Position1.operator-(Gauss_Position)).Dot(OutField_Y),
					(Position1.operator-(Gauss_Position)).Dot(OutField_Z));
				Unwrapped_Gauss_Phase[i][j] = Gauss_Phase_Circular(Frequency, Target_W0, Position2);
			}
		}


	//������In_Field���������ϵ���λ��������н��۵�
	vector <vector<double>> Mirror_Amplitude(N_Mirror, vector<double>(N_Mirror, 0));
	vector <vector<double>> Wrapped_Mirror_Phase(N_Mirror, vector<double>(N_Mirror, 0));
	for (int i = 0; i < N_Mirror; i++)
		for (int j = 0; j < N_Mirror; j++)
		{
			{
				Vector3 Position1 = Mirror_Position
					+ Mirror_X.operator*(Initial_Mirror[i][j].x)
					+ Mirror_Y.operator*(Initial_Mirror[i][j].y)
					+ Mirror_Z.operator*(Initial_Mirror[i][j].z);
				Vector3 Position2((Position1.operator-(InField_Position)).Dot(InField_X1),
					(Position1.operator-(InField_Position)).Dot(InField_Y1),
					(Position1.operator-(InField_Position)).Dot(InField_Z1));
				//20180322����
				Mirror_Amplitude[i][j]= abs(Calculate_SinglePoint(Frequency, InField_ds, N_InField,
					Infield_splittimes, Infield_splitinfo, Position2));
				Wrapped_Mirror_Phase[i][j] = arg(Calculate_SinglePoint(Frequency, InField_ds, N_InField,
					Infield_splittimes, Infield_splitinfo, Position2));
			}
		}

	vector <vector<double>> Unwrapped_Mirror_Phase(N_Mirror, vector<double>(N_Mirror, 0));
	//�����λ����Unwrap_2D����
	ofstream OutFile0;
	OutFile0.open("Original_Distribution.txt");//��������ļ�
	for (int i = 0; i < N_Mirror; i++)
		for (int j = 0; j < N_Mirror; j++)
		{
		  OutFile0 << Mirror_Amplitude[i][j] << " " << Wrapped_Mirror_Phase[i][j] << endl;
		}
	OutFile0.close();//�����ϣ��ر��ĵ�
	

	Unwrapped_Mirror_Phase= Unwrap_2D(N_Mirror, Wrapped_Mirror_Phase);
	//�����λ�Ľ��۵�

	ofstream OutFile1;
	OutFile1.open("Unwrapped_Mirror_Phase.txt");//��������ļ�
	for (int i = 0; i < N_Mirror; i++)
	for (int j = 0; j < N_Mirror; j++)
   {
	OutFile1 << Unwrapped_Mirror_Phase[i][j] << endl;
	}
	OutFile1.close();//�����ϣ��ر��ĵ�
	

	//���ڳ���������ֻ����һ�Σ�����ʹ�ù̶������Theta��������
	vector <vector<double>> Delta(N_Mirror, vector<double>(N_Mirror, 0));
	for (int i = 0; i < N_Mirror; i++)
		for (int j = 0; j < N_Mirror; j++)
		{
			{
				Delta[i][j] = (Unwrapped_Gauss_Phase[i][j] - Unwrapped_Mirror_Phase[i][j])
					*Lamda / 4.0 / Pi / cos(Theta);
			}
		}
	
	double Delta_Middle = Delta[floor(N_Mirror / 2.0)][floor(N_Mirror / 2.0)];
	vector <vector<Vector3>> Corrected_Mirror(N_Mirror, vector<Vector3>(N_Mirror, 0));
	for (int i = 0; i < N_Mirror; i++)
		for (int j = 0; j < N_Mirror; j++)
		{
			{
				Delta[i][j] = Delta[i][j]- Delta_Middle;//ʹ�м��Deltaֵ����
				Corrected_Mirror[i][j] = Initial_Mirror[i][j] + Vector3(0, 0, Delta[i][j]);
			}
		}
	return Corrected_Mirror;
}






