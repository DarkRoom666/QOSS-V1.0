#include "PhaseUnwrap.h"

//һά��λ���۵�����ά��λ���۵�������Ҫ�õ�����
//�����Phase�������۵����ģ���arg����atan2�����õ�����Χ�ǣ�-Pi��Pi]
vector <double> Unwrap_1D(int N0, vector <double> Phase)
{
	vector <double> Phase1(N0);
	for (int k = 0; k < N0; k++)
	{
		Phase1[k] = Phase[k];
	}

	for (int k = 0; k < N0-1; k++)
	{
		if ((abs(Phase1[k] - Phase1[k + 1]) >= 0.6 * 2 * Pi) & (Phase1[k] > Phase1[k+1]))
		{
			for (int k1 = k+1; k1 < N0; k1++)
			{
				Phase1[k1] = Phase1[k1] + 2 * Pi;
			}
		}

		if ((abs(Phase1[k] - Phase1[k + 1]) >= 0.6 * 2 * Pi) & (Phase1[k] < Phase1[k + 1]))
		{
			for (int k1 = k + 1; k1 < N0; k1++)
			{
				Phase1[k1] = Phase1[k1] - 2 * Pi;
			}
		}
	}
	return Phase1;
}

//��ά��λ���۵�����,��N0*N0��λ�����������
vector <vector<double>> Unwrap_2D(int N0, vector <vector <double>> Phase)
{
	double Middle = N0 / 2;
	int Middle_Point = floor(Middle);//�ҵ����ĵ�

	//�����ӵ�һ�飬ѡ����۵��Ĺ��м��Middle_Point���м���ʱ��ѡ����λ����������
	vector <vector<double>> Phase1(N0, (vector <double>(N0, 0)));
	int RotationIndex = 0;//���ѡ���о���0��ѡ���о���1
	int Index_i = 0; int Index_j = 0;//��ʼ��
	for(int i = 0; i < N0 - 1; i++)
	{
		if (abs(Phase[Middle_Point][i + 1] - Phase[Middle_Point][i]) >= (2 * 0.6*Pi))
			Index_i = Index_i + 1;//���ĵ������еĲ����������

		if (abs(Phase[i + 1][Middle_Point] - Phase[i][Middle_Point]) >= (2 * 0.6*Pi))
			Index_j = Index_j + 1;//���ĵ������еĲ����������
	}

	if (Index_i >= Index_j)//ѡ�����ĵ�������
	{
		for (int i = 0; i < N0; i++)
			for (int j = 0; j < N0; j++)
			{
				{
					Phase1[i][j] = Phase[i][j];//ֱ�Ӹ���
				}
			}
	}
	else
	{
		RotationIndex = 1;//ѡ���У���Ϊ1
		for (int i = 0; i < N0; i++)
			for (int j = 0; j < N0; j++)
			{
				{
					Phase1[i][j] = Phase[j][i];//��ԭ����ת�Ⱥ���
				}
			}
	}
	
	//���ڶ�άVector�ľ������޷�ֱ�ӷ���ĳһ�У���˶���Phase2ΪPhase1��ת��
	//Phase2��ĳһ�м�ΪPhase1��ĳһ��
	vector <vector<double>> Phase2(N0, (vector <double>(N0, 0)));
	for (int i = 0; i < N0; i++)
		for (int j = 0; j < N0; j++)
	{
		{
			Phase2[i][j] = Phase1[j][i];//���Ʋ�ת��
		}
	}

	//�Ƚ���Middle_Point�����еĽ��۵�
	vector <double> UnwrappedMiddleRow(N0);
	UnwrappedMiddleRow = Unwrap_1D(N0, Phase1[Middle_Point]);
	//��Phase[Middle_Point][Middle_Point]Ϊ��׼ֵУ׼
	double delta1= UnwrappedMiddleRow[Middle_Point] - Phase1[Middle_Point][Middle_Point];
	for (int k = 0; k < N0; k++)
	{
		UnwrappedMiddleRow[k] = UnwrappedMiddleRow[k] - delta1;
	}

	//�����Phase�ĸ����н��н��۵�
	vector <vector<double>> Phase3(N0, (vector <double>(N0, 0)));
	vector <vector<double>> Phase4(N0, (vector <double>(N0, 0)));
	for (int k = 0; k < N0; k++)
	{
		Phase3[k] = Unwrap_1D(N0, Phase2[k]);
	}
	for (int i = 0; i < N0; i++)
		for (int j = 0; j < N0; j++)
		{
			{
				Phase4[i][j] = Phase3[j][i];//Phase3ת�ȵõ�Phase4,Phase4��Ϊ��Phase���۵��Ľ��
			}
		}

	vector <vector<double>> Phase5(N0, (vector <double>(N0, 0)));
	for (int i = 0; i < N0; i++)
		for (int j = 0; j < N0; j++)
		{
			{
				double delta2 = Phase4[Middle_Point][j] - UnwrappedMiddleRow[j];
				Phase5[i][j] = Phase4[i][j] - delta2;
				//������У�����м�����
			}
		}
	//�����������Ķ�ά��λ���۵����õ��ԴֲڵĽ��۵���Ķ�ά��λ�ֲ�Phase5
	
	//���濪ʼ���к�����ƽ��
	//QualityMap��Ѱ��ÿһ���ϵ���λ��������
	vector <vector<int>> QualityMap(N0, (vector <int>(N0, 0)));
	vector <vector<int>> QualityMap_Expanded(N0, (vector <int>(N0, 0)));
	//���������ǵ�1�к�N0-1����λ�����������������ϵ���λ������UnwrappedMiddleRow��
	for (int i = 0; i < N0; i++)
		for (int j = 0; j < N0-1; j++)
		{
			{
				QualityMap[i][j] = 0;//û�в������ͱ�0
				QualityMap_Expanded[i][j] = 0;//û�в������ͱ�0

				if(abs(Phase5[i][j + 1] - Phase5[i][j])>=(2*0.6*Pi))
				QualityMap[i][j] = 1;//�в������ͱ�1
			}
		}

	//ȥ���㷶Χ����չ
	int BasicRange = 7;//���ڲ������������ߵ���չ
	int LittleRange = 5;//���ڲ������������ܵ�С��չ
	for (int i = 0; i < N0; i++)
		for (int j = 0; j < N0; j++)
		{
			{

				if (QualityMap[i][j] == 1)
				{
					//�������Ե
					if (j < BasicRange)//�趨��ֵ
					{
						for (int k = 0; k < j + BasicRange; k++)
						{
							QualityMap_Expanded[i][k] = 1;//����������Χ
						}
					}
					//�����ұ�Ե
					if (j > ((N0 - 1) - BasicRange))//�趨��ֵ
					{
						for (int k = j - BasicRange; k < N0; k++)
						{
							QualityMap_Expanded[i][k] = 1;//����������Χ
						}
					}
					//���м�
					if ((j >= BasicRange) & (j <= ((N0 - 1) - BasicRange)))//�趨��ֵ
					{
						for (int k1 = j - BasicRange; k1 <= j + BasicRange; k1++)
						{
							QualityMap_Expanded[i][k1] = 1;//����������Χ

							if ((i >= LittleRange) & (i <= ((N0 - 1) - LittleRange)))//�����ܵ�С��չ
							{
								for (int k2 = i - LittleRange; k2 <= i + LittleRange; k2++)
								{
									QualityMap_Expanded[k2][k1] = 1;//����������������Χ
								}
							}
						}

					}
				}
					//if������
			}
		}

	for (int i = 0; i < N0; i++)
		{
			QualityMap_Expanded[i][0] = 0;
			QualityMap_Expanded[i][N0-1] = 0;//�����Ե�����ڲ���������Ȼ��Ӱ�������ֵ
			QualityMap_Expanded[Middle_Point][i] = 0;//�м�һ�в����ڲ�����
		}

	//��ʼ��ֵ����Ϊ���������ֲ�ֵ
	//�Ȳ�0��Middle_Point-1��
	for (int i = 0; i < Middle_Point; i++)
		for (int j = 0; j < N0; j++)
		{
			{
				if (QualityMap_Expanded[i][j] == 1)
				{
					double Dleft = 1;//�������ߵĵ�QualityMap_Expanded�϶�����1
					double Ddown = 1;//���������ĵ�QualityMap_Expanded�϶�����1�������������������Ddown������
					double Dright = 1; double Dup = 1; //�������ǳ�ʼ��

					for (int k = j + 1; k < N0; k++)
					{
						if (QualityMap_Expanded[i][k] == 0)
						{
							Dright = k - j;//�����Ҳ��������
							break;//�ҵ�������ѭ��
						}
					}

					for (int k = i + 1; k <= Middle_Point; k++)
					{
						if (QualityMap_Expanded[k][j] == 0)
						{
							Dup = k - i;//�����Ҳ��������
							break;//�ҵ�������ѭ��
						}
					}

					//������������������棬��ʱphase5(i,j)���������������Ծ��뷴�Ȳ�ֵ�õ�
					if (i == 0)
						Phase5[i][j] = 1/Dleft / (1 / Dleft + 1 / Dright + 1 / Dup)*Phase5[i][j - Dleft]
						             + 1/Dright / (1 / Dleft + 1 / Dright + 1 / Dup)*Phase5[i][j + Dright]
						             + 1/Dup / (1 / Dleft + 1 / Dright + 1 / Dup)*Phase5[i + Dup][j];
					else 
						Phase5[i][j] = 1/Dleft / (1 / Dleft + 1 / Dright + 1 / Dup + 1 / Ddown)*Phase5[i][j - Dleft]
							         + 1/Dright / (1 / Dleft + 1 / Dright + 1 / Dup + 1 / Ddown)*Phase5[i][j + Dright]
							         + 1/Dup / (1 / Dleft + 1 / Dright + 1 / Dup + 1 / Ddown)*Phase5[i + Dup][j]
						             + 1/Ddown / (1 / Dleft + 1 / Dright + 1 / Dup + 1 / Ddown)*Phase5[i - Ddown][j];

					//�������ϲ�ֵ�󣬸õ�Ĳ������Ա�����
					QualityMap_Expanded[i][j] = 0;
				}

			}
		}

	//�ٲ�N0-1��Middle_Point+1��
	for (int i = N0-1; i > Middle_Point; i--)
		for (int j = 0; j < N0; j++)
		{
			{
				if (QualityMap_Expanded[i][j] == 1)
				{
					double Dleft = 1;//�������ߵĵ�QualityMap_Expanded�϶�����1
					double Dup = 1;//���������ĵ�QualityMap_Expanded�϶�����1�������������������Dup������
					double Dright = 1; double Ddown = 1; //�������ǳ�ʼ��

					for (int k = j + 1; k < N0; k++)
					{
						if (QualityMap_Expanded[i][k] == 0)
						{
							Dright = k - j;//�����Ҳ��������
							break;//�ҵ�������ѭ��
						}
					}

					for (int k = i - 1; k >= Middle_Point; k--)
					{
						if (QualityMap_Expanded[k][j] == 0)
						{
							Ddown = i - k;//�����Ҳ��������
							break;//�ҵ�������ѭ��
						}
					}

					//������������������棬��ʱphase5[i][j]���������������Ծ��뷴�Ȳ�ֵ�õ�
					if (i == N0 - 1)
						Phase5[i][j] = 1 / Dleft / (1 / Dleft + 1 / Dright + 1 / Ddown)*Phase5[i][j - Dleft]
						             + 1 / Dright / (1 / Dleft + 1 / Dright + 1 / Ddown)*Phase5[i][j + Dright]
						             + 1 / Ddown / (1 / Dleft + 1 / Dright + 1 / Ddown)*Phase5[i - Ddown][j];
					else//���������������ĸ����ֵ�õ�
						Phase5[i][j] = 1 / Dleft / (1 / Dleft + 1 / Dright + 1 / Dup + 1 / Ddown)*Phase5[i][j - Dleft]
						             + 1 / Dright / (1 / Dleft + 1 / Dright + 1 / Dup + 1 / Ddown)*Phase5[i][j + Dright]
						             + 1 / Dup / (1 / Dleft + 1 / Dright + 1 / Dup + 1 / Ddown)*Phase5[i + Dup][j]
						             + 1 / Ddown / (1 / Dleft + 1 / Dright + 1 / Dup + 1 / Ddown)*Phase5[i - Ddown][j];
					//�������ϲ�ֵ�󣬸õ�Ĳ������Ա�����
					QualityMap_Expanded[i][j] = 0;
				}

			}
		}

//���濴�������֮ǰ��û�б�ת�ȹ�����־��RotationIndex��0����û�б�ת�ȹ���1����ת�ȹ�
	if (RotationIndex == 0)
		return Phase5;
	else 
	{
		vector <vector<double>> Phase6(N0, (vector <double>(N0, 0)));
		for (int i = 0; i < N0; i++)
			for (int j = 0; j < N0; j++)
			{
				{
					Phase6[i][j] = Phase5[j][i];//������ƺ�ľ���ת�Ⱥ󷵻�
				}
			}
		return Phase6;
	}

}
