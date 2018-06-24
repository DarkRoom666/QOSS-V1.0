#include "..\include\WaveguideRadiator.h"
#include <fstream>
#include "../util/Constant_Var.h"

WaveguideRadiator::WaveguideRadiator()
{
	N_width = 0;
	M_depth = 0;
	fileAddress = "";
    isAmPhs = true;
	beginY = 0;
	beginX = 0;
}

WaveguideRadiator::~WaveguideRadiator()
{
}

void WaveguideRadiator::readData()
{
	ifstream file(fileAddress);
	if (!file.is_open())
		return;
	string temp;
	while (beginY--)
		getline(file, temp);
	for (int i = 0; i < N_width; i++)
	{
		for (int j = 0; j < M_depth; j++)
		{
			getline(file, temp);
			istringstream tempLine(temp);
			int tempBegin = beginX;
			double temp1;
			while (tempBegin--)
				tempLine >> temp1;
			double a1, p1, a2, p2;
			tempLine >> a1 >> p1 >> a2 >> p2;
			if (isAmPhs)
			{
				p1 += data[3];
				a1 *= data[2];
				p2 += data[3];
				a2 *= data[2];
				Ex_unit[i][j] = complex<double>(a1*cos(p1 / 180 * Pi),
					a1*sin(p1 / 180 * Pi));
				Ey_unit[i][j] = complex<double>(a2*cos(p2 / 180 * Pi),
					a2*sin(p2 / 180 * Pi));
			}
			else
			{
				Ex_unit[i][j] = data[2] * complex<double>(
					a1 * cos(data[3] / 180 * Pi), p1 * sin(data[3] / 180 * Pi));
				Ey_unit[i][j] = data[2] * complex<double>(
					a2 * cos(data[3] / 180 * Pi), p2 * sin(data[3] / 180 * Pi));
			}

		}
	}
	file.close();
}

void WaveguideRadiator::setBasePara(const std::string& _fileAddress, int _N_width,
	int _M_depth, int _beginY, int _beginX, double _ds,
	bool _isAmPhs)
{
	fileAddress = _fileAddress;
	M_depth = _M_depth;
	N_width = _N_width;
	beginY = _beginY;
	beginX = _beginX;
	isAmPhs = _isAmPhs;
	ds = _ds;
}

void WaveguideRadiator::setLayout(int _gapCol, int _gapRow, int _numNN)
{
	gapCol = _gapCol;
	gapRow = _gapRow;
	numNN = _numNN;
}

void WaveguideRadiator::setEachUnit(const std::vector<std::vector<double>>& _amVec,
	const std::vector<std::vector<double>>& _phsVec)
{
	int n = _amVec.size();
	amVec.resize(n);
	phsVec.resize(n);
	for (int i = 0; i < n; ++i)
	{
		amVec[i].resize(n);
		phsVec[i].resize(n);
	}

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			amVec[i][j] = _amVec[i][j];
			phsVec[i][j] = _phsVec[i][j];
		}
}

void WaveguideRadiator::genField(Field* field)
{
	vector<vector<complex<double>>> Ex, Ey;
	int N = N_width * numNN + gapRow * (numNN - 1);
	int M = M_depth * numNN + gapCol * (numNN - 1);

	Ex.resize(N);
	Ey.resize(N);
	for (int i = 0; i < N; ++i)
	{
		Ex[i].resize(M);
		Ey[i].resize(M);
	}

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < M; ++j)
		{
			Ex[i][j] = 0;
			Ey[i][j] = 0;
		}

	int indexI = 0;
	int indexJ = 0;
	for (int x = 0; x < N; ++x)
		for (int y = 0; y < M; ++y)
		{
			for (int i = 0; i < N; ++i)
				for (int j = 0; j < M; ++j)
				{
					indexI = i + indexI * (N_width + gapRow);
					indexJ = j + indexJ * (M_depth + gapCol);

					Ex[indexI][indexJ] = amVec[x][y] * complex<double>(
						Ex_unit[i][j].real() * cos(phsVec[x][y] / 180 * Pi),
						Ex_unit[i][j].imag() * sin(phsVec[x][y] / 180 * Pi));
					Ey[indexI][indexJ] = amVec[x][y] * complex<double>(
						Ey_unit[i][j].real() * cos(phsVec[x][y] / 180 * Pi),
						Ey_unit[i][j].imag() * sin(phsVec[x][y] / 180 * Pi));
				}
		}
	


	field->setNM(N, M);
	field->setDs(ds);
	field->setIsSource(true);
	field->setField(Ex, Ey);
}

void WaveguideRadiator::updateData()
{
}
