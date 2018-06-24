/*
*	created by liyun 2018/6/24
*   function 波导辐射器
*   version 1.0
*/
#ifndef WAVEGUIDERADIATOR_H
#define WAVEGUIDERADIATOR_H

#include "BasicParameters.h"
#include <complex>
#include <string>
#include "Field.h"

class WaveguideRadiator : public BasicParameters
{
public:
	WaveguideRadiator();
	~WaveguideRadiator();

	void readData();

	void setBasePara(const std::string& _fileAddress, int _N_width,
		int _M_depth, int _beginY, int _beginX, double _ds,
		bool _isAmPhs = true);

	void setLayout(int _gapCol, int _gapRow, int _numNN);

	void setEachUnit(const std::vector<std::vector<double>>& _amVec,
		const std::vector<std::vector<double>>& _phsVec);

	// 根据输入生产N*N的场
	void genField(Field* field);

	void updateData();

private:

	int N_width, M_depth;
	std::string fileAddress;
	double ds;
	bool isAmPhs;
	int beginY;
	int beginX;

	int gapCol;
	int gapRow;
	int numNN;  // N * N

	std::vector<std::vector<double>> amVec;
	std::vector<std::vector<double>> phsVec;

	vector<vector<complex<double>>> Ex_unit, Ey_unit;


};

#endif //WAVEGUIDERADIATOR_H