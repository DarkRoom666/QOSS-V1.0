#pragma once
#ifndef HUYGENSBOXDATA_H
#define HUYGENSBOXDATA_H

//five-face HuygensBoxData
#include <complex>
#include <vector>
#include <string>
using namespace std;
class HuygensBox {
public:
	HuygensBox();
	~HuygensBox();


	int SetHuygensParas(int _Nx, int _Ny, int _Nz, float _dx, float _dy, float _dz, float _cx, float _cy, float _cz, float _frequency);
	int SetHuygensData(complex<float>* _Huygensin, int _index);
	int ReadHuygensFile(const std::string & file, int _index);
	int GetCurrentList(vector<complex<float>> &_Jx, vector<complex<float>> &_Jy, vector<complex<float>> &_Jz, vector<complex<float>> &_Jmx, vector<complex<float>> &_Jmy, vector<complex<float>> &_Jmz);
	int GetPosList(vector<float> & _posx, vector<float> & _posy, vector<float> & _posz, vector<float> &_dss);
	float GetFreq(void);
	float CountPower(void);
private:
		int Nx;
		int Ny;
		int Nz; 
		float dx;
		float dy;
		float dz;
		float cx;
		float cy;
		float cz;
		float frequency;
		int Size;
		
		vector<complex<float>> HuygensData;

};
#endif