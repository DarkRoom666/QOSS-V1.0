#include "HuygensBoxData.h"

HuygensBox::HuygensBox() {
	Nx = 1;
}

HuygensBox::~HuygensBox() {

}

int HuygensBox::SetHuygensParas(int _Nx, int _Ny, int _Nz, float _dx, float _dy, float _dz, float _cx, float _cy, float _cz, float _frequency) {
	Nx = _Nx;
	Ny = _Ny;
	Nz = _Nz;
	Size = (Nx*Ny + Nx*Nz * 2 + Ny*Nz * 2)*4;
	dx = _dx;
	dy = _dy;
	dz = _dz;
	cx = _cx;
	cy = _cy;
	cz = _cz;
	frequency = _frequency;
	return 0;
}

int HuygensBox::SetHuygensData(complex<float>* _Huygensin, int _index) {
	int shift = Size*_index;
	
	HuygensData.resize(Size);
	for (int i = 0; i < Size; i++) {
		HuygensData[i] = _Huygensin[i + shift];
	}
	return 0;
}

int HuygensBox::ReadHuygensFile(const std::string & file, int _index) {
	int Nfreq;
	float tempf;
	float tempr;
	float tempi;
	complex<float> tempc;
	
	FILE* readfile;
	readfile = fopen(file.c_str(), "rb");
	//大小
	fread(&tempi, sizeof(int), 1, readfile);	Nx = tempi;
	fread(&tempi, sizeof(int), 1, readfile);	Ny = tempi;
	fread(&tempi, sizeof(int), 1, readfile);	Nz = tempi;
	Size = (Nx*Ny + Nx*Nz * 2 + Ny*Nz * 2) * 4;
	//离散间隔
	fread(&tempf, sizeof(float), 1, readfile);	dx = tempf;
	fread(&tempf, sizeof(float), 1, readfile);	dy = tempf;
	fread(&tempf, sizeof(float), 1, readfile);	dz = tempf;
	//中间位置
	fread(&tempf, sizeof(float), 1, readfile);	cx = tempf;
	fread(&tempf, sizeof(float), 1, readfile);	cy = tempf;
	fread(&tempf, sizeof(float), 1, readfile);	cz = tempf;
	//频率先读个数，再读场值
	fread(&tempi, sizeof(int), 1, readfile);	Nfreq = tempi;
	if (Nfreq - 1 < _index) return -1;	//要是序列号超出范围就报错返回

	for (int i = 0; i < Nfreq; i++) {
		fread(&tempf, sizeof(float), 1, readfile);
		if (i == _index) frequency = tempf;
	}
	
	HuygensData.resize(Size);
	//读取数据
	for (int i = 0; i < Nfreq; i++) {
		if (i == _index) {//要读的
			for (int ii = 0; ii < Size; ii++) {
				fread(&tempc, sizeof(complex<float>), 1, readfile);
				HuygensData[ii] = tempc;
			}
		}
		else {
			for (int ii = 0; ii < Size; ii++) {
				fread(&tempc, sizeof(complex<float>), 1, readfile);
			}
		}
	}
	return 0;
}

int HuygensBox::GetCurrentList(vector<complex<float>> &_Jx, vector<complex<float>> &_Jy, vector<complex<float>> &_Jz, vector<complex<float>> &_Jmx, vector<complex<float>> &_Jmy, vector<complex<float>> &_Jmz) {
	_Jx.resize(Size/4);
	_Jy.resize(Size/4);
	_Jz.resize(Size/4);
	_Jmx.resize(Size/4);
	_Jmy.resize(Size/4);
	_Jmz.resize(Size/4);

	int ShiftFace;
	int ShiftEx;	int ShiftEy;	int ShiftEz;
	int ShiftHx;	int ShiftHy;	int ShiftHz;
	int SizeXOY = Nx*Ny;
	int SizeXOZ = Nx*Nz;
	int SizeYOZ = Ny*Nz;

	//第一个惠更斯面	XOY_1	建立标准流程
	float nx = 0.0;	 float ny = 0.0;	float nz = 1.0;	//这是表面法向向量
	//PowerSum = complex<float>(0.0, 0.0);	//这是功率积分
	//Shift 是在HuygensData中的位移
	complex<float> zero = complex<float>(0.0, 0.0);
	ShiftFace = 0;
	ShiftEx = 0;
	ShiftEy = SizeXOY;
	ShiftEz = 0;	//Empty
	ShiftHx = SizeXOY * 2;
	ShiftHy = SizeXOY * 3;
	ShiftHz = 0;	//Empty
	//从HuygensBoxData中读取场值
	for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			_Jx[i + j*Nx] = -nz*HuygensData[i + j*Nx + ShiftHy + ShiftFace];	//ny*Hz - nz*Hy		-nz*Hy
			_Jy[i + j*Nx] = nz*HuygensData[i + j*Nx + ShiftHx + ShiftFace];			//nz*Hx - nx*Hz		nz*Hx
			_Jz[i + j*Nx] = zero;																								//nx*Hy	- ny*Hx		0
			_Jmx[i + j*Nx] = nz*HuygensData[i + j*Nx + ShiftEy + ShiftFace];		//nz*Ey - ny*Ez		nz*Ey
			_Jmy[i + j*Nx] = -nz*HuygensData[i + j*Nx + ShiftEx + ShiftFace];		//nx*Ez - nz*Ex		-nz*Ex
			_Jmz[i + j*Nx] = zero;																							//ny*Ex - nx*Ey		0
																																				//功率积分
		}
	}

	//第二个惠更斯面	XOZ_0	
	nx = 0.0;	ny = -1.0;	nz = 0.0;	//表面法向
	//HuygensBoxData中的偏移量
	ShiftFace = (SizeXOY) * 4;
	ShiftEx = 0;
	ShiftEy = 0;//Empty
	ShiftEz = SizeXOZ;
	ShiftHx = SizeXOZ * 2;
	ShiftHy = 0;//Empty
	ShiftHz = SizeXOZ * 3;
	for (int k = 0; k < Nz; k++) {
		for (int i = 0; i < Nx; i++) {
			_Jx[i + k*Nx + SizeXOY] = ny*HuygensData[i + k*Nx + ShiftHz + ShiftFace];	//ny*Hz - nz*Hy		ny*Hz
			_Jy[i + k*Nx + SizeXOY] = zero;												//nz*Hx - nx*Hz		0
			_Jz[i + k*Nx + SizeXOY] = -ny*HuygensData[i + k*Nx + ShiftHx + ShiftFace];	//nx*Hy	- ny*Hx		-ny*Hx
			_Jmx[i + k*Nx + SizeXOY] = -ny*HuygensData[i + k*Nx + ShiftEz + ShiftFace];	//nz*Ey - ny*Ez		-ny*Ez
			_Jmy[i + k*Nx + SizeXOY] = zero;																		//nx*Ez - nz*Ex		0
			_Jmz[i + k*Nx + SizeXOY] = ny*HuygensData[i + k*Nx + ShiftEx + ShiftFace];	//ny*Ex - nx*Ey		ny*Ex

		}
	}

	//第三个惠更斯面	XOZ_1	
	nx = 0.0;	ny = 1;	nz = 0.0;	//表面法向
	//设置HuygensBoxData中的偏移量
	ShiftFace = (SizeXOY + SizeXOZ) * 4;
	ShiftEx = 0;
	ShiftEy = 0;//Empty
	ShiftEz = SizeXOZ;
	ShiftHx = SizeXOZ * 2;
	ShiftHy = 0;//Empty
	ShiftHz = SizeXOZ * 3;

	//设置表面电磁流
	for (int k = 0; k < Nz; k++) {
		for (int i = 0; i < Nx; i++) {
			_Jx[i + k*Nx + SizeXOY + SizeXOZ] = ny*HuygensData[i + k*Nx + ShiftHz + ShiftFace];		//ny*Hz - nz*Hy		ny*Hz
			_Jy[i + k*Nx + SizeXOY + SizeXOZ] = zero;																		//nz*Hx - nx*Hz		0
			_Jz[i + k*Nx + SizeXOY + SizeXOZ] = -ny*HuygensData[i + k*Nx + ShiftHx + ShiftFace];		//nx*Hy	- ny*Hx		-ny*Hx
			_Jmx[i + k*Nx + SizeXOY + SizeXOZ] = -ny*HuygensData[i + k*Nx + ShiftEz + ShiftFace];		//nz*Ey - ny*Ez		-ny*Ez
			_Jmy[i + k*Nx + SizeXOY + SizeXOZ] = zero;																		//nx*Ez - nz*Ex		0
			_Jmz[i + k*Nx + SizeXOY + SizeXOZ] = ny*HuygensData[i + k*Nx + ShiftEx + ShiftFace];		//ny*Ex - nx*Ey		ny*Ex

		}
	}

	//第四个惠更斯面	YOZ_0	
	nx = -1.0;	ny = 0;	nz = 0.0;	//表面法向
	//设置HuygensBoxData中的偏移
	ShiftFace = (SizeXOY + SizeXOZ * 2) * 4;
	ShiftEx = 0;//Empty
	ShiftEy = 0;
	ShiftEz = SizeYOZ;
	ShiftHx = 0;//Empty
	ShiftHy = SizeYOZ * 2;
	ShiftHz = SizeYOZ * 3;
	//设置表面电磁流
	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
			_Jx[j + k*Ny + SizeXOY + SizeXOZ * 2] = zero;																		//ny*Hz - nz*Hy		0
			_Jy[j + k*Ny + SizeXOY + SizeXOZ * 2] = -nx*HuygensData[j + k*Ny + ShiftHz + ShiftFace];	//nz*Hx - nx*Hz		-nx*Hz
			_Jz[j + k*Ny + SizeXOY + SizeXOZ * 2] = nx*HuygensData[j + k*Ny + ShiftHy + ShiftFace];		//nx*Hy	- ny*Hx		nx*Hy
			_Jmx[j + k*Ny + SizeXOY + SizeXOZ * 2] = zero;																		//nz*Ey - ny*Ez		0
			_Jmy[j + k*Ny + SizeXOY + SizeXOZ * 2] = nx*HuygensData[j + k*Ny + ShiftEz + ShiftFace];	//nx*Ez - nz*Ex		nx*Ez
			_Jmz[j + k*Ny + SizeXOY + SizeXOZ * 2] = -nx*HuygensData[j + k*Ny + ShiftEy + ShiftFace];	//ny*Ex - nx*Ey		-nx*Ey

		}
	}

	//第五个惠更斯面	YOZ_1	
	nx = 1.0;	ny = 0;	nz = 0.0;	//表面法向向量
	//HuygensBoxData中的偏移
	ShiftFace = (SizeXOY + SizeXOZ * 2 + SizeYOZ) * 4;
	ShiftEx = 0;//Empty
	ShiftEy = 0;
	ShiftEz = SizeYOZ;
	ShiftHx = 0;//Empty
	ShiftHy = SizeYOZ * 2;
	ShiftHz = SizeYOZ * 3;
	//设置表面电磁流
	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
			_Jx[j + k*Ny + SizeXOY + SizeXOZ * 2 + SizeYOZ] = zero;												//ny*Hz - nz*Hy		0
			_Jy[j + k*Ny + SizeXOY + SizeXOZ * 2 + SizeYOZ] = -nx*HuygensData[j + k*Ny + ShiftHz + ShiftFace];	//nz*Hx - nx*Hz		-nx*Hz
			_Jz[j + k*Ny + SizeXOY + SizeXOZ * 2 + SizeYOZ] = nx*HuygensData[j + k*Ny + ShiftHy + ShiftFace];	//nx*Hy	- ny*Hx		nx*Hy
			_Jmx[j + k*Ny + SizeXOY + SizeXOZ * 2 + SizeYOZ] = zero;											//nz*Ey - ny*Ez		0
			_Jmy[j + k*Ny + SizeXOY + SizeXOZ * 2 + SizeYOZ] = nx*HuygensData[j + k*Ny + ShiftEz + ShiftFace];	//nx*Ez - nz*Ex		nx*Ez
			_Jmz[j + k*Ny + SizeXOY + SizeXOZ * 2 + SizeYOZ] = -nx*HuygensData[j + k*Ny + ShiftEy + ShiftFace];	//ny*Ex - nx*Ey		-nx*Ey

		}
	}

	return 0;
}

int HuygensBox::GetPosList(vector<float> & _posx, vector<float> & _posy, vector<float> & _posz, vector<float> &_dss) {
	//惠更斯面的位置信息
	vector<float> px;	px.resize(Nx);
	for (int i = 0; i < Nx; i++) { px[i] = (i + 0.5 - Nx / 2.0)*dx + cx; }

	vector<float> py;	py.resize(Ny);
	for (int j = 0; j < Ny; j++) { py[j] = (j + 0.5 - Ny / 2.0)*dy + cy; }

	vector<float> pz;	pz.resize(Nz);
	for (int k = 0; k < Nz; k++) { pz[k] = (k + 0.5 - Nz / 2.0)*dz + cz; }

	_posx.resize(Size/4);
	_posy.resize(Size/4);
	_posz.resize(Size/4);
	_dss.resize(Size/4);
	float xx, yy, zz;
	//各个面的位置
	float z_xoy1 = cz + Nz*dz / 2.0;
	float y_xoz0 = cy - Ny*dy / 2.0;
	float y_xoz1 = cy + Ny*dy / 2.0;
	float x_yoz0 = cx - Nx*dx / 2.0;
	float x_yoz1 = cx + Nx*dx / 2.0;
	int Shift;	int index;
	//XOY
	Shift = 0; zz = z_xoy1;
	for (int j = 0; j < Ny; j++) {
		yy = py[j];
		for (int i = 0; i < Nx; i++) {
			xx = px[i];
			index = i + j*Nx + Shift;
			_posx[index] = xx;
			_posy[index] = yy;
			_posz[index] = zz;
			_dss[index] = dx*dy;
		}
	}
	//XOZ0
	Shift = Nx*Ny; yy = y_xoz0;
	for (int k = 0; k < Nz; k++) {
		zz = pz[k];
		for (int i = 0; i < Nx; i++) {
			xx = px[i];
			index = i + k*Nx + Shift;
			_posx[index] = xx;
			_posy[index] = yy;
			_posz[index] = zz;
			_dss[index] = dx*dz;
		}
	}
	//XOZ1
	Shift = Nx*Ny + Nx*Nz; yy = y_xoz1;
	for (int k = 0; k < Nz; k++) {
		zz = pz[k];
		for (int i = 0; i < Nx; i++) {
			xx = px[i];
			index = i + k*Nx + Shift;
			_posx[index] = xx;
			_posy[index] = yy;
			_posz[index] = zz;
			_dss[index] = dx*dz;
		}
	}
	//YOZ0 
	Shift = Nx*Ny + Nx*Nz * 2;	xx = x_yoz0;
	for (int k = 0; k < Nz; k++) {
		zz = pz[k];
		for (int j = 0; j < Ny; j++) {
			yy = py[j];
			index = j + k*Ny + Shift;
			_posx[index] = xx;
			_posy[index] = yy;
			_posz[index] = zz;
			_dss[index] = dy*dz;
		}
	}
	//YOZ1 
	Shift = Nx*Ny + Nx*Nz * 2 + Ny*Nz;	xx = x_yoz1;
	for (int k = 0; k < Nz; k++) {
		zz = pz[k];
		for (int j = 0; j < Ny; j++) {
			yy = py[j];
			index = j + k*Ny + Shift;
			_posx[index] = xx;
			_posy[index] = yy;
			_posz[index] = zz;
			_dss[index] = dy*dz;
		}
	}
	return 0;
}

float HuygensBox::GetFreq(void) {
	float result = frequency;
	return result;
}

float HuygensBox::CountPower(void) {
	float result;
	return result;
}