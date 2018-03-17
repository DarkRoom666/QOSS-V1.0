#include "ApertureRadiation.h"



ApertureRadiation::ApertureRadiation()
{

}

ApertureRadiation::~ApertureRadiation()
{
	Free_1D(Freqs);	Freqs = NULL;//��̬������ڴ棬�����
}

void ApertureRadiation::SetReturnFloat(void(*returnFloat)(float, void*), void *_user)
{
	this->returnFloat = returnFloat;
	this->user = _user;
}

void ApertureRadiation::SetLogfile(fstream* logAdd) {
	this->Logfile = logAdd;
}

void ApertureRadiation::SetUpHuygens(Position3D _SourceCenter, int _Nx, int _Ny, int _Nz, float _dx, float _dy, float _dz) {
	SourceCenter = _SourceCenter;
	Nx = _Nx;
	Ny = _Ny;
	Nz = _Nz;
	dx = _dx;
	dy = _dy;
	dz = _dz;
}

void ApertureRadiation::SetUpFrequency(int _NumFreq, float _CenterFreq, float _Bandwidth){
	NumFreq = _NumFreq;
	Freqs = Allocate_1D(Freqs, NumFreq);
	if (NumFreq > 1) {//��ֹһ��Ƶ��
		for (int i = 0; i < NumFreq; i++) {
			Freqs[i] = _Bandwidth*float(i)/(NumFreq - 1) + _CenterFreq - _Bandwidth*0.5;
		}
	}
	else { //����һ��Ƶ��
		Freqs[0] = _CenterFreq;
	}
}

//����
void ApertureRadiation::SetUpPropagatedAperture(Position3D _ApertureCenter, Vector3D _ApertureDirection, Vector3D _UDirection, Vector3D _VDirection, float _Lu, float _Lv, int _Nu, int _Nv) {
	ApertureCenter = _ApertureCenter;	//�����������	
	ApertureDirection = _ApertureDirection;		//�������ָ��
	UDirection = _UDirection;	//U direction
	VDirection = _VDirection;	//V Direction
	//Vector3D key;	//ֻ��Ϊ�˵���crossProduct
	//VDirection = key.crossProduct(ApertureDirection, UDirection);	//V direction
	Lu = _Lu;							//�������U����ߴ�			
	Lv = _Lv;							//�������V����ߴ�
	Nu = _Nu;							//�������U�����������
	Nv = _Nv;							//�������V�����������
	du = Lu / (Nu - 1);					//ע�⣺ NuΪ������
	dv = Lv / (Nv - 1);					//ע�⣺ NvΪ������
	//cout << "discrete interval of Aperture: " << du << " " << dv << endl;
} 

//���㴫����ĳ� - �ӻݸ�˹���� ��Ƶ
void ApertureRadiation::Propagation5FaceBox(vector<vector<complex<float>>>& _Eu, vector<vector<complex<float>>>& _Ev, complex<float>* _HuygensBoxData, int _OmpNum, int _Freqindex) {
	//���������ڴ���ⲿ���룡 ����ֻ�Ǽ���
	//����������ķ�ʽ���к�������
	if(returnFloat)
	{
		returnFloat(0.0, user);
	}


	int i, j, k;
	int ii, jj;
	//�ݸ�˹���λ����Ϣ
	float* px;	px = Allocate_1D(px, Nx);	for (i = 0; i < Nx; i++) { px[i] = (i + 0.5 - Nx / 2.0)*dx + SourceCenter.X(); }
	float* py;	py = Allocate_1D(py, Ny);	for (j = 0; j < Ny; j++) { py[j] = (j + 0.5 - Ny / 2.0)*dy + SourceCenter.Y(); }
	float* pz;	pz = Allocate_1D(pz, Nz);	for (k = 0; k < Nz; k++) { pz[k] = (k + 0.5 - Nz / 2.0)*dz + SourceCenter.Z(); }

	float x_yoz_0, x_yoz_1, y_xoz_0, y_xoz_1, z_xoy_1;//�������λ��
	x_yoz_0 = SourceCenter.X() - Nx*dx / 2.0;	//yoz-���λ�� X
	x_yoz_1 = SourceCenter.X() + Nx*dx / 2.0;	//yoz+���λ�� X
	y_xoz_0 = SourceCenter.Y() - Ny*dy / 2.0;	//xoz-���λ�� Y
	y_xoz_1 = SourceCenter.Y() + Ny*dy / 2.0;	//xoz+���λ�� Y
	z_xoy_1 = SourceCenter.Z() + Nz*dz / 2.0;	//xoy+���λ�� Z
	//�ݸ�˹������Ƶ����Ϣ
	float computeFreq;	computeFreq = Freqs[_Freqindex];

	float nx;	float ny;  float nz;
	int ShiftFreq;	int ShiftFace;
	int ShiftEx;	int ShiftEy;	int ShiftEz;
	int ShiftHx;	int ShiftHy;	int ShiftHz;
	
	int N1, N2;

	//Eigen::Aligned128;
	Eigen::setCpuCacheSizes(32768,262144,2621440);
	
	//Array������
	ArrayXcf Jx;	ArrayXcf Jy;	ArrayXcf Jz;	//�������	����ų�
	ArrayXcf Jmx;	ArrayXcf Jmy;	ArrayXcf Jmz;	//�������	����糡
	ArrayXf Px;		ArrayXf Py;		ArrayXf Pz;				//�洢Դ��λ��
	
	//omp������м�洢
	ArrayXXcf Eu_omp;	Eu_omp.resize(Nu*Nv, _OmpNum);
	ArrayXXcf Ev_omp;	Ev_omp.resize(Nu*Nv, _OmpNum);

	complex<float> zero(0.0, 0.0);
	float dsds;
	//���ݸ�˹��ĳߴ�
	int SizeXOY = Nx*Ny;
	int SizeXOZ = Nx*Nz;
	int SizeYOZ = Ny*Nz;
	//�ݸ�˹����ܲ�������
	int TotalSize = (SizeXOY + SizeXOZ * 2 + SizeYOZ * 2);

	//��������ߴ�	���Ǵ洢�������������������ֵ����λ��
	Jx.resize(TotalSize);	Jy.resize(TotalSize);	Jz.resize(TotalSize);
	Jmx.resize(TotalSize);	Jmy.resize(TotalSize);	Jmz.resize(TotalSize);
	Px.resize(TotalSize);	Py.resize(TotalSize);	Pz.resize(TotalSize);


	//���ñ������ ��ȡ�ݸ�˹����
	//��һ���ݸ�˹��	XOY_1	������׼����
	nx = 0.0;	ny = 0.0;	nz = 1.0;	//���Ǳ��淨������
	//Shift ����_HuygensBoxData�е�λ��
	ShiftFreq = _Freqindex*(SizeXOY + SizeXOZ * 2 + SizeYOZ * 2) * 4;	//
	ShiftFace = 0;
	ShiftEx = 0;
	ShiftEy = Nx*Ny * 1;
	ShiftEz = 0;	//Empty
	ShiftHx = Nx*Ny * 2;
	ShiftHy = Nx*Ny * 3;
	ShiftHz = 0;	//Empty
	dsds = dx*dy;
	//��HuygensBoxData�ж�ȡ��ֵ
	for (j = 0; j < Ny; j++) {
		for (i = 0; i < Nx; i++) {
			Jx[i + j*Nx] = -nz*_HuygensBoxData[i + j*Nx + ShiftHy + ShiftFace + ShiftFreq] * dsds;		//ny*Hz - nz*Hy		-nz*Hy
			Jy[i + j*Nx] = nz*_HuygensBoxData[i + j*Nx + ShiftHx + ShiftFace + ShiftFreq] * dsds;		//nz*Hx - nx*Hz		nz*Hx
			Jz[i + j*Nx] = zero;																		//nx*Hy	- ny*Hx		0
			Jmx[i + j*Nx] = nz*_HuygensBoxData[i + j*Nx + ShiftEy + ShiftFace + ShiftFreq] * dsds;		//nz*Ey - ny*Ez		nz*Ey
			Jmy[i + j*Nx] = -nz*_HuygensBoxData[i + j*Nx + ShiftEx + ShiftFace + ShiftFreq] * dsds;		//nx*Ez - nz*Ex		-nz*Ex
			Jmz[i + j*Nx] = zero;																		//ny*Ex - nx*Ey		0
		}
	} 
	//����λ��
	for (j = 0; j < Ny; j++) {
		for (i = 0; i < Nx; i++) {
			Px[i + j*Nx] = px[i];
			Py[i + j*Nx] = py[j];
			Pz[i + j*Nx] = z_xoy_1;
		}
	}
	//�ڶ����ݸ�˹��	XOZ_0	
	nx = 0.0;	ny = -1.0;	nz = 0.0;	//���淨��
	//HuygensBoxData�е�ƫ����
	ShiftFace = (SizeXOY)*4;
	ShiftEx = 0;
	ShiftEy = 0;//Empty
	ShiftEz = SizeXOZ;	/
	ShiftHx = SizeXOZ * 2;
	ShiftHy = 0;//Empty
	ShiftHz = SizeXOZ * 3;
	dsds = dx*dz;
	for (k = 0; k < Nz; k++) {
		for (i = 0; i < Nx; i++) {
			Jx[i + k*Nx + SizeXOY] = ny*_HuygensBoxData[i + k*Nx + ShiftHz + ShiftFace + ShiftFreq] * dsds;		//ny*Hz - nz*Hy		ny*Hz
			Jy[i + k*Nx + SizeXOY] = zero;																		//nz*Hx - nx*Hz		0
			Jz[i + k*Nx + SizeXOY] = -ny*_HuygensBoxData[i + k*Nx + ShiftHx + ShiftFace + ShiftFreq] * dsds;	//nx*Hy	- ny*Hx		-ny*Hx
			Jmx[i + k*Nx + SizeXOY] = -ny*_HuygensBoxData[i + k*Nx + ShiftEz + ShiftFace + ShiftFreq] * dsds;	//nz*Ey - ny*Ez		-ny*Ez
			Jmy[i + k*Nx + SizeXOY] = zero;																		//nx*Ez - nz*Ex		0
			Jmz[i + k*Nx + SizeXOY] = ny*_HuygensBoxData[i + k*Nx + ShiftEx + ShiftFace + ShiftFreq] * dsds;	//ny*Ex - nx*Ey		ny*Ex
		}
	}
	//����λ��
	for (k = 0; k < Nz; k++) {
		for (i = 0; i < Nx; i++) {
			Px[i + k*Nx + SizeXOY] = px[i];
			Py[i + k*Nx + SizeXOY] = y_xoz_0;
			Pz[i + k*Nx + SizeXOY] = pz[k];
		}
	}
	//�������ݸ�˹��	XOZ_1	
	nx = 0.0;	ny = 1;	nz = 0.0;	//���淨��
	//����HuygensBoxData�е�ƫ����
	ShiftFace = (SizeXOY + SizeXOZ) * 4;
	ShiftEx = 0;
	ShiftEy = 0;//Empty
	ShiftEz = SizeXOZ;
	ShiftHx = SizeXOZ * 2;
	ShiftHy = 0;//Empty
	ShiftHz = SizeXOZ * 3;
	dsds = dx*dz;
	//���ñ�������
	for (k = 0; k < Nz; k++) {
		for (i = 0; i < Nx; i++) {
			Jx[i + k*Nx + SizeXOY + SizeXOZ] = ny*_HuygensBoxData[i + k*Nx + ShiftHz + ShiftFace + ShiftFreq] * dsds;		//ny*Hz - nz*Hy		ny*Hz
			Jy[i + k*Nx + SizeXOY + SizeXOZ] = zero;																		//nz*Hx - nx*Hz		0
			Jz[i + k*Nx + SizeXOY + SizeXOZ] = -ny*_HuygensBoxData[i + k*Nx + ShiftHx + ShiftFace + ShiftFreq] * dsds;		//nx*Hy	- ny*Hx		-ny*Hx
			Jmx[i + k*Nx + SizeXOY + SizeXOZ] = -ny*_HuygensBoxData[i + k*Nx + ShiftEz + ShiftFace + ShiftFreq] * dsds;		//nz*Ey - ny*Ez		-ny*Ez
			Jmy[i + k*Nx + SizeXOY + SizeXOZ] = zero;																		//nx*Ez - nz*Ex		0
			Jmz[i + k*Nx + SizeXOY + SizeXOZ] = ny*_HuygensBoxData[i + k*Nx + ShiftEx + ShiftFace + ShiftFreq] * dsds;		//ny*Ex - nx*Ey		ny*Ex
		}
	}
	//���õ����λ��
	for (k = 0; k < Nz; k++) {
		for (i = 0; i < Nx; i++) {
			Px[i + k*Nx + SizeXOY + SizeXOZ] = px[i];
			Py[i + k*Nx + SizeXOY + SizeXOZ] = y_xoz_1;
			Pz[i + k*Nx + SizeXOY + SizeXOZ] = pz[k];
		}
	}

	//���ĸ��ݸ�˹��	YOZ_0	
	nx = -1.0;	ny = 0;	nz = 0.0;	//���淨��
	//����HuygensBoxData�е�ƫ��
	ShiftFace = (SizeXOY + SizeXOZ*2) * 4;
	ShiftEx = 0;//Empty
	ShiftEy = 0;
	ShiftEz = SizeYOZ;
	ShiftHx = 0;//Empty
	ShiftHy = SizeYOZ * 2;
	ShiftHz = SizeYOZ * 3;
	dsds = dy*dz;
	//���ñ�������
	for (k = 0; k < Nz; k++) {
		for (j = 0; j < Ny; j++) {
			Jx[j + k*Ny + SizeXOY + SizeXOZ*2] = zero;																		//ny*Hz - nz*Hy		0
			Jy[j + k*Ny + SizeXOY + SizeXOZ*2] = -nx*_HuygensBoxData[j + k*Ny + ShiftHz + ShiftFace + ShiftFreq] * dsds;	//nz*Hx - nx*Hz		-nx*Hz
			Jz[j + k*Ny + SizeXOY + SizeXOZ*2] = nx*_HuygensBoxData[j + k*Ny + ShiftHy + ShiftFace + ShiftFreq] * dsds;		//nx*Hy	- ny*Hx		nx*Hy
			Jmx[j + k*Ny + SizeXOY + SizeXOZ*2] = zero;																		//nz*Ey - ny*Ez		0
			Jmy[j + k*Ny + SizeXOY + SizeXOZ*2] = nx*_HuygensBoxData[j + k*Ny + ShiftEz + ShiftFace + ShiftFreq] * dsds;	//nx*Ez - nz*Ex		nx*Ez
			Jmz[j + k*Ny + SizeXOY + SizeXOZ*2] = -nx*_HuygensBoxData[j + k*Ny + ShiftEy + ShiftFace + ShiftFreq] * dsds;	//ny*Ex - nx*Ey		-nx*Ey
		}
	}
	//����λ��
	for (k = 0; k < Nz; k++) {
		for (j = 0;j < Ny; j++) {
			Px[j + k*Ny + SizeXOY + SizeXOZ*2] = x_yoz_0;
			Py[j + k*Ny + SizeXOY + SizeXOZ*2] = py[j];
			Pz[j + k*Ny + SizeXOY + SizeXOZ*2] = pz[k];
		}
	}

	//������ݸ�˹��	YOZ_1	
	nx = 1.0;	ny = 0;	nz = 0.0;	//���淨������
	//HuygensBoxData�е�ƫ��
	ShiftFace = (SizeXOY + SizeXOZ * 2 + SizeYOZ) * 4;
	ShiftEx = 0;//Empty
	ShiftEy = 0;
	ShiftEz = SizeYOZ;
	ShiftHx = 0;//Empty
	ShiftHy = SizeYOZ * 2;
	ShiftHz = SizeYOZ * 3;
	dsds = dy*dz;
	//���ñ�������
	for (k = 0; k < Nz; k++) {
		for (j = 0; j < Ny; j++) {
			Jx[j + k*Ny + SizeXOY + SizeXOZ*2 + SizeYOZ] = zero;																	//ny*Hz - nz*Hy		0
			Jy[j + k*Ny + SizeXOY + SizeXOZ*2 + SizeYOZ] = -nx*_HuygensBoxData[j + k*Ny + ShiftHz + ShiftFace + ShiftFreq] * dsds;	//nz*Hx - nx*Hz		-nx*Hz
			Jz[j + k*Ny + SizeXOY + SizeXOZ*2 + SizeYOZ] = nx*_HuygensBoxData[j + k*Ny + ShiftHy + ShiftFace + ShiftFreq] * dsds;	//nx*Hy	- ny*Hx		nx*Hy
			Jmx[j + k*Ny + SizeXOY + SizeXOZ*2 + SizeYOZ] = zero;																	//nz*Ey - ny*Ez		0
			Jmy[j + k*Ny + SizeXOY + SizeXOZ*2 + SizeYOZ] = nx*_HuygensBoxData[j + k*Ny + ShiftEz + ShiftFace + ShiftFreq] * dsds;	//nx*Ez - nz*Ex		nx*Ez
			Jmz[j + k*Ny + SizeXOY + SizeXOZ*2 + SizeYOZ] = -nx*_HuygensBoxData[j + k*Ny + ShiftEy + ShiftFace + ShiftFreq] * dsds;	//ny*Ex - nx*Ey		-nx*Ey
		}
	}
	//���õ����λ��
	for (k = 0; k < Nz; k++) {
		for (j = 0; j < Ny; j++) {
			Px[j + k*Ny + SizeXOY + SizeXOZ*2 + SizeYOZ] = x_yoz_1;
			Py[j + k*Ny + SizeXOY + SizeXOZ*2 + SizeYOZ] = py[j];
			Pz[j + k*Ny + SizeXOY + SizeXOZ*2 + SizeYOZ] = pz[k];
		}
	}

	//�����Դ�����ã�
	//���ж�Դ��OpenMp����

	//���ԣ�
	//TotalSize = SizeXOY;

	int* Size;	Size = Allocate_1D(Size, _OmpNum);
	int* Start;	Start = Allocate_1D(Start, _OmpNum);
	(*Logfile) << "Omp Division for Huygens Propagation: Total Size (number): " << TotalSize << endl;
	for (int id = 0; id < _OmpNum; id++) {
		Size[id] = TotalSize / _OmpNum;
		int rr = TotalSize % _OmpNum;
		if (id < rr && rr != 0) Size[id] += 1;
		if (id == 0) {
			Start[id] = 0;
		}
		else
		{
			Start[id] = Start[id-1] + Size[id-1];
		}
		(*Logfile) << "StartLoc: " << Start[id]<<", Size: "<<Size[id]<<", thread "<<id << endl;
	}


	/*
	Nyn_cal[i] = (Ny_cal) / threadNum;
	int rr = (Ny_cal) % threadNum;
	if (i<rr && rr != 0)  Nyn_cal[i] += 1;
	if (i == 0) {
		Nys_cal[i] = 1;
	}
	else {
		Nys_cal[i] = Nys_cal[i - 1] + Nyn_cal[i - 1];
	}
	*/
	//��ʼ����
	omp_set_num_threads(_OmpNum);
	#pragma omp parallel
	{
		//float dsds; //���СƬ
		float Tx, Ty, Tz;
		float ww = 2 * PIf * computeFreq;
		float kk = 2 * PIf * computeFreq / C_Speedf;	//wavenumber
		complex<float> unit(1.0, 0.0);
		complex<float> uniti(0.0, 1.0);
		complex<float> TEx(0.0, 0.0);
		complex<float> TEy(0.0, 0.0);
		complex<float> TEz(0.0, 0.0);
		complex<float> item1(0.0, 0.0); complex<float> item2(0.0, 0.0); complex<float> item3(0.0, 0.0);	//temp

		int id = omp_get_thread_num();
		if (id == 0)
		{
			//cout << "Start OMP Division for Huygens Propagation." << endl;
			//cout << "Start OMP Division for Huygens Propagation @ Freq(GHz): " << computeFreq / 1.0e9 << endl;
			(*Logfile) << endl;
			(*Logfile) << "Start OMP Division for Huygens Propagation @ Freq(GHz): " << computeFreq / 1.0e9 << endl;

		}
		//Array������
		ArrayXcf Jx_omp;	ArrayXcf Jy_omp;	ArrayXcf Jz_omp;	//�������	����ų�
		ArrayXcf Jmx_omp;	ArrayXcf Jmy_omp;	ArrayXcf Jmz_omp;	//�������	����糡
		ArrayXf Px_omp;		ArrayXf Py_omp;		ArrayXf Pz_omp;		//�洢Դ��λ��

		ArrayXf Rx_omp;		ArrayXf Ry_omp;		ArrayXf Rz_omp;		ArrayXf R_omp;	//Դ���۲���λ��ʸ��
		ArrayXcf Coe1_omp;	ArrayXcf Coe2_omp;	ArrayXcf Coe3_omp;	//�洢ϵ��
		ArrayXcf Valx_omp;	ArrayXcf Valy_omp;	ArrayXcf Valz_omp;	//�洢3��ļ�����
		ArrayXcf Temp1_omp;	ArrayXcf Temp2_omp;	ArrayXcf Temp3_omp;	//�洢�м����
		ArrayXcf expR_omp;	//�洢ָ����

		//��������ߴ�
		Jx_omp.resize(Size[id]);	Jy_omp.resize(Size[id]);	Jz_omp.resize(Size[id]);
		Jmx_omp.resize(Size[id]);	Jmy_omp.resize(Size[id]);	Jmz_omp.resize(Size[id]);
		Rx_omp.resize(Size[id]);	Ry_omp.resize(Size[id]);	Rz_omp.resize(Size[id]);
		Px_omp.resize(Size[id]);	Py_omp.resize(Size[id]);	Pz_omp.resize(Size[id]);
		Coe1_omp.resize(Size[id]);	Coe2_omp.resize(Size[id]);	Coe3_omp.resize(Size[id]);
		Valx_omp.resize(Size[id]);	Valy_omp.resize(Size[id]);	Valz_omp.resize(Size[id]);
		Temp1_omp.resize(Size[id]);	Temp2_omp.resize(Size[id]);	Temp3_omp.resize(Size[id]);
		expR_omp.resize(Size[id]);	
		
		//LoadFromSource
		for (int i = 0; i < Size[id]; i++) {
			Jx_omp[i] = Jx[i + Start[id]];		Jy_omp[i] = Jy[i + Start[id]];		Jz_omp[i] = Jz[i + Start[id]];
			Jmx_omp[i] = Jmx[i + Start[id]];	Jmy_omp[i] = Jmy[i + Start[id]];	Jmz_omp[i] = Jmz[i + Start[id]];
			Px_omp[i] = Px[i + Start[id]];		Py_omp[i] = Py[i + Start[id]];		Pz_omp[i] = Pz[i + Start[id]];
		}
		//Start Omp Running
		//��������
		for (int jj = 0; jj < Nv; jj++) {
			for (int ii = 0; ii < Nu; ii++) {
				//�۲���λ��
				Tx = ApertureCenter.X() + (jj - (Nv-1) / 2.0)*dv*VDirection.X() + (ii - (Nu-1) / 2.0)*du*UDirection.X();
				Ty = ApertureCenter.Y() + (jj - (Nv-1) / 2.0)*dv*VDirection.Y() + (ii - (Nu-1) / 2.0)*du*UDirection.Y();
				Tz = ApertureCenter.Z() + (jj - (Nv-1) / 2.0)*dv*VDirection.Z() + (ii - (Nu-1) / 2.0)*du*UDirection.Z();
				//if(id == 0) cout << Tx << " " << Ty << " " << Tz << "    ";
				//Set R and ExpR
				//ȫ�� �����ĸ���ֱ�ӼӼ��˳��������ǳ����

				Rx_omp = Tx - Px_omp;	
				Ry_omp = Ty - Py_omp;	
				Rz_omp = Tz - Pz_omp;
				R_omp = Rx_omp.square() + Ry_omp.square() + Rz_omp.square();
				R_omp = R_omp.sqrt();

				expR_omp = (-kk*uniti)*R_omp;		
				expR_omp = expR_omp.exp();

				item1 = unit / (4 * PIf * ww * Eps0f * uniti); // 1 / (4*PI*j*w*eps)
				item1 = item1;
				Coe1_omp = (float(3.0) - (kk*kk)*(R_omp*R_omp) + (uniti*float(3.0*kk)*R_omp)) / (R_omp*R_omp*R_omp*R_omp*R_omp);
				Coe1_omp = item1 * Coe1_omp;	Coe1_omp = Coe1_omp * expR_omp;

				Coe2_omp = float(2.0) * (float(1.0) + (uniti*kk)*R_omp) / (R_omp*R_omp*R_omp);
				Coe2_omp = item1 * Coe2_omp;	Coe2_omp = Coe2_omp * expR_omp;

				item2 = float(1.0) / (float(4.0)*PIf);
				item2 = item2;

				Coe3_omp = (float(-1.0) - (uniti*kk)*R_omp) / (R_omp*R_omp*R_omp);
				Coe3_omp = item2 * Coe3_omp;	Coe3_omp = Coe3_omp * expR_omp;

				//Ex
				Valx_omp = Coe1_omp * (Ry_omp*Rx_omp*Jy_omp - Ry_omp*Ry_omp*Jx_omp - Rz_omp*Rz_omp*Jx_omp + Rz_omp*Rx_omp*Jz_omp)
					+ Coe2_omp * (Jx_omp)
					+Coe3_omp * (Rz_omp*Jmy_omp - Ry_omp*Jmz_omp);
				//Ey 
				Valy_omp = Coe1_omp * (Rz_omp*Ry_omp*Jz_omp - Rz_omp*Rz_omp*Jy_omp - Rx_omp*Rx_omp*Jy_omp + Rx_omp*Ry_omp*Jx_omp)
					+ Coe2_omp * (Jy_omp)
					+Coe3_omp * (Rx_omp*Jmz_omp - Rz_omp*Jmx_omp);
				//Ez 
				Valz_omp = Coe1_omp * (Rx_omp*Rz_omp*Jx_omp - Rx_omp*Rx_omp*Jz_omp - Ry_omp*Ry_omp*Jz_omp + Ry_omp*Rz_omp*Jy_omp)
					+ Coe2_omp * (Jz_omp)
					+Coe3_omp * (Ry_omp*Jmx_omp - Rx_omp*Jmy_omp);

				TEx = Valx_omp.sum();
				TEy = Valy_omp.sum();
				TEz = Valz_omp.sum();

				Eu_omp(ii+jj*Nu,id) = TEx * UDirection.X() + TEy * UDirection.Y() + TEz * UDirection.Z();
				//E dot u = Eu
				Ev_omp(ii+jj*Nu,id) = TEx * VDirection.X() + TEy * VDirection.Y() + TEz * VDirection.Z();
				//E dot v = Ev
			}//ii
			if (id == 0) {
				//cout << ".";
				(*Logfile) << ".";
				if (returnFloat)
				{
					returnFloat(float(jj*1.0 / Nv*100.0), user);
				}			
			}
		}//jj
		//������ɣ�����
		Jx_omp.resize(0);		Jy_omp.resize(0);		Jz_omp.resize(0);
		Jmx_omp.resize(0);		Jmy_omp.resize(0);		Jmz_omp.resize(0);
		Rx_omp.resize(0);		Ry_omp.resize(0);		Rz_omp.resize(0);
		Px_omp.resize(0);		Py_omp.resize(0);		Pz_omp.resize(0);
		Coe1_omp.resize(0);		Coe2_omp.resize(0);		Coe3_omp.resize(0);
		Valx_omp.resize(0);		Valy_omp.resize(0);		Valz_omp.resize(0);
		Temp1_omp.resize(0);	Temp2_omp.resize(0);	Temp3_omp.resize(0);
		expR_omp.resize(0);		
	}//omp
	//cout << "Huygens Propagation Done!" << endl;
	(*Logfile) << "Huygens Propagation Done!" << endl;
	for (int v = 0; v < Nv; v++) {
		for (int u = 0; u < Nu; u++) {
			for (int id = 0; id < _OmpNum; id++) {
				_Eu[u][v] = _Eu[u][v] + Eu_omp(u + v*Nu, id);
				_Ev[u][v] = _Ev[u][v] + Ev_omp(u + v*Nu, id);
			}
		}
	}

	//��д�ɾ����������ʽ

	//������� �����̬���ɵ��ڴ�
	Free_1D(px);	px = NULL;	Free_1D(py);	py = NULL;	Free_1D(pz);	pz = NULL;
	//Eigen Part
	Jx.resize(0);		Jy.resize(0);		Jz.resize(0);
	Jmx.resize(0);		Jmy.resize(0);		Jmz.resize(0);
}
