#include "FDTDRadiator.h"
#include "Radiator.h"	//Setup Radiator Model
#include "FDTD.h"		//Setup FDTD Computation
#include "ApertureRadiation.h"	//SetupHuygensRadiation
#include "SourceModeGenerationT.h"	//Setup Source Model Aperture
#include "Global_Vars.h"
#include "Constant_Val.h"

using namespace calculation;
#define FDTDRADIATOR  _declspec(dllexport)
FDTDRadiator::FDTDRadiator() {
	FDTDRadiator::returnInt = NULL;
	FDTDRadiator::returnFloat = NULL;
	FDTDRadiator::user = NULL;
	FDTDRadiator::FILENAME = NULL;//20180320

	Ex_Port = NULL;
	Ey_Port = NULL;
	Hx_Port = NULL;
	Hy_Port = NULL;

	Eps = NULL;
}
FDTDRadiator::~FDTDRadiator() {
	Free_2D(Ex_Port);	Ex_Port = NULL;
	Free_2D(Ey_Port);	Ey_Port = NULL;
	Free_2D(Hx_Port);	Hx_Port = NULL;
	Free_2D(Hy_Port);	Hy_Port = NULL;

	Free_3D(Eps, Nz_model); Eps = NULL;
}

void FDTDRadiator::SetInt(int input) 
{
	FDTDRadiator::number = input;
}

void FDTDRadiator::SetpFun(void(*pFun)(int)) 
{
	this->pFun = pFun;
}
void FDTDRadiator::SetReturnInt(void(*returnInt)(int, void*), void *_user)
{
	this->returnInt = returnInt;
	this->user = _user;
}
void FDTDRadiator::SetReturnFloat(void(*returnFloat)(float, void*), void *_user)
{
	this->returnFloat = returnFloat;
	this->user = _user;
}
//获取电场分布
void FDTDRadiator::GetProFieldE(vector<vector<complex<double>>>& _Eu, vector<vector<complex<double>>>& _Ev, int _Nu, int _Nv)
{
	for (int u = 0; u < _Nu; u++) {
		for (int v = 0; v < _Nv; v++) {
			complex<double> temp;
			temp = complex<double>(double(FDTDRadiator::Eu[u][v].real()), double(FDTDRadiator::Eu[u][v].imag()));
			_Eu[u][v] = temp;
			temp = complex<double>(double(FDTDRadiator::Ev[u][v].real()), double(FDTDRadiator::Ev[u][v].imag()));
			_Ev[u][v] = temp;
		}
	}
}
//获取磁场分布
void FDTDRadiator::GetProFieldH(vector<vector<complex<double>>>& _Hu, vector<vector<complex<double>>>& _Hv, int _Nu, int _Nv)
{
	for (int u = 0; u < _Nu; u++) {
		for (int v = 0; v < _Nv; v++) {
			complex<double> temp;
			temp = complex<double>(double(FDTDRadiator::Hu[u][v].real()), double(FDTDRadiator::Hu[u][v].imag()));
			_Hu[u][v] = temp;
			temp = complex<double>(double(FDTDRadiator::Hv[u][v].real()), double(FDTDRadiator::Hv[u][v].imag()));
			_Hv[u][v] = temp;
		}
	}
}

void FDTDRadiator::GetProPowerRatio(double& _PowerRatio) {
	_PowerRatio = FDTDRadiator::PowerRatio;
}

void FDTDRadiator::LoadProFieldE(const char* _filename,vector<vector<complex<double>>>& _Eu, vector<vector<complex<double>>>& _Ev, int _Nu, int _Nv)
{
	FILE* Binread;
	Binread = fopen(_filename, "rb");
	int tempNu, tempNv;
	fread(&tempNu, sizeof(int), 1, Binread);	//Nu
	fread(&tempNv, sizeof(int), 1, Binread);	//Nv
	float temp;
	fread(&temp, sizeof(float), 1, Binread);	//du
	fread(&temp, sizeof(float), 1, Binread);	//dv
	fread(&temp, sizeof(float), 1, Binread);	//freq
	//fread(&temp, sizeof(float), 1, Binread);	//PowerRatio
	//FDTDRadiator::PowerRatio = double(temp);
	//注意 文件里是complex<float>的，目标位置是complex<double>的
	complex<float> readbuf;
	for (int v = 0; v<_Nv; v++) {
		for (int u = 0; u<_Nu; u++) {
			fread(&readbuf, sizeof(complex<float>), 1, Binread);
			_Eu[u][v] = complex<double>(double(readbuf.real()), double(readbuf.imag()));

		}
	}
	for (int v = 0; v<_Nv; v++) {
		for (int u = 0; u<_Nu; u++) {
			fread(&readbuf, sizeof(complex<float>), 1, Binread);
			_Ev[u][v] = complex<double>(double(readbuf.real()), double(readbuf.imag()));
		}
	}
	fclose(Binread);
}

void FDTDRadiator::LoadProFieldH(const char* _filename, vector<vector<complex<double>>>& _Hu, vector<vector<complex<double>>>& _Hv, int _Nu, int _Nv)
{
	FILE* Binread;
	Binread = fopen(_filename, "rb");
	int tempNu, tempNv;
	fread(&tempNu, sizeof(int), 1, Binread);	//Nu
	fread(&tempNv, sizeof(int), 1, Binread);	//Nv
	float temp;
	fread(&temp, sizeof(float), 1, Binread);	//du
	fread(&temp, sizeof(float), 1, Binread);	//dv
	fread(&temp, sizeof(float), 1, Binread);	//freq
	//fread(&temp, sizeof(float), 1, Binread);	//PowerRatio
	//FDTDRadiator::PowerRatio = double(temp);
	
	//注意 文件里是complex<float>的，目标位置是complex<double>的
	complex<float> readbuf;
	//笨笨的shift
	for (int v = 0; v<_Nv; v++) {
		for (int u = 0; u<_Nu; u++) {
			fread(&readbuf, sizeof(complex<float>), 1, Binread);
			//_Eu[u][v] = complex<double>(double(readbuf.real()), double(readbuf.imag()));

		}
	}
	for (int v = 0; v<_Nv; v++) {
		for (int u = 0; u<_Nu; u++) {
			fread(&readbuf, sizeof(complex<float>), 1, Binread);
			//_Ev[u][v] = complex<double>(double(readbuf.real()), double(readbuf.imag()));
		}
	}
	//磁场
	for (int v = 0; v<_Nv; v++) {
		for (int u = 0; u<_Nu; u++) {
			fread(&readbuf, sizeof(complex<float>), 1, Binread);
			_Hu[u][v] = complex<double>(double(readbuf.real()), double(readbuf.imag()));

		}
	}
	for (int v = 0; v<_Nv; v++) {
		for (int u = 0; u<_Nu; u++) {
			fread(&readbuf, sizeof(complex<float>), 1, Binread);
			_Hv[u][v] = complex<double>(double(readbuf.real()), double(readbuf.imag()));
		}
	}
	fclose(Binread);
}

void FDTDRadiator::SelectDT(float _dx, float _dy, float _dz){
	
	FDTDRadiator::NN = 5;
	float tempdt;
	//三个方向的间距中取最小的dt
	FDTDRadiator::Requireddt = _dx / C_Speedf * float(0.5);
	tempdt = _dy / C_Speedf * float(0.5); if (FDTDRadiator::Requireddt > tempdt) FDTDRadiator::Requireddt = tempdt;
	tempdt = _dz / C_Speedf * float(0.5); if (FDTDRadiator::Requireddt > tempdt) FDTDRadiator::Requireddt = tempdt;
	FDTDRadiator::cdt = float(1 / FDTDRadiator::Frequency / (NN*4.0));
	while (FDTDRadiator::cdt > FDTDRadiator::Requireddt) {
		FDTDRadiator::NN++;
		FDTDRadiator::cdt = float(1.0 / Frequency / (NN*4.0));
	}

}

void FDTDRadiator::SetUpLowOrderVlasovRadiator(int _WG_m, int _WG_n, double _Frequency, double _Radius, double _F, int _Ns, int _OmpNum) {
	//Set I Load Paras
	FDTDRadiator::WG_m = _WG_m;
	FDTDRadiator::WG_n = _WG_n;
	FDTDRadiator::Frequency = _Frequency;
	FDTDRadiator::Radius = _Radius;
	FDTDRadiator::F = _F;
	FDTDRadiator::Rotation = 0;
	FDTDRadiator::SourceKind = 1;//低阶模式
	FDTDRadiator::SourceType = 1;//TE模式
	FDTDRadiator::Ns = _Ns;
	FDTDRadiator::Nx_exc = _Ns;
	FDTDRadiator::Ny_exc = _Ns;
	FDTDRadiator::OmpNum = _OmpNum;
	//logfile.open("./FDTDRadiator.log", ios::out);
}

void FDTDRadiator::SetUpAperturePlane(Position3D _AperturePosition, Vector3D _ApertureDirection, Vector3D _UDirection, Vector3D _VDirection, int _Nu, int _Nv, double _Lu, double _Lv) {
	FDTDRadiator::AperturePosition = _AperturePosition;
	FDTDRadiator::ApertureDirection = _ApertureDirection;
	FDTDRadiator::UDirection = _UDirection;
	FDTDRadiator::VDirection = _VDirection;
	FDTDRadiator::Nu = _Nu;	
	FDTDRadiator::Nv = _Nv;
	FDTDRadiator::Lu = float(_Lu);	
	FDTDRadiator::Lv = float(_Lv);
	//再申请个内存吧 先U再V便于理解
	FDTDRadiator::Eu.resize(Nu);	
	FDTDRadiator::Ev.resize(Nu);
	for (int u = 0; u < Nu; u++) {
		FDTDRadiator::Eu[u].resize(Nv);
		FDTDRadiator::Ev[u].resize(Nv);
	}
	//归零初始化
	for (int u = 0; u < Nu; u++) {
		for (int v = 0; v < Nv; v++) {
			FDTDRadiator::Eu[u][v] = complex<float>(0.0, 0.0);
			FDTDRadiator::Ev[u][v] = complex<float>(0.0, 0.0);
		}
	}
	//磁场
	FDTDRadiator::Hu.resize(Nu);
	FDTDRadiator::Hv.resize(Nu);
	for (int u = 0; u < Nu; u++) {
		FDTDRadiator::Hu[u].resize(Nv);
		FDTDRadiator::Hv[u].resize(Nv);
	}
	//归零初始化
	for (int u = 0; u < Nu; u++) {
		for (int v = 0; v < Nv; v++) {
			FDTDRadiator::Hu[u][v] = complex<float>(0.0, 0.0);
			FDTDRadiator::Hv[u][v] = complex<float>(0.0, 0.0);
		}
	}
}

void FDTDRadiator::WriteApertureDataToFile(const char* _filename) {
	FILE* Binwrite;
	Binwrite = fopen(_filename, "wb");
	float temp;
	fwrite(&(FDTDRadiator::Nu), sizeof(int), 1, Binwrite);
	fwrite(&(FDTDRadiator::Nv), sizeof(int), 1, Binwrite);
	temp = Lu / (Nu-1);	//du
	fwrite(&temp, sizeof(float), 1, Binwrite);
	temp = Lv / (Nv-1); //dv
	fwrite(&temp, sizeof(float), 1, Binwrite);
	temp = float(FDTDRadiator::Frequency);
	fwrite(&temp, sizeof(float), 1, Binwrite);	//还得改哦 多频的情况
	//temp = float(FDTDRadiator::PowerRatio);
	//fwrite(&temp, sizeof(float), 1, Binwrite);	//这个是功率传输系数，还得改哦， 多频的情况
	//先写电场
	for(int v=0; v<Nv; v++){
		for (int u=0; u<Nu; u++) {
			fwrite(&(FDTDRadiator::Eu[u][v]), sizeof(complex<float>), 1, Binwrite);
		}
	}
	for (int v = 0; v<Nv; v++) {
		for (int u = 0; u<Nu; u++) {
			fwrite(&(FDTDRadiator::Ev[u][v]), sizeof(complex<float>), 1, Binwrite);
		}
	}
	//再写磁场
	for (int v = 0; v < Nv; v++) {
		for (int u = 0; u < Nu; u++) {
			fwrite(&(FDTDRadiator::Hu[u][v]), sizeof(complex<float>), 1, Binwrite);
		}
	}
	//再写磁场
	for (int v = 0; v < Nv; v++) {
		for (int u = 0; u < Nu; u++) {
			fwrite(&(FDTDRadiator::Hv[u][v]), sizeof(complex<float>), 1, Binwrite);
		}
	}
	fclose(Binwrite);
}

//读取端口激励场分布
void FDTDRadiator::SetUpExcRadiator(string _filename) {
	//读文件
	fstream file;
	file.open(_filename, ios::in);
	string temp;
	getline(file, temp);
	istringstream tempLine(temp);
	double tx, ty, tz, rx, ry, rz, rth;
	tempLine >> tx >> ty >> tz >> rx >> ry >> rz >> rth;
	double ds;
	tempLine >> Nx_exc >> Ny_exc >> ds >> ds;
	//注意，大软件的面文件顺序是先Nx 再Ny, 而FDTD计算中是先Ny再Nx，在此进行转换
	Ex_Port = Allocate_2D(Ex_Port, Ny_exc, Nx_exc);
	Ey_Port = Allocate_2D(Ey_Port, Ny_exc, Nx_exc);
	Hx_Port = Allocate_2D(Hx_Port, Ny_exc, Nx_exc);
	Hy_Port = Allocate_2D(Hy_Port, Ny_exc, Nx_exc);
	//挨行读场分布
	float exa, exd, eya, eyd, hxa, hxd, hya, hyd;
	for (int i = 0; i < Nx_exc; i++) {
		for (int j = 0; j < Ny_exc; j++) {
			
			getline(file, temp);
			istringstream perline(temp);
			perline >> exa >> exd >> eya >> eyd >> hxa >> hxd >> hya >> hyd;

			Ex_Port[j][i] = complex<float>(exa*cos(exd*PIf / 180.0), exa*sin(exd*PIf / 180.0));
			Ey_Port[j][i] = complex<float>(eya*cos(eyd*PIf / 180.0), eya*sin(eyd*PIf / 180.0));
			Hx_Port[j][i] = complex<float>(hxa*cos(hxd*PIf / 180.0), hxa*sin(hxd*PIf / 180.0));
			Hy_Port[j][i] = complex<float>(hya*cos(hyd*PIf / 180.0), hya*sin(hyd*PIf / 180.0));
		}
	}
	file.close();
}
//读取模型分布
void FDTDRadiator::SetUpModelRadiator(string _filename) {
	FILE* modelin;
	modelin = fopen(_filename.c_str(), "rb");
	int readi; float readf;
	if (Eps != NULL) {
		Free_3D(Eps, Nz_model);
	}

	fread(&readi, sizeof(int), 1, modelin);	FDTDRadiator::Nx_model = readi;
	fread(&readi, sizeof(int), 1, modelin);	FDTDRadiator::Ny_model = readi;
	fread(&readi, sizeof(int), 1, modelin);	FDTDRadiator::Nz_model = readi;

	FDTDRadiator::Eps = Allocate_3D(Eps, Nz_model, Ny_model, Nx_model);

	fread(&readi, sizeof(int), 1, modelin); FDTDRadiator::Shiftx = readi;
	fread(&readi, sizeof(int), 1, modelin); FDTDRadiator::Shifty = readi;

	fread(&readf, sizeof(float), 1, modelin); FDTDRadiator::dx = readf;
	fread(&readf, sizeof(float), 1, modelin); FDTDRadiator::dy = readf;
	fread(&readf, sizeof(float), 1, modelin); FDTDRadiator::dz = readf;

	fread(&readf, sizeof(float), 1, modelin); FDTDRadiator::cx = readf;
	fread(&readf, sizeof(float), 1, modelin); FDTDRadiator::cy = readf;
	fread(&readf, sizeof(float), 1, modelin); FDTDRadiator::cz = readf;

	for (int k = 0; k < Nz_model; k++) {
		for (int j = 0; j < Ny_model; j++) {
			for (int i = 0; i < Nx_model; i++) {
				fread(&readf,sizeof(float),1,modelin);
				Eps[k][j][i] = readf;
			}
		}
	}
	fclose(modelin);
}

void FDTDRadiator::SetUpCommonFDTD(double _freq, double _ompNum, int _N_spa, int _timemode, int _huygensmode) {
	Frequency = _freq;
	OmpNum = _ompNum;
	N_spa = _N_spa;
	timemode = _timemode;
	huygensmode = _huygensmode;
}

void FDTDRadiator::runCommonFDTD(void) {
	fstream logfile;
	logfile.open("./FDTDRadiator.log", ios::out);
	logfile << "***This is the log file for Common FDTD Radiator Running.***" << endl;
	if (returnInt)
	{
		returnInt(0, user);
	}
	CFDTD FDTD;
	FDTD.SetLogfile(&logfile);
	FDTD.SetReturnFloat(FDTDRadiator::returnFloat, FDTDRadiator::user);
	logfile << "   Model Domain Size, Nx: " << Nx_model << ", Ny: " << Ny_model << ", Nz: " << Nz_model << endl;
	logfile << "   Leads to " << Nx_model*Ny_model*Nz_model / 1.0e6 << " Mcells" << endl;
	if (FDTDRadiator::timemode == 0) {//点频形式
		SelectDT(dx, dy, dz);
		logfile << "  Discrete interval in X: " << dx << "m, in Y: " << dy << "m, in Z: " << dz << "m." << endl;
		logfile << "  Discrete interval in lambda, X: " << dx*FDTDRadiator::Frequency / C_Speedf << ", Y: " << dy*FDTDRadiator::Frequency / C_Speedf << ", Z: " << dz*FDTDRadiator::Frequency / C_Speedf << "." << endl;
		logfile << "  Required max dt: " << FDTDRadiator::Requireddt << ", Used dt: " << FDTDRadiator::cdt << ", T contains" << 1.0 / FDTDRadiator::cdt / FDTDRadiator::Frequency << "dt. Ref:" << FDTDRadiator::NN * 4 << endl;
		
		//初始化
		FDTD.Initial(OmpNum, NN, Nx_model, Ny_model, Nz_model, N_spa, 8, cdt, dx, dy, dz, Frequency);
		//申请内存
		FDTD.MemoryAllocate(timemode, 1, 0);
		//导入模型
		FDTD.SetupModel(Eps,cx,cy,cz);
		//设置激励
		FDTD.SetExcPort(N_spa, Nx_exc, Ny_exc, Shiftx+N_spa, Shifty+N_spa, Ex_Port, Ey_Port, Hx_Port, Hy_Port);
		if (returnInt)
		{
			returnInt(1, user);//开始FDTD 计算
		}
		FDTD.Update();
		if (returnInt)
		{
			returnInt(2, user);//完成FDTD 计算
		}

	}
	logfile.close();
}

void FDTDRadiator::run() {
	//低阶TE模式
	fstream logfile;
	logfile.open("./FDTDRadiator.log", ios::out);
	logfile << "***This is the log file for FDTD Radiator Running.***" << endl;

	if (returnInt) // 如果没有注册则不回调
	{
		returnInt(0, user);
	}

	//不同的频点-高阶模的交散半径不一样 最好逐点频计算
	//不同的频点-低阶模的变化要小一些 可以带宽计算
	CFDTD FDTD;
	FDTD.SetLogfile(&logfile);			//设置logfile
	FDTD.SetReturnFloat(FDTDRadiator::returnFloat, FDTDRadiator::user);	//设置回调函数
	ApertureRadiation HuygensPro;
	HuygensPro.SetLogfile(&logfile);	//设置logfile
	HuygensPro.SetReturnFloat(FDTDRadiator::returnFloat, FDTDRadiator::user);	//设置回调函数
	Position3D SourceCenter;


	//低阶辐射器运行过程
	if (SourceKind == 1 && SourceType == 1) {
		//Part I Set Port Paras
		SourceModeGenerationT Port;
		Port.SetSource_Circular(SourceKind, SourceType, Rotation, WG_m, WG_n, Frequency, Radius);
		Port.SetOutputProperty(Ns);
		Port.FieldCalculation_CircularT();

		//Buffer for Read Port Fields 
		vector<vector<complex<double>>> Ex0;
		vector<vector<complex<double>>> Ey0;
		vector<vector<complex<double>>> Hx0;
		vector<vector<complex<double>>> Hy0;
		Port.GetEX(Ex0);	Port.GetEY(Ey0);
		Port.GetHX(Hx0);	Port.GetHY(Hy0);
		//Write WaveGuide Mode Field Distribution into FDTD Buffer
		//Truncate Double into Float

		Ex_Port = Allocate_2D(Ex_Port, Ns, Ns);	Ey_Port = Allocate_2D(Ey_Port, Ns, Ns);
		Hx_Port = Allocate_2D(Hx_Port, Ns, Ns);	Hy_Port = Allocate_2D(Hy_Port, Ns, Ns);
		for (int j = 0; j < Ns; j++) {
			for (int i = 0; i < Ns; i++) {
				Ex_Port[j][i] = complex<float>(float(Ex0[i][j].real()), float(Ex0[i][j].imag()));
				Ey_Port[j][i] = complex<float>(float(Ey0[i][j].real()), float(Ey0[i][j].imag()));
				Hx_Port[j][i] = complex<float>(float(Hx0[i][j].real()), float(Hx0[i][j].imag()));
				Hy_Port[j][i] = complex<float>(float(Hy0[i][j].real()), float(Hy0[i][j].imag()));
			}
		}
		//cout << " Excitation Port Field Generated" << endl;
		logfile << " Excitation Port Field Generated" << endl;
		/*
		FILE *WGField;
		WGField = fopen("./WaveguideField.dat", "wb");

		fwrite(&Ns, sizeof(int), 1, WGField);
		for (int j = 0; j < Ns; j++) {
			for (int i = 0; i < Ns; i++) {
				fwrite(&Ex_Port[j][i], sizeof(complex<float>), 1, WGField);	//Ex.real
				fwrite(&Ey_Port[j][i], sizeof(complex<float>), 1, WGField);	//1
				fwrite(&Hx_Port[j][i], sizeof(complex<float>), 1, WGField);	//2
				fwrite(&Hy_Port[j][i], sizeof(complex<float>), 1, WGField);	//3
			}
		}
		fclose(WGField);
		//cout << " Output Excitation Port Field done!" << endl;
		logfile << " Output Excitation Port Field done!" << endl;
		*/
		Ex0.~vector(); Ey0.~vector(); Hx0.~vector(); Hy0.~vector();
		//Ready to Feed FDTD
		Port.GetCircularWaveguideProperty(PHId, THETAd, Rcd, Lcd);
		//Check Parameters OK
		FDTDRadiator::Rc = float(Rcd);
		FDTDRadiator::Lc = float(Lcd);
		FDTDRadiator::Lp = Lc / float(4.0);
		//cout << " Cut Length Lc: " << Lc*1.0e3<<"mm " << endl;
		logfile<< " Cut Length Lc: " << FDTDRadiator::Lc*float(1.0e3) << "mm " << endl;

		Port.~SourceModeGenerationT();

		//Part II Setup Model
		//建立辐射器的剖分
		//SetUp Radiator Mesh
		CRadiator Radiator(FDTDRadiator::SourceKind, FDTDRadiator::SourceType, 0, FDTDRadiator::WG_m, FDTDRadiator::WG_n, FDTDRadiator::Ns, FDTDRadiator::Ns, C_Speedf / float(Frequency), float(FDTDRadiator::Radius), FDTDRadiator::Lc, FDTDRadiator::Lp, 0, 0, 0, 0, 0, 0);
		Radiator.SetLogfile(&logfile);
		if (FDTDRadiator::SourceKind == 1 && FDTDRadiator::SourceType == 1) {//低阶TE圆电辐射器
			Radiator.ResetRadiator_RoundL(float(C_Speed / FDTDRadiator::Frequency), float(FDTDRadiator::Radius), FDTDRadiator::Ns, FDTDRadiator::Ns, FDTDRadiator::N_spa, FDTDRadiator::Lc, FDTDRadiator::Lp);
			Radiator.GenerateCellArray_RoundL();
			if (FDTDRadiator::F > 0) {
				Radiator.SetFirstMirror_RoundL(FDTDRadiator::F); //待查验
			}
			//cout << " LowOrder Radiator and Mirror Parameters:" << endl;
			//cout << "   Radius:" << Radius << "m, Cut Propagation Section Length: " << Lc << "m, Mirror Height: " << Radiator.MirrorHeight << "m, Mirror Focus Length:" << F << "m." << endl;
			//cout << "   Model Domain Size, Nx: " << Radiator.Nx_model << ", Ny: " << Radiator.Ny_model << ", Nz: " << Radiator.Nz_model << endl;
			//cout << "   Leads to " << Radiator.Nx_model*Radiator.Ny_model*Radiator.Nz_model / 1.0e6 << " Mcells" << endl;
			logfile << " LowOrder Radiator and Mirror Parameters:" << endl;
			logfile << "   Radius:" << FDTDRadiator::Radius << "m, Cut Propagation Section Length: " << FDTDRadiator::Lc << "m, Mirror Height: " << Radiator.MirrorHeight << "m, Mirror Focus Length:" << FDTDRadiator::F << "m." << endl;
			logfile << "   Model Domain Size, Nx: " << Radiator.Nx_model << ", Ny: " << Radiator.Ny_model << ", Nz: " << Radiator.Nz_model << endl;
			logfile << "   Leads to " << Radiator.Nx_model*Radiator.Ny_model*Radiator.Nz_model / 1.0e6 << " Mcells" << endl;

			Nx_model = Radiator.Nx_model;
			Ny_model = Radiator.Ny_model;
			Nz_model = Radiator.Nz_model;
			dx = Radiator.dx;
			dy = Radiator.dy;
			dz = Radiator.dz;


			//int _timemode, int _Nfreq, float _BW
			//选取DT 使其恰好是中心频率的4倍整数分之一
			SelectDT(dx, dy, dz);

			logfile << "  Discrete interval in X: " << dx << "m, in Y: " << dy << "m, in Z: " << dz << "m." << endl;
			logfile << "  Discrete interval in lambda, X: " << dx*FDTDRadiator::Frequency / C_Speedf << ", Y: " << dy*FDTDRadiator::Frequency / C_Speedf << ", Z: " << dz*FDTDRadiator::Frequency / C_Speedf << "." << endl;
			logfile << "  Required max dt: " << FDTDRadiator::Requireddt << ", Used dt: " << FDTDRadiator::cdt << ", T contains" << 1.0 / FDTDRadiator::cdt / FDTDRadiator::Frequency << "dt. Ref:" << FDTDRadiator::NN * 4 << endl;
			//ParameterPass
			FDTD.Initial(FDTDRadiator::OmpNum, FDTDRadiator::NN, Nx_model, Ny_model, Nz_model, FDTDRadiator::N_spa, 8, FDTDRadiator::cdt, dx, dy, dz, FDTDRadiator::Frequency);
			FDTD.MemoryAllocate(0, 1, 0);	//这个是单频的形式		
			//int _threadNum, int _NN, int _Nx_model, int _Ny_model, int _Nz_model, int _Nspa, int _num_pml, float _dt, float _dx, float _dy, float _dz, float _Frequency
			if (FDTDRadiator::F>0) {
				Nx_exc = Ns;	Ny_exc = Ns;
				Nz_exc = Radiator.Pz_model;
				Shiftx = Radiator.Px_model + N_spa;
				Shifty = Radiator.Py_model + N_spa;
			}
			else {
				Nx_exc = Ns;	Ny_exc = Ns;
				Nz_exc = Radiator.Pz_model;
				Shiftx = Radiator.Px_model + 2 * N_spa;
				Shifty = Radiator.Py_model + 2 * N_spa;
			}

			//Load PEC Structure From Radiator
			FDTD.SetupModel(Radiator.EpsMap, Radiator.Esig);
			//int*** _EpsMap, float*** _Esig
			FDTD.SetExcPort(Nz_exc, Nx_exc, Ny_exc, Shiftx, Shifty, Ex_Port, Ey_Port, Hx_Port, Hy_Port);

			Radiator.~CRadiator();
		}
	}
	else {
	

	}


	
	//SetUpHuygens - Freq;
	HuygensPro.SetUpFrequency(FDTD.Nfreq, FDTD.freq, FDTD.BW);
	HuygensPro.SetUpPropagatedAperture(FDTDRadiator::AperturePosition, FDTDRadiator::ApertureDirection, FDTDRadiator::UDirection, FDTDRadiator::VDirection, FDTDRadiator::Lu, FDTDRadiator::Lv, FDTDRadiator::Nu, FDTDRadiator::Nv);
	
	//int _Nz_exc, int _Nx_exc, int _Ny_exc, int _Px_exc, int Py_exc, complex<float>** _Ex_Port, complex<float>** _Ey_Port, complex<float>** _Hx_Port, complex<float>** _Hy_Port

	//SetUpHuygens - Position 注意：低阶辐射器本身朝-X方向，一级镜朝+X方向
	SourceCenter.setX(-(Shiftx - N_spa*2)*dx*0.5);
	SourceCenter.setY(0.0);
	SourceCenter.setZ(FDTD.Nz_DFT*FDTD.dz*0.5 - FDTDRadiator::Lp - FDTD.N_spa*FDTD.dz);
	//cout << "SourceCenter, X, Y, Z, m: " << SourceCenter.X() << " " << SourceCenter.Y() << " " << SourceCenter.Z() << endl;
	//cout << "Huygens Box Span, X, Y, Z, mm: " << FDTD.Nx_DFT*FDTD.dx << " " << FDTD.Ny_DFT*FDTD.dy << " " << FDTD.Nz_DFT*FDTD.dz << endl;
	logfile << "SourceCenter, X, Y, Z, m: " << SourceCenter.X() << " " << SourceCenter.Y() << " " << SourceCenter.Z() << endl;
	logfile << "Huygens Box Span, X, Y, Z, mm: " << FDTD.Nx_DFT*FDTD.dx << " " << FDTD.Ny_DFT*FDTD.dy << " " << FDTD.Nz_DFT*FDTD.dz << endl;
	HuygensPro.SetUpHuygens(SourceCenter, FDTD.Nx_DFT, FDTD.Ny_DFT, FDTD.Nz_DFT, FDTD.dx, FDTD.dy, FDTD.dz);
		
	//getchar();
	if (returnInt)
	{
		returnInt(1, user);//开始FDTD 计算
	}
	FDTD.Update();
	if (returnInt)
	{
		returnInt(2, user);//开始Huygens外推计算
	}
	HuygensPro.Propagation5FaceBox(FDTDRadiator::Eu, FDTDRadiator::Ev,FDTDRadiator::Hu,FDTDRadiator::Hv,FDTDRadiator::PowerRatio, FDTD.HuygensBoxData, FDTDRadiator::OmpNum, 0);
	FDTDRadiator::WriteApertureDataToFile("./PropagatedEField.dat");

	//cout << "Calculation Completed. ENTER to exit. " << endl;
	logfile << "Calculation Complete." << endl;
	logfile << "Clean Memories" << endl;
	if (returnInt)
	{
		returnInt(3, user);//完成计算

	}
	FDTD.~FDTD();
	HuygensPro.~ApertureRadiation();
	logfile.close();	//写文件


	//写回调位置！！！！！！！！！！！！！！！！！

}