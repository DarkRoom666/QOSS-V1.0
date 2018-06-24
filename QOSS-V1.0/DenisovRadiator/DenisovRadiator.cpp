#include "../DenisovRadiator/DenisovRadiator.h"
#include "../DenisovRadiator/CodeJin/OMTfunctions.h"
#include "../DenisovRadiator/CodeJin/Mathematical_Functions_Jin.h"
#include "../DenisovRadiator/CodeJin/SourceModeGenerationD.h"	//modified class for Denisov Design
#include "../DenisovRadiator/CodeJin/SourceModeGenerationT.h"	//modified class for FDTD output
#include <fstream>
using namespace calculation;

DenisovRadiator::DenisovRadiator()
{
	frequency = 140e9;
	radius = 16.8e-3;
	m = 22;
	n = 6;
	Denisovlog.open("./Denisov.log",ios::out);
}

DenisovRadiator::~DenisovRadiator()
{
	Denisovlog.close();
}

void DenisovRadiator::SetupDenisov(double _freq, double _radius, int _m, int _n) 
{	
	frequency = _freq;
	radius = _radius;
	m = _m;
	n = _n;
	//m = 6; n = 2;
	Denisovlog << "Denisov Basic Parameters Set" << endl;
	QString message = "Denisov Basic Paremeters Set";
	emit SendText(message);
	emit SendValue(100.0);
	
	Denisovlog << "Test for Root Finding" << endl;
	Denisovlog << m << "th order, " << n << "th zero points is " << rootbessel(m, n) << endl;
	Denisovlog << m << "th order, " << n << "th zero points of derivate is " << rootdbessel(m, n) << endl;

}

void DenisovRadiator::SetTangentialEField(int _Nx, int _Ny) {
	DenisovRadiator::Nx = _Nx;
	DenisovRadiator::Ny = _Ny;
	//Allocate Memory
	DenisovRadiator::Ex.resize(Nx);
	DenisovRadiator::Ey.resize(Nx);
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			Ex[i].resize(Ny);
			Ey[i].resize(Ny);
		}
	}
}

void DenisovRadiator::SetSurfaceCurrent(int _Nphi, int _Sparse) {
	DenisovRadiator::Nphi = _Nphi;
	DenisovRadiator::Nsparse = _Sparse;
	DenisovRadiator::Nheight = DenisovRadiator::Nz / _Sparse;
	DenisovRadiator::J.resize(DenisovRadiator::Nphi);
	for (int p = 0; p < DenisovRadiator::Nphi; p++) {
		DenisovRadiator::J[p].resize(DenisovRadiator::Nheight);
	}
}

void DenisovRadiator::GetTangentialEField(vector<vector<complex<double>>> &_Ex, vector<vector<complex<double>>> &_Ey, int _Nx, int _Ny) {
	for (int i = 0; i < _Nx; i++) {
		for (int j = 0; j < _Ny; j++) {
			_Ex[i][j] = Ex[i][j];
			_Ey[i][j] = Ey[i][j];
		}
	}
}

void DenisovRadiator::GetSurfaceCurrent(vector<vector<double>> &_J, int _Nphi, int _Nheight) {
	for (int p = 0; p < _Nphi; p++) {
		for (int h = 0; h < _Nheight; h++) {
			_J[p][h] = J[p][h];
		}
	}
}

void DenisovRadiator::SetupDenisovParas(double _delBeta1, double _delBeta2, int _l1, int _l2, double _z0_1, double _z0_2, double _ls_1, double _ls_2, double _lc_1, double _lc_2, double _Hz, double _dz, int _Nz, double _mag1, double _mag2) {
	DenisovRadiator::delBeta1 = _delBeta1;
	DenisovRadiator::delBeta2 = _delBeta2;
	DenisovRadiator::l1 = _l1;
	DenisovRadiator::l2 = _l2;
	DenisovRadiator::z0_1 = _z0_1;
	DenisovRadiator::z0_2 = _z0_2;
	DenisovRadiator::ls_1 = _ls_1;
	DenisovRadiator::ls_2 = _ls_2;
	DenisovRadiator::lc_1 = _lc_1;
	DenisovRadiator::lc_2 = _lc_2;
	DenisovRadiator::Hz = _Hz;
	DenisovRadiator::dz = _dz;
	DenisovRadiator::Nz = _Nz;
	DenisovRadiator::mag1 = _mag1;
	DenisovRadiator::mag2 = _mag2;
	DenisovRadiator::Zz.resize(DenisovRadiator::Nz);	for (int i = 0; i < DenisovRadiator::Nz; i++) { DenisovRadiator::Zz[i] = (i + 0.5)*dz; }
	DenisovRadiator::sig1.resize(DenisovRadiator::Nz);	for (int i = 0; i < DenisovRadiator::Nz; i++) { DenisovRadiator::sig1[i] = sig(mag1, z0_1, lc_1, ls_1, Zz[i]); }
	DenisovRadiator::sig2.resize(DenisovRadiator::Nz);	for (int i = 0; i < DenisovRadiator::Nz; i++) { DenisovRadiator::sig2[i] = sig(mag2, z0_2, lc_2, ls_2, Zz[i]); }
	DenisovRadiator::dsig1.resize(DenisovRadiator::Nz);	for (int i = 0; i < DenisovRadiator::Nz; i++) { DenisovRadiator::dsig1[i] = dsig(mag1, z0_1, lc_1, ls_1, Zz[i]); }
	DenisovRadiator::dsig2.resize(DenisovRadiator::Nz);	for (int i = 0; i < DenisovRadiator::Nz; i++) { DenisovRadiator::dsig2[i] = dsig(mag2, z0_2, lc_2, ls_2, Zz[i]); }


}

void DenisovRadiator::GetDenisovTurbulence(vector<double>& _sig1, vector<double>& _sig2) {
	for (int i = 0; i < Nz; i++) {
		_sig1[i] = DenisovRadiator::sig1[i];
		_sig2[i] = DenisovRadiator::sig2[i];
	}
}

void DenisovRadiator::SetupModeField(int _Nx, int _Ny) {

}

void DenisovRadiator::run(){
	QString message = "Denisov Simulation is Runing";
	emit SendText(message);
	double ii = 0.0;
	emit SendValue(ii);
	emit SendBtnClose();

	double lambda;  lambda = C_Speed / DenisovRadiator::frequency;
	double k0;		k0 = 2 * Pi / lambda;
	double ww;		ww = 2 * Pi * DenisovRadiator::frequency;

	Coes9.resize(0);	Coes9.resize(10);

	//����Denisov����
	//�����ģʽ�ڴ�
	for (int i = 0; i < DenisovRadiator::Nx; i++) {
		for (int j = 0; j < DenisovRadiator::Ny; j++) {
			DenisovRadiator::Ex[i][j] = 0;
			DenisovRadiator::Ey[i][j] = 0;
		}
	}
	for (int p = 0; p < DenisovRadiator::Nphi; p++) {
		for (int h = 0; h < DenisovRadiator::Nheight; h++) {
			DenisovRadiator::J[p][h] = 0.0;
		}
	}

	//����ģʽ����
	//Eigen::VectorXi MVec;
	//Eigen::VectorXi NVec;
	//Chi root of dbessel
	//Eigen::VectorXd ChiVec;
	//kp ������  
	//Eigen::VectorXd KpVec;
	//kz ������
	//Eigen::VectorXd KzVec;

	//����ģʽ����
	//vector<int> MVec;
	//vector<int> NVec;
	//Chi root of dbessel
	vector<double> ChiVec;
	//kp ������  
	vector<double> KpVec;
	//kz ������
	vector<double> KzVec;
	//ģʽϵ������
	Eigen::VectorXcd CoeVec;
	//ת������
	Eigen::MatrixXcd TransArray;
	//Store
	//Eigen::ArrayXXcd StoreCoeVec;

	if (DenisovRadiator::m == 22 && DenisovRadiator::n == 6) {
		mspan = 15;
		nspan = 9;
	}
	else if (DenisovRadiator::m == 26 && DenisovRadiator::n == 6) {
		mspan = 15;
		nspan = 9;
	}
	else if (DenisovRadiator::m == 28 && DenisovRadiator::n == 8) {
		mspan = 15;
		nspan = 9;
	}
	else if (DenisovRadiator::m == 6 && DenisovRadiator::n == 2) {
		mspan = 9;
		nspan = 3;
	}
	DenisovRadiator::mspan = mspan;
	DenisovRadiator::nspan = nspan;

	CoeVec.resize(mspan*nspan);
	ChiVec.resize(mspan*nspan);
	KpVec.resize(mspan*nspan);
	KzVec.resize(mspan*nspan);
	TransArray.resize(mspan*nspan, mspan*nspan);
	StoreCoeVec.resize(mspan*nspan, DenisovRadiator::Nz+1);

	//Set Up M.N Vector
	MVec.resize(mspan*nspan);
	NVec.resize(mspan*nspan);
	for (int mm = 0; mm < mspan; mm++) {
		for (int nn = 0; nn < nspan; nn++) {
			MVec[nn + mm*nspan] = DenisovRadiator::m + mm - (mspan - 1) / 2;
			NVec[nn + mm*nspan] = DenisovRadiator::n + nn - (nspan - 1) / 2;
		}
	}
	for (int i = 0; i < mspan*nspan; i++) {
		//Initialize CoeVec0
		if (MVec[i] == DenisovRadiator::m && NVec[i] == DenisovRadiator::n) {
			CoeVec[i] = 1.0;
		}
		else {
			CoeVec[i] = 0.0;
		}
		StoreCoeVec(i, 0) = CoeVec(i);
		//Initialize ChiVec
		ChiVec[i] = rootdbessel(MVec[i],NVec[i]);
		//Initialize KpVec
		KpVec[i] = ChiVec[i] / DenisovRadiator::radius;
		//Initialize KzVec
		if (k0*k0 - KpVec[i] * KpVec[i] > 0) {
			KzVec[i] = sqrt(k0*k0 - KpVec[i] * KpVec[i]);
		}
		else {
			KzVec[i] = 0.0;
		}
	}
	fstream Matrixout;
	//�����������OpenMp
	//Open MP
	int ompNum = 16;
	vector<int> r_s;	r_s.resize(ompNum);
	vector<int> r_n;	r_n.resize(ompNum);
	for (int i = 0; i<ompNum; i++) {
		r_n[i] = (mspan*nspan) / ompNum;
		int rest = (mspan*nspan) % ompNum;
		if (i<rest && rest != 0)  r_n[i] += 1;
		if (i == 0) {
			r_s[i] = 0;
		}
		else {
			r_s[i] = r_s[i - 1] + r_n[i - 1];
		}
	}

	for (int nn = 0; nn < Nz; nn++) {//Propagation Distance
		float zp = (nn + 0.5)*DenisovRadiator::dz;
		int rr, cc;
		TransArray.setZero();

		omp_set_num_threads(ompNum);
		#pragma omp parallel
		{	
			int id = omp_get_thread_num();
			//Fill TransArray
			for (int rr_omp = r_s[id]; rr_omp < r_s[id]+r_n[id]; rr_omp++) {	//row
				int mr = MVec[rr_omp];	//m0 in Matlab
				int nr = NVec[rr_omp];	//n0 in Matlab
				double dz_omp = DenisovRadiator::dz;
				double radius_omp = DenisovRadiator::radius;
									//TransArray(ir, ir) = -j*kkzmn_alfa(m0, n0)*exp(-j*kkzmn_alfa(m0, n0)*dz / 2);
				double middle = (-1.0*0.5)*KzVec[rr_omp] * dz_omp;
				TransArray(rr_omp, rr_omp) = (complex<double>(0.0, -1.0)) * KzVec[rr_omp] * exp(complex<double>(0.0, middle));

				complex<double> see = TransArray(rr_omp, rr_omp)*dz_omp;
				see = see;

				double sig1a = DenisovRadiator::sig1[nn] / radius_omp;
				double dsig1a = DenisovRadiator::dsig1[nn] / radius_omp;
				double sig2a = DenisovRadiator::sig2[nn] / radius_omp;
				double dsig2a = DenisovRadiator::dsig2[nn] / radius_omp;
				double freq = DenisovRadiator::frequency;
				double ra = radius_omp;
				double db1 = DenisovRadiator::delBeta1;
				double db2 = DenisovRadiator::delBeta2;
				int ll1 = DenisovRadiator::l1;
				int ll2 = DenisovRadiator::l2;
				double kzr = KzVec[rr_omp];
				double Chir = ChiVec[rr_omp];

				for (int cc_omp = 0; cc_omp < mspan*nspan; cc_omp++) {	//column
					int mc = MVec[cc_omp];	//m1 in Matlab
					int nc = NVec[cc_omp];	//m1 in Matlab
					complex < double > Temp;


					double kzc = KzVec[cc_omp];

					double Chic = ChiVec[cc_omp];

					Temp = CQR(freq, ra, db1, db2, ll1, ll2, zp, sig1a, dsig2a, sig2a, dsig2a, kzr, kzc, Chir, Chic, mr, mc);
					//Temp = CQR(freq,ra,db1,db2,ll1,ll2,)
					Temp = Temp * exp(complex<double>(0.0, -1.0*kzc * dz_omp*0.5));
					TransArray(rr_omp, cc_omp) = TransArray(rr_omp, cc_omp) + Temp;
					Temp = Temp;
				}//cc
				see = see;
			}//rr
		}//omp

		CoeVec = CoeVec + (TransArray*CoeVec)*DenisovRadiator::dz;

		for (int rr = 0; rr < mspan*nspan; rr++) {
			StoreCoeVec(rr, nn + 1) = CoeVec(rr);
		}
	/*
	for (int nn = 0; nn < Nz; nn++) {//Propagation Distance
		float zp = (nn + 0.5)*DenisovRadiator::dz;
		int rr, cc;
		TransArray.setZero();

		//Fill TransArray
		for (rr = 0; rr < mspan*nspan; rr++) {	//row
			int mr = MVec[rr];	//m0 in Matlab
			int nr = NVec[rr];	//n0 in Matlab
			//TransArray(ir, ir) = -j*kkzmn_alfa(m0, n0)*exp(-j*kkzmn_alfa(m0, n0)*dz / 2);
			double middle = (-1.0*0.5)*KzVec[rr] * DenisovRadiator::dz;
			TransArray(rr, rr) = (complex<double>(0.0, -1.0)) * KzVec[rr] * exp(complex<double>(0.0, middle));
			
			complex<double> see = TransArray(rr, rr)*DenisovRadiator::dz;
			see = see;
			
			double sig1a = DenisovRadiator::sig1[nn] / DenisovRadiator::radius;
			double dsig1a = DenisovRadiator::dsig1[nn] / DenisovRadiator::radius;
			double sig2a = DenisovRadiator::sig2[nn] / DenisovRadiator::radius;
			double dsig2a = DenisovRadiator::dsig2[nn] / DenisovRadiator::radius;
			double freq = DenisovRadiator::frequency;
			double ra = DenisovRadiator::radius;
			double db1 = DenisovRadiator::delBeta1;
			double db2 = DenisovRadiator::delBeta2;
			int ll1 = DenisovRadiator::l1;
			int ll2 = DenisovRadiator::l2;
			double kzr = KzVec[rr];
			double Chir = ChiVec[rr];

			for (cc = 0; cc < mspan*nspan; cc++) {	//column
				int mc = MVec[cc];	//m1 in Matlab
				int nc = NVec[cc];	//m1 in Matlab
				complex < double > Temp;


				double kzc = KzVec[cc];

				double Chic = ChiVec[cc];

				Temp = CQR(freq, ra, db1, db2, ll1, ll2, zp, sig1a, dsig2a, sig2a, dsig2a, kzr,kzc,Chir,Chic,mr,mc);
				//Temp = CQR(freq,ra,db1,db2,ll1,ll2,)
				Temp = Temp * exp(complex<double>(0.0, -1.0*kzc * DenisovRadiator::dz*0.5));
				TransArray(rr, cc) = TransArray(rr, cc) + Temp;
				Temp = Temp;
			}//cc
			see = see;
		}//rr

	
		CoeVec = CoeVec + (TransArray*CoeVec)*DenisovRadiator::dz;
		for (int rr = 0; rr < mspan*nspan; rr++) {
			StoreCoeVec(rr, nn + 1) = CoeVec(rr);
		}
		*/

		double temp;
		//Take Coefficients
		int midN = (nspan - 1) / 2;
		int midM = (mspan - 1) / 2;
		if (DenisovRadiator::m > 6 ) {
			//11 TE m-2 n+1 22 6 - 20 7	%TE m - 2 n + 1  22 6 - 20 7
			DenisovRadiator::Coes9[0] = pow(abs(CoeVec(midN + 1 + (midM - 2)*nspan)), 2);
			//12 TE m+1 n 22 6 - 23 6	%TE m + 1 n  22 6 - 23 6
			DenisovRadiator::Coes9[1] = pow(abs(CoeVec(midN + 0 + (midM + 1)*nspan)), 2);
			//13 TE m+4 n-1 22 6 - 26 5	%TE m + 4 n - 1 22 6 - 26 5
			DenisovRadiator::Coes9[2] = pow(abs(CoeVec(midN - 1 + (midM + 4)*nspan)), 2);
			//21 TE m-3 n+1 22 6 - 19 7	%TE m - 3 n + 1  22 6 - 19 7
			DenisovRadiator::Coes9[3] = pow(abs(CoeVec(midN + 1 + (midM - 3)*nspan)), 2);
			//22 TE m n 22 6 - 22 6		%TE m n    22 6 - 22 6
			DenisovRadiator::Coes9[4] = pow(abs(CoeVec(midN + 0 + (midM + 0)*nspan)), 2);
			//temp = abs(CoeVec[midN + 0 + (midM)*nspan]); temp = temp*temp;
			//DenisovRadiator::Coes9[4] = temp;
			//23 TE m+3 n-1 22 6 - 25 5	%TE m + 3 n - 1   22 6 - 25 5
			DenisovRadiator::Coes9[5] = pow(abs(CoeVec(midN - 1 + (midM + 3)*nspan)), 2);
			//31 TE m-4 n+1 22 6 - 18 7	%TE m - 4 n + 1  22 6 - 18 7
			DenisovRadiator::Coes9[6] = pow(abs(CoeVec(midN + 1 + (midM - 4)*nspan)), 2);
			//32 TE m-1 n 22 6 - 21 6	%TE m - 1 n  22 6 - 21 6
			DenisovRadiator::Coes9[7] = pow(abs(CoeVec(midN + 0 + (midM - 1)*nspan)), 2);
			//33 TE m+2 n-1 22 6 - 24 5	%TE m + 2 n - 1 22 6 - 24 5
			DenisovRadiator::Coes9[8] = pow(abs(CoeVec(midN - 1 + (midM + 2)*nspan)), 2);

			DenisovRadiator::Coes9[9] = CoeVec.cwiseAbs2().sum();

		}


		//�Ѽ���õ��ĸ�����ģʽ�Ĺ��ʱ��ʷ��ؽ���
		//total main neighbor, corner
		double CoeTotal = Coes9[9];
		double CoeMain = Coes9[4];
		double CoeNeighbor = Coes9[1] + Coes9[3] + Coes9[5] + Coes9[7];
		double CoeCorner = Coes9[0] + Coes9[2] + Coes9[6] + Coes9[8];
		emit SendCoefficients(CoeTotal, CoeMain, CoeNeighbor, CoeCorner, nn);


		if ((nn + 1) % (Nz/5) == 0) { 
			omp_set_num_threads(ompNum);
			#pragma omp parallel
			{
				int id = omp_get_thread_num();
				for (int rr_omp = r_s[id]; rr_omp < r_s[id]+r_n[id]; rr_omp++) {
					if (abs(CoeVec(rr_omp))*abs(CoeVec(rr_omp)) > 0.01) {
						int mr = MVec[rr_omp];
						int nr = NVec[rr_omp];
						SourceModeGenerationD Source(2, 1, 2, mr, nr, DenisovRadiator::frequency, DenisovRadiator::radius, 0, 0, DenisovRadiator::Nx);
						//2, 1, 2, m, n, frequency, radius, 0, 0, Nx
						//Source.SetSource_Circular(2,1,2,mr,nr,DenisovRadiator::frequency,DenisovRadiator::radius);
						//Source.SetOutputProperty(Nx);
						Source.FieldCalculation_Circular();

						vector<vector<complex<double>>> TEx;
						vector<vector<complex<double>>> TEy;
						TEx.resize(Nx);	TEy.resize(Ny);
						for (int i = 0; i < Nx; i++) {
							TEx[i].resize(Ny);
							TEy[i].resize(Nz);
						}

						Source.GetEX(TEx); Source.GetEY(TEy);

						for (int i = 0; i < Nx; i++) {
							for (int j = 0; j < Ny; j++) {
								Ex[i][j] += TEx[i][j] * CoeVec[rr_omp];
								Ey[i][j] += TEy[i][j] * CoeVec[rr_omp];
							}
						}
					}//if
				}//rr
			}//omp
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					//Ex[i][j] = sqrt(Ex[i][j].real()*Ex[i][j].real()+ Ex[i][j].imag()*Ex[i][j].imag());
					//Ey[i][j] = sqrt(Ey[i][j].real()*Ey[i][j].real() + Ey[i][j].imag()*Ey[i][j].imag());
				}
			}
			
			//finished add
			//emit SendTangentialEField(Ex, Ey, Nx, Ny);
			emit SendFieldUpdate();

			emit SendValue(double(nn*1.0 / Nz * 100));
		}
	}//nn

	//emit SendValue(double(Max));
	message = "Denisov Simulation is Done";
	emit SendText(message);
	//for (int p = 0; p < DenisovRadiator::Nphi; p++) {
		//for (int h = 0; h < DenisovRadiator::Nheight; h++) {
			//DenisovRadiator::J[p][h] = h - p*2;
		//}
	//}
	message = "Generating Surface Currents";
	emit SendText(message);
	
	//����������� _omp
	for(int h = 0; h<DenisovRadiator::Nheight;h++){//Height
		vector<complex<double>> HPhi_r;	HPhi_r.resize(DenisovRadiator::Nphi);
		vector<complex<double>> HZ_r;	HZ_r.resize(DenisovRadiator::Nphi);
		for (int p = 0; p < DenisovRadiator::Nphi; p++) {
			HPhi_r[p] = complex<double>(0.0, 0.0);
			HZ_r[p] = complex<double>(0.0, 0.0);
		}
		omp_set_num_threads(ompNum);
		#pragma omp parallel
		{
			int id = omp_get_thread_num();
			int hh = h*DenisovRadiator::Nsparse + Nsparse/2;
			for (int rr_omp = r_s[id]; rr_omp < r_s[id]+r_n[id]; rr_omp++) {//Mode
				complex<double> Coe = StoreCoeVec(rr_omp, hh);
				if (abs(CoeVec(rr_omp))*abs(CoeVec(rr_omp)) > 0.01) {//Judge
					int mr = MVec[rr_omp];
					int nr = NVec[rr_omp];
					SourceModeGenerationD Source(2, 1, 2, mr, nr, DenisovRadiator::frequency, DenisovRadiator::radius, 0, 0, DenisovRadiator::Nx);

					vector<complex<double>> THPhi;
					vector<complex<double>> THz;
					THPhi.resize(DenisovRadiator::Nphi);	THz.resize(DenisovRadiator::Nphi);

					Source.GetJ_R(THPhi, THz, DenisovRadiator::Nphi);

					for (int p = 0; p < DenisovRadiator::Nphi; p++) {
						HPhi_r[p] += THPhi[p] * Coe;
						HZ_r[p] += THz[p] * Coe;
					}
				}//if
			}//rr
		}
		for (int p = 0; p < Nphi; p++) {
			DenisovRadiator::J[p][h] = sqrt(HPhi_r[p].real()*HPhi_r[p].real() + HPhi_r[p].imag()*HPhi_r[p].imag() + HZ_r[p].real()*HZ_r[p].real() + HZ_r[p].imag()*HZ_r[p].imag());
		}
		emit SendValue(double(h*1.0 / DenisovRadiator::Nheight * 100));
	}//for h
	//Send SurfaceCurrent
	emit SendCurrentUpdate();

	emit SendBtnOpen();
	emit SendValue(ii);
}

void DenisovRadiator::GetExcitationField(vector<vector<complex<double>>> &_Ex, vector<vector<complex<double>>> &_Ey, vector<vector<complex<double>>> &_Hx, vector<vector<complex<double>>> &_Hy, int _index, int _N) {
	vector<vector<complex<double>>> tempEx;
	vector<vector<complex<double>>> tempEy;
	vector<vector<complex<double>>> tempHx;
	vector<vector<complex<double>>> tempHy;
	//vector<complex<double>> CoeVec;
	//�м�������洢���ɵĳ����
	tempEx.resize(_N); for (int i = 0; i < _N; i++) { tempEx[i].resize(_N); }
	tempEy.resize(_N); for (int i = 0; i < _N; i++) { tempEy[i].resize(_N); }
	tempHx.resize(_N); for (int i = 0; i < _N; i++) { tempHx[i].resize(_N); }
	tempHy.resize(_N); for (int i = 0; i < _N; i++) { tempHy[i].resize(_N); }

	Eigen::VectorXcd CoeVec;
	CoeVec.resize(mspan*nspan);
	//��ȡ���ֵ
	for (int i = 0; i < mspan*nspan; i++) {
		CoeVec(i) = StoreCoeVec(i,_index);
	}
	double max = 0;
	vector<bool> cal; cal.resize(mspan*nspan);
	//�����ֵ
	for (int i = 0; i < mspan*nspan; i++) {
		if (abs(CoeVec(i)) > max) max = abs(CoeVec(i));
	}
	//����Ҫ�����ֵ(0.005*Max)
	for (int i = 0; i < mspan*nspan; i++) {
		if (abs(CoeVec(i)) > max*5e-3) cal[i] = true;
		else cal[i] = false;
	}
	//��ʼ���ɳ�
	for (int i = 0; i < mspan*nspan; i++) {
		if (cal[i]) {//�������
			int mr = MVec[i];
			int nr = NVec[i];

			SourceModeGenerationD Source(2, 1, 2, mr, nr, DenisovRadiator::frequency, DenisovRadiator::radius, 0, 0, _N);
			Source.FieldCalculation_CircularT();//���㳡
			Source.GetEX(tempEx);
			Source.GetEY(tempEy);
			Source.GetHX(tempHx);
			Source.GetHY(tempHy);
			
			for (int ii = 0; ii < _N; ii++) {
				for (int jj = 0; jj < _N; jj++) {
					_Ex[ii][jj] = _Ex[ii][jj] + tempEx[ii][jj];
					_Ey[ii][jj] = _Ey[ii][jj] + tempEy[ii][jj];
					_Hx[ii][jj] = _Hx[ii][jj] + tempHx[ii][jj];
					_Hy[ii][jj] = _Hy[ii][jj] + tempHy[ii][jj];
				}
			}
		}
	}


}
#include "moc_DenisovRadiator.cpp"