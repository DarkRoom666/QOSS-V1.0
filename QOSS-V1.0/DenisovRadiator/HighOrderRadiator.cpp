#include "../DenisovRadiator/HighOrderRadiator.h"

HighOrderRadiator::HighOrderRadiator(void) {

}

HighOrderRadiator::~HighOrderRadiator(void) {

}
//���÷������Ļ�������
void HighOrderRadiator::SetRadiatorBasicParas(double _Ra, double _zcut, double _phic) {
	Ra = _Ra;
	Zcut = _zcut;
	phic = _phic*Pi/180;
}
//���÷������Ŷ�����
void HighOrderRadiator::SetRadiatorTurbulenceParas(double _delBeta1, double _delBeta2, double _z0_1, double _z0_2, double _ls_1, double _ls_2, double _lc_1, double _lc_2, double _mag1, double _mag2) {
	l1 = 1;
	l2 = 3;
	delBeta1 = _delBeta1;
	z0_1 = _z0_1;
	ls_1 = _ls_1;
	lc_1 = _lc_1;
	mag1 = _mag1;
	delBeta2 = _delBeta2;
	z0_2 = _z0_2;
	ls_2 = _ls_2;
	lc_2 = _lc_2;
	mag2 = _mag2;
}
//����һ�����Ĳ���
void HighOrderRadiator::SetFirstMirrorParas(double _F1, double _F2, double _Rc, double _Phi, double _Lcut) {
	F1 = _F1;
	F2 = _F2;
	Rc = _Rc;
	Phi = _Phi;
	Lcut = _Lcut;
}
//����һ�����ĵ���
int HighOrderRadiator::GetFirstMirrorLattice(vector<vector<Vector3>> &_lattice, double _ds) {
	//������Ƶ����Ķ���
	int Nz;		//��ֱά��
	int Nphi;	//ˮƽά��
	//�жϵ�������
	Nz = int(Lcut*1.8 / _ds);	double dz = Lcut*1.8 / Nz;
	double Rr; Rr = (F1 + F2) / 2;
	Nphi = int(Pi*Rr / _ds);	double dphi = Pi / Nphi;
	//������������С
	_lattice.resize(Nz+1);
	for (int i = 0; i < Nz+1; i++) {
		_lattice[i].resize(Nphi+1);
	}
	//��д����
	double x,y,z,l;
	double t;
	double lg = F1 + sqrt((F1 + F2)*(F1 + F2) + Rc*Rc);
	for (int i = 0; i < Nz+1; i++) {
		z = Lcut*0.5 + F1 / tan(Phi) + i*dz - Lcut*0.9;
		for (int j = 0; j < Nphi+1; j++) {
			t = j*dphi - 0.03*Pi;
			l = 4 * Rc*Rc + 4 * F2*F2 - 4 * lg*lg + 8 * Rc*F2*cos(t) - 4 * t*t*Rc*Rc + 4 * t*Rc*Rc*Pi - 8 * t*Rc*lg - Pi*Pi*Rc*Rc + 4 * Pi*Rc*lg;
			l = l / (-8 * F2*sin(t) - 8 * t*Rc + 4 * Pi*Rc - 8 * lg);
			x = -Rc*cos(t) - l*sin(t);
			y = Rc*sin(t) - l*cos(t);
			_lattice[i][j].set(x,y,z);
		}
	}
	
	return 0;
}
//���ط������ĵ���-1 ע�⣬������1ά������ 
int HighOrderRadiator::GetRadiatorLattice1(vector<Vector3> &_lattice, double _ds) {
	//������Ƶ�Ķ���
	int Nphi;
	int Nz;
	//�ж���Ҫ���ٵ� (Բ����)
	Nz = int(Zcut / _ds);	double dz = Zcut / Nz;
	Nphi = int(Pi*Ra / _ds); double dphi = Pi / Nphi;
	//�ж���Ҫ���ٵ� (�пڶ�)
	int Nz1;
	Nz1 = int(Lcut*0.5 / _ds);	double dz1 = Lcut*0.5 / Nz1;
	vector<int> Nphi_cut; vector<double> dphi_cut;
	Nphi_cut.resize(Nz1);	dphi_cut.resize(Nz1);
	for (int i = 0; i < Nz1; i++) {
		double Lphi = Pi*Ra* ( (Nz1 - (i + 1))*dz1 / (Lcut*0.5) );
		Nphi_cut[i] = int(Lphi / _ds);	
		if (Nphi_cut[i] != 0) { dphi_cut[i] = Lphi / Ra / Nphi_cut[i]; }
		else { dphi_cut[i] = 0; }
	}
	//�пڶ�����Ҫ�ĵ���
	int NCut; NCut = 0;
	for (int i = 0; i < Nz1; i++) {
		NCut = NCut + Nphi_cut[i];
	}
	NCut = NCut + Nz1;//ÿһ����Phi������Nphi+1��
	//ȫ���ĵ���
	int NTotal = (Nz + 1)*(Nphi + 1) + NCut;
	_lattice.resize(NTotal);
	//�������Բ���Σ��������Ŷ���
	double x, y, z;
	double r;
	double t;
	for (int i = 0; i < Nz+1; i++) {
		z = -Zcut + i*dz;
		for (int j = 0; j < Nphi + 1; j++) {
			t = phic + j*dphi;
			r = Ra + sig(mag1, z0_1, lc_1, ls_1, z + Zcut)*cos(delBeta1*(z + Zcut) - l1*t)
				+ sig(mag2, z0_2, lc_2, ls_2, z + Zcut)*cos(delBeta2*(z + Zcut) - l2*t);
			x = r*cos(t);
			y = r*sin(t);
			//��ѭ��������phi
			_lattice[j + i*(Nphi + 1)].set(x, y, z);
		}
	}
	//Ȼ������пڶ�(�������Ŷ�)
	int Shift = (Nphi + 1)*(Nz + 1);
	int index = 0 + Shift;
	for (int i = 0; i < Nz1; i++) {//i�Ǹ߶�
		z = (i + 1)*dz1;
		//phi ����
		for (int j = 0; j < Nphi_cut[i]+1; j++) {//j��phi
			t = phic + Pi* ( (i + 1)*dz1 / (Lcut*0.5)) + j*dphi_cut[i];
			x = cos(t)*Ra;
			y = sin(t)*Ra;
			_lattice[index].set(x, y, z);
			index++;
		}
	}

	return 0;
}
//���ط������ĵ���-2 
int HighOrderRadiator::GetRadiatorLattice2(vector<Vector3> &_lattice, double _ds) {
	//������Ƶ�Ķ���
	int Nphi;
	int Nz;
	//�ж���Ҫ���ٵ� (Բ����)
	Nz = int(Zcut / _ds);	double dz = Zcut / Nz;
	Nphi = int(Pi*Ra / _ds); double dphi = Pi / Nphi;
	//�ж���Ҫ���ٵ� (�пڶ��°��)
	int Nz1 = int(Lcut*0.5 / _ds); double dz1 = Lcut*0.5 / Nz1;
	int Nphi1 = Nphi;			   double dphi1 = dphi;
	//�ж���Ҫ���ٵ� (�пڶ��ϰ��)

	int Nz2 = int(Lcut*0.5 / _ds);	double dz2 = Lcut*0.5 / Nz2;
	vector<int> Nphi_cut; vector<double> dphi_cut;
	Nphi_cut.resize(Nz2);	dphi_cut.resize(Nz2);
	for (int i = 0; i < Nz2; i++) {
		double Lphi = Pi*Ra* ((Nz2 - (i + 1))*dz2 / (Lcut*0.5));
		Nphi_cut[i] = int(Lphi / _ds);
		if (Nphi_cut[i] != 0) { dphi_cut[i] = Lphi / Ra / Nphi_cut[i]; }
		else { dphi_cut[i] = 0; }
	}
	//�пڶ��ϰ������Ҫ�ĵ���
	int NCut; NCut = 0;
	for (int i = 0; i < Nz2; i++) {
		NCut = NCut + Nphi_cut[i];
	}
	NCut = NCut + Nz2;//ÿһ����Phi������Nphi+1��
					  //ȫ���ĵ���
	int NTotal = (Nz + 1)*(Nphi + 1) + Nz1*(Nphi1+1) + NCut;
	_lattice.resize(NTotal);
	//�������Բ���Σ��������Ŷ���
	double x, y, z;
	double r;
	double t;
	for (int i = 0; i < Nz + 1; i++) {
		z = -Zcut + i*dz;
		for (int j = 0; j < Nphi + 1; j++) {
			t = phic + Pi + j*dphi;
			r = Ra + sig(mag1, z0_1, lc_1, ls_1, z + Zcut)*cos(delBeta1*(z + Zcut) - l1*t)
				+ sig(mag2, z0_2, lc_2, ls_2, z + Zcut)*cos(delBeta2*(z + Zcut) - l2*t);
			x = r*cos(t);
			y = r*sin(t);
			//��ѭ��������phi
			_lattice[j + i*(Nphi + 1)].set(x, y, z);
		}
	}
	//Ȼ������пڶ�(�������Ŷ�) �°��
	int Shift = (Nphi + 1)*(Nz + 1);
	for (int i = 0; i < Nz1; i++) {
		z =  (i+1)*dz1;
		for (int j = 0; j < Nphi1 + 1; j++) {
			t = phic + Pi + j*dphi1;
			r = Ra;
			x = r*cos(t);
			y = r*sin(t);
			//��ѭ��������phi
			_lattice[j + i*(Nphi1 + 1)+Shift].set(x, y, z);
		}
	}
	//Ȼ������пڶ�(�������Ŷ�) �ϰ��
	Shift = (Nphi + 1)*(Nz + 1) + (Nphi1 + 1)*Nz1;
	int index = 0 + Shift;
	for (int i = 0; i < Nz2; i++) {//i�Ǹ߶�
		z = (i + 1)*dz2 + Lcut*0.5;
		//phi ����
		for (int j = 0; j < Nphi_cut[i] + 1; j++) {//j��phi
			t = phic + Pi + Pi* ((i + 1)*dz2 / (Lcut*0.5)) + j*dphi_cut[i];
			x = cos(t)*Ra;
			y = sin(t)*Ra;
			_lattice[index].set(x, y, z);
			index++;
		}
	}

	return 0;

}
//���ؾֲ�λ�õİ뾶
inline double HighOrderRadiator::R(double _phi, double _t) {
	double result;
	_t = _t + Zcut;
	double sig1 = sig(mag1, z0_1, lc_1, ls_1, _t);
	double sig2 = sig(mag2, z0_2, lc_2, ls_2, _t);
	result = Ra + sig1*cos(delBeta1*_t - l1*_phi) + sig2*cos(delBeta2*_t - l2*_phi);
	return result;
}
//����λ��
inline Vector3 HighOrderRadiator::Position(double _phi, double _t) {
	Vector3 result;
	double Rr = R(_phi,_t);
	result.set(Rr*cos(_phi), Rr*sin(_phi), _t);
	return result;
}
//����_phi�ݶȷ���
inline Vector3 HighOrderRadiator::DelPhiVec(double _phi, double _t) {
	Vector3 result;
	double Rr = R(_phi, _t);
	_t = _t + Zcut;
	double sig1 = sig(mag1, z0_1, lc_1, ls_1, _t);
	double sig2 = sig(mag2, z0_2, lc_2, ls_2, _t);
	double dRp = sig1*sin(delBeta1*_t - l1*_phi)*l1 + sig2*sin(delBeta2*_t - l2*_phi)*l2;
	result.set( dRp*cos(_phi) + Rr*(-sin(_phi)), dRp*sin(_phi) + Rr*cos(_phi),0);
	result.Normalization();
	return result;
}
//����_t�ݶȷ���
inline Vector3 HighOrderRadiator::DelTVec(double _phi, double _t) {
	Vector3 result; 
	double Rr = R(_phi, _t);
	_t = _t + Zcut;
	double sig1 = sig(mag1, z0_1, lc_1, ls_1, _t);
	double sig2 = sig(mag2, z0_2, lc_2, ls_2, _t);
	double dsig1 = dsig(mag1, z0_1, lc_1, ls_1, _t);
	double dsig2 = dsig(mag2, z0_2, lc_2, ls_2, _t);
	double dRt = dsig1*cos(delBeta1*_t - l1*_phi) + sig1*delBeta1*(-sin(delBeta1*_t - l1*_t))
		       + dsig2*cos(delBeta2*_t - l2*_phi) + sig2*delBeta2*(-sin(delBeta2*_t - l2*_t));
	result.set( dRt*cos(_phi),
				dRt*sin(_phi),
				1);
	result.Normalization();
	return result;
}
//���ط�����
inline Vector3 HighOrderRadiator::Normal(double _phi, double _t) {
	Vector3 dPhi;
	Vector3 dT;
	Vector3 result;
	double Rr = R(_phi, _t);
	_t = _t + Zcut;
	double sig1 = sig(mag1, z0_1, lc_1, ls_1, _t);
	double sig2 = sig(mag2, z0_2, lc_2, ls_2, _t);
	double dsig1 = dsig(mag1, z0_1, lc_1, ls_1, _t);
	double dsig2 = dsig(mag2, z0_2, lc_2, ls_2, _t);
	double dRp = sig1*sin(delBeta1*_t - l1*_phi)*l1 + sig2*sin(delBeta2*_t - l2*_phi)*l2;
	double dRt = dsig1*cos(delBeta1*_t - l1*_phi) + sig1*delBeta1*(-sin(delBeta1*_t - l1*_t))
		+ dsig2*cos(delBeta2*_t - l2*_phi) + sig2*delBeta2*(-sin(delBeta2*_t - l2*_t));
	
	dPhi.set(dRp*cos(_phi) + Rr*(-sin(_phi)),  dRp*sin(_phi) + Rr*cos(_phi),  0);
	dPhi.Normalization();
	
	dT.set(dRt*cos(_phi), dRt*sin(_phi), 1);
	dT.Normalization();
	
	result = dPhi.Cross(dT);
	result.Normalization();
	return result;
}
//�ж��Ƿ��ڷ����������� _phi���뷶Χ 0~2Pi
inline bool HighOrderRadiator::OnRadiatorSurface(double _phi, double _t) {
	bool on;
	_t = _t + Zcut;
	//����������
	if (_t<0 || _t>Zcut + Lcut){on = false;}
	//�Ŷ�+������
	else if (_t <= Zcut) { on = true; }
	else {//�п�
		if (_phi > phic) {
			if (_t <= (_phi - phic) / (2 * Pi)*Lcut + Zcut) on = true;
			else on = false;
		}
		else {
			if (_t <= (_phi + 2 * Pi - phic) / (2 * Pi)*Lcut + Zcut) on = true;
			else on = false;
		}
	}
	return on;
}

