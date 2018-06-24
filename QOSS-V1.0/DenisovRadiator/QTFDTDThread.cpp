#include "QTFDTDThread.h"
//#include "FDTDRadiator.h"
#include <cmath>

using namespace userInterface;
void QTFDTDThread::run()
{
	if(!newmode){
		FDTDradiator->SetReturnInt(QTFDTDThread::setMainValue, this);	//����״̬��ѯ-�����Ļص�����
		FDTDradiator->SetReturnFloat(QTFDTDThread::setSlaverValue, this);	//����״̬����-�������Ļص�����

		if (this->m == 0) {
			FDTDradiator->SetUpLowOrderVlasovRadiator(this->m, this->n, this->frequency, this->radius, 20e-3, this->Ns, 14);
		}
		else {
		
		}


		//SourceCenter and Aperture Parameters;
		//���ó������ -Ĭ��(��λm)

		// ��ʱд�� �Ժ��е�һ�����侵���ĵ����²���
		//��Lc�õ�phi
		//Lc = 2 * Radius*tan(Pi * 0.5 - phi1);	//2*a*cotphi
		double phi; phi = atan(this->radius*2.0 / this->lcut);
		Position3D AperturePosition;
		AperturePosition.setX(-this->prodis*sin(phi));
		AperturePosition.setY(0);
		AperturePosition.setZ(this->prodis*cos(phi) + this->lcut*0.5);
		Vector3D UDirection;
		UDirection.setX(cos(phi));
		UDirection.setY(0);
		UDirection.setZ(sin(phi));
		UDirection.normalize();
		Vector3D VDirection;
		VDirection.setX(0);
		VDirection.setY(1.0);
		VDirection.setZ(0.0);
		VDirection.normalize();

		Vector3D ApertureDirection;
		ApertureDirection = ApertureDirection.crossProduct(UDirection, VDirection);

		int Nu;	int Nv;	Nu = this->Na; Nv = this->Na;
		double Lu;	Lu = this->aperlen;
		double Lv;	Lv = this->aperlen;

		//���ó�������λ��	�ھ�����		 ����ָ��ʸ��       U����ʸ��  V����ʸ�� ��������� ����
		FDTDradiator->SetUpAperturePlane(AperturePosition, ApertureDirection, UDirection, VDirection, Nu, Nv, Lu, Lv);
		//����

		//������� ע�⣺�ڴ��ڵ��ô����룬����DLL�����룡
		vector<vector<complex<double>>> Eu;
		vector<vector<complex<double>>> Ev;
		Eu.resize(Nu);	Ev.resize(Nu);
		for (int u = 0; u < Nu; u++) {
			Eu[u].resize(Nv);
			Ev[u].resize(Nv);
		}
		FDTDradiator->run();
		//Get Surface Field
		FDTDradiator->GetProFieldE(Eu, Ev, Nu, Nv);	//�Ѻ����ӿڸĳ���double ��������Ӧ
													//FDTDradiator->GetProPowerRatio(QTFDTDThread::PowerRatio);

													//���ǰ���Ѽ����꣬��ֱ�Ӷ�ȡ�ļ�	���->run ��->getProFieldE;	--���� 20180325
													//FDTDradiator->LoadProFieldE("./PropagatedEField.dat",Eu,Ev,Nu,Nv);

													//VTK��ͼ��
													/*
													GraphTrans graphTrans;
													graphTrans.setGraphTransPar(AperturePosition.X(), AperturePosition.Y(), AperturePosition.Z(),
													0, 1, 0, 15.06);
													double ds = Lu / (Nu - 1);
													field->setNM(Nu, Nu);
													field->setPlane(graphTrans, ds);
													field->setField(Eu, Ev);
													field->setShowPara(1, 1, 0);
													//field->updateData();
													*/

	}
	else {
		FDTDradiator->SetReturnInt(QTFDTDThread::setMainValue, this);	//����״̬��ѯ-�����Ļص�����
		FDTDradiator->SetReturnFloat(QTFDTDThread::setSlaverValue, this);	//����״̬����-�������Ļص�����
		FDTDradiator->SetUpModelRadiator(modelfile);
		FDTDradiator->SetUpExcRadiator(excfile);
		FDTDradiator->SetUpCommonFDTD(frequency, OmpNum, N_spa, timemode, huygensmode);
		FDTDradiator->runCommonFDTD();
	}
}

void QTFDTDThread::setMainValue(int val, void *user)
{
	((QTFDTDThread*)user)->sendMainValue(val);
}

void QTFDTDThread::setSlaverValue(float val, void *user)
{
	((QTFDTDThread*)user)->sendSlaverValue(val);
}

void QTFDTDThread::killFDTD()
{
	// free FDTD
	deleteLater();
}

void QTFDTDThread::GetProPowerRatio(double& _PowerRatio) {
	_PowerRatio = QTFDTDThread::PowerRatio;
}

void QTFDTDThread::setModelFile(string _filename) {
	modelfile = _filename;
}

void QTFDTDThread::setExcFile(string _filename) {
	excfile = _filename;

}

void QTFDTDThread::setComputation(double _freq, int _ompNum, int _N_spa, int _timemode, int _huygensmode) {
	frequency = _freq;
	OmpNum = _ompNum;
	N_spa = _N_spa;
	timemode = _timemode;
	huygensmode = _huygensmode;
	newmode = true;
	//set FDTDRadiator
}

void QTFDTDThread::setFDTDcalculated(bool _in) {
	FDTDcalculated = _in;
}

void QTFDTDThread::setNewMode(bool _in) {
	newmode = _in;
}