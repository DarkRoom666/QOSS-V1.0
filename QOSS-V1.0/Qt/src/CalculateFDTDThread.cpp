#include "../include/CalculateFDTDThread.h"
#include "../../FDTDRadiator/FDTDRadiator.h"

using namespace userInterface;
void CalculateFDTDThread::run()
{
	FDTDradiator->SetReturnInt(CalculateFDTDThread::setMainValue, this);	//����״̬��ѯ-�����Ļص�����
	FDTDradiator->SetReturnFloat(CalculateFDTDThread::setSlaverValue, this);	//����״̬����-�������Ļص�����

	FDTDradiator->SetUpLowOrderVlasovRadiator(0, 1, 44e9, 16e-3, 20e-3, 80, 4);

	//SourceCenter and Aperture Parameters;
	//���ó������ -Ĭ��(��λm)

	// ��ʱд�� �Ժ��е�һ�����侵���ĵ����²���
	Position3D AperturePosition;	AperturePosition.setX(0.02);			AperturePosition.setY(0);	AperturePosition.setZ(0.2824);
	Vector3D UDirection;			UDirection.setX(cos(15.06 * 3.1415926 / 180));	UDirection.setY(0);	UDirection.setZ(-sin(15.06 * 3.1415926 / 180));	UDirection.normalize();
	Vector3D VDirection;			VDirection.setX(0);						VDirection.setY(1.0);		VDirection.setZ(0.0);					VDirection.normalize();
	Vector3D ApertureDirection;		ApertureDirection = ApertureDirection.crossProduct(UDirection, VDirection);
	int Nu;	int Nv;	Nu = 201; Nv = 201;
	double Lu;	Lu = 0.3;	//100 lambda
	double Lv;	Lv = 0.3;	//100 lambda

											//���ó�������λ��	�ھ�����		 ����ָ��ʸ��       U����ʸ��  V����ʸ�� ��������� ����
	FDTDradiator->SetUpAperturePlane(AperturePosition, ApertureDirection, UDirection, VDirection, Nu, Nv, Lu, Lv);
	//����
	//FDTDradiator->SetProFieldFileName("./test1.dat");	//����������ֲ�������

	//������� ע�⣺�ڴ��ڵ��ô����룬����DLL�����룡
	vector<vector<complex<double>>> Eu;
	vector<vector<complex<double>>> Ev;
	Eu.resize(Nu);	Ev.resize(Nu);
	for (int u = 0; u < Nu; u++) {
		Eu[u].resize(Nv);
		Ev[u].resize(Nv);
	}
	//FDTDradiator->run();
	//Get Surface Field
	//FDTDradiator->GetProFieldE(Eu, Ev, Nu, Nv);	//�Ѻ����ӿڸĳ���double ��������Ӧ

	//���ǰ���Ѽ����꣬��ֱ�Ӷ�ȡ�ļ�	���->run ��->getProFieldE;	--���� 20180325
	FDTDradiator->LoadProFieldE("./PropagatedEField.dat",Eu,Ev,Nu,Nv);

	GraphTrans graphTrans;
	graphTrans.setGraphTransPar(AperturePosition.X(), AperturePosition.Y(), AperturePosition.Z(),
		0, 1, 0, 15.06);
	double ds = Lu / (Nu - 1); 
	field->setNM(Nu, Nu);
	field->setPlane(graphTrans, ds); 
	field->setField(Eu, Ev);
	field->setShowPara(1, 1, 0);
	//field->updateData();

}

void CalculateFDTDThread::setMainValue(int val, void *user)
{
	((CalculateFDTDThread*)user)->sendMainValue(val);
}

void CalculateFDTDThread::setSlaverValue(float val, void *user)
{
	((CalculateFDTDThread*)user)->sendSlaverValue(val);
}

void CalculateFDTDThread::killFDTD()
{
	// free FDTD
	deleteLater();
}
