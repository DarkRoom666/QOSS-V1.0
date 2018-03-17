#include "../include/CalculateFDTDThread.h"
#include <../Calculation/FDTDRadiator.h>

using namespace userInterface;
void CalculateFDTDThread::run()
{
	FDTDradiator->SetReturnInt(CalculateFDTDThread::setMainValue, this);	//����״̬��ѯ-�����Ļص�����
	FDTDradiator->SetReturnFloat(CalculateFDTDThread::setSlaverValue, this);	//����״̬����-�������Ļص�����

	FDTDradiator->SetUpLowOrderVlasovRadiator(0, 1, 44e9, 16e-3, 20e-3, 80, 20);

	//SourceCenter and Aperture Parameters;
	//���ó������ -Ĭ��(��λm)
	Position3D AperturePosition;	AperturePosition.setX(0.04);			AperturePosition.setY(0);	AperturePosition.setZ(0.35 - 0.0205);
	Vector3D UDirection;			UDirection.setX(cos(15 * 3.1415926 / 180));	UDirection.setY(0);			UDirection.setZ(-sin(15 * 3.1415926 / 180));	UDirection.normalize();
	Vector3D VDirection;			VDirection.setX(0);						VDirection.setY(1.0);		VDirection.setZ(0.0);					VDirection.normalize();
	Vector3D ApertureDirection;		ApertureDirection = ApertureDirection.crossProduct(UDirection, VDirection);
	int Nu;	int Nv;	Nu = 121; Nv = 121;
	double Lu;	Lu = 2.998e8 / 44e9 * 60;	//100 lambda
	double Lv;	Lv = 2.998e8 / 44e9 * 60;	//100 lambda

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
	FDTDradiator->run();

	double ds = 0.001; // ��
	field->setNM(Nu, Nu);
	field->setPlane(GraphTrans(), ds); // ��
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
