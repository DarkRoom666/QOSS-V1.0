#include "CUDAPhysicalOptics.h"
#include "CUDAPOCalculation.h"
#include "GraphTrans.h"
#include "Field.h"
#include "STLMirror.h"
#include "Position3D.h"//λ��
#include "Matrix4D.h"	//ת�ƾ���
#include "Vector3D.h"	//����
#include "Raytracing.h"
#include "../QOSS-V1.0/util/Constant_Var.h"
#include "GEMS_Memory.h"
#include "HuygensBoxData.h"
#include "Normal2GraphTrans.h"
#include <vector>

using namespace calculation;
using namespace std;


CUDAPhysicalOptics::CUDAPhysicalOptics()
{
	this->returnFloat = NULL;
	this->user = NULL;
	//FieldName = NULL;
	//STLName = NULL;
	//this->FieldName.resize(100);
	//this->STLName.resize(100);
}

CUDAPhysicalOptics::~CUDAPhysicalOptics()
{
	//if (FieldName) delete[] FieldName;
	//if (STLName) delete[] STLName;
}

void CUDAPhysicalOptics::setInputFile(const std::string & file)
{	
	
	InFieldName = file;
	InFieldfile = true;
	/*
	int size = file.size();
	FieldName = new char[size];
	for (int i = 0; i < size; i++) {
		FieldName[i] = file.c_str()[i];
	}*/
}

void CUDAPhysicalOptics::setOutputFile(const std::string & file)
{

	InFieldName = file;
	OutFieldfile = true;
	/*
	int size = file.size();
	FieldName = new char[size];
	for (int i = 0; i < size; i++) {
	FieldName[i] = file.c_str()[i];
	}*/
}

void CUDAPhysicalOptics::setModelFile(const std::string & file)
{
	STLName = file;
	STLfile = true;
	/*
	int size = file.size();
	STLName = new char[size];
	for (int i = 0; i < size; i++) {
		STLName[i] = file.c_str()[i];
	}*/

}

void CUDAPhysicalOptics::setHuygensFile(const std::string & file) {
	HuygensName = file;
	Huygensfile = true;
}

void CUDAPhysicalOptics::setModelPtr(void *_polydata) {
	polydataptr = _polydata;
	STLfile = false;
}

void CUDAPhysicalOptics::setHuygensPtr(void *_huygens) {
	huygensptr = _huygens;
	Huygensfile = false;
}

void CUDAPhysicalOptics::setInField(void *_infield) {
	infieldptr = _infield;
	InFieldfile = false;
}

void CUDAPhysicalOptics::setOutField(void *_outfield) {
	outfieldptr = _outfield;
	OutFieldfile = false;
}

int CUDAPhysicalOptics::calculateHuygensPro(int _index) {
	//�������������ϵ �����ǻݸ�˹���� ����ǿ���
	HuygensBox* HBox;
	if (Huygensfile) {
		HBox = new HuygensBox;
		HBox->ReadHuygensFile(HuygensName, _index);
	}
	else {
		HBox = (HuygensBox*) huygensptr;
	}
	Field* outField;
	if (OutFieldfile) {
		outField = new Field;
		outField->setFileAddress(OutFieldName);
		outField->readData();
	}
	else {
		outField = (Field *)outfieldptr;
	}
	//Դ
	vector<complex<float>> Jmx;
	vector<complex<float>> Jmy;
	vector<complex<float>> Jmz;
	vector<complex<float>> Jx;
	vector<complex<float>> Jy;
	vector<complex<float>> Jz;
	vector<float> pxi;
	vector<float> pyi;
	vector<float> pzi;
	vector<float> dss;
	
	HBox->GetCurrentList(Jx,Jy,Jz,Jmx,Jmy,Jmz);
	HBox->GetPosList(pxi, pyi, pzi, dss);

	CUDAPOCalculation cudapo;
	//���ûݸ�˹Դ
	cudapo.setFrequency(HBox->GetFreq());
	cudapo.setHuygensCurrentInput(Jx, Jy, Jz, Jmx, Jmy, Jmz);
	cudapo.setHuygensPosInput(pxi, pyi, pzi, dss);

	GraphTrans gt_out;//�������λ�úͷ�����Ϣ
	int M, N;
	double ds;//��ɢ���
	Vector3 u_out;	u_out.set(1.0, 0.0, 0.0);
	Vector3 v_out;	v_out.set(0.0, 1.0, 0.0);
	Vector3 n_out;	n_out.set(0.0, 0.0, 1.0);
	Vector3 p_cen_out;
	Vector3 p_center;
	outField->getSourcePara(gt_out, M, N, ds);
	//��������λ��
	p_cen_out.set(gt_out.getTrans_x(), gt_out.getTrans_y(), gt_out.getTrans_z());

	//������ת��
	Vector3D  RotateAxis_out(gt_out.getRotate_x(), gt_out.getRotate_y(), gt_out.getRotate_z());
	Matrix4D  RotateMatrix_out = Matrix4D::getRotateMatrix(gt_out.getRotate_theta(), RotateAxis_out);

	u_out = RotateMatrix_out*u_out;
	v_out = RotateMatrix_out*v_out;
	n_out = RotateMatrix_out*n_out;

	//���ó��䳡λ��
	cudapo.setOutputAperture(u_out, v_out, p_cen_out, float(ds), M, N);
	//��������������
	cudapo.calculateHuygens2E();
	
	//�ļ�������
	complex<double>** Eu_out = NULL;	Eu_out = Allocate_2D(Eu_out, M, N);
	complex<double>** Ev_out = NULL;	Ev_out = Allocate_2D(Ev_out, M, N);
	cudapo.getOutApertureE(Eu_out, Ev_out);
	//�ļ������
	ofstream outfile("./outPO_H.txt");
	outfile << gt_out.getTrans_x() << " "
		<< gt_out.getTrans_y() << " "
		<< gt_out.getTrans_z() << " "
		<< gt_out.getRotate_x() << " "
		<< gt_out.getRotate_y() << " "
		<< gt_out.getRotate_z() << " "
		<< gt_out.getRotate_theta() << " "
		<< M << " " << N << " " << ds << endl;
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
		{
			outfile
				<< abs(Eu_out[i][j]) << " " << arg(Eu_out[i][j]) * 180 / Pi << " "
				<< abs(Ev_out[i][j]) << " " << arg(Ev_out[i][j]) * 180 / Pi << " " << endl;
		}
	//д�ļ�����������
	if (returnFloat) // ���û��ע���򲻻ص�
	{
		returnFloat(100, user);
	}

	Free_2D(Eu_out);
	Free_2D(Ev_out);

	if(Huygensfile) delete HBox;
	if(OutFieldfile) delete outField;
	return 0;
}

int CUDAPhysicalOptics::calculateHuygensSTLPro(int _index) {
	//�������������ϵ �����ǻݸ�˹���� ����ǿ��� �м仹�и�������
	//����
	HuygensBox* HBox;
	if (Huygensfile) {
		HBox = new HuygensBox;	//˭new ˭delete
		HBox->ReadHuygensFile(HuygensName, _index);
	}
	else {
		HBox = (HuygensBox*)huygensptr;
	}
	//��ȡ��������-�ļ�
	STLMirror stlMirror;
	vtkPolyData* poly;
	if (STLfile) {
		stlMirror.setNameFile(STLName);
		stlMirror.readData();
	}
	else {
		//ָ�봫��STL
		 poly = (vtkPolyData*)polydataptr;
	}
	//���
	Field* outField;
	if (OutFieldfile) {
		outField = new Field;
		outField->setFileAddress(OutFieldName);
		outField->readData();
	}
	else {
		outField = (Field *)outfieldptr;
	}
	//Դ
	vector<complex<float>> Jmx;	vector<complex<float>> Jmy;	vector<complex<float>> Jmz;
	vector<complex<float>> Jx;	vector<complex<float>> Jy;	vector<complex<float>> Jz;
	vector<float> pxi;			vector<float> pyi;			vector<float> pzi;
	vector<float> dss;

	HBox->GetCurrentList(Jx, Jy, Jz, Jmx, Jmy, Jmz);
	HBox->GetPosList(pxi, pyi, pzi, dss);

	CUDAPOCalculation cudapo;

	//���ûݸ�˹Դ
	cudapo.setFrequency(HBox->GetFreq());
	cudapo.setHuygensCurrentInput(Jx, Jy, Jz, Jmx, Jmy, Jmz);
	cudapo.setHuygensPosInput(pxi, pyi, pzi, dss);
	//���÷�����
	if (STLfile) cudapo.setReflectSTL(stlMirror.getPolyData());
	else cudapo.setReflectSTL(polydataptr);
	//��һ������
	cudapo.calculateHuygens2H();

	//��ȡ����ı������
	vector<complex<double>> Hx_Current;	vector<complex<double>> Hy_Current;	vector<complex<double>> Hz_Current;

	int NumCurrent;
	if(STLfile) NumCurrent = stlMirror.getPolyData()->GetNumberOfCells();
	else NumCurrent = poly->GetNumberOfCells();
	
	Hx_Current.resize(NumCurrent);	Hy_Current.resize(NumCurrent);	Hz_Current.resize(NumCurrent);
	//��ȡ��һ�������� �������
	cudapo.getSTLlistHfield(Hx_Current, Hy_Current, Hz_Current);
	cudapo.cleaninput();
	cudapo.cleanoutput();

	//���õڶ�����������-�������
	if (STLfile)	cudapo.setSTLCurrentSourceZeroOrder(stlMirror.getPolyData(), Hx_Current, Hy_Current, Hz_Current);
	else cudapo.setSTLCurrentSourceZeroOrder(poly, Hx_Current, Hy_Current, Hz_Current);

	//���õڶ����������-���򳡷ֲ�
	GraphTrans gt_out;//�������λ�úͷ�����Ϣ
	int M, N;
	double ds;//��ɢ���
	Vector3 u_out;	u_out.set(1.0, 0.0, 0.0);
	Vector3 v_out;	v_out.set(0.0, 1.0, 0.0);
	Vector3 n_out;	n_out.set(0.0, 0.0, 1.0);
	Vector3 p_cen_out;
	Vector3 p_center;
	outField->getSourcePara(gt_out, M, N, ds);
	//��������λ��
	p_cen_out.set(gt_out.getTrans_x(), gt_out.getTrans_y(), gt_out.getTrans_z());

	//������ת��
	Vector3D  RotateAxis_out(gt_out.getRotate_x(), gt_out.getRotate_y(), gt_out.getRotate_z());
	Matrix4D  RotateMatrix_out = Matrix4D::getRotateMatrix(gt_out.getRotate_theta(), RotateAxis_out);

	u_out = RotateMatrix_out*u_out;
	v_out = RotateMatrix_out*v_out;
	n_out = RotateMatrix_out*n_out;

	//���ó��䳡λ��
	cudapo.setOutputAperture(u_out, v_out, p_cen_out, float(ds), M, N);
	//�ڶ�������
	cudapo.calculateS2F();

	//�ļ�������
	complex<double>** Eu_out = NULL;	Eu_out = Allocate_2D(Eu_out, M, N);
	complex<double>** Ev_out = NULL;	Ev_out = Allocate_2D(Ev_out, M, N);
	cudapo.getOutApertureE(Eu_out,Ev_out);
	//�ļ������
	ofstream outfile("./outPO_HSTL.txt");
	outfile << gt_out.getTrans_x() << " "
		<< gt_out.getTrans_y() << " "
		<< gt_out.getTrans_z() << " "
		<< gt_out.getRotate_x() << " "
		<< gt_out.getRotate_y() << " "
		<< gt_out.getRotate_z() << " "
		<< gt_out.getRotate_theta() << " "
		<< M << " " << N << " " << ds << endl;

	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
		{
			outfile
				<< abs(Eu_out[i][j]) << " " << arg(Eu_out[i][j]) * 180 / Pi << " "
				<< abs(Ev_out[i][j]) << " " << arg(Ev_out[i][j]) * 180 / Pi << " " << endl;
		}


	if(Huygensfile) delete HBox;
	if (OutFieldfile) delete outField;
	return 0;
}

void CUDAPhysicalOptics::calculate(double fre, double dis)
{	
	//��ȡ���볡�ֲ�-�ļ�
	Field* inputField;
	if(InFieldfile){
		inputField = new Field;
		inputField->setFileAddress(InFieldName);
		inputField->readData();
	}
	else {
		inputField = (Field*) infieldptr;
	}

	//��ȡ��������-�ļ�
	STLMirror stlMirror;
	vtkPolyData* poly;
	if (STLfile) {	
		stlMirror.setNameFile(STLName);
		stlMirror.readData();
	}
	else {
		poly = (vtkPolyData*) polydataptr;
	}

	GraphTrans gt_inc;//�������λ�úͷ�����Ϣ
	int M, N;
	double ds;//��ɢ���
	Vector3 u_input;	u_input.set(1.0, 0.0, 0.0);
	Vector3 v_input;	v_input.set(0.0, 1.0, 0.0);
	Vector3 n_input;	n_input.set(0.0, 0.0, 1.0);
	Vector3 p_cen_in;
	Vector3 p_center;
	inputField->getSourcePara(gt_inc, M, N, ds);
	//��������λ��
	p_cen_in.set(gt_inc.getTrans_x(), gt_inc.getTrans_y(), gt_inc.getTrans_z());

	//������ת��
	Vector3D  RotateAxis_in(gt_inc.getRotate_x(), gt_inc.getRotate_y(), gt_inc.getRotate_z());
	Matrix4D  RotateMatrix_in = Matrix4D::getRotateMatrix(gt_inc.getRotate_theta(), RotateAxis_in);

	u_input = RotateMatrix_in*u_input;
	v_input = RotateMatrix_in*v_input;
	n_input = RotateMatrix_in*n_input;

	vector<vector<complex<double>>> Eu;
	vector<vector<complex<double>>> Ev;

	Eu = inputField->getEx();	Ev = inputField->getEy();
	
	
	CUDAPOCalculation cudapo;
	// ���õ�λ
	cudapo.SetReturnFloat(returnFloat, user);
	// ����Ƶ��
	cudapo.setFrequency(fre);
	//�������䳡
	cudapo.setPlaneApertureEField_D(Eu,Ev,u_input,v_input,p_cen_in,float(ds),M,N);
	//���÷�������
	if(STLfile) cudapo.setReflectSTL(stlMirror.getPolyData());
	else cudapo.setReflectSTL(poly);

	if (returnFloat) // ���û��ע���򲻻ص�
	{
		returnFloat(30, user);
	}

	//���㳡������
	cudapo.calculateF2S();
	//��ȡ����ı������
	vector<complex<double>> Hx_Current;
	vector<complex<double>> Hy_Current;
	vector<complex<double>> Hz_Current;

	int NumCurrent;
	if (STLfile) NumCurrent = stlMirror.getPolyData()->GetNumberOfCells();
	else NumCurrent = poly->GetNumberOfCells();

	Hx_Current.resize(NumCurrent);
	Hy_Current.resize(NumCurrent);
	Hz_Current.resize(NumCurrent);
	//��ȡ��һ�������� �������
	cudapo.getSTLlistHfield(Hx_Current,Hy_Current,Hz_Current);
	cudapo.cleaninput();
	cudapo.cleanoutput();

	//���õڶ�����������-�������
	if(STLfile)	cudapo.setSTLCurrentSourceZeroOrder(stlMirror.getPolyData(), Hx_Current, Hy_Current, Hz_Current);
	else cudapo.setSTLCurrentSourceZeroOrder(poly, Hx_Current, Hy_Current, Hz_Current);

	if (returnFloat) // ���û��ע���򲻻ص�
	{
		returnFloat(60, user);
	}

	//���³��䳡��λ��
	p_center = p_cen_in;

	Vector3 tempReflect;	//�������
	Vector3 InterPoint;		//����
	GraphTrans gt_out; gt_out = gt_inc;
	RayTracing RT(&stlMirror);
	bool tempIsIntersect = false;
	RT.calcReflect(p_center, n_input, tempReflect, InterPoint, tempIsIntersect);
	tempReflect.Normalization();
	double z0 = CalDistance(p_center,InterPoint);
	
	//����ƽ�� �������浽�����潻��
	p_center = p_center + n_input*z0;
	gt_out.updateTranslate(p_center);
	//������ת
	updateSource_n(tempReflect, gt_out);
	//���·����ƽ�� �ӷ����潻�㵽������
	p_center = p_center + tempReflect*dis;
	gt_out.updateTranslate(p_center);
	Vector3 u_out;	u_out.set(1.0, 0.0, 0.0);
	Vector3 v_out;	v_out.set(0.0, 1.0, 0.0);
	Vector3 n_out;	n_out.set(0.0, 0.0, 1.0);
	Vector3 p_cen_out;	p_cen_out = p_center;

	//������ת��
	Vector3D  RotateAxis_out(gt_out.getRotate_x(), gt_out.getRotate_y(), gt_out.getRotate_z());
	Matrix4D  RotateMatrix_out = Matrix4D::getRotateMatrix(gt_out.getRotate_theta(), RotateAxis_out);

	u_out = RotateMatrix_out*u_out;
	v_out = RotateMatrix_out*v_out;
	n_out = RotateMatrix_out*n_out;
	//���ó��䳡λ��
	cudapo.setOutputAperture(u_out,v_out,p_cen_out,float(ds),M,N);
	//CUDAPO����
	cudapo.calculateS2F();

	if (returnFloat) // ���û��ע���򲻻ص�
	{
		returnFloat(90, user);
	}
	//�ļ�������
	complex<double>** Eu_out = NULL;	Eu_out = Allocate_2D(Eu_out, M, N);
	complex<double>** Ev_out = NULL;	Ev_out = Allocate_2D(Ev_out, M, N);
	cudapo.getOutApertureE(Eu_out,Ev_out);
	//�ļ������
	ofstream outfile("./outPO.txt");
	outfile << gt_out.getTrans_x() << " "
		<< gt_out.getTrans_y() << " "
		<< gt_out.getTrans_z() << " "
		<< gt_out.getRotate_x() << " "
		<< gt_out.getRotate_y() << " "
		<< gt_out.getRotate_z() << " "
		<< gt_out.getRotate_theta() << " "
		<< M << " " << N << " " << ds << endl;
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
		{
			outfile
				<< abs(Eu_out[i][j]) << " " << arg(Eu_out[i][j]) * 180 / Pi << " "
				<< abs(Ev_out[i][j]) << " " << arg(Ev_out[i][j]) * 180 / Pi << " " << endl;
		}
	//д�ļ�����������
	if (returnFloat) // ���û��ע���򲻻ص�
	{
		returnFloat(100, user);
	}

	Free_2D(Eu_out);
	Free_2D(Ev_out);

	if(InFieldfile) delete inputField;
}

void CUDAPhysicalOptics::SetReturnFloat(void(*returnFloat)(float, void *), void * _user)
{
	this->returnFloat = returnFloat;
	this->user = _user;
}
