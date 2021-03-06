#include "CUDAPhysicalOptics.h"
#include "CUDAPOCalculation.h"
#include "GraphTrans.h"
#include "Field.h"
#include "STLMirror.h"
#include "Position3D.h"//位置
#include "Matrix4D.h"	//转移矩阵
#include "Vector3D.h"	//向量
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
	//设置输入输出关系 输入是惠更斯盒子 输出是口面
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
	//源
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
	//设置惠更斯源
	cudapo.setFrequency(HBox->GetFreq());
	cudapo.setHuygensCurrentInput(Jx, Jy, Jz, Jmx, Jmy, Jmz);
	cudapo.setHuygensPosInput(pxi, pyi, pzi, dss);

	GraphTrans gt_out;//出射面的位置和方向信息
	int M, N;
	double ds;//离散间隔
	Vector3 u_out;	u_out.set(1.0, 0.0, 0.0);
	Vector3 v_out;	v_out.set(0.0, 1.0, 0.0);
	Vector3 n_out;	n_out.set(0.0, 0.0, 1.0);
	Vector3 p_cen_out;
	Vector3 p_center;
	outField->getSourcePara(gt_out, M, N, ds);
	//口面中心位置
	p_cen_out.set(gt_out.getTrans_x(), gt_out.getTrans_y(), gt_out.getTrans_z());

	//设置旋转：
	Vector3D  RotateAxis_out(gt_out.getRotate_x(), gt_out.getRotate_y(), gt_out.getRotate_z());
	Matrix4D  RotateMatrix_out = Matrix4D::getRotateMatrix(gt_out.getRotate_theta(), RotateAxis_out);

	u_out = RotateMatrix_out*u_out;
	v_out = RotateMatrix_out*v_out;
	n_out = RotateMatrix_out*n_out;

	//设置出射场位置
	cudapo.setOutputAperture(u_out, v_out, p_cen_out, float(ds), M, N);
	//外推至近场口面
	cudapo.calculateHuygens2E();
	
	//文件输出输出
	complex<double>** Eu_out = NULL;	Eu_out = Allocate_2D(Eu_out, M, N);
	complex<double>** Ev_out = NULL;	Ev_out = Allocate_2D(Ev_out, M, N);
	cudapo.getOutApertureE(Eu_out, Ev_out);
	//文件输出：
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
	//写文件结束，结束
	if (returnFloat) // 如果没有注册则不回调
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
	//设置输入输出关系 输入是惠更斯盒子 输出是口面 中间还有个反射面
	//输入
	HuygensBox* HBox;
	if (Huygensfile) {
		HBox = new HuygensBox;	//谁new 谁delete
		HBox->ReadHuygensFile(HuygensName, _index);
	}
	else {
		HBox = (HuygensBox*)huygensptr;
	}
	//读取计算网格-文件
	STLMirror stlMirror;
	vtkPolyData* poly;
	if (STLfile) {
		stlMirror.setNameFile(STLName);
		stlMirror.readData();
	}
	else {
		//指针传递STL
		 poly = (vtkPolyData*)polydataptr;
	}
	//输出
	Field* outField;
	if (OutFieldfile) {
		outField = new Field;
		outField->setFileAddress(OutFieldName);
		outField->readData();
	}
	else {
		outField = (Field *)outfieldptr;
	}
	//源
	vector<complex<float>> Jmx;	vector<complex<float>> Jmy;	vector<complex<float>> Jmz;
	vector<complex<float>> Jx;	vector<complex<float>> Jy;	vector<complex<float>> Jz;
	vector<float> pxi;			vector<float> pyi;			vector<float> pzi;
	vector<float> dss;

	HBox->GetCurrentList(Jx, Jy, Jz, Jmx, Jmy, Jmz);
	HBox->GetPosList(pxi, pyi, pzi, dss);

	CUDAPOCalculation cudapo;

	//设置惠更斯源
	cudapo.setFrequency(HBox->GetFreq());
	cudapo.setHuygensCurrentInput(Jx, Jy, Jz, Jmx, Jmy, Jmz);
	cudapo.setHuygensPosInput(pxi, pyi, pzi, dss);
	//设置反射面
	if (STLfile) cudapo.setReflectSTL(stlMirror.getPolyData());
	else cudapo.setReflectSTL(polydataptr);
	//第一步计算
	cudapo.calculateHuygens2H();

	//提取计算的表面电流
	vector<complex<double>> Hx_Current;	vector<complex<double>> Hy_Current;	vector<complex<double>> Hz_Current;

	int NumCurrent;
	if(STLfile) NumCurrent = stlMirror.getPolyData()->GetNumberOfCells();
	else NumCurrent = poly->GetNumberOfCells();
	
	Hx_Current.resize(NumCurrent);	Hy_Current.resize(NumCurrent);	Hz_Current.resize(NumCurrent);
	//读取第一步计算结果 表面电流
	cudapo.getSTLlistHfield(Hx_Current, Hy_Current, Hz_Current);
	cudapo.cleaninput();
	cudapo.cleanoutput();

	//设置第二部计算输入-表面电流
	if (STLfile)	cudapo.setSTLCurrentSourceZeroOrder(stlMirror.getPolyData(), Hx_Current, Hy_Current, Hz_Current);
	else cudapo.setSTLCurrentSourceZeroOrder(poly, Hx_Current, Hy_Current, Hz_Current);

	//设置第二部计算输出-切向场分布
	GraphTrans gt_out;//出射面的位置和方向信息
	int M, N;
	double ds;//离散间隔
	Vector3 u_out;	u_out.set(1.0, 0.0, 0.0);
	Vector3 v_out;	v_out.set(0.0, 1.0, 0.0);
	Vector3 n_out;	n_out.set(0.0, 0.0, 1.0);
	Vector3 p_cen_out;
	Vector3 p_center;
	outField->getSourcePara(gt_out, M, N, ds);
	//口面中心位置
	p_cen_out.set(gt_out.getTrans_x(), gt_out.getTrans_y(), gt_out.getTrans_z());

	//设置旋转：
	Vector3D  RotateAxis_out(gt_out.getRotate_x(), gt_out.getRotate_y(), gt_out.getRotate_z());
	Matrix4D  RotateMatrix_out = Matrix4D::getRotateMatrix(gt_out.getRotate_theta(), RotateAxis_out);

	u_out = RotateMatrix_out*u_out;
	v_out = RotateMatrix_out*v_out;
	n_out = RotateMatrix_out*n_out;

	//设置出射场位置
	cudapo.setOutputAperture(u_out, v_out, p_cen_out, float(ds), M, N);
	//第二不计算
	cudapo.calculateS2F();

	//文件输出输出
	complex<double>** Eu_out = NULL;	Eu_out = Allocate_2D(Eu_out, M, N);
	complex<double>** Ev_out = NULL;	Ev_out = Allocate_2D(Ev_out, M, N);
	cudapo.getOutApertureE(Eu_out,Ev_out);
	//文件输出：
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
	//读取输入场分布-文件
	Field* inputField;
	if(InFieldfile){
		inputField = new Field;
		inputField->setFileAddress(InFieldName);
		inputField->readData();
	}
	else {
		inputField = (Field*) infieldptr;
	}

	//读取计算网格-文件
	STLMirror stlMirror;
	vtkPolyData* poly;
	if (STLfile) {	
		stlMirror.setNameFile(STLName);
		stlMirror.readData();
	}
	else {
		poly = (vtkPolyData*) polydataptr;
	}

	GraphTrans gt_inc;//入射面的位置和方向信息
	int M, N;
	double ds;//离散间隔
	Vector3 u_input;	u_input.set(1.0, 0.0, 0.0);
	Vector3 v_input;	v_input.set(0.0, 1.0, 0.0);
	Vector3 n_input;	n_input.set(0.0, 0.0, 1.0);
	Vector3 p_cen_in;
	Vector3 p_center;
	inputField->getSourcePara(gt_inc, M, N, ds);
	//口面中心位置
	p_cen_in.set(gt_inc.getTrans_x(), gt_inc.getTrans_y(), gt_inc.getTrans_z());

	//设置旋转：
	Vector3D  RotateAxis_in(gt_inc.getRotate_x(), gt_inc.getRotate_y(), gt_inc.getRotate_z());
	Matrix4D  RotateMatrix_in = Matrix4D::getRotateMatrix(gt_inc.getRotate_theta(), RotateAxis_in);

	u_input = RotateMatrix_in*u_input;
	v_input = RotateMatrix_in*v_input;
	n_input = RotateMatrix_in*n_input;

	vector<vector<complex<double>>> Eu;
	vector<vector<complex<double>>> Ev;

	Eu = inputField->getEx();	Ev = inputField->getEy();
	
	
	CUDAPOCalculation cudapo;
	// 设置单位
	cudapo.SetReturnFloat(returnFloat, user);
	// 设置频率
	cudapo.setFrequency(fre);
	//设置入射场
	cudapo.setPlaneApertureEField_D(Eu,Ev,u_input,v_input,p_cen_in,float(ds),M,N);
	//设置反射网格
	if(STLfile) cudapo.setReflectSTL(stlMirror.getPolyData());
	else cudapo.setReflectSTL(poly);

	if (returnFloat) // 如果没有注册则不回调
	{
		returnFloat(30, user);
	}

	//计算场到电流
	cudapo.calculateF2S();
	//提取计算的表面电流
	vector<complex<double>> Hx_Current;
	vector<complex<double>> Hy_Current;
	vector<complex<double>> Hz_Current;

	int NumCurrent;
	if (STLfile) NumCurrent = stlMirror.getPolyData()->GetNumberOfCells();
	else NumCurrent = poly->GetNumberOfCells();

	Hx_Current.resize(NumCurrent);
	Hy_Current.resize(NumCurrent);
	Hz_Current.resize(NumCurrent);
	//读取第一步计算结果 表面电流
	cudapo.getSTLlistHfield(Hx_Current,Hy_Current,Hz_Current);
	cudapo.cleaninput();
	cudapo.cleanoutput();

	//设置第二部计算输入-表面电流
	if(STLfile)	cudapo.setSTLCurrentSourceZeroOrder(stlMirror.getPolyData(), Hx_Current, Hy_Current, Hz_Current);
	else cudapo.setSTLCurrentSourceZeroOrder(poly, Hx_Current, Hy_Current, Hz_Current);

	if (returnFloat) // 如果没有注册则不回调
	{
		returnFloat(60, user);
	}

	//更新出射场的位置
	p_center = p_cen_in;

	Vector3 tempReflect;	//反射光线
	Vector3 InterPoint;		//交点
	GraphTrans gt_out; gt_out = gt_inc;
	RayTracing RT(&stlMirror);
	bool tempIsIntersect = false;
	RT.calcReflect(p_center, n_input, tempReflect, InterPoint, tempIsIntersect);
	tempReflect.Normalization();
	double z0 = CalDistance(p_center,InterPoint);
	
	//更新平移 从入射面到反射面交点
	p_center = p_center + n_input*z0;
	gt_out.updateTranslate(p_center);
	//更新旋转
	updateSource_n(tempReflect, gt_out);
	//更新反射后平移 从反射面交点到出射面
	p_center = p_center + tempReflect*dis;
	gt_out.updateTranslate(p_center);
	Vector3 u_out;	u_out.set(1.0, 0.0, 0.0);
	Vector3 v_out;	v_out.set(0.0, 1.0, 0.0);
	Vector3 n_out;	n_out.set(0.0, 0.0, 1.0);
	Vector3 p_cen_out;	p_cen_out = p_center;

	//设置旋转：
	Vector3D  RotateAxis_out(gt_out.getRotate_x(), gt_out.getRotate_y(), gt_out.getRotate_z());
	Matrix4D  RotateMatrix_out = Matrix4D::getRotateMatrix(gt_out.getRotate_theta(), RotateAxis_out);

	u_out = RotateMatrix_out*u_out;
	v_out = RotateMatrix_out*v_out;
	n_out = RotateMatrix_out*n_out;
	//设置出射场位置
	cudapo.setOutputAperture(u_out,v_out,p_cen_out,float(ds),M,N);
	//CUDAPO计算
	cudapo.calculateS2F();

	if (returnFloat) // 如果没有注册则不回调
	{
		returnFloat(90, user);
	}
	//文件输出输出
	complex<double>** Eu_out = NULL;	Eu_out = Allocate_2D(Eu_out, M, N);
	complex<double>** Ev_out = NULL;	Ev_out = Allocate_2D(Ev_out, M, N);
	cudapo.getOutApertureE(Eu_out,Ev_out);
	//文件输出：
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
	//写文件结束，结束
	if (returnFloat) // 如果没有注册则不回调
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
