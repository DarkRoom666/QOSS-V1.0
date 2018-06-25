//Written by Ming Jin, 2018
//version 0.0
//Physical Optics (Zero Order Current Description on STL Triangles) by CUDA Computation
//Requires dense mesh!  1/6~1/8 lambda

//����Դ���� ��ת���ɵ���
//�������

#pragma once

#ifndef CUDAPHYSICALOPTICS_H
#define CUDAPHYSICALOPTICS_H

#include <string>

class _declspec(dllexport) CUDAPhysicalOptics
{
public:
	CUDAPhysicalOptics();
	~CUDAPhysicalOptics();

	void setInputFile(const std::string & file);
	void setModelFile(const std::string & file);
	void setMode(int _mode);
	

	void calculate(double fre, double dis);
	void SetReturnFloat(void(*returnFloat)(float, void*), void*_user);// ע��ص�����
	
private:


	int Mode;	//0 Huygens Box(J,H/JM,E) to AnyAperture(E)
				//1 Huygens Box(J,H/JM,E) to Mirror(J,H) Then to AnyAperture(E)
				//2 Aperture(JM,E) to Mirror(J,H) Then to AnyAperture(E)
				//3 Mirrors(J,H) to Aperture(E)

	void(*returnFloat)(float, void*);
	void *user; // �ص���������ָ��

protected:
	std::string FieldName;
	std::string STLName;
};


#endif // CUDAPHYSICALOPTICS