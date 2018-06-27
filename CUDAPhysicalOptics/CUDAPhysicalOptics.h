//Written by Ming Jin, 2018
//version 1.0
//Physical Optics (Zero Order Current Description on STL Triangles) by CUDA Computation
//Requires dense mesh!  1/6~1/8 lambda

//输入源：场 会转化成电流
//输出：场

#pragma once

#ifndef CUDAPHYSICALOPTICS_H
#define CUDAPHYSICALOPTICS_H

#include <string>

class _declspec(dllexport) CUDAPhysicalOptics
{
public:
	CUDAPhysicalOptics();
	~CUDAPhysicalOptics();
	//文件接口
	void setInputFile(const std::string & file);
	void setOutputFile(const std::string & file);
	void setModelFile(const std::string & file);
	void setHuygensFile(const std::string & file);
	//指针接口
	void setModelPtr(void *_polydata);
	void setHuygensPtr(void * _huygens);
	void setInField(void *_infield);
	void setOutField(void * _outfield);
	
	int calculateHuygensPro(int _index);//惠更斯外推
	int calculateHuygensSTLPro(int _index);//惠更斯经镜面反射后的场
	void calculate(double fre, double dis);//同PVVA-STL
	void SetReturnFloat(void(*returnFloat)(float, void*), void*_user);// 注册回调函数
	
private:


	int Mode;	//0 Huygens Box(J,H/JM,E) to AnyAperture(E)
				//1 Huygens Box(J,H/JM,E) to Mirror(J,H) Then to AnyAperture(E)
				//2 Aperture(JM,E) to Mirror(J,H) Then to AnyAperture(E)
				//3 Mirrors(J,H) to Aperture(E)

	void(*returnFloat)(float, void*);
	void *user; // 回调函数的类指针
	void *huygensptr;
	void *polydataptr;
	void *infieldptr;
	void *outfieldptr;

protected:
	std::string InFieldName;
	std::string OutFieldName;
	std::string STLName;
	std::string HuygensName;
	bool STLfile;
	bool InFieldfile;
	bool OutFieldfile;
	bool Huygensfile;
};


#endif // CUDAPHYSICALOPTICS