/*
*	created by liyun 2017/12/4
*   function �������� ������
*   y^2 = 4Fx  // fΪ����
*   R Ϊy�����ֵ
*   version 1.0
*/

#ifndef PARABOLICCYLINDER_H
#define PARABOLICCYLINDER_H

#include "QuadricSurfaceMirror.h"
class ParabolicCylinder : public  QuadricSurfaceMirror
{
public:
	// ����Ĭ�ϵ�ƽ�澵 û�д������
	ParabolicCylinder(const GraphTrans & _graphTrans);
	// ����ƽ�澵 �������
	ParabolicCylinder(const GraphTrans & _graphTrans, const std::vector<double>& parameter);

	void setParameter(double focus, double yMax, double zMin, double zMax);
private:
	double focus;
	double yMax;
	double zMin;
	double zMax;
};

#endif // PARABOLICCYLINDER_H