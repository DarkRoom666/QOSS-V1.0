/*
*	created by liyun 2017/12/7
*   function �������� ������
*   x^2 / a^2 + y^2 / b^2 + z^2 / c^2 = 1 
*   R Ϊx,y�����ֵ
*   version 1.0
*/
#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "QuadricSurfaceMirror.h"
class Ellipsoid : public  QuadricSurfaceMirror
{
public:
	// ����Ĭ�ϵ�ƽ�澵 û�д������
	Ellipsoid(const GraphTrans & _graphTrans);
	// ����ƽ�澵 �������
	Ellipsoid(const GraphTrans & _graphTrans, const std::vector<double>& parameter);

	void setParameter(double a, double b, double c, double theta);

	virtual QTreeWidgetItem* getTree();
private:
	double a;
	double b;
	double c;
	double theta;
};

#endif // ELLIPSOID_H