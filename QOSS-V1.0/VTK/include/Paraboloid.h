/*
*	created by liyun 2017/12/6
*   function �������� ������
*   x^2 + y^2 = 4Fx  // fΪ����
*   R Ϊx,y�����ֵ
*   version 1.0
*/

#ifndef PARABOLOID_H
#define PARABOLOID_H

#include "QuadricSurfaceMirror.h"
class Paraboloid : public  QuadricSurfaceMirror
{
public:
	// ����Ĭ�ϵ�ƽ�澵 û�д������
	Paraboloid(const GraphTrans & _graphTrans);
	// ����ƽ�澵 �������
	Paraboloid(const GraphTrans & _graphTrans, const std::vector<double>& parameter);

	void setParameter(double focus, double radius);

	void setFocus(double);
	double getFocus() const { return focus; }

	void setRadius(double);
	double getRadius() const { return radius; }

	virtual QTreeWidgetItem* getTree();
private:
	double focus;
	double radius;
};

#endif // PARABOLOID_H