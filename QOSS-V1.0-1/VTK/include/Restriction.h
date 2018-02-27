/*
*	created by liyun 2017/12/25  Merry Christmas  
*   function ���ӵ��������� ��Բ��������������ѡ��
*   version 1.0
*/

#ifndef RESTRICTION_H
#define RESTRICTION_H

#include "BasicParameters.h"
#include <vector>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <QTreeWidgetItem>
#include <QString>

class Restriction : public BasicParameters
{
public:
	Restriction(); // ���̳���
	Restriction(const Restriction&); // ���̳���

	// ����Ĭ�ϵ�ƽ�澵 û�д������
	Restriction(const GraphTrans & _graphTrans);
	// ����ƽ�澵 �������
	Restriction(const GraphTrans & _graphTrans, const std::vector<double> parameter);
	virtual ~Restriction();

	void calPolyData(double ds = 0);
	vtkSmartPointer<vtkPolyData> getPolyData() const;

	virtual void updateData();

	// ���������߼��㷴����ߺͽ���
	virtual void calcReflectedRay(const Vector3&, Vector3&, Vector3&) {};

	QTreeWidgetItem* getTree();

	vtkSmartPointer<vtkActor> getActor() const;
	void calActor();

	void setDataByNum(int, double);

private:

	//����ÿ��ģ�͵���ʾ���ʷ�����
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkPolyData> polyData;

	bool isTransparent; // �Ƿ���ʾ
	vtkSmartPointer<vtkProperty> property;

};

#endif // RESTRICTION_H