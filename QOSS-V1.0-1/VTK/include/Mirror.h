/*
*	created by liyun 2017/10/23
*   function ���ӵ������ 
*   version 1.0
*/

#ifndef MIRROR_H  
#define MIRROR_H

#include "util/GraphTrans.h"
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <QTreeWidgetItem>
#include <QString>

#include <vector>
#include <../util/Vector3.h>

#include "BasicParameters.h"
#include "Restriction.h"

enum MirrorsType
{
	PLANEMIRROR = 0,
	QUADRICSURFACE,
	PARABOLICCYLINDER,
	PARABOLOID,
	ELLIPSOID,
};

class actor;
class Mirror : public BasicParameters
{
public:
	Mirror();
	virtual ~Mirror();

	MirrorsType getMirrorsType() const;

	// cal(����) and get actor 
	vtkSmartPointer<vtkActor> getActor() const;
	virtual void calActor();

	// cal(����) and get polyData 
	vtkSmartPointer<vtkPolyData> getPolyData() const;
	virtual void calPolyData(double ds = 0) = 0;

	virtual void updateData() = 0;

	// ���������߼��㽻��ͷ���
	virtual void calcReflectedRay(const Vector3&, Vector3&, Vector3&) = 0;

	void setSelected(bool);

	virtual QTreeWidgetItem* getTree() { return nullptr; };

	void addRestriction(Restriction*);
	void setRestriction(int,Restriction*);
	void removeRestriction(int);

	// ������������������
	void removeRestrictionAll();

	// ����������������
	void clearRestrictionAll();

	Restriction* getRestriction(int) const;
	// �ƶ�mirror�е�Restriction���ƶ���ԭmirror��RestrictionΪ��
	void moveRestriction(Mirror*);
	const vector<Restriction*>& getRestrictionAll() const { return restrictions; }

	void switchIsTransparent();
	bool getIsTransparent() const { return isTransparent; }

	void switchIsShow();
	bool getIsShow() const { return isShow; }

protected:

	

	MirrorsType type;

	//����ÿ��ģ�͵���ʾ���ʷ�����
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkPolyData> polyData;

	bool isTransparent; // �Ƿ�͸��
	bool isShow; // �Ƿ���ʾ
	vtkSmartPointer<vtkProperty> property;

	vector<Restriction*> restrictions;
private:

};



#endif // MIRROR_H