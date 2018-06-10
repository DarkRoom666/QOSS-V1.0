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
#include <vtkAxesActor.h>
#include <vector>
#include <../util/Vector3.h>
#include <vtkjsoncpp/json/json.h>

#include "BasicParameters.h"
#include "Restriction.h"

enum MirrorsType
{
	PLANEMIRROR = 0,
	QUADRICSURFACE,
	PARABOLICCYLINDER,
	PARABOLOID,
	ELLIPSOID,
	STLMIRROR,
	PHSCORMIRROR
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

	virtual void calcRestriction();

	void switchIsTransparent();
	bool getIsTransparent() const { return isTransparent; }

	void switchIsShow();
	bool getIsShow() const { return isShow; }

	vtkSmartPointer<vtkAxesActor> getActorAxes() const;
	void calcActorAxes();

	// ���STL ��ʽ�ļ�
	void saveSTL();

	virtual Json::Value getDataJson(const string& dir, int index) const;

protected:

	MirrorsType type;

	//����ÿ��ģ�͵���ʾ���ʷ�����
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkPolyData> polyData;

	vtkSmartPointer<vtkAxesActor> actorAxes;

	vector<Restriction*> restrictions;

	bool isTransparent; // �Ƿ�͸��
	bool isShow; // �Ƿ���ʾ
	vtkSmartPointer<vtkProperty> property;




private:

};



#endif // MIRROR_H