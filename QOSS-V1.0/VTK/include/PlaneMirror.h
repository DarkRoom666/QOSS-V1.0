/*
*	created by liyun 2017/10/24
*   function ƽ�澵(Ĭ�ϵľ���)
*   version 1.0
*/

#ifndef PLANEMIRROR_H  
#define PLANEMIRROR_H

#include "Mirror.h"
#include <vector>

class PlaneMirror: public Mirror
{
public:
	// ����Ĭ�ϵ�ƽ�澵 û�д������
	PlaneMirror(const GraphTrans & _graphTrans);
	// ����ƽ�澵 �������
	PlaneMirror(const GraphTrans & _graphTrans, const std::vector<double> parameter);
	virtual ~PlaneMirror();

	virtual void calPolyData(double ds = 0);

	virtual void updateData();

	// ���������߼��㷴����ߺͽ���
	virtual void calcReflectedRay(const Vector3&, Vector3&, Vector3&);

	virtual QTreeWidgetItem* getTree();

private:


};



#endif // PLANEMIRROR_H