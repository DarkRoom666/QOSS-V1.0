/*
*	created by liyun 2017/11/22
*   function ���澵��2Dλ���Լ���������3Dλ�� �������Ƶĳߴ��С
*   version 1.0
*/

#ifndef MIRRORPOSITION_H
#define MIRRORPOSITION_H

#include <vector>
#include "../util/Vector3.h"
#include "../util/GraphTrans.h"
namespace calculation
{
	class MirrorPosition
	{
	public:
		MirrorPosition();
		~MirrorPosition();

		void setMirrorNum(int num);
		int getMirrorNum() const { return mirrorNum; }

		void setLightDirection(int num, const Vector3& light);
		Vector3 getLightDirection(int num) const;

		void setLightLength(int num, double light);
		double getLightLength(int num) const;

		void setStartingPoint(const Vector3& light);
		Vector3 getStartingPoint() const;

		// ֻ���ò�����
		void setLightPhi(int num, double phi);
		// ͨ���޸ļн� Ӱ����һ������
		void calcLightPhi(int num, double phi);

		void getBoundaryByDefault(double &_length, double &_width);


		// �������ɾ��ӵĳ�ʼλ��
		void getInitialPosition(vector<GraphTrans> &position);

	private:
		void calcLightPhi(int num);

		std::vector<Vector3> lightDirection;
		std::vector<double> lightLength;
		std::vector<double> lightPhi;

		Vector3 org;
		int mirrorNum;
		double scale;
	};
}


#endif // MIRRORPOSITION_H
