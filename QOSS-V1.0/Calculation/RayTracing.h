/*
*	created by liyun 2017/11/27
*   function �������׷�� ���뾵�ӵ�ָ�� ���������ߺͽ����Լ��Ƿ��ཻ
*   version 1.0
*/

#ifndef RAYTRACING_H  
#define RAYTRACING_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include "../util/Vector3.h"
#include "../VTK/include/Mirror.h"

#include "../Calculation/Matrix4D.h"
#include "../Calculation/Vector3D.h"

namespace calculation
{
	class RayTracing
	{
	public:
		RayTracing(Mirror* _mirror);

		~RayTracing();

		void setMirror(Mirror* mirror);

		// ���������ߺͽ����Լ��Ƿ��ཻ
		void calcReflect(const Vector3& startPiont,
			const Vector3& direction,
			Vector3 &reflex, Vector3 &intersection,
			bool &isIntersect);

		void calcReflectBatch(const vector<vector<Vector3>>& startPiont,
			const vector<vector<Vector3>>& direction,
			vector<vector<Vector3>> &reflex, vector<vector<Vector3>> &intersection,
			vector<vector<bool>> &isIntersect);

		// ����������ߺͷ����������
		static Vector3 reflectLight(const Vector3& a, const Vector3& n);

		// ����ǽ�����Ҫ ��ǰ����ñ任����
		void calcMatrix();
	private:

		// ����ģ���ʷ����ݼ��㷴��
		void calcReflectByPolyDataBatch(const vector<vector<Vector3>>& startPiont,
			const vector<vector<Vector3>>& direction,
			vector<vector<Vector3>> &reflex, vector<vector<Vector3>> &intersection,
			vector<vector<bool>> &isIntersect);

		void calcReflectByPolyData(const Vector3& startPiont,
			const Vector3& direction, Vector3 &reflex,
			Vector3 &intersection,
			bool &isIntersect);

		void calcReflectByQuadricSurface(const Vector3& startPiont,
			const Vector3& direction, Vector3 &reflex,
			Vector3 &intersection,
			bool &isIntersect);

		//calcReflectByPolyData

		// ��������ֱ���ཻ�ж�
		bool isIntersect(const Vector3 &orig, const Vector3 &dir,
			const Vector3 &v0, const Vector3 &v1, const Vector3 &v2,
			Vector3 &intersection, double &t);

		bool ray_CurvedSurface(const vector<double> & a, Vector3 n, Vector3 org,
			double &t, Vector3 &interPoint);

		Mirror* mirror;

		bool isCalcMatrix;

		// ��������ϵת��ģ�͵��������ϵ�������������ת��ƽ��
		Matrix4D R_rotatMatrix;
		Matrix4D R_translateMatrix;

		// ģ�͵��������ϵת�������������
		Matrix4D rotatMatrix; 
		Matrix4D translateMatrix; 

	};

}

#endif // RAYTRACING_H

