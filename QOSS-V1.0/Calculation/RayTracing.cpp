#include "RayTracing.h"

#include "../util/Definition.h"
#include "../VTK/include/Restriction.h"


calculation::RayTracing::RayTracing(Mirror* _mirror)
	:mirror(_mirror)
{
	isCalcMatrix = false;
}

calculation::RayTracing::~RayTracing()
{
}

void calculation::RayTracing::setMirror(Mirror * mirror)
{
	this->mirror = mirror;
}

void calculation::RayTracing::calcNormalOfLine_Mirror(const Vector3 & startPiont,
	const Vector3 & direction, Vector3 & normal, Vector3 & intersection,
	bool & isIntersect, double & t)
{
	switch (mirror->getMirrorsType())
	{
	case PLANEMIRROR:
		calcNormalOfLine_MirrorByPolyData(startPiont, direction, normal,
			intersection, isIntersect, t);
		break;
	case QUADRICSURFACE:
	case PARABOLICCYLINDER:
	case PARABOLOID:
	case ELLIPSOID:
		calcNormalOfLine_MirrorByQuadricSurface(startPiont, direction, normal, 
			intersection, isIntersect, t);
		break;
	default:
		break;
	}
}

void calculation::RayTracing::calcReflectBatch(const vector<vector<Vector3>>& startPiont,
	const vector<vector<Vector3>>& direction,
	vector<vector<Vector3>> &reflex, vector<vector<Vector3>> &intersection,
	vector<vector<bool>> &isIntersect)
{
	switch (mirror->getMirrorsType())
	{
	case PLANEMIRROR:
		calcReflectByPolyDataBatch(startPiont, direction, reflex, intersection, isIntersect);
		break;
	default:
		break;
	}
}

void calculation::RayTracing::calcReflect(const Vector3 & startPiont, const Vector3 & direction, 
	Vector3 & reflex, Vector3 & intersection, bool & isIntersect)
{
	switch (mirror->getMirrorsType())
	{
	case PLANEMIRROR:
		calcReflectByPolyData(startPiont, direction, reflex, intersection, isIntersect);
		break;
	case QUADRICSURFACE:
	case PARABOLICCYLINDER:
	case PARABOLOID:
	case ELLIPSOID:
		calcReflectByQuadricSurface(startPiont, direction, reflex, intersection, isIntersect);
		break;
	default:
		break;
	}
}


void calculation::RayTracing::calcReflectByPolyDataBatch(const vector<vector<Vector3>>& startPiont,
	const vector<vector<Vector3>>& direction,
	vector<vector<Vector3>> &reflex, vector<vector<Vector3>> &intersection,
	vector<vector<bool>> &isIntersect)
{
	for (int i = 0; i < startPiont.size(); i++)
		for (int j = 0; j < startPiont[i].size(); j++)
		{
			bool isTmep = false;
			calcReflectByPolyData(startPiont[i][j], direction[i][j], reflex[i][j],
				intersection[i][j], isTmep);
			isIntersect[i][j] = isTmep;
		}
}

void calculation::RayTracing::calcReflectByPolyData(const Vector3 & startPiont, 
	const Vector3 & direction, Vector3 & reflex, 
	Vector3 & intersection, bool & isIntersect)
{
	vtkSmartPointer<vtkPolyData> polyData = mirror->getPolyData();
	int EleNum = polyData->GetNumberOfCells();
	double t;
	for (int i = 0; i < EleNum; i++)  //���뷴����Ľ���
	{
		vtkIdList * p;
		p = polyData->GetCell(i)->GetPointIds();
		double * point;
		point = polyData->GetPoint(p->GetId(0));
		Vector3 NodesXYZ1(point[0], point[1], point[2]);
		point = polyData->GetPoint(p->GetId(1));
		Vector3 NodesXYZ2(point[0], point[1], point[2]);
		point = polyData->GetPoint(p->GetId(2));
		Vector3 NodesXYZ3(point[0], point[1], point[2]);
		if (this->isIntersect(startPiont, direction, NodesXYZ1,
			NodesXYZ2, NodesXYZ3, intersection, t))
		{
			if (t >= 0)
			{
				Vector3 tempa = NodesXYZ1 - NodesXYZ2;
				Vector3 tempb = NodesXYZ1 - NodesXYZ3;
				Vector3 n_light = tempa.Cross(tempb);  //������

				isIntersect = true;
				reflex = reflectLight(direction, n_light);
				return;
			}
		}
	}
	isIntersect = false;
}

void calculation::RayTracing::calcReflectByQuadricSurface(const Vector3 & startPiont, 
	const Vector3 & direction, Vector3 & reflex, Vector3 & intersection, bool & isIntersect)
{
	double t;
	const vector<double>& tempData = mirror->getData();
	
	if (!isCalcMatrix)
	{
		this->calcMatrix();
	}

	// ����������ϵת��ģ�͵��������ϵ
	Vector3 tempStartPiont = R_translateMatrix[0] * R_rotatMatrix[0] * startPiont;
	Vector3 tempDirection = R_rotatMatrix[0] * direction; // ����ֻ������ת

	if (ray_CurvedSurface(tempData, tempDirection, tempStartPiont, t, intersection))
	{
		const vector<Restriction*> restrictions = mirror->getRestrictionAll();
		// ��ģ�͵��������ϵת����������ϵ
		Vector3 intersectionGlobal = translateMatrix[0] * rotatMatrix[0] * intersection;
		bool isInRestriction = true;
		for (int i = 0; i < restrictions.size(); i++)
		{
			const vector<double>& tempRestrictionData = restrictions[i]->getData();
			// ����������ϵת�������������������ϵ
			Vector3 intersectionLocal = R_translateMatrix[i + 1] *
				R_rotatMatrix[i + 1] * intersectionGlobal;
			if (intersectionLocal.z > tempRestrictionData[1] ||
				intersectionLocal.z < 0.0)
			{
				isInRestriction = false;
				break;
			}
			double tempRadius = intersectionLocal.x * intersectionLocal.x +
				intersectionLocal.y * intersectionLocal.y;

			if (tempRadius > tempRestrictionData[0] * tempRestrictionData[0])
			{
				isInRestriction = false;
				break;
			}
		}
		
		if (!isInRestriction) //��������������
		{
			intersection = startPiont; // �ý��������� ���򲻱� ������޽���ʱ���⴦��
			isIntersect = false;
			return;
		}

		double x = 2 * tempData[0] * intersection.x + tempData[3] * intersection.y +
			tempData[5] * intersection.z + tempData[6];
		double y = 2 * tempData[1] * intersection.y + tempData[3] * intersection.x +
			tempData[4] * intersection.z + tempData[7];
		double z = 2 * tempData[2] * intersection.z + tempData[4] * intersection.y +
			tempData[5] * intersection.x + tempData[8];
		Vector3 tempn(x, y, z);
		if (tempn.Dot(tempDirection) > 0.0)
			tempn.set(-x, -y, -z);

		reflex = reflectLight(tempDirection, tempn);
		isIntersect = true;

		// ��ģ�͵��������ϵת����������ϵ
		intersection = intersectionGlobal;
		reflex = rotatMatrix[0] * reflex; // ���·���

	}
	else
	{
		intersection = startPiont; // �ý��������� ���򲻱� ������޽���ʱ���⴦��
		isIntersect = false;
	}
			

}

void calculation::RayTracing::calcNormalOfLine_MirrorByQuadricSurface(
	const Vector3 & startPiont, const Vector3 & direction, 
	Vector3 & normal, Vector3 & intersection, bool & isIntersect, double & t)
{
	const vector<double>& tempData = mirror->getData();

	if (!isCalcMatrix)
	{
		this->calcMatrix();
	}

	// ����������ϵת��ģ�͵��������ϵ
	Vector3 tempStartPiont = R_translateMatrix[0] * R_rotatMatrix[0] * startPiont;
	Vector3 tempDirection = R_rotatMatrix[0] * direction; // ����ֻ������ת

	double t1, t2;
	bool isOK1 = false;
	bool isOK2 = false;
	Vector3 intersection1, intersection2;
	line_CurvedSurface(tempData, tempDirection, tempStartPiont, t1, t2, isOK1,
		isOK2, intersection1, intersection2);
	if (isOK1&&isOK2)
	{
		Vector3 intersectionGlobal1 = translateMatrix[0] * 
			rotatMatrix[0] * intersection1;
		Vector3 intersectionGlobal2 = translateMatrix[0] *
			rotatMatrix[0] * intersection2;
		bool isInRestriction1 = isInRestriction(intersectionGlobal1);
		bool isInRestriction2 = isInRestriction(intersectionGlobal2);
		if (isInRestriction1&&isInRestriction2) // �����ⶼ���� ȡ����ֵС��
		{
			if (abs(t1) > abs(t2))
			{
				t = t2;
				intersection = intersection2;
			}
			else
			{
				t = t1;
				intersection = intersection1;
			}
		}
		else if (isInRestriction1)
		{
			t = t1;
			intersection = intersection1;
		}
		else if (isInRestriction2)
		{
			t = t2;
			intersection = intersection2;
		}
		else
		{
			intersection = startPiont; // �ý��������� ���򲻱� ������޽���ʱ���⴦��
			isIntersect = false;
			return;
		}
	}
	else if (isOK1)
	{
		Vector3 intersectionGlobal1 = translateMatrix[0] *
			rotatMatrix[0] * intersection1;
		if (isInRestriction(intersectionGlobal1))
		{
			t = t1;
			intersection = intersection1;
		}
		else
		{
			intersection = startPiont; // �ý��������� ���򲻱� ������޽���ʱ���⴦��
			isIntersect = false;
			return;
		}
	}
	else if (isOK2)
	{
		Vector3 intersectionGlobal2 = translateMatrix[0] *
			rotatMatrix[0] * intersection2;
		if (isInRestriction(intersectionGlobal2))
		{
			t = t2;
			intersection = intersection2;
		}
		else
		{
			intersection = startPiont; // �ý��������� ���򲻱� ������޽���ʱ���⴦��
			isIntersect = false;
			return;
		}
	}
	else
	{
		intersection = startPiont; // �ý��������� ���򲻱� ������޽���ʱ���⴦��
		isIntersect = false;
		return;
	}
	double x = 2 * tempData[0] * intersection.x + tempData[3] * intersection.y +
		tempData[5] * intersection.z + tempData[6];
	double y = 2 * tempData[1] * intersection.y + tempData[3] * intersection.x +
		tempData[4] * intersection.z + tempData[7];
	double z = 2 * tempData[2] * intersection.z + tempData[4] * intersection.y +
		tempData[5] * intersection.x + tempData[8];
	Vector3 tempn(x, y, z);
	if (tempn.Dot(tempDirection) > 0.0)
		tempn.set(-x, -y, -z);

	intersection = translateMatrix[0] *
		rotatMatrix[0] * intersection;
	normal = rotatMatrix[0] * tempn; // ���·���
	isIntersect = true;
	
}

void calculation::RayTracing::calcNormalOfLine_MirrorByPolyData(
	const Vector3 & startPiont,const Vector3 & direction, Vector3 & normal, 
	Vector3 & intersection, bool & isIntersect, double & t)
{
	vtkSmartPointer<vtkPolyData> polyData = mirror->getPolyData();
	int EleNum = polyData->GetNumberOfCells();
	for (int i = 0; i < EleNum; i++)  //���뷴����Ľ���
	{
		vtkIdList * p;
		p = polyData->GetCell(i)->GetPointIds();
		double * point;
		point = polyData->GetPoint(p->GetId(0));
		Vector3 NodesXYZ1(point[0], point[1], point[2]);
		point = polyData->GetPoint(p->GetId(1));
		Vector3 NodesXYZ2(point[0], point[1], point[2]);
		point = polyData->GetPoint(p->GetId(2));
		Vector3 NodesXYZ3(point[0], point[1], point[2]);
		if (this->isIntersect(startPiont, direction, NodesXYZ1,
			NodesXYZ2, NodesXYZ3, intersection, t))
		{			
			Vector3 tempa = NodesXYZ1 - NodesXYZ2;
			Vector3 tempb = NodesXYZ1 - NodesXYZ3;
			normal = tempa.Cross(tempb);  //������

			isIntersect = true;
			return;			
		}
	}
	isIntersect = false;
}

bool calculation::RayTracing::isIntersect(const Vector3 & orig, const Vector3 & dir,
	const Vector3 & v0, const Vector3 & v1, const Vector3 & v2, 
	Vector3 & intersection, double & t)
{
	double u, v;
	// E1
	Vector3 E1 = v1 - v0;

	// E2
	Vector3 E2 = v2 - v0;

	// P
	Vector3 P = dir.Cross(E2);

	// determinant
	double det = E1.Dot(P);

	Vector3 T;
	T = orig - v0;

	// If determinant is near zero, ray lies in plane of triangle
	//if (det < 0.00000001 && det > -0.00000001)
	//	return false;

	// Calculate u and make sure u <= 1
	u = T.Dot(P);
	double fInvDet = 1.0f / det;
	u *= fInvDet;

	if (u < 0.0 || u > 1)
		return false;

	// Q
	Vector3 Q = T.Cross(E1);

	// Calculate v and make sure u + v <= 1
	v = dir.Dot(Q);
	v *= fInvDet;
	if (v < 0.0 || u + v > 1)
		return false;

	// Calculate t, scale parameters, ray intersects triangle
	t = E2.Dot(Q);
	t *= fInvDet;

	intersection = orig + dir * t;

	return true;
}

bool calculation::RayTracing::ray_CurvedSurface(const vector<double> &a, 
	Vector3 n, Vector3 org, double & t, Vector3 & interPoint)
{
	double x0 = org.x, y0 = org.y, z0 = org.z;
	double x1 = n.x, y1 = n.y, z1 = n.z;

	double A = a[0] * x1 * x1 + a[1] * y1 * y1 + a[2] * z1 * z1 + a[3] * x1 * y1 +
		a[4] * z1 * y1 + a[5] * x1 * z1;
	double B = 2 * a[0] * x1 * x0 + 2 * a[1] * y1 * y0 + 2 * a[2] * z1 * z0 +
		a[3] * (x0 * y1 + x1 * y0) + a[4] * (z0 * y1 + z1 * y0) + a[5] * (z0 * x1 + z1 * x0) +
		a[6] * x1 + a[7] * y1 + a[8] * z1;
	double C = a[0] * x0 * x0 + a[1] * y0 * y0 + a[2] * z0 * z0 +
		a[3] * x0 * y0 + a[4] * z0 * y0 + a[5] * x0 * z0 +
		a[6] * x0 + a[7] * y0 + a[8] * z0 + a[9];

	if (A < -THRESHOLD || A > THRESHOLD)
	{
		double temp = B * B - 4 * A * C;
		if (temp >= 0)
			temp = pow(temp, 0.5);
		else
			return false;

		double tempt1, tempt2;
		tempt1 = (-B + temp) / 2.0 / A;
		tempt2 = (-B - temp) / 2.0 / A; // �����ʽ��������

		if (tempt1 >= 0.0 && tempt2 >= 0.0) // �����ڵ���0 ��Ҫ�ȱȽ��Ƿ��ڶ�������������
		{
			bool isIn1 = false;
			Vector3 interPoint1(x0 + x1 * tempt1, y0 + y1 * tempt1, z0 + z1 * tempt1);	
			// �ж��Ƿ��ڸ�����������
			if (a[10] - THRESHOLD < interPoint1.x  &&
				interPoint1.x < a[11] + THRESHOLD &&
				a[12] - THRESHOLD < interPoint1.y &&
				interPoint1.y < a[13] + THRESHOLD &&
				a[14] - THRESHOLD < interPoint1.z &&
				interPoint1.z < a[15] + THRESHOLD)
			{
				isIn1 = true;
			}
			bool isIn2 = false;
			Vector3 interPoint2(x0 + x1 * tempt2, y0 + y1 * tempt2, z0 + z1 * tempt2);
			// �ж��Ƿ��ڸ�����������
			if (a[10] - THRESHOLD < interPoint2.x  &&
				interPoint2.x < a[11] + THRESHOLD &&
				a[12] - THRESHOLD < interPoint2.y &&
				interPoint2.y < a[13] + THRESHOLD &&
				a[14] - THRESHOLD < interPoint2.z &&
				interPoint2.z < a[15] + THRESHOLD)
			{
				isIn2 = true;
			}
			if (isIn1 && isIn2)
			{
				if (tempt1 > tempt2)
				{
					interPoint = interPoint1;
					return true;
				}
				else
				{
					interPoint = interPoint2;
					return true;
				}
					
			}
			else if (isIn1 && !isIn2)
			{
				interPoint = interPoint1;
				return true;
			}
			else if (!isIn1 && isIn2)
			{
				interPoint = interPoint2;
				return true;
			}
			else 
			{
				return false;
			}
			
		}
		else if (tempt1 < 0.0 && tempt2 < 0.0) // ��С��0 �޽�
		{
			return false;
		}
		else                           // ȡ��ֵ
		{
			if (tempt1 < tempt2)
				t = tempt2;
			else
				t = tempt1;
		}
	}
	else                          // ֻ��һ�����㣬�뷨��ƽ��
		t = -C / B;

	if (t < 0.0)
		return false;
	else
	{
		// �ж��Ƿ��ڸ�����������
		interPoint.set(x0 + x1 * t, y0 + y1 * t, z0 + z1 * t);
		if (a[10] - THRESHOLD < interPoint.x  &&
			interPoint.x < a[11] + THRESHOLD &&
			a[12] - THRESHOLD < interPoint.y &&
			interPoint.y < a[13] + THRESHOLD &&
			a[14] - THRESHOLD < interPoint.z &&
			interPoint.z < a[15] + THRESHOLD)
		{
			return true;
		}
		else
			return false;
	}
}

void calculation::RayTracing::line_CurvedSurface(const vector<double>& a,
	Vector3 n, Vector3 org, double & t1, double & t2, bool & isOk1, bool & isOk2,
	Vector3 &interPoint1, Vector3 &interPoint2)
{
	double x0 = org.x, y0 = org.y, z0 = org.z;
	double x1 = n.x, y1 = n.y, z1 = n.z;

	double A = a[0] * x1 * x1 + a[1] * y1 * y1 + a[2] * z1 * z1 + a[3] * x1 * y1 +
		a[4] * z1 * y1 + a[5] * x1 * z1;
	double B = 2 * a[0] * x1 * x0 + 2 * a[1] * y1 * y0 + 2 * a[2] * z1 * z0 +
		a[3] * (x0 * y1 + x1 * y0) + a[4] * (z0 * y1 + z1 * y0) + a[5] * (z0 * x1 + z1 * x0) +
		a[6] * x1 + a[7] * y1 + a[8] * z1;
	double C = a[0] * x0 * x0 + a[1] * y0 * y0 + a[2] * z0 * z0 +
		a[3] * x0 * y0 + a[4] * z0 * y0 + a[5] * x0 * z0 +
		a[6] * x0 + a[7] * y0 + a[8] * z0 + a[9];

	if (A < -THRESHOLD || A > THRESHOLD)
	{
		double temp = B * B - 4 * A * C;
		if (temp >= 0)
			temp = pow(temp, 0.5);
		else
		{
			isOk1 = false;
			isOk2 = false;
			return;
		}
		double tempt1, tempt2;
		tempt1 = (-B + temp) / 2.0 / A;
		tempt2 = (-B - temp) / 2.0 / A; // �����ʽ��������

		interPoint1.set(x0 + x1 * tempt1, y0 + y1 * tempt1, z0 + z1 * tempt1);
		// �ж��Ƿ��ڸ�����������
		if (a[10] - THRESHOLD < interPoint1.x  &&
			interPoint1.x < a[11] + THRESHOLD &&
			a[12] - THRESHOLD < interPoint1.y &&
			interPoint1.y < a[13] + THRESHOLD &&
			a[14] - THRESHOLD < interPoint1.z &&
			interPoint1.z < a[15] + THRESHOLD)
		{
			t1 = tempt1;
			isOk1 = true;
		}
		interPoint2.set(x0 + x1 * tempt2, y0 + y1 * tempt2, z0 + z1 * tempt2);
		// �ж��Ƿ��ڸ�����������
		if (a[10] - THRESHOLD < interPoint2.x  &&
			interPoint2.x < a[11] + THRESHOLD &&
			a[12] - THRESHOLD < interPoint2.y &&
			interPoint2.y < a[13] + THRESHOLD &&
			a[14] - THRESHOLD < interPoint2.z &&
			interPoint2.z < a[15] + THRESHOLD)
		{
			t2 = tempt2;
			isOk2 = true;
		}
	}
	else // ֻ��һ�����㣬�뷨��ƽ��
	{
		isOk2 = false;
		double tempt1 = -C / B;
		interPoint1.set(x0 + x1 * tempt1, y0 + y1 * tempt1, z0 + z1 * tempt1);
		if (a[10] - THRESHOLD < interPoint1.x  &&
			interPoint1.x < a[11] + THRESHOLD &&
			a[12] - THRESHOLD < interPoint1.y &&
			interPoint1.y < a[13] + THRESHOLD &&
			a[14] - THRESHOLD < interPoint1.z &&
			interPoint1.z < a[15] + THRESHOLD)
		{
			t1 = tempt1;
			isOk1 = true;
		}
		else
		{
			isOk1 = false;
		}
	}
}

bool calculation::RayTracing::isInRestriction(const Vector3 & intersectionGlobal)
{
	const vector<Restriction*> restrictions = mirror->getRestrictionAll();
	for (int i = 0; i < restrictions.size(); i++)
	{
		const vector<double>& tempRestrictionData = restrictions[i]->getData();
		// ����������ϵת�������������������ϵ
		Vector3 intersectionLocal = R_translateMatrix[i + 1] *
			R_rotatMatrix[i + 1] * intersectionGlobal;
		if (intersectionLocal.z > tempRestrictionData[1] ||
			intersectionLocal.z < 0.0)
		{
			return false;
		}
		double tempRadius = intersectionLocal.x * intersectionLocal.x +
			intersectionLocal.y * intersectionLocal.y;

		if (tempRadius > tempRestrictionData[0] * tempRestrictionData[0])
		{
			return false;
			break;
		}
	}
	return true;
}

Vector3 calculation::RayTracing::reflectLight(const Vector3 & a, const Vector3 & n)
{
	//�ȵ�λ��
	double absa = pow(a.Dot(a), 0.5);
	double absn = pow(n.Dot(n), 0.5);
	Vector3 tempa = a * (1 / absa);
	Vector3 tempn = n * (1 / absn);
	double I = 2 * tempn.Dot(tempa);
	if (I < 0)
		I = -I;
	else
		tempa = Vector3(0.0, 0.0, 0.0) - tempa;

	return tempn * I + tempa;
}

void calculation::RayTracing::calcMatrix()
{
	const vector<Restriction*> restrictions = mirror->getRestrictionAll();
	int num = restrictions.size() + 1;
	R_rotatMatrix.resize(num);
	R_translateMatrix.resize(num);
	rotatMatrix.resize(num);
	translateMatrix.resize(num);

	GraphTrans tempGraphTrans;

	for (int i = 0; i < restrictions.size()+1; ++i)
	{
		if(0 == i)
			tempGraphTrans = mirror->getGraphTrans();
		else
			tempGraphTrans = restrictions[i-1]->getGraphTrans();
		// ��������ϵת��ģ�͵��������ϵ�������������ת��ƽ��
		Vector3D RotateAsix(tempGraphTrans.getRotate_x(),
			tempGraphTrans.getRotate_y(), tempGraphTrans.getRotate_z());
		R_rotatMatrix[i] = Matrix4D::getRotateMatrix(-tempGraphTrans.getRotate_theta(),
			RotateAsix);
		Vector3D rotatTranslate(tempGraphTrans.getTrans_x(), tempGraphTrans.getTrans_y(),
			tempGraphTrans.getTrans_z());
		rotatTranslate = R_rotatMatrix[i] * rotatTranslate; // ����ת��ƽ�ƣ���ƽ�Ƶ�����һ����ת��
		R_translateMatrix[i] = Matrix4D::getTranslateMatrix(rotatTranslate * (-1));

		// ģ�͵��������ϵת�������������
		rotatMatrix[i] = Matrix4D::getRotateMatrix(tempGraphTrans.getRotate_theta(),
			RotateAsix);
		translateMatrix[i] = Matrix4D::getTranslateMatrix(tempGraphTrans.getTrans_x(), 
			tempGraphTrans.getTrans_y(),tempGraphTrans.getTrans_z());
	}

	isCalcMatrix = true;
}
