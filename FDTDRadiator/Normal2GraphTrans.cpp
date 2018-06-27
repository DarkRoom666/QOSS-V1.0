#include "Normal2GraphTrans.h"
double CalDistance(const Vector3 &a, const Vector3 &b)
{
	return pow(pow((a.x - b.x), 2) + pow((a.y - b.y), 2) + pow((a.z - b.z), 2), 0.5);
}

void updateSource_n(const Vector3& new_n, GraphTrans & _gt)
{
	//new_n.Normalization();
	if (new_n.x != 0 || new_n.y != 0 || new_n.z != 1)
	{
		Vector3 rotate_axis = Vector3(0, 0, 1).Cross(new_n); // Ðý×ªÖá
		double rotate_theta = acos(Vector3(0, 0, 1).Dot(new_n));
		rotate_theta = rotate_theta / Pi * 180;
		_gt.updateRotate(rotate_axis, rotate_theta);
	}
	else
	{
		_gt.updateRotate(Vector3(0, 0, 1), 0);
	}
}