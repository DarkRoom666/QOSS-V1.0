#pragma once
//��������ά������������-�������˵�˺Ͳ��

#ifndef VECTOR3F_H
#define VECTOR3F_H

#include<iostream>

class  Vector3f
{
public:
	Vector3f(float x = 0, float y = 0, float z = 0)
		:x(x), y(y), z(z)
	{

	};

	~Vector3f()
	{
	};

	//����"-"
	Vector3f operator - (const Vector3f &v) const
	{
		return Vector3f(x - v.x, y - v.y, z - v.z);
	};

	//����"+"
	Vector3f operator + (const Vector3f &v) const
	{
		return Vector3f(x + v.x, y + v.y, z + v.z);
	};

	//������
	Vector3f operator * (const float &v) const
	{
		return Vector3f(x * v, y * v, z * v);
	};

	//���
	double Dot(const Vector3f &v) const
	{
		return x*v.x + y*v.y + z*v.z;
	};

	//���
	Vector3f Cross(const Vector3f &v) const
	{
		return Vector3f(
			y * v.z - z * v.y,
			z * v.x - x * v.z,
			x * v.y - y * v.x);
	}

	void Normalization()
	{
		double temp = pow((x*x + y*y + z*z), 0.5);
		x /= temp;
		y /= temp;
		z /= temp;
	}

	double Length() const
	{
		return pow((x*x + y*y + z*z),0.5);
	}

	void display()
	{
		std::cout << x << " " << y << " " << z << std::endl;
	}

	void set(double x1, double y1, double z1)
	{
		x = x1; y = y1; z = z1;
	}

	double x, y, z;
};



#endif // VECTOR3F_H