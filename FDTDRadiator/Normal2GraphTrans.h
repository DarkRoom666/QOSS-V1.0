#pragma once
#include "vector3.h"
#include "Matrix4D.h"
#include "GraphTrans.h"
#include "Constant_Val.h"

double CalDistance(const Vector3 &a, const Vector3 &b);

void updateSource_n(const Vector3& new_n, GraphTrans & _gt);
