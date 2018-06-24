#pragma once
//这是一个计算阶乘以及的专用函数

#ifndef MATHEMATICAL_FUNCTIONS_JIN_H
#define MATHEMATICAL_FUNCTIONS_JIN_H

#include <iostream>
#include <cmath>
#include <complex>//复数类，C++自带
#include <vector>//vector用于存放数组与矩阵数据  
#include <string>
#include "../util/Constant_Var.h"

using namespace std;

//Add by Jin Ming
double rootdbessel(int m, int n);	//Find the nth root of mth order derivate bessel function of the first kind

double rootbessel(int m, int n);	//Find the nth root of mth order bessel function of the first kind

double kzmnTE(int m, int n, double lambda, double radius);

#endif // 
//20171026已经验证