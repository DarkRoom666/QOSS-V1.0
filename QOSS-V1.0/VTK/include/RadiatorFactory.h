/*
*	created by liyun 2017/11/29
*   function �������ַ���������
*   version 1.0
*/

#ifndef RADIATORFACTORY_H  
#define RADIATORFACTORY_H

#include "Radiator.h"

#include <vector>
#include <memory>

class RadiatorFactory
{
public:

	static std::shared_ptr<Radiator> getRadiator(Radiator::RadiatorType type,
		shared_ptr<calculation::SourceModeGeneration>);

private:
	RadiatorFactory() {};
	~RadiatorFactory() {};

};



#endif // RADIATORFACTORY_H
