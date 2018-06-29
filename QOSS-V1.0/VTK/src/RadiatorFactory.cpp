#include "..\include\RadiatorFactory.h"
#include "..\include\LowOrderRadiator.h"
#include "..\include\HighRadiatorShow.h"

std::shared_ptr<Radiator> RadiatorFactory::getRadiator(Radiator::RadiatorType type,
	shared_ptr<calculation::SourceModeGeneration> source)
{
	std::shared_ptr<Radiator> res;
	switch (type)
	{
	case Radiator::LOWORDER:
		res = std::make_shared<LowOrderRadiator>(source);
		break;
	case Radiator::HIGHT:
		//res = std::make_shared<HighRadiatorShow>(source);
		break;
	default:
		break;
	}
	return res;
}

