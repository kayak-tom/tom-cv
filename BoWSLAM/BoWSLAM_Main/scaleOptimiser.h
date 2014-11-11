#pragma once

#include "util/exception.h"

//relscale d,g + a triple of nodes

namespace NSLAMMap
{
	class CSLAMMap;
	class CCycle;
}

class CScaleOptimiser
{
public:
	virtual void run(const bool bVerbose) = 0;
	virtual void addCycle(NSLAMMap::CCycle & c, const bool bVerbose) = 0;
	virtual ~CScaleOptimiser() {}
};
