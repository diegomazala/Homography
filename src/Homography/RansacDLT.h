#ifndef __RANSAC_DLT_H__
#define __RANSAC_DLT_H__


#include <vector>
#include "DLT.h"
#include "Points.h"


class RansacDLT
{
public:

	static DLT solve(const Points& points);
	static std::vector<int> randomIndices(int count, int min, int max);
};



#endif // __RANSAC_DLT_H__s