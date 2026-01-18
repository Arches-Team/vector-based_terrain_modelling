#pragma once

#include <QtCore/qpoint.h>

#include "libs/evector.h"


namespace utils
{
	constexpr float eps = 1e-8f;
	constexpr float pi = 3.141592653589793f;
	Vector2 toVec2(QPointF p);

	constexpr int sizeKernel()
	{
		//TODO : cyclic reference
		int max = 13;//GaussianKernel::size;
		// max = std::max(max, DetailsKernel::size);
		return max;
	}
}
