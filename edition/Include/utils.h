#pragma once

#include <iostream>
#include <QtCore/qpoint.h>
#include <QtGui/qpolygon.h>
#include "libs/evector.h"

//#include "Kernels/GaussianKernel.h"

namespace utils
{
	constexpr float eps = 1e-8f;
	constexpr float pi = 3.141592653589793f;
	//QPointF getCentroid(QPolygonF poly);
	Vector2 toVec2(QPointF p);
	// std::vector<Kernel> loadCSV(const std::string& filename);

	constexpr int sizeKernel()
	{
		//TODO : cyclic reference
		int max = 13;//GaussianKernel::size;
		// max = std::max(max, DetailsKernel::size);
		return max;
	}
}
