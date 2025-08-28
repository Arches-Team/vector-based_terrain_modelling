#include "Kernels.h"

#include <QTextStream>
#include <QtCore/qfile.h>

#include "npy.hpp"

#include "Kernels/GaussianKernel.h"
#include "Kernels/FactoryKernel.h"
#include "utils.h"

void Kernels::loadCSVFile(const QString& filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly))
	{
		std::cout << file.errorString().toStdString() << std::endl;
		return;
	}

	QString data = file.readAll();
	auto lines = data.split('\n');

	auto primitivesStr = lines[0].split(',');

	file.close();

	clear();

	int nbParam = utils::sizeKernel();
	int i = 0;
	std::vector<float> kernel;
	int id = -1;
	for (const auto& pStr : primitivesStr)
	{
		float cell = std::stof(pStr.toStdString());
		if (i == 0)
			id = static_cast<int>(cell);
		else
			kernel.push_back(cell);
		i++;

		if (i == nbParam)
		{
			m_kernels.emplace_back(kernel::create(static_cast<KernelType>(id), kernel));
			i = 0;
			kernel.clear();
		}
	}
}

void Kernels::loadNPYFile(const QString& filename)
{
	clear();

	auto [data, shape, fortran_order] = npy::read_npy<float>(filename.toStdString());

	for (int i = 0; i < data.size(); i += GaussianKernel::nbParamInit)
	{
		std::vector<float> k(GaussianKernel::nbParamInit);
		for (int j = 0; j < GaussianKernel::nbParamInit; ++j)
		{
			k[j] = data[i + j];
		}

		m_kernels.emplace_back(std::make_unique<GaussianKernel>(k, true));
	}

	sort();
}

void Kernels::saveCSVFile(const QString& filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
	{
		std::cout << file.errorString().toStdString() << std::endl;
		return;
	}
	QTextStream out(&file);

	bool isBeginning = true;
	for (auto& k : m_kernels)
	{
		auto data = k->get();
		for (const float& d : data)
		{
			if (!isBeginning)
				out << ",";
			isBeginning = false;
			out << QString::number(d);
		}
	}
	out << "\n";

	file.close();
}

std::vector<float> Kernels::getArray()
{
	std::vector<float> ret(utils::sizeKernel() * size());
	int i{0};
	for (auto& k : m_kernels)
	{
		auto a = k->get();
		std::ranges::copy(a, ret.begin() + i);
		i += utils::sizeKernel();
	}
	return ret;
}

void Kernels::sort()
{
	std::ranges::sort(m_kernels, [](const std::unique_ptr<Kernel>& a, const std::unique_ptr<Kernel>& b)
								 {
									return *a > *b;
								 }
	);
}

float Kernels::getMaxAbsAmplitude()
{
	float maxAmplitude = 0.;
	for (auto& k : m_kernels)
		maxAmplitude = std::max(std::abs(k->amplitude()), maxAmplitude);
	return maxAmplitude;
}

void Kernels::clear(const std::vector<Kernel*>& kernels)
{
	std::vector<int> idx{};
	for (auto& kRef : kernels)
	{
		for (int i = 0; i < m_kernels.size(); ++i)
		{
			if (kRef == m_kernels[i].get())
			{
				idx.push_back(i - static_cast<int>(idx.size()));
				break;
			}
		}
	}
	for (const int i : idx)
		m_kernels.erase(m_kernels.begin() + i);
}

void Kernels::clear()
{
	m_kernels.clear();
}

void Kernels::removeIf(const std::function<bool(const std::unique_ptr<Kernel>&)>& predicate)
{
	const auto newEnd = std::ranges::remove_if(m_kernels, 
		[&predicate](const std::unique_ptr<Kernel>& k) {return predicate(k); }
	).begin();

	m_kernels.erase(newEnd, m_kernels.end());
}

Kernel& Kernels::operator[](const int i)
{
	return *m_kernels[i];
}
