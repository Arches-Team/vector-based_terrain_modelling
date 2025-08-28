#pragma once

#include <vector>

#include "libs/gabornoise.h"
#include "Kernels/Kernel.h"

class Kernels
{
private:
	std::vector<std::unique_ptr<Kernel>> m_kernels;

public:
	Kernels() = default;
	Kernels(const Kernels& k)
	{
		for (const auto& kernel : k.m_kernels)
			m_kernels.push_back(kernel->clone());
	}
	Kernels& operator=(const Kernels& k)
	{
		m_kernels.clear();
		for (const auto& kernel : k.m_kernels)
			m_kernels.push_back(kernel->clone());
		return *this;
	}

	void loadCSVFile(const QString& filename);
	void loadNPYFile(const QString& filename);
	void saveCSVFile(const QString& filename);
	int size() { return static_cast<int>(m_kernels.size()); }

	std::vector<float> getArray();
	void sort();

	template <typename T>
	Kernel* add(const std::vector<float>& v)
	{
		return m_kernels.emplace_back(std::make_unique<T>(v)).get();
	}

	Kernel* add(std::unique_ptr<Kernel> k)
	{
		return m_kernels.emplace_back(std::move(k)).get();
	}

	// Statistics
	float getMaxAbsAmplitude();

	void clear(const std::vector<Kernel*>& kernels);
	void clear();
	void removeIf(const std::function<bool(const std::unique_ptr<Kernel>&)>& predicate);
	bool empty() const { return m_kernels.empty(); }

	// Probably not the best idea to expose private parameters.
	std::vector<std::unique_ptr<Kernel>>& getKernels() { return m_kernels; }

	Kernel& operator[](int i);

	// Iterator methods
	std::vector<std::unique_ptr<Kernel>>::iterator begin() { return m_kernels.begin(); }
	std::vector<std::unique_ptr<Kernel>>::iterator end() { return m_kernels.end(); }
	std::vector<std::unique_ptr<Kernel>>::const_iterator begin() const { return m_kernels.begin(); }
	std::vector<std::unique_ptr<Kernel>>::const_iterator end() const { return m_kernels.end(); }
	std::vector<std::unique_ptr<Kernel>>::const_iterator cbegin() const { return m_kernels.cbegin(); }
	std::vector<std::unique_ptr<Kernel>>::const_iterator cend() const { return m_kernels.cend(); }
};

