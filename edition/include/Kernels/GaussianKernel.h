#pragma once

#include <vector>

#include "libs/ellipse.h"

#include "Kernel.h"

class GaussianKernel : public Kernel
{
private:
	float m_beta;
	void print(std::ostream& os) const override;

	std::vector<float> getArray() const override
	{
		return std::vector{
			static_cast<float>(type()), m_scaleX, m_scaleY, m_theta,
				m_amplitude* m_modAmplitude, m_posX, m_posY, m_beta
		};
	}

public:
	GaussianKernel(const float scaleX, const float scaleY, const float modScale, const float theta, const float amplitude,
	       const float posX, const float posY, const float beta):
		Kernel(scaleX * modScale, scaleY * modScale, theta, amplitude, posX, posY), m_beta(beta)
	{
	}

	explicit GaussianKernel(const std::vector<float>& v, bool fromOptim = false): Kernel(v[0] * v[2], v[1] * v[2], v[3], v[4], v[5], v[6]), m_beta(v[7])
	{
		if (fromOptim)
		{
			m_scaleX *= 2.f;
			m_scaleY *= 2.f;
			m_theta = 3 * M_PI / 2 - m_theta;
		}
	}

	std::unique_ptr<Kernel> clone() const
	{
		return std::make_unique<GaussianKernel>(*this);
	}

	static constexpr int nbParamInit = 8;
	static constexpr int size = 8;
	KernelType type() const override { return KernelType::GAUSSIAN; }

	float& beta() { return m_beta; }
};

