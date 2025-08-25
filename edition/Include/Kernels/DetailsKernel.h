#pragma once

#include <vector>

#include "libs/ellipse.h"
#include "libs/evector.h"

#include "Kernel.h"

class DetailsKernel : public Kernel
{
private:
	float m_scaleXOrig, m_scaleYOrig, m_thetaOrig, m_thetaWODeform;
	float m_textureX, m_textureY;
	void print(std::ostream& os) const override;

	std::vector<float> getArray() const override
	{
		return std::vector{
			static_cast<float>(type()), m_scaleX, m_scaleY, m_theta,
				m_amplitude* m_modAmplitude, m_posX, m_posY, m_textureX, m_textureY, m_scaleXOrig, m_scaleYOrig, m_thetaOrig, m_thetaWODeform
		};
	}

public:
	DetailsKernel(const float scaleX, const float scaleY, const float theta, const float amplitude,
		const float posX, const float posY, const float textureX, const float textureY) :
		Kernel(scaleX, scaleY, theta, amplitude, posX, posY), m_textureX(textureX), m_textureY(textureY), m_scaleXOrig(scaleX), m_scaleYOrig(scaleY), m_thetaOrig(theta), m_thetaWODeform(theta)
	{
	}

	DetailsKernel(const float scaleX, const float scaleY, const float theta, const float amplitude,
		const float posX, const float posY, const float textureX, const float textureY, const float scaleXOrig, const float scaleYOrig, const float thetaOrig, const float thetaWODeform) :
		Kernel(scaleX, scaleY, theta, amplitude, posX, posY), m_textureX(textureX), m_textureY(textureY), m_scaleXOrig(scaleXOrig), m_scaleYOrig(scaleYOrig), m_thetaOrig(thetaOrig), m_thetaWODeform(thetaWODeform)
	{
	}

	explicit DetailsKernel(const std::vector<float>& v) : Kernel(v[0], v[1], v[2], v[3], v[4], v[5]), m_textureX(v[6]), m_textureY(v[7]), m_scaleXOrig(v[0]), m_scaleYOrig(v[1]), m_thetaOrig(v[2]), m_thetaWODeform(v[2])
	{
	}

	std::unique_ptr<Kernel> clone() const
	{
		return std::make_unique<DetailsKernel>(*this);
	}

	static constexpr int nbParamInit = 8;
	static constexpr int size = 13;
	KernelType type() const override { return KernelType::DETAILS; }

	float& textureX() { return m_textureX; }
	float& textureY() { return m_textureY; }

	void rotate(float theta, const Vector2& o)
	{
		Kernel::rotate(theta, o);
		m_thetaWODeform -= theta;
	}
};

