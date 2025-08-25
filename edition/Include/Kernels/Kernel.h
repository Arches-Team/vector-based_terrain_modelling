#pragma once

#include <vector>

#include "libs/ellipse.h"
#include "libs/evector.h"

#include "utils.h"

enum class KernelType
{
	GAUSSIAN = 0,
	DETAILS = 1,
};

class Kernel
{
protected:
	float m_scaleX, m_scaleY, m_theta, m_amplitude, m_posX, m_posY;
	float m_modAmplitude{1.};

	virtual void print(std::ostream& os) const = 0;

	virtual std::vector<float> getArray() const = 0;

public:
	Kernel(const float scaleX, const float scaleY, const float theta, const float amplitude,
	       const float posX, const float posY):
		m_scaleX(scaleX), m_scaleY(scaleY), m_theta(theta), m_amplitude(amplitude), m_posX(posX), m_posY(posY)
	{
	}

	explicit Kernel(const std::vector<float>& v): m_scaleX(v[0]), m_scaleY(v[1]), m_theta(v[2]),
	                                              m_amplitude(v[3]), m_posX(v[4]), m_posY(v[5])
	{
	}

	virtual std::unique_ptr<Kernel> clone() const = 0;

	virtual KernelType type() const = 0;

	std::vector<float> get() const
	{
		auto ret = getArray();
		// Padding
		while (ret.size() < utils::sizeKernel())
			ret.emplace_back(0.);
		return ret;
	}

	void setModAmplitude(const float val) { m_modAmplitude = val; }
	void setAmplitude(const float val) { m_amplitude = val; };
	bool isShown() const { return (m_modAmplitude - utils::eps) > 0.; }

	Vector2 pos() const { return Vector2(m_posX, m_posY); }
	float& scaleX() { return m_scaleX; }
	float& scaleY() { return m_scaleY; }
	float& posX() { return m_posX; }
	float& posY() { return m_posY; }
	float& theta() { return m_theta; }
	float& amplitude() { return m_amplitude; }

	double distCenter(const Vector2& p) const;
	Ellipse2 getEllipse(const float sigma = 3.f) const;
	bool isInside(const Vector2& p) const;

	void scale(float scale, const Vector2& o);
	void scale(Vector2 axis, float factor);
	void translate(const Vector2& offset);
	void setPos(const Vector2& pos);
	virtual void rotate(float theta, const Vector2& o);

	// Normalizing the kernel by setting scaleX as the major axis and theta \in [0 ; pi].
	// Does not change visually the kernel.
	void normalize();

	float volume() const { return 2.f * utils::pi * m_scaleX * m_scaleY * std::abs(m_amplitude); }
	bool operator<(const Kernel& k) const { return volume() < k.volume(); }
	bool operator>(const Kernel& k) const { return !(*this < k); }

	friend std::ostream& operator<<(std::ostream& os, const Kernel& k)
	{
		k.print(os);
		return os;
	}
};

