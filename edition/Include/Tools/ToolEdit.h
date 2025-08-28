#pragma once

#include "libs/heightfield.h"

#include "Tool.h"

class ToolEdit : public Tool
{
public:
	enum class Mode
	{
		ERASE,
		AMPLITUDE,
		AMPLITUDELR,
		AMPLITUDEHR,
		WARP,
		MOVE,
	};

private:
	bool m_canEdit{false};
	double m_radius{1250 * 0.25};
	Vector m_intersectionPoint;
	Vector m_clickedIntersectionPoint;
	Vector2 m_oldPos;

	Mode m_mode{Mode::ERASE};

	Vector2 m_oldMouseScreenPos;
	std::vector<Kernel*> m_selectedPrimitives;

public:
	ToolEdit(VectorTerrainRaytracingWidget* parent, Kernels& kernels): Tool(parent, kernels)
	{
	}

	void mouseMoveEvent(QMouseEvent* e) override;
	void mousePressEvent(QMouseEvent* e) override;
	void mouseReleaseEvent(QMouseEvent* e) override;
	void mouseWheelEvent(QWheelEvent* e) override;

	void erasePrimitives(const Vector& intersectionPoint, const Vector2& boxSize) const;
	void changeAmplitude(const Vector2& mouseScreenPos, double size_min = 0., double size_max = 1.0);
	void warp(const Vector2& mouseScreenPos);
	void move(const Vector2& intersection);
	void defineRenderer(const ScalarField2* h, const Vector& p);

	ToolType type() { return ToolType::EDIT; }

	void setRadius(const double radius) { m_radius = radius; }
	double getRadius() const { return m_radius; }
	void setMode(const Mode& mode) { m_mode = mode; }
    
};

