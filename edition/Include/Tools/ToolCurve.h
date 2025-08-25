#pragma once
#include "curve.h"
#include "Tool.h"

class ToolCurve : public Tool
{
private:
	enum class State
	{
		FIRST_CURVE,
		LASSO,
		SECOND_CURVE,
	};

	State m_state{State::FIRST_CURVE};
	QVector<Vector2> m_curvePointsFrom{};
	QVector<Vector2> m_curvePointsTo{};
	QuadricCurve2 m_curveFrom{};
	QuadricCurve2 m_curveTo{};

	QPolygonF m_lasso{};
	bool m_editLasso{false};

	std::unique_ptr<MayaSimpleRenderer> m_rendererCurveTo{nullptr};
	std::unique_ptr<MayaSimpleRenderer> m_rendererLasso{nullptr};

	void moveGaussians(double distMax) const;
	void deleteRenderer(std::unique_ptr<MayaSimpleRenderer>& renderer);

public:
	ToolCurve(GaussianTerrainRaytracingWidget* parent, Kernels& kernels): Tool(parent, kernels)
	{
	}

	void mouseMoveEvent(QMouseEvent* e) override;
	void mousePressEvent(QMouseEvent* e) override;

	~ToolCurve() override
	{
		if (m_rendererCurveTo) m_rendererCurveTo->DeleteBuffers();
		if (m_rendererLasso) m_rendererLasso->DeleteBuffers();
	}

	void render() const override
	{
		if (m_rendererCurveTo) m_rendererCurveTo->Draw();
		if (m_rendererLasso) m_rendererLasso->Draw();
		this->Tool::render();
	}

	ToolType type() override { return ToolType::CURVE; }
};

