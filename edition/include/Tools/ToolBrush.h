#pragma once

#include "libs/polygon.h"

#include "Tool.h"
#include "Kernels/Kernel.h"

class ToolBrush : public Tool
{
private:
	enum class State
	{
		START,
		LASSO,
		MOVE,
	};

	Polygon2 m_lasso;

	State m_state{State::START};
	bool m_canMove{false};

	std::vector<Kernel*> m_selectedKernels;

	void movePrimitives(const Vector2& dest);
	void addToSelection();

	void updateRenderer();
	void updateShowKernels() const;

	void reset();

	QPointF m_oldPos;

public:
	ToolBrush(VectorTerrainRaytracingWidget* parent, Kernels& kernels): Tool(parent, kernels)
	{
	}

	void mouseMoveEvent(QMouseEvent* e) override;
	void mousePressEvent(QMouseEvent* e) override;
	void mouseReleaseEvent(QMouseEvent* e) override;
	void mouseWheelEvent(QWheelEvent* e) override;
	void keyPressedEvent(QKeyEvent* e) override;

	void apply() override;
	void update() override;

	void saveBrush(const QString& filename);
	void loadBrush(const QString& filename);

	ToolType type() { return ToolType::MOVE; }
};

