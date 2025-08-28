#pragma once

#include <QtGui/QMouseEvent>

#include "libs/mayarender.h"
#include "Kernels.h"

class VectorTerrainRaytracingWidget;

enum class ToolType
{
	NONE,
	HAND,
	MOVE,
	EDIT,
	GRAPH,
};

class Tool
{
protected:
	std::unique_ptr<MayaSimpleRenderer> m_renderer{nullptr};
	VectorTerrainRaytracingWidget* m_parent{nullptr};
	Kernels& m_kernels;

	virtual void deleteRenderer()
	{
		if (m_renderer)
		{
			m_renderer->DeleteBuffers();
			m_renderer = nullptr;
		}
	}

public:
	explicit Tool(VectorTerrainRaytracingWidget* parent, Kernels& kernels): m_parent(parent),
		m_kernels(kernels)
	{
	}

	virtual ~Tool()
	{
		// When the tool is destroyed, apply the effects (if any).
		apply();
		deleteRenderer();
	}

	virtual void render() const
	{
		if (m_renderer)
			m_renderer->Draw();
	}

	virtual ToolType type() { return ToolType::NONE; }

	virtual void mouseMoveEvent(QMouseEvent* e)
	{
	}

	virtual void mousePressEvent(QMouseEvent* e)
	{
	}

	virtual void mouseReleaseEvent(QMouseEvent* e)
	{
	}

	virtual void mouseWheelEvent(QWheelEvent* e)
	{
	}

	virtual void keyPressedEvent(QKeyEvent* e)
	{
	}

	virtual void apply()
	{
	}

	virtual void update()
	{
	}
};

