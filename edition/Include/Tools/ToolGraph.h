#pragma once

#include "VectorTerrainRaytracingWidget.h"
#include "Kernels/Kernel.h"
#include "graph.h"
#include "Tool.h"

class ToolGraph : public Tool
{
private:
	struct KernelData
	{
		Kernel* kernel;
		double curvilignAbscissa;
		double dist;
	};

	struct KernelLink
	{
		Vector2 origOldPos;
		Vector2 destOldPos;
		// Kernel, curvilign abscissa
		std::vector<KernelData> kernels;
	};

	struct Spring
	{
		Node* orig;
		Node* dest;
		double l0;
		double k;
	};

	// https://stackoverflow.com/questions/15160889/how-can-i-make-an-unordered-set-of-pairs-of-integers-in-c
	struct pair_hash
	{
		inline std::size_t operator()(const std::pair<int, int>& v) const
		{
			return v.first * 31 + v.second;
		}
	};

	Graph m_graphCrest;
	Graph m_graphRiver;
	std::unordered_map<const Edge*, KernelLink> m_edgeLinksKernel;
	std::unordered_map<Kernel*, const Edge*> m_kernelLinksEdge;
	std::vector<Spring> m_springs;
	Vector2 m_oldPos;
	bool m_canMove{false};
	Node* m_nodeToMove{nullptr};
	std::unordered_set<const Node*> m_nodesSelection;
	int m_depth{0};
	bool m_showGraph{true};
	bool m_translateOnly{true};
	bool m_scaleKernel{ true };
	// Spring stiffness
	double m_k{1.};

	GLBuffer m_hfBuffer;

	bool m_influenceEnabled{ false };
	double m_influenceRadius{ 0.2 };
	std::unordered_map<Node*, Vector2> m_movedNodes; // Node, offset

	Vector2 m_lastIntersection;
	float m_blendThreshold{ 1.f };

	static constexpr int m_nodeCircleSegment = 50;
	int m_nbLinesRenderer{0};
	std::unique_ptr<MayaSimpleRendererColors> m_graphRenderer{nullptr};

	void associateKernelsToEdges();
	void createSprings();
	void updateSpringsForce(Node* fixedNode=nullptr) const;
	void moveNode(Node* node, const Vector2& offset);
	void updateKernels();
	void fillHoles();
	void initGraph();

	void fillSelectionNodes(Node* firstNode);

	void recordGraph();

protected:
	void deleteRenderer() override
	{
		if (m_graphRenderer)
		{
			m_graphRenderer->DeleteBuffers();
			m_graphRenderer = nullptr;
		}
		Tool::deleteRenderer();
	}

public:
	ToolGraph(VectorTerrainRaytracingWidget* parent, Kernels& kernels);

	void render() const override
	{
		if (m_graphRenderer) m_graphRenderer->Draw();
		this->Tool::render();
	}

	void mouseMoveEvent(QMouseEvent* e) override;
	Node* getClosestNode(const Vector2& pos);
	void mousePressEvent(QMouseEvent* e) override;
	void mouseReleaseEvent(QMouseEvent* e) override;
	void mouseWheelEvent(QWheelEvent* e) override;
	void keyPressedEvent(QKeyEvent* e) override;

	void setDepth(int val) { m_depth = std::max(0, val); }
	void setShowGraph(bool show);
	void setTranslateOnly(bool translate) { m_translateOnly = translate; }
	void setStiffness(double k) { m_k = k; }
	void setBlendThreshold(double threshold) { m_blendThreshold = threshold; }
	void setScaleKernel(bool scale) { m_scaleKernel = scale; }

	void setInfluenceRegionEnabled(bool enable) { m_influenceEnabled = enable; }

	void updateRenderer();

	ToolType type() override { return ToolType::GRAPH; }
};

