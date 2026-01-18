#pragma once

#include "libs/realtime.h"

#include "Tools/Tool.h"
#include "Tools/ToolEdit.h"
#include "Kernels.h"
#include "MayaSimpleRendererColors.h"


class VectorTerrainRaytracingWidget : public TerrainRaytracingWidget
{
	Q_OBJECT

public:
	VectorTerrainRaytracingWidget();
	~VectorTerrainRaytracingWidget() override;
	void initializeGL() override;
	void reloadShaders() override;
	void setNbPrimitives(int val);
	void showInfluence(const bool& show);
	int getNbPrimitives();

	// Accept ext = {"npy", "cvs"}
	void loadPrimitivesFile(const QString& filename, const QString& ext);
	void saveCSVFile(const QString& filename);

	void setRenderResolution(int resolution);
	int getRenderResolution() const { return m_hfSize; }
	
	void addDetailsKernel(const ScalarField2& gt);

	QString getRecordName(const QString& suffix = "");
	void recordHF(const QString& suffix = "", const ScalarField2* sf = nullptr);
	bool saveLogs() const { return m_saveLogs; }
	void exportToRes(int res);

	void paintGL() override;

	void setAlbedo(const QImage&);
	void updateInfluenceRenderers();
	void rasterizePrimitives();

	const ScalarField2* getHF() const;
	void setTool(const ToolType& type);
	void applyTool() const { m_currentTool->apply(); }
	void updatePrimitivesBuffer();
	GLBuffer& getHFBuffer() { return hfBuffer; };

	float getBrushThreshold() const { return m_brushThreshold; }

private:
	void loadShader();

	void resetCam();
	void initHF(int size = 256);
	void computeAccelerationGrid();

	static bool isControlShiftAltPressed(QMouseEvent* e);
	static bool isShiftAltPressed(QMouseEvent* e);

	GLBuffer m_ssbo_primitives;
	static constexpr int m_gridSize = 25;
	int m_hfSize{ 256 };
	// TODO: compute dynamically maxPerCell
	static constexpr int m_maxPerCell = 2000;
	int m_maxRange{ 200 };
	GLBuffer m_gridCellCountsBuffer;
	GLBuffer m_gridCellMappingsBuffer;

	ScalarField2 m_details;
	ScalarField2 m_OriginalDetails;
	GLBuffer m_detailsBuffer;
	void setDetails(const ScalarField2& gt);
	Kernels m_originalKernels;

	std::unique_ptr<MayaSimpleRendererColors> m_influenceRenderer{nullptr};

	GLuint m_rasterizerShader;
	GLuint m_accelerationGridShader;

	QImage m_texture;

	Kernels m_kernels;

	std::unique_ptr<Tool> m_currentTool{};

	bool m_showInfluence{false};
	int m_nbPrimitivesToShow{1};
	float m_brushThreshold{0.2f};

	static constexpr int m_influenceCircleSegment = 50;

	float m_noiseLevel{1.};

	QString m_logFolder;

	bool m_saveLogs{ true };

signals:
	void nbPrimitivesChanged(int newVal);
	void updateDepthGraph(int val);

public slots:
	void mouseMoveEvent(QMouseEvent*) override;
	void mousePressEvent(QMouseEvent*) override;
	void mouseReleaseEvent(QMouseEvent* e) override;
	void wheelEvent(QWheelEvent* e) override;
	void keyPressEvent(QKeyEvent*) override;

	void saveBrush();
	void openBrush(QString filename="");
	void updateBrushThreshold(int val);
	void clear();

	void updateDepthGraphTool(int val) const;
	void updateShowGraphTool(bool show) const;
	void updateTranslateOnlyGraphTool(bool translate) const;
	void updateStiffnessGraphTool(int val) const;
	void updateBlendGraphTool(int val) const;
	void setScaleGraphTool(bool scale) const;
	void setInfluenceRegionGraphTool(bool enable) const;

	void setEditMode(const ToolEdit::Mode& mode) const;

	void setSaveLogs(bool save) { m_saveLogs = save; }

	void setNoiseLevel(int val);
	void setMaxRange(int maxRange)
	{
	    m_maxRange = maxRange;
	    rasterizePrimitives();
	}
	void nbPrimitivesChanged();
};

