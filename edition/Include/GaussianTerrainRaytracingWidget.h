#pragma once
#include "libs/realtime.h"
#include "libs/mayarender.h"

#include "Tools/Tool.h"
#include "Tools/ToolEdit.h"
#include "Kernels.h"
#include "MayaSimpleRendererColors.h"


class GaussianTerrainRaytracingWidget : public TerrainRaytracingWidget
{
	Q_OBJECT

public:
	GaussianTerrainRaytracingWidget();
	~GaussianTerrainRaytracingWidget();
	void initializeGL() override;
	void ReloadShaders() override;
	void SetNbGaussians(int val);
	void useGaussians(const bool& use);
	void showInfluence(const bool& show);
	int getNbGaussians();
	void openGaussiansCSVFile(const QString& filename);
	void openGaussiansNPYFile(const QString& filename);
	void saveGaussiansCSVFile(const QString& filename);

	void setRenderResolution(int resolution);
	int getRenderResolution() const { return m_hfSize; }
	
	void addDetailsKernel(const ScalarField2& gt);

	QString getRecordName(const QString& suffix = "");
	void recordHF(const QString& suffix = "", const ScalarField2* sf = nullptr);
	void exportToRes(int res);

	void paintGL() override;

	void SetAlbedo(const QImage&);
	void updateInfluenceRenderers();
	void rasterizeGaussians();

	const ScalarField2* getHF() const;
	void setTool(const ToolType& type);
	void applyTool() const { m_currentTool->apply(); }
	void UpdateGaussiansBuffer();
	GLBuffer& getHFBuffer() { return hfBuffer; };

	// TODO: move in ToolBrush
	float getBrushThreshold() const { return m_brushThreshold; }

	//Debug
	void normalizeKernels();
	bool createDummyTerrain();	
	void timing(int frames);


private:
	void loadShader();

	void resetCam();
	void initHF(int size = 256);
	void computeAccelerationGrid();

	static bool isControlShiftAltPressed(QMouseEvent* e);
	static bool isShiftAltPressed(QMouseEvent* e);

	GLBuffer m_ssbo_gaussians;
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
	int m_nbGaussiansToShow{1};
	float m_brushThreshold{0.2f};

	static constexpr int m_influenceCircleSegment = 50;

	float m_noiseLevel{1.};

	QString m_logFolder;

	bool m_saveLogs{ true };

signals:
	void nbGaussiansChanged(int newVal);
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
	    rasterizeGaussians();
	}
	void nbGaussiansChanged();
};

