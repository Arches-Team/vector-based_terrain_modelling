#pragma once

#include "libs/heightfield.h"
#include "libs/maya.h"

#include "ui_main.h"
#include "VectorTerrainRaytracingWidget.h"

class QteWindow : public QMainWindow
{
  Q_OBJECT

private:
  Ui::Assets m_uiw; //!< Interface : QtDesigner.
  std::unique_ptr<VectorTerrainRaytracingWidget> m_raywidget;
  HeightField m_hf;

  std::map<QString, QString> m_templateBrushes;

public:
  QteWindow();
  void displayHeightfield(bool setCamera = false);

private:
  void createActions();
  void enableAllTools();

public slots:
  void openHeightfield();
  void openPrimitivesFile();
  void savePrimitivesFile();
  void reloadShader();
  void openGroundTruth();
  void exportHighRes();

  void dragEnterEvent(QDragEnterEvent*) override;
  void dropEvent(QDropEvent*) override;

  void updateRayStep(int val);
  void updateRayEps(int val);
  void updateNbPrimitives(int val);
  void setMaxPrimitives(int val);
  void updateShowInfluence(bool show);

  void handTool();
  void moveTool();
  void eraseTool();
  void graphTool();
  void applyTool();

  void loadTemplateBrush();
  void refreshTemplateBrush();

  // Tool Edit
  void setEditErase();
  void setEditAmplitude();
  void setEditAmplitudeLR();
  void setEditAmplitudeHR();
  void setEditWarp();
  void setEditMove();

  void setRenderResolution(int resolution);
};
