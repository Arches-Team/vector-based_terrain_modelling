#ifndef __Qte__
#define __Qte__

#include "libs/heightfield.h"
#include "libs/maya.h"

#include "ui_main.h"
#include "GaussianTerrainRaytracingWidget.h"

class QteWindow : public QMainWindow
{
  Q_OBJECT

private:
  Ui::Assets uiw; //!< Interface : QtDesigner.
  std::unique_ptr<GaussianTerrainRaytracingWidget> raywidget;
  HeightField m_hf;

  std::map<QString, QString> m_templateBrushes;

public:
  QteWindow();
  void DisplayHeightfield(bool setCamera = false);

private:
  void createActions();
  void enableAllTools();

public slots:
  void OpenHeightfield();
  void OpenGaussianFile();
  void SaveGaussianFile();
  void reloadShader();
  void OpenGroundTruth();
  void ExportHighRes();

  void dragEnterEvent(QDragEnterEvent*) override;
  void dropEvent(QDropEvent*) override;

  void updateRayStep(int val);
  void updateRayEps(int val);
  void updateNbGaussians(int val);
  void setMaxGaussians(int val);
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

#endif

