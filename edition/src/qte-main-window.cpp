#include "qte.h"

#include <direct.h>
#include <QtWidgets/qmessagebox.h>
#include <QtWidgets/QFileDialog>

#include "Eigen/Dense"
#include "libs/heightfieldshader.h"

#include "graph.h"
#include "Tools/ToolEdit.h"




/*!
\class QteWindow qte.h
\brief This class implements a main window.
*/

/*!
\brief Create the main window.
*/
QteWindow::QteWindow()
{
    // Loading interface
    uiw.setupUi(this);

    uiw.btn_handTool->setEnabled(false);
    uiw.toolOptionsWidget->setCurrentIndex(0);

    // Ray tracing widget
    {
        const Camera camera(Vector(-10.0, -10.0, 10.0));
        raywidget = std::make_unique<GaussianTerrainRaytracingWidget>();
        auto* GLlayout = new QGridLayout;
        GLlayout->addWidget(raywidget.get(), 0, 0);
        GLlayout->setContentsMargins(0, 0, 0, 0);
        uiw.raytracingwidget->setLayout(GLlayout);
        raywidget->SetCamera(camera);
    }

    uiw.slider_nbgaussians->setMaximum(raywidget->getNbGaussians());

    createActions();
    setAcceptDrops(true);

    refreshTemplateBrush();
}

void QteWindow::dragEnterEvent(QDragEnterEvent* e)
{
    if (e->mimeData()->hasUrls())
    {
        e->acceptProposedAction();
    }
}

void QteWindow::DisplayHeightfield(bool setCamera)
{
    MayaGeometry mg("Heightfield", m_hf.CreateMesh());
    HeightFieldShader shader(m_hf);
    auto texture = shader.ShadedRelief();

    raywidget->SetHeightField(&m_hf);
    raywidget->UseElevationShading(true);

    if (setCamera)
    {
        auto cam = Camera::View(mg.GetBox());
        raywidget->SetCamera(cam);
    }
}

void QteWindow::dropEvent(QDropEvent* e)
{
    foreach(const QUrl & url, e->mimeData()->urls())
    {
        QString fileName = url.toLocalFile();
        m_hf = ScalarField2(Box2(1250.0), QImage(fileName));
        m_hf.Scale(Vector(1.0, 1.0, 0.005));
    }
}

void QteWindow::updateRayStep(const int val)
{
    raywidget->SetSteps(val);
}

void QteWindow::updateRayEps(const int val)
{
    raywidget->SetEpsilon(static_cast<float>(val) / 1000.f);
}

void QteWindow::updateNbGaussians(const int val)
{
    raywidget->SetNbGaussians(val);
    uiw.label_nbgaussians->setText(QString("%1 primitives").arg(val));
}

void QteWindow::setMaxGaussians(const int val)
{
    raywidget->SetNbGaussians(val);
    uiw.label_nbgaussians->setText(QString("%1 primitives").arg(val));
    uiw.slider_nbgaussians->setRange(1, val);
    uiw.slider_nbgaussians->setValue(val);
}

void QteWindow::updateShowInfluence(const bool show)
{
    raywidget->showInfluence(show);
}

void QteWindow::handTool()
{
    enableAllTools();
    raywidget->setTool(ToolType::HAND);
    uiw.btn_handTool->setEnabled(false);
    uiw.toolOptionsWidget->setCurrentIndex(0);
}

void QteWindow::moveTool()
{
    enableAllTools();
    raywidget->setTool(ToolType::MOVE);
    uiw.btn_moveTool->setEnabled(false);
    uiw.toolOptionsWidget->setCurrentIndex(1);
}

void QteWindow::eraseTool()
{
    enableAllTools();
    raywidget->setTool(ToolType::EDIT);
    uiw.btn_eraseTool->setEnabled(false);
    uiw.toolOptionsWidget->setCurrentIndex(3);
    uiw.radio_erase->setChecked(true);
    uiw.radio_amplitude->setChecked(false);
}

void QteWindow::graphTool()
{
    enableAllTools();
    raywidget->setTool(ToolType::GRAPH);
    uiw.btn_graphTool->setEnabled(false);
    uiw.toolOptionsWidget->setCurrentIndex(2);
}

void QteWindow::applyTool()
{
    raywidget->applyTool();
}

void QteWindow::loadTemplateBrush()
{
    auto selectedBrush = uiw.list_brushes->currentItem()->text();
    if (selectedBrush.isEmpty())
        return;

    raywidget->openBrush(m_templateBrushes[selectedBrush]);
}

void QteWindow::refreshTemplateBrush()
{
    m_templateBrushes.clear();
    uiw.list_brushes->clear();

    auto brushFolder = std::string(SOLUTION_DIR) + "/data/brushes/templates/";

    _mkdir(brushFolder.c_str());

    // Add all csv files in the template folder
    for (const auto& entry : std::filesystem::directory_iterator(brushFolder)) {
        if (entry.is_regular_file() && entry.path().extension() == ".csv") {
            std::string fileName = entry.path().stem().string();
            QString filePath = QString::fromStdString(entry.path().string());
            m_templateBrushes[QString::fromStdString(fileName)] = filePath;
        }
    }

    for (const auto& kv : m_templateBrushes)
        uiw.list_brushes->addItem(kv.first);

}

void QteWindow::setEditErase()
{
    raywidget->setEditMode(ToolEdit::Mode::ERASE);
}

void QteWindow::setEditAmplitude()
{
    raywidget->setEditMode(ToolEdit::Mode::AMPLITUDE);
}

void QteWindow::setEditAmplitudeLR()
{
    raywidget->setEditMode(ToolEdit::Mode::AMPLITUDELR);
}

void QteWindow::setEditAmplitudeHR()
{
    raywidget->setEditMode(ToolEdit::Mode::AMPLITUDEHR);
}

void QteWindow::setEditWarp()
{
    raywidget->setEditMode(ToolEdit::Mode::WARP);
}

void QteWindow::setEditMove()
{
    raywidget->setEditMode(ToolEdit::Mode::MOVE);
}

void QteWindow::setRenderResolution(int resolution)
{
    raywidget->setRenderResolution(resolution);
}


/*!
\brief Create callbacks between member slots and user interface.
*/
void QteWindow::createActions()
{
    //File
    connect(uiw.actionExit, SIGNAL(triggered()), this, SLOT(close()));

    // Button connections
    connect(uiw.openGaussianFile, SIGNAL(clicked()), this, SLOT(OpenGaussianFile()));
    connect(uiw.btn_saveGaussian, SIGNAL(clicked()), this, SLOT(SaveGaussianFile()));
    connect(uiw.btn_loadShader, SIGNAL(clicked()), this, SLOT(reloadShader()));
    connect(uiw.btn_openTerrain, SIGNAL(clicked()), this, SLOT(OpenHeightfield()));
    connect(uiw.btn_addDetails, SIGNAL(clicked()), this, SLOT(OpenGroundTruth()));
    connect(uiw.btn_exporthighres, SIGNAL(clicked()), this, SLOT(ExportHighRes()));
    connect(uiw.btn_clear, SIGNAL(clicked()), raywidget.get(), SLOT(clear()));

    // Tools
    connect(uiw.btn_handTool, SIGNAL(clicked()), this, SLOT(handTool()));
    connect(uiw.btn_moveTool, SIGNAL(clicked()), this, SLOT(moveTool()));
    connect(uiw.btn_eraseTool, SIGNAL(clicked()), this, SLOT(eraseTool()));
    connect(uiw.btn_graphTool, SIGNAL(clicked()), this, SLOT(graphTool()));
    connect(uiw.btn_applyTool, SIGNAL(clicked()), this, SLOT(applyTool()));

    connect(raywidget.get(), SIGNAL(nbGaussiansChanged(int)), this, SLOT(setMaxGaussians(int)));

    connect(uiw.slider_nbgaussians, SIGNAL(valueChanged(int)), this, SLOT(updateNbGaussians(int)));
    connect(uiw.check_showInfluence, SIGNAL(clicked(bool)), this, SLOT(updateShowInfluence(bool)));

    connect(uiw.slider_noiseLevel, SIGNAL(valueChanged(int)), raywidget.get(), SLOT(setNoiseLevel(int)));
    connect(uiw.slider_amplitude, SIGNAL(valueChanged(int)), raywidget.get(), SLOT(setMaxRange(int)));

    // Brush
    connect(uiw.btn_saveBrush, SIGNAL(clicked()), raywidget.get(), SLOT(saveBrush()));
    connect(uiw.btn_saveBrush, SIGNAL(clicked()), this, SLOT(refreshTemplateBrush()));
    connect(uiw.btn_openBrush, SIGNAL(clicked()), raywidget.get(), SLOT(openBrush()));
    connect(uiw.slider_brushthreshold, SIGNAL(valueChanged(int)), raywidget.get(), SLOT(updateBrushThreshold(int)));
    connect(uiw.list_brushes, SIGNAL(itemDoubleClicked(QListWidgetItem*)), this, SLOT(loadTemplateBrush()));

    // Graph
    connect(uiw.slider_depthGraph, SIGNAL(valueChanged(int)), raywidget.get(), SLOT(updateDepthGraphTool(int)));
		connect(raywidget.get(), SIGNAL(updateDepthGraph(int)), uiw.slider_depthGraph, SLOT(setValue(int)));

    connect(uiw.check_influence, SIGNAL(clicked(bool)), raywidget.get(), SLOT(setInfluenceRegionGraphTool(bool)));
    connect(uiw.slider_blendThresholdGraph, SIGNAL(valueChanged(int)), raywidget.get(), SLOT(updateBlendGraphTool(int)));

    // Edit
    connect(uiw.radio_erase, SIGNAL(clicked()), this, SLOT(setEditErase()));
    connect(uiw.radio_amplitude, SIGNAL(clicked()), this, SLOT(setEditAmplitude()));
    connect(uiw.radio_amplitude_lr, SIGNAL(clicked()), this, SLOT(setEditAmplitudeLR()));
    connect(uiw.radio_amplitude_hr, SIGNAL(clicked()), this, SLOT(setEditAmplitudeHR()));
    connect(uiw.radio_warp, SIGNAL(clicked()), this, SLOT(setEditWarp()));
    connect(uiw.radio_move, SIGNAL(clicked()), this, SLOT(setEditMove()));

    // Debug
    connect(uiw.check_saveLogs, SIGNAL(clicked(bool)), raywidget.get(), SLOT(setSaveLogs(bool)));

    // Render resolution
    connect(uiw.action128x128, &QAction::triggered, this, [this]{setRenderResolution(128); });
    connect(uiw.action256x256, &QAction::triggered, this, [this]{setRenderResolution(256); });
    connect(uiw.action512x512, &QAction::triggered, this, [this]{setRenderResolution(512); });
    connect(uiw.action1024x1024, &QAction::triggered, this, [this]{setRenderResolution(1024); });
    connect(uiw.action2048x2048, &QAction::triggered, this, [this]{setRenderResolution(2048); });
}

void QteWindow::enableAllTools()
{
    uiw.btn_handTool->setEnabled(true);
    uiw.btn_moveTool->setEnabled(true);
    uiw.btn_eraseTool->setEnabled(true);
    uiw.btn_graphTool->setEnabled(true);
}

void QteWindow::OpenHeightfield()
{
    const QString filename = QFileDialog::getOpenFileName(this, tr("Open Heightfield"), QString(),
                                                          tr("Image files(*.jpg, *.png)"));

    m_hf = ScalarField2(Box2(800.0), QImage(filename));
    m_hf.Scale(Vector(1.0, 1.0, 0.001));

    DisplayHeightfield(true);
}

void QteWindow::OpenGaussianFile()
{
    const QString filename = QFileDialog::getOpenFileName(this, tr("Open Gaussian file"),
                                                          QString::fromStdString(
                                                              std::string(SOLUTION_DIR) + "/data/"),
                                                          tr("Gaussians files(*.csv *.npy)"));
    const QFileInfo fileInfo{filename};
    const auto ext = fileInfo.suffix();
    if (ext == "csv")
        raywidget->openGaussiansCSVFile(filename);
    else if (ext == "npy")
        raywidget->openGaussiansNPYFile(filename);

    raywidget->UseGreenBrownYellowShading(true);
    emit raywidget->nbGaussiansChanged(raywidget->getNbGaussians());

    raywidget->recordHF("hf");
}

void QteWindow::SaveGaussianFile()
{
    const QString filename = QFileDialog::getSaveFileName(this, tr("Save Gaussian CSV file"),
        QString::fromStdString(
            std::string(SOLUTION_DIR) + "/data/"),
        tr("Gaussians files(*.csv)"));

    raywidget->saveGaussiansCSVFile(filename);
}

void QteWindow::reloadShader()
{
    raywidget->ReloadShaders();
    updateNbGaussians(uiw.slider_nbgaussians->value());
}

void QteWindow::OpenGroundTruth()
{
    const QString filename = QFileDialog::getOpenFileName(this, tr("Open Heightfield"), QString(),
        tr("Image files(*.jpg, *.png)"));

    if (filename.isEmpty())
        return;

    QImage image(filename);
    const int res = raywidget->getRenderResolution();

    if (image.width() != res && image.height() != res)
    {
        switch (QMessageBox::question(
            this,
            tr("Vector terrain"),
            tr("The ground truth you selected is not the same size as the current resolution. Would you like to scale the ground truth?"),

            QMessageBox::Yes |
            QMessageBox::Cancel,

            QMessageBox::Cancel))
        {
        case QMessageBox::Yes:
            image = image.scaled(res, res, Qt::KeepAspectRatio, Qt::SmoothTransformation);
            break;
        case QMessageBox::Cancel:
            return;
            break;
        default:
            return;
            break;
        }
    }

    auto gt = ScalarField2(Box2(800.0), image);
    gt.SetRange(0, 1);

    raywidget->addDetailsKernel(gt);
}

void QteWindow::ExportHighRes()
{
    int id = uiw.combo_exporthighres->currentIndex();
    int res = std::pow(2, (id + 7));
    raywidget->exportToRes(res);
}
