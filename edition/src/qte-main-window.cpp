#include "qte.h"


#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <unistd.h>
#endif

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
    m_uiw.setupUi(this);

    m_uiw.btn_handTool->setEnabled(false);
    m_uiw.toolOptionsWidget->setCurrentIndex(0);

    // Ray tracing widget
    {
        const Camera camera(Vector(-10.0, -10.0, 10.0));
        m_raywidget = std::make_unique<VectorTerrainRaytracingWidget>();
        auto* GLlayout = new QGridLayout;
        GLlayout->addWidget(m_raywidget.get(), 0, 0);
        GLlayout->setContentsMargins(0, 0, 0, 0);
        m_uiw.raytracingwidget->setLayout(GLlayout);
        m_raywidget->SetCamera(camera);
    }

    m_uiw.slider_nbprimitives->setMaximum(m_raywidget->getNbPrimitives());

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

void QteWindow::displayHeightfield(bool setCamera)
{
    MayaGeometry mg("Heightfield", m_hf.CreateMesh());
    HeightFieldShader shader(m_hf);
    auto texture = shader.ShadedRelief();

    m_raywidget->SetHeightField(&m_hf);
    m_raywidget->UseElevationShading(true);

    if (setCamera)
    {
        auto cam = Camera::View(mg.GetBox());
        m_raywidget->SetCamera(cam);
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
    m_raywidget->SetSteps(val);
}

void QteWindow::updateRayEps(const int val)
{
    m_raywidget->SetEpsilon(static_cast<float>(val) / 1000.f);
}

void QteWindow::updateNbPrimitives(const int val)
{
    m_raywidget->setNbPrimitives(val);
    m_uiw.label_nbprimitives->setText(QString("%1 primitives").arg(val));
}

void QteWindow::setMaxPrimitives(const int val)
{
    m_raywidget->setNbPrimitives(val);
    m_uiw.label_nbprimitives->setText(QString("%1 primitives").arg(val));
    m_uiw.slider_nbprimitives->setRange(1, val);
    m_uiw.slider_nbprimitives->setValue(val);
}

void QteWindow::updateShowInfluence(const bool show)
{
    m_raywidget->showInfluence(show);
}

void QteWindow::handTool()
{
    enableAllTools();
    m_raywidget->setTool(ToolType::HAND);
    m_uiw.btn_handTool->setEnabled(false);
    m_uiw.toolOptionsWidget->setCurrentIndex(0);
}

void QteWindow::moveTool()
{
    enableAllTools();
    m_raywidget->setTool(ToolType::MOVE);
    m_uiw.btn_moveTool->setEnabled(false);
    m_uiw.toolOptionsWidget->setCurrentIndex(1);
}

void QteWindow::eraseTool()
{
    enableAllTools();
    m_raywidget->setTool(ToolType::EDIT);
    m_uiw.btn_eraseTool->setEnabled(false);
    m_uiw.toolOptionsWidget->setCurrentIndex(3);
    m_uiw.radio_erase->setChecked(true);
    m_uiw.radio_amplitude->setChecked(false);
}

void QteWindow::graphTool()
{
    enableAllTools();
    m_raywidget->setTool(ToolType::GRAPH);
    m_uiw.btn_graphTool->setEnabled(false);
    m_uiw.toolOptionsWidget->setCurrentIndex(2);
}

void QteWindow::applyTool()
{
    m_raywidget->applyTool();
}

void QteWindow::loadTemplateBrush()
{
    auto selectedBrush = m_uiw.list_brushes->currentItem()->text();
    if (selectedBrush.isEmpty())
        return;

    m_raywidget->openBrush(m_templateBrushes[selectedBrush]);
}

void QteWindow::refreshTemplateBrush()
{
    m_templateBrushes.clear();
    m_uiw.list_brushes->clear();

    auto brushFolder = std::string(SOLUTION_DIR) + "/data/brushes/templates/";


    #ifdef _WIN32
    _mkdir(brushFolder.c_str());
    #else
    mkdir(brushFolder.c_str(), 0755);
    #endif

    // Add all csv files in the template folder
    for (const auto& entry : std::filesystem::directory_iterator(brushFolder)) {
        if (entry.is_regular_file() && entry.path().extension() == ".csv") {
            std::string fileName = entry.path().stem().string();
            QString filePath = QString::fromStdString(entry.path().string());
            m_templateBrushes[QString::fromStdString(fileName)] = filePath;
        }
    }

    for (const auto& kv : m_templateBrushes)
        m_uiw.list_brushes->addItem(kv.first);

}

void QteWindow::setEditErase()
{
    m_raywidget->setEditMode(ToolEdit::Mode::ERASE);
}

void QteWindow::setEditAmplitude()
{
    m_raywidget->setEditMode(ToolEdit::Mode::AMPLITUDE);
}

void QteWindow::setEditAmplitudeLR()
{
    m_raywidget->setEditMode(ToolEdit::Mode::AMPLITUDELR);
}

void QteWindow::setEditAmplitudeHR()
{
    m_raywidget->setEditMode(ToolEdit::Mode::AMPLITUDEHR);
}

void QteWindow::setEditWarp()
{
    m_raywidget->setEditMode(ToolEdit::Mode::WARP);
}

void QteWindow::setEditMove()
{
    m_raywidget->setEditMode(ToolEdit::Mode::MOVE);
}

void QteWindow::setRenderResolution(int resolution)
{
    m_raywidget->setRenderResolution(resolution);
}


/*!
\brief Create callbacks between member slots and user interface.
*/
void QteWindow::createActions()
{
    //File
    connect(m_uiw.actionExit, SIGNAL(triggered()), this, SLOT(close()));

    // Button connections
    connect(m_uiw.openPrimitivesFile, SIGNAL(clicked()), this, SLOT(openPrimitivesFile()));
    connect(m_uiw.btn_savePrimitives, SIGNAL(clicked()), this, SLOT(savePrimitivesFile()));
    connect(m_uiw.btn_loadShader, SIGNAL(clicked()), this, SLOT(reloadShader()));
    connect(m_uiw.btn_openTerrain, SIGNAL(clicked()), this, SLOT(openHeightfield()));
    connect(m_uiw.btn_addDetails, SIGNAL(clicked()), this, SLOT(openGroundTruth()));
    connect(m_uiw.btn_exporthighres, SIGNAL(clicked()), this, SLOT(exportHighRes()));
    connect(m_uiw.btn_clear, SIGNAL(clicked()), m_raywidget.get(), SLOT(clear()));

    // Tools
    connect(m_uiw.btn_handTool, SIGNAL(clicked()), this, SLOT(handTool()));
    connect(m_uiw.btn_moveTool, SIGNAL(clicked()), this, SLOT(moveTool()));
    connect(m_uiw.btn_eraseTool, SIGNAL(clicked()), this, SLOT(eraseTool()));
    connect(m_uiw.btn_graphTool, SIGNAL(clicked()), this, SLOT(graphTool()));
    connect(m_uiw.btn_applyTool, SIGNAL(clicked()), this, SLOT(applyTool()));

    connect(m_raywidget.get(), SIGNAL(nbPrimitivesChanged(int)), this, SLOT(setMaxPrimitives(int)));

    connect(m_uiw.slider_nbprimitives, SIGNAL(valueChanged(int)), this, SLOT(updateNbPrimitives(int)));
    connect(m_uiw.check_showInfluence, SIGNAL(clicked(bool)), this, SLOT(updateShowInfluence(bool)));

    connect(m_uiw.slider_noiseLevel, SIGNAL(valueChanged(int)), m_raywidget.get(), SLOT(setNoiseLevel(int)));
    connect(m_uiw.slider_amplitude, SIGNAL(valueChanged(int)), m_raywidget.get(), SLOT(setMaxRange(int)));

    // Brush
    connect(m_uiw.btn_saveBrush, SIGNAL(clicked()), m_raywidget.get(), SLOT(saveBrush()));
    connect(m_uiw.btn_saveBrush, SIGNAL(clicked()), this, SLOT(refreshTemplateBrush()));
    connect(m_uiw.btn_openBrush, SIGNAL(clicked()), m_raywidget.get(), SLOT(openBrush()));
    connect(m_uiw.slider_brushthreshold, SIGNAL(valueChanged(int)), m_raywidget.get(), SLOT(updateBrushThreshold(int)));
    connect(m_uiw.list_brushes, SIGNAL(itemDoubleClicked(QListWidgetItem*)), this, SLOT(loadTemplateBrush()));

    // Graph
    connect(m_uiw.slider_depthGraph, SIGNAL(valueChanged(int)), m_raywidget.get(), SLOT(updateDepthGraphTool(int)));
		connect(m_raywidget.get(), SIGNAL(updateDepthGraph(int)), m_uiw.slider_depthGraph, SLOT(setValue(int)));

    connect(m_uiw.check_influence, SIGNAL(clicked(bool)), m_raywidget.get(), SLOT(setInfluenceRegionGraphTool(bool)));
    connect(m_uiw.slider_blendThresholdGraph, SIGNAL(valueChanged(int)), m_raywidget.get(), SLOT(updateBlendGraphTool(int)));

    // Edit
    connect(m_uiw.radio_erase, SIGNAL(clicked()), this, SLOT(setEditErase()));
    connect(m_uiw.radio_amplitude, SIGNAL(clicked()), this, SLOT(setEditAmplitude()));
    connect(m_uiw.radio_amplitude_lr, SIGNAL(clicked()), this, SLOT(setEditAmplitudeLR()));
    connect(m_uiw.radio_amplitude_hr, SIGNAL(clicked()), this, SLOT(setEditAmplitudeHR()));
    connect(m_uiw.radio_warp, SIGNAL(clicked()), this, SLOT(setEditWarp()));
    connect(m_uiw.radio_move, SIGNAL(clicked()), this, SLOT(setEditMove()));

    // Debug
    connect(m_uiw.check_saveLogs, SIGNAL(clicked(bool)), m_raywidget.get(), SLOT(setSaveLogs(bool)));

    // Render resolution
    connect(m_uiw.action128x128, &QAction::triggered, this, [this]{setRenderResolution(128); });
    connect(m_uiw.action256x256, &QAction::triggered, this, [this]{setRenderResolution(256); });
    connect(m_uiw.action512x512, &QAction::triggered, this, [this]{setRenderResolution(512); });
    connect(m_uiw.action1024x1024, &QAction::triggered, this, [this]{setRenderResolution(1024); });
    connect(m_uiw.action2048x2048, &QAction::triggered, this, [this]{setRenderResolution(2048); });
}

void QteWindow::enableAllTools()
{
    m_uiw.btn_handTool->setEnabled(true);
    m_uiw.btn_moveTool->setEnabled(true);
    m_uiw.btn_eraseTool->setEnabled(true);
    m_uiw.btn_graphTool->setEnabled(true);
}

void QteWindow::openHeightfield()
{
    const QString filename = QFileDialog::getOpenFileName(this, tr("Open Heightfield"), QString(),
                                                          tr("Image files(*.jpg, *.png)"));

    m_hf = ScalarField2(Box2(800.0), QImage(filename));
    m_hf.Scale(Vector(1.0, 1.0, 0.001));

    displayHeightfield(true);
}

void QteWindow::openPrimitivesFile()
{
    const QString filename = QFileDialog::getOpenFileName(this, tr("Open Primitives file"),
                                                          QString::fromStdString(
                                                              std::string(SOLUTION_DIR) + "/data/"),
                                                          tr("Primitives files(*.csv *.npy)"));
    const QFileInfo fileInfo{filename};
    const auto ext = fileInfo.suffix();
    m_raywidget->loadPrimitivesFile(filename, ext);

    m_raywidget->UseGreenBrownYellowShading(true);
    emit m_raywidget->nbPrimitivesChanged(m_raywidget->getNbPrimitives());

    m_raywidget->recordHF("hf");
}

void QteWindow::savePrimitivesFile()
{
    const QString filename = QFileDialog::getSaveFileName(this, tr("Save Primitives CSV file"),
        QString::fromStdString(
            std::string(SOLUTION_DIR) + "/data/"),
        tr("Primitives files(*.csv)"));

    m_raywidget->saveCSVFile(filename);
}

void QteWindow::reloadShader()
{
    m_raywidget->reloadShaders();
    updateNbPrimitives(m_uiw.slider_nbprimitives->value());
}

void QteWindow::openGroundTruth()
{
    const QString filename = QFileDialog::getOpenFileName(this, tr("Open Heightfield"), QString(),
        tr("Image files(*.jpg, *.png)"));

    if (filename.isEmpty())
        return;

    QImage image(filename);
    const int res = m_raywidget->getRenderResolution();

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

    m_raywidget->addDetailsKernel(gt);
}

void QteWindow::exportHighRes()
{
    int id = m_uiw.combo_exporthighres->currentIndex();
    int res = std::pow(2, (id + 7));
    m_raywidget->exportToRes(res);
}
