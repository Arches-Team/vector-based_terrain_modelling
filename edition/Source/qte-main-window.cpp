
#include "qte.h"
#include <QtWidgets/QFileDialog>
#include "heightfieldshader.h"
#include "graph.h"
#include "Tools/ToolEdit.h"
#include "Eigen/Dense"
#include <math.h>
#include <tuple>
#include <Kernels/GaussianKernel.h>
#include <QtWidgets/qmessagebox.h>
#include <direct.h>

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
    uiw.check_saveLogs->setVisible(false);
    uiw.btn_curveTool->setVisible(false);
    uiw.btn_AdrienGraph->setVisible(false);
    uiw.btn_test->setVisible(false);


    // Connecting QT Actions
    createActions();

    setAcceptDrops(true);

    // Fill brush template
    refreshTemplateBrush();
}

/*!
\brief Filter drag event.
*/
void QteWindow::dragEnterEvent(QDragEnterEvent* e)
{
    if (e->mimeData()->hasUrls())
    {
        e->acceptProposedAction();
    }
}

/*!
\brief Display the heightfield.
*/
void QteWindow::DisplayHeightfield(bool setCamera)
{
    MayaGeometry mg("Heightfield", m_hf.CreateMesh());
    HeightFieldShader shader(m_hf);
    auto texture = shader.ShadedRelief();

    raywidget->SetHeightField(&m_hf);
    //raywidget->UseGreenBrownYellowShading(true);
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

/*!
* \brief Move tool
*/
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

void QteWindow::curveTool()
{
    enableAllTools();
    raywidget->setTool(ToolType::CURVE);
    uiw.btn_curveTool->setEnabled(false);
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

void QteWindow::setRender128()
{
    raywidget->setRenderResolution(128);
}

void QteWindow::setRender256()
{
    raywidget->setRenderResolution(256);
}

void QteWindow::setRender512()
{
    raywidget->setRenderResolution(512);
}

void QteWindow::setRender1024()
{
    raywidget->setRenderResolution(1024);
}

void QteWindow::setRender2048()
{
    raywidget->setRenderResolution(2048);
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
    connect(uiw.btn_curveTool, SIGNAL(clicked()), this, SLOT(curveTool()));
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
    //connect(uiw.btn_loadTemplateBrush, SIGNAL(clicked()), this, SLOT(loadTemplateBrush()));
    //connect(uiw.btn_refreshTemplate, SIGNAL(clicked()), this, SLOT(refreshTemplateBrush()));

    // Graph
    connect(uiw.slider_depthGraph, SIGNAL(valueChanged(int)), raywidget.get(), SLOT(updateDepthGraphTool(int)));
		connect(raywidget.get(), SIGNAL(updateDepthGraph(int)), uiw.slider_depthGraph, SLOT(setValue(int)));

    //connect(uiw.check_showGraph, SIGNAL(clicked(bool)), raywidget.get(), SLOT(updateShowGraphTool(bool)));
    //connect(uiw.check_fillHolesGraph, SIGNAL(clicked(bool)), raywidget.get(), SLOT(setFillHolesGraphTool(bool)));
    //connect(uiw.check_scale, SIGNAL(clicked(bool)), raywidget.get(), SLOT(setScaleGraphTool(bool)));
    connect(uiw.check_influence, SIGNAL(clicked(bool)), raywidget.get(), SLOT(setInfluenceRegionGraphTool(bool)));
    //connect(uiw.check_translateOnlyGraph, SIGNAL(clicked(bool)), raywidget.get(), SLOT(updateTranslateOnlyGraphTool(bool)));
    //connect(uiw.slider_stiffnessGraph, SIGNAL(valueChanged(int)), raywidget.get(), SLOT(updateStiffnessGraphTool(int)));
    connect(uiw.slider_blendThresholdGraph, SIGNAL(valueChanged(int)), raywidget.get(), SLOT(updateBlendGraphTool(int)));

    // Edit
    connect(uiw.radio_erase, SIGNAL(clicked()), this, SLOT(setEditErase()));
    connect(uiw.radio_amplitude, SIGNAL(clicked()), this, SLOT(setEditAmplitude()));
    connect(uiw.radio_amplitude_lr, SIGNAL(clicked()), this, SLOT(setEditAmplitudeLR()));
    connect(uiw.radio_amplitude_hr, SIGNAL(clicked()), this, SLOT(setEditAmplitudeHR()));
    connect(uiw.radio_warp, SIGNAL(clicked()), this, SLOT(setEditWarp()));
    connect(uiw.radio_move, SIGNAL(clicked()), this, SLOT(setEditMove()));

    // Debug
    connect(uiw.btn_AdrienGraph, SIGNAL(clicked()), this, SLOT(adrienComputeGraph()));
    connect(uiw.btn_test, SIGNAL(clicked()), this, SLOT(testButton()));
    connect(uiw.check_saveLogs, SIGNAL(clicked(bool)), raywidget.get(), SLOT(setSaveLogs(bool)));

    // Render resolution
    connect(uiw.action128x128, SIGNAL(triggered()), this, SLOT(setRender128()));
    connect(uiw.action256x256, SIGNAL(triggered()), this, SLOT(setRender256()));
    connect(uiw.action512x512, SIGNAL(triggered()), this, SLOT(setRender512()));
    connect(uiw.action1024x1024, SIGNAL(triggered()), this, SLOT(setRender1024()));
    connect(uiw.action2048x2048, SIGNAL(triggered()), this, SLOT(setRender2048()));
}

/*!
* \brief Enable all tools
*/
void QteWindow::enableAllTools()
{
    uiw.btn_handTool->setEnabled(true);
    uiw.btn_moveTool->setEnabled(true);
    uiw.btn_eraseTool->setEnabled(true);
    uiw.btn_curveTool->setEnabled(true);
    uiw.btn_graphTool->setEnabled(true);

}

/*!
\brief Action binded to the button
*/
void QteWindow::OpenHeightfield()
{
    const QString filename = QFileDialog::getOpenFileName(this, tr("Open Heightfield"), QString(),
                                                          tr("Image files(*.jpg, *.png)"));

    m_hf = ScalarField2(Box2(800.0), QImage(filename));
    m_hf.Scale(Vector(1.0, 1.0, 0.001));

    DisplayHeightfield(true);
}

/*!
* \brief Open a gaussian file
*/
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
    //raywidget->UseElevationShading(true);
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
    //gt.Scale(Vector(1.0, 1.0, 0.001));
    gt.SetRange(0, 1);

    raywidget->addDetailsKernel(gt);
}

void QteWindow::ExportHighRes()
{
    int id = uiw.combo_exporthighres->currentIndex();
    int res = std::pow(2, (id + 7));
    raywidget->exportToRes(res);
}

void QteWindow::adrienComputeGraph()
{
    const ScalarField2* sf = raywidget->getHF();

    if (sf == nullptr)
    {
        raywidget->openGaussiansCSVFile(
            QString::fromStdString(std::string(SOLUTION_DIR) + "/data/cuda/alaska_clean.csv"));
        sf = raywidget->getHF();
        raywidget->UseElevationShading(true);
    }

    const auto hf = HeightField(*sf);

    Graph crestGraph, riverGraph;

    constexpr double distance = 0.1;
    GraphControl::createCrestRiverGraph(crestGraph, riverGraph, hf);
    riverGraph.print(QString::fromStdString(std::string(SOLUTION_DIR) + "/tmp/testRiver.pdf"));
    riverGraph.reduceGraph(distance);
    riverGraph.print(QString::fromStdString(std::string(SOLUTION_DIR) + "/tmp/testReduceRiver.pdf"));
    crestGraph.print(QString::fromStdString(std::string(SOLUTION_DIR) + "/tmp/testCrest.pdf"));
    crestGraph.reduceGraph(distance);
    crestGraph.print(QString::fromStdString(std::string(SOLUTION_DIR) + "/tmp/testReduceCrest.pdf"));

    auto points = crestGraph.getNodesPosition();
    const auto mesh = Mesh2::Delaunay(points);
    std::cout << "Graph size " << crestGraph.size() << " mesh points " << mesh.VertexSize() << " mesh indexes " << mesh.
        IndexSize() << " mesh triangle " << mesh.TriangleSize() << std::endl;
}

void QteWindow::testButton()
{
    std::cout << "Test button" << std::endl;
    //Kernel kernel(1, 0.5, 1., 0., 1., 0., 0., 1.);
    //// Points to test
    //const std::vector testPoints = {
    //    Vector2(-1, -1),
    //    Vector2(-0.5, 0),
    //    Vector2(0, -0.5),
    //    Vector2(0.5, 0.5),
    //    Vector2(1, 1),
    //    Vector2(-1, 1),
    //    Vector2(1, -1),
    //    Vector2(0, 0),
    //    Vector2(0.5, 0),
    //    Vector2(0, 0.25)
    //};

    //// Run the tests
    //for (const auto& point : testPoints)
    //{
    //    std::cout << "Point(" << point[0] << "," << point[1] << ") inside ellipse: "
    //        << kernel.isInside(point) << std::endl;
    //}

    // Draw the graph in a QImage with a qgraphicscene
    //Ellipse2 ellipse(5, 10);

    //QPen pen;
    //pen.setWidth(1);
    //QBrush brush;
    //QGraphicsScene scene;
    //ellipse.Draw(scene, pen, brush);

    //ellipse = ellipse.Scaled(1.);
    ////ellipse = ellipse.Rotated(45);
    //pen.setColor(Qt::blue);
    //ellipse.Draw(scene, pen, brush);

    //

    //pen.setColor(Qt::red);
    ////ellipse.Draw(scene, pen, brush);

    //Draw::CreateVector(scene, QString::fromStdString(std::string(SOLUTION_DIR) + "/tmp/ellipse.pdf"));

    //float a = 2.f;
    //float b = 1.f;
    //float theta = utils::pi / 4.f;
    //Eigen::Vector2f axis(0, 1);

    //auto angleVector = [](Eigen::Vector2f x, Eigen::Vector2f y) -> float
    //    {
    //        return atan2(y[1] * x[0] - y[0] * x[1], y[0] * x[0] + y[1] * x[1]);
    //    };

    //auto x = [](float a, float b, float theta, float t) -> float
    //    {
    //        return a * cos(t) * cos(theta) - b * sin(t) * sin(theta);
    //    };
    //auto y = [](float a, float b, float theta, float t) -> float
    //    {
    //        return a * cos(t) * sin(theta) + b * sin(t) * cos(theta);
    //    };

    //auto ellipse_parameters = [](double A, double B, double C, double D, double E, double F) -> std::tuple<double, double, double>
    //    {
    //        double theta = 0.5 * atan2(-B, C - A);
    //        double denominator = B * B - 4 * A * C;

    //        double sqrt_part = sqrt((A - C) * (A - C) + B * B);

    //        double a = -sqrt(2 * ((denominator * F) * ((A + C) - sqrt_part))) / denominator;
    //        double b = -sqrt(2 * ((denominator * F) * ((A + C) + sqrt_part))) / denominator;

    //        return std::make_tuple(a, b, theta);
    //    };

    //
    //Eigen::Vector2f xAxis(0, 1);
    //float factor = 1.05f;
    //float angle = angleVector(xAxis, axis);
    //theta -= angle;

    //Eigen::MatrixXf A(5, 6);

    //constexpr int n = 5;
    //for (int i = 0; i < n; i++)
    //{
    //    float t = i * 2 * utils::pi / n;
    //    float current_x = x(a, b, theta, t);
    //    float current_y = y(a, b, theta, t);
    //    current_x *= factor;
    //    A(i, 0) = current_x * current_x;
    //    A(i, 1) = current_x * current_y;
    //    A(i, 2) = current_y * current_y;
    //    A(i, 3) = current_x;
    //    A(i, 4) = current_y;
    //    A(i, 5) = 1.0;
    //}

    //Eigen::MatrixXf ATA = A.transpose() * A;
    //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver(ATA);

    //// Extract the eigenvector associated with the minimum eigenvalue
    //Eigen::VectorXf z = solver.eigenvectors().col(0);

    //// Calculate ellipse parameters
    //auto [recoved_a, recoved_b, recoved_theta] = ellipse_parameters(z(0), z(1), z(2), z(3), z(4), z(5));

    //recoved_theta += angle;

    //std::cout << "Recovered param: " << recoved_a << ", " << recoved_b << ", " << recoved_theta << std::endl;
    //std::cout << "Real param: " << a << ", " << b << ", " << theta << std::endl;

    //auto k = GaussianKernel({ 1.f, 3.f, 1.f, float(M_PI)/2.f + float(M_PI)*0.01f, 1.f, 0.f, 0.f, 1.f});
    /*auto k = GaussianKernel({ 0.005f, 0.1f, 1.f, 0.f, 0.1f, -0.4f, 0.f, 1.f });
    std::cout << "Before " << k << std::endl;
    Vector2 axis(1, 0);
    float factor{ 0.5 };
    k.scale(axis, factor);
    std::cout << "Scaling of " << factor << " on axis " << axis << std::endl;
    std::cout << "After " << k << std::endl;
    */

    //auto k = GaussianKernel({ 1.f, 3.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f });
    //auto e = k.getEllipse(1.0);

    //std::cout << "Initial Kernel: " << k << "\nInitial Ellipse: " << e << std::endl;

    //// Define the scaling axes and factors to test.
    //std::vector<Vector2> axes = { Vector2(1, 0), Vector2(0, 1), Vector2(1, 1), Vector2(-1, 0) };
    //std::vector<float> factors = { 0.5, 1.0, 1.5, 2.0 };

    //std::cout << e.Scaled(Vector2(1, 1), false) << std::endl;
    //std::cout << e.Scaled(Vector2(0.5, 1), false) << std::endl;
    //std::cout << e.Scaled(Vector2(0.5, 0.5), false) << std::endl;

    //auto rotate = [](Vector2 x, double angle) -> Vector2
    //{
    //    float cosTheta = std::cos(angle);
    //    float sinTheta = std::sin(angle);
    //    return Vector2(x[0] * cosTheta - x[1] * sinTheta, x[0] * sinTheta + x[1] * cosTheta);
    //};

    //auto angleVector = [](Vector2 x, Vector2 y) -> double
    //{
    //    return atan2(y[1] * x[0] - y[0] * x[1], y[0] * x[0] + y[1] * x[1]);
    //};

    //// Loop over each axis and factor.
    //for (size_t i = 0; i < axes.size(); ++i) {
    //    for (size_t j = 0; j < factors.size(); ++j) {
    //        Vector2 axis = axes[i];
    //        float factor = factors[j];

    //        // Scale the kernel and the ellipse.
    //        k.scale(axis, factor);
    //        axis = Normalized(axis);
    //        auto angle = angleVector(Vector2::X, axis);
    //        axis = rotate(axis, -angle);
    //        axis *= factor;
    //        axis = rotate(axis, angle);

    //        float ux = 1.0 + (factor - 1.0) * std::abs(axis[0]);
    //        float uy = 1.0 + (factor - 1.0) * std::abs(axis[1]);

    //        auto scaledEllipse = e.Scaled(Vector2(1 / ux, 1 / uy), false);

    //        // Output the results.
    //        std::cout << "\nTest " << (i * factors.size() + j + 1) << ": " << std::endl;
    //        std::cout << "Scaling of " << factor << " on axis " << axis << std::endl;
    //        std::cout << "U " << Vector2(ux, uy) << std::endl;
    //        std::cout << "After scaling - Kernel: " << k << "\nEllipse: " << scaledEllipse << std::endl;

    //        // Reset kernel and ellipse for the next test.
    //        k = GaussianKernel({ 1.f, 3.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f });
    //        e = k.getEllipse(1.0);
    //    }
    //}

    /*if (!raywidget->createDummyTerrain())
    {
        raywidget->normalizeKernels();    
        raywidget->updateInfluenceRenderers();
        raywidget->rasterizeGaussians();
    }*/
    //raywidget->export2K();
    const std::vector<int> resolutions{ 256, 512, 1024, 2048, 4096 };
    for(const auto res: resolutions)
    {
        raywidget->setRenderResolution(res);
        const int frames = 100;
        std::cout << "Running for " << raywidget->getRenderResolution() << "x" << raywidget->getRenderResolution() << " resolution." << std::endl;
        std::vector<std::string> files{ "500.npy", "1k.npy", "2k.npy", "4k.npy", "7k.npy", "14k.npy", "20k.npy", "30k.npy", "60k.npy" };
        for (const auto& f : files)
        {
            raywidget->openGaussiansNPYFile(
                QString::fromStdString(std::string(SOLUTION_DIR) + "/Data/timing/" + f));

            std::cout << f << " primitives" << std::endl;
            raywidget->timing(frames);
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    
    
}
