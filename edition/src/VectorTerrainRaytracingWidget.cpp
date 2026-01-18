#include "VectorTerrainRaytracingWidget.h"

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <unistd.h>
#endif

#include <QtGui/QMouseEvent>
#include <QtGui/qpainter.h>
#include <QtWidgets/qfiledialog.h>
#include <QtWidgets/qmessagebox.h>

#include "libs/scalarfield.h"
#include "libs/mayageometry.h"
#include "libs/noise.h"

#include "Tools/ToolEdit.h"
#include "Tools/ToolBrush.h"

#include "Tools/ToolGraph.h"
#include <Kernels/DetailsKernel.h>
#include <Kernels/GaussianKernel.h>


VectorTerrainRaytracingWidget::VectorTerrainRaytracingWidget(): m_accelerationGridShader(), m_rasterizerShader()
{
	// Create the log folder
	std::time_t now = std::time(nullptr);
	std::tm localTime;


	#ifdef _WIN32
    localtime_s(&localTime, &now);
	#else
    localtime_r(&now, &localTime);
	#endif
	char buffer[80];
	std::strftime(buffer, sizeof(buffer), "%d_%m_%Y_%H_%M_%S", &localTime);

	std::string mainName = std::string(SOLUTION_DIR) + "/logs/";


	#ifdef _WIN32
		if (_mkdir(mainName.c_str()) != 0)
	#else
		if (mkdir(mainName.c_str(), 0755) != 0)
	#endif

	{
		std::cerr << "Failed to create the main log folder: " << mainName << ". The folder probably already exists." << std::endl;
	}

	std::string folderName = mainName + buffer;


	#ifdef _WIN32
		if (_mkdir(mainName.c_str()) != 0)
	#else
		if (mkdir(mainName.c_str(), 0755) != 0)
	#endif

	{
		std::cerr << "Failed to create log folder: " << folderName << std::endl;
		folderName = "";
	}

	m_logFolder = QString::fromStdString(folderName);
}

VectorTerrainRaytracingWidget::~VectorTerrainRaytracingWidget()
{
	if (m_influenceRenderer)
	{
		m_influenceRenderer->DeleteBuffers();
	}
}

void VectorTerrainRaytracingWidget::initializeGL()
{
	TerrainRaytracingWidget::initializeGL();

	m_ssbo_primitives.Generate();
	m_gridCellCountsBuffer.Generate();
	m_gridCellMappingsBuffer.Generate();
	m_detailsBuffer.Generate();

	m_gridCellCountsBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(int) * m_gridSize * m_gridSize, nullptr,
	                               GL_STREAM_READ);
	m_gridCellMappingsBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(int) * m_gridSize * m_gridSize * m_maxPerCell,
	                                 nullptr, GL_STREAM_READ);
	loadShader();

	connect(this, SIGNAL(nbPrimitivesChanged(int)), this, SLOT(nbPrimitivesChanged()));

	clear();
}

void VectorTerrainRaytracingWidget::reloadShaders()
{
	TerrainRaytracingWidget::reloadShaders();
	loadShader();
}

void VectorTerrainRaytracingWidget::setNbPrimitives(int val)
{
	m_nbPrimitivesToShow = val;
	glUseProgram(shaderProgram);
	glUniform1i(32, val);

	glUseProgram(m_rasterizerShader);
	glUniform1i(glGetUniformLocation(m_rasterizerShader, "showNbPrimitives"), val);
	glUseProgram(0);

	updateInfluenceRenderers();
	computeAccelerationGrid();
	rasterizePrimitives();

	recordHF("showNbPrimitives");
}

void VectorTerrainRaytracingWidget::showInfluence(const bool& show)
{
	m_showInfluence = show;
	updateInfluenceRenderers();
}

int VectorTerrainRaytracingWidget::getNbPrimitives()
{
	return m_kernels.size();
}

void VectorTerrainRaytracingWidget::loadPrimitivesFile(const QString& filename, const QString& ext)
{
	initHF();
	if (ext == "csv")
		m_kernels.loadCSVFile(filename);
	else if (ext == "npy")
		m_kernels.loadNPYFile(filename);

	emit nbPrimitivesChanged(m_kernels.size());
	rasterizePrimitives();
	resetCam();
	
	getHF();
	m_originalKernels = Kernels(m_kernels);
}

void VectorTerrainRaytracingWidget::saveCSVFile(const QString& filename)
{
	m_kernels.saveCSVFile(filename);
}

void VectorTerrainRaytracingWidget::setRenderResolution(int resolution)
{
	m_hfSize = resolution;
	delete hf;
	hf = nullptr;
	initHF(m_hfSize);
	updatePrimitivesBuffer();
	emit nbPrimitivesChanged(m_kernels.size());
	rasterizePrimitives();

	std::cout << "Set render resolution to " << resolution << std::endl;
	recordHF("new_res");
}

void VectorTerrainRaytracingWidget::addDetailsKernel(const ScalarField2& gt)
{
	setDetails(gt);

	auto box = Box2(Vector2(0), Vector2(m_details.GetSizeX(), m_details.GetSizeY()));
	float detailsKernelSize = 0.04f;
	auto points = box.Poisson((detailsKernelSize / 3.f) * box[1][0], 2000);

	for (const auto& point : points)
	{
		auto pointNorm = point;
		pointNorm -= box[0];
		pointNorm /= (box[1] - box[0]);
		pointNorm = (pointNorm * 2.);
		pointNorm[0] -= 1.;
		pointNorm[1] -= 1.;
		
		m_kernels.add<DetailsKernel>({ detailsKernelSize, detailsKernelSize, 0., 0.5f, static_cast<float>(pointNorm[1]), static_cast<float>(pointNorm[0]),
			(static_cast<float>(pointNorm[0])+1.f)/2.f, (static_cast<float>(pointNorm[1]) + 1.f) / 2.f });
	}

	m_originalKernels = Kernels(m_kernels);

	emit nbPrimitivesChanged(getNbPrimitives());
	rasterizePrimitives();
}

/*
* \brief record the current height field in the log folder
* \param suffix: the suffix to add to the filename
* \param sf: the scalar field to record, if nullptr, the current height field is recorded
*/
void VectorTerrainRaytracingWidget::recordHF(const QString& suffix, const ScalarField2* sf)
{
	if (!m_saveLogs)
		return;

	getHF();

	const auto filename = getRecordName(suffix);

	if (filename.isEmpty())
		return;
	if(sf)
		sf->CreateImage(true).save(filename);
	else
	{
		//hf->CreateImage(true).save(filename);
		hf->CreateImage(0., bbox[1][1], true).save(filename);
		if (m_showInfluence)
		{
			const auto filenameInfluence = getRecordName(suffix + "_influence");
			QGraphicsScene scene;

			const Color colorPos(255, 0, 0, 255);
			const Color colorPosl(2050, 0, 150, 150);
			const Color colorNegl(150, 0, 200, 150);
			const Color colorNeg(0, 0, 255, 255);

			const float maxAmplitude = m_kernels.getMaxAbsAmplitude();

			for (int i = 0; i < m_nbPrimitivesToShow; i++)
			{
				// Alpha value
				float alpha = std::abs(m_kernels[i].amplitude() / maxAmplitude);
				auto color =  ( m_kernels[i].amplitude() > 0. ? colorPos * alpha + (1.-alpha)* colorPosl : colorNeg * alpha + (1. - alpha) * colorNegl);

				auto ellipse = m_kernels[i].getEllipse(1.5f);
				auto pos = ellipse.Center();
				pos[0] += 1.;
				pos[1] += 1.;
				pos = (pos / 2.) * 255.;
				ellipse = Ellipse2(pos, ellipse.A()*128, ellipse.B()*128, ellipse.Axis());
				
				ellipse.Draw(scene, QPen(color.GetQt(), 0.5), QBrush());
			}

			const auto image = Draw::CreateImage(scene, 2048, QRectF(0, 0, 255, 255));
			image.scaled(QSize(2048, 2048)).save(filenameInfluence);
		}
	}
}

/*
* \brief Initialize the height field with a dummy flat scalar field
*/
void VectorTerrainRaytracingWidget::initHF(int size)
{
	if (hf == nullptr)
	{
		// Create a dummy scalar field with hard coded constant to render primitives
		hf = new ScalarField2(Box2(1250.0), size, size, 0.);
		(*hf)[1] = 800;
		UpdateInternal();
		K = 1.;
	}
}

bool VectorTerrainRaytracingWidget::isControlShiftAltPressed(QMouseEvent* e)
{
	return (e->modifiers() & Qt::ControlModifier) || isShiftAltPressed(e);
}

bool VectorTerrainRaytracingWidget::isShiftAltPressed(QMouseEvent* e)
{
	return (e->modifiers() & Qt::ShiftModifier) || (e->modifiers() & Qt::AltModifier);
}

void VectorTerrainRaytracingWidget::resetCam()
{
	auto box = hf->GetBox();
	double a, b;
	hf->GetRange(a, b);
	const auto cam = Camera::View(Box(Vector(box[0][0], box[0][1], a), Vector(box[1][0], box[1][1], b)));
	SetCamera(cam);
}

void VectorTerrainRaytracingWidget::paintGL()
{
	glUseProgram(shaderProgram);

	glEnable(GL_DEPTH_TEST);
	TerrainRaytracingWidget::paintGL();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(camera.Eye()[0], camera.Eye()[1], camera.Eye()[2], camera.At()[0], camera.At()[1], camera.At()[2],
	          camera.Up()[0], camera.Up()[1], camera.Up()[2]);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(Math::RadianToDegree(camera.GetAngleOfViewV(width(), height())),
	               static_cast<GLdouble>(width()) / static_cast<GLdouble>(height()), camera.GetNear(), camera.GetFar());

	if (m_influenceRenderer && m_showInfluence)
		m_influenceRenderer->Draw();
	if (m_currentTool)
		m_currentTool->render();
}

void VectorTerrainRaytracingWidget::setAlbedo(const QImage& img)
{
	m_texture = img;
	SetAlbedo(img);
}

void VectorTerrainRaytracingWidget::updateInfluenceRenderers()
{
	m_influenceRenderer = nullptr;

	if (!m_showInfluence)
		return;

	if (!m_kernels.empty())
	{
		const auto h = HeightField(*hf);
		QVector<VectorFloat> points(m_influenceCircleSegment * m_kernels.size() * 2);
		QVector<ColorFloat> colors(m_influenceCircleSegment * m_kernels.size() * 2);
		auto boxSize = Vector2(h.GetBox().Size());

		const ColorFloat colorNeg(0.2f, 0.2f, 0.8f);
		const ColorFloat colorPos(0.9f, 0.8f, 0.8f);
		const float maxAmplitude = m_kernels.getMaxAbsAmplitude();

		int cpt = 0;


		for (int i = 0; i < m_kernels.size(); i++)
		{
			auto color = m_kernels[i].amplitude() > 0. ? colorPos : colorNeg;

			// Alpha value
			color[3] = std::clamp((std::abs(m_kernels[i].amplitude()) / maxAmplitude)+0.f, 0.4f, 1.f);

			if (i >= m_nbPrimitivesToShow)
				color[3] = 0.;

			auto ellipse = m_kernels[i].getEllipse(1.f);


			for (int j = 1; j < m_influenceCircleSegment; j++)
			{
				Vector2 previousPos = ellipse.Vertex(2. * Math::Pi * static_cast<double>(j - 1) / m_influenceCircleSegment) *
					1250.;
				Vector2 pos = ellipse.Vertex(2. * Math::Pi * static_cast<double>(j) / m_influenceCircleSegment) * 1250.;
				points[cpt] = VectorFloat(h.Vertex(previousPos));
				colors[cpt] = color;
				cpt++;
				points[cpt] = VectorFloat(h.Vertex(pos));
			    colors[cpt] = color;
				cpt++;
			}
			points[cpt] = points[cpt - m_influenceCircleSegment * 2 + 2];
			colors[cpt] = color;
			cpt++;
			points[cpt] = VectorFloat(h.Vertex(ellipse.Vertex(0.) * 1250.));
			colors[cpt] = color;
			cpt++;
		}

		m_influenceRenderer = std::make_unique<MayaSimpleRendererColors>(points, colors, GL_LINES);
		m_influenceRenderer->setDepthTest(false);
		m_influenceRenderer->setLineWidth(1.);
	}
}

/*
* \brief Set the current tool type
*/
void VectorTerrainRaytracingWidget::setTool(const ToolType& type)
{
	switch (type)
	{
	case ToolType::HAND:
		m_currentTool = nullptr;
		break;
	case ToolType::MOVE:
		m_currentTool = std::make_unique<ToolBrush>(this, m_kernels);
		break;
	case ToolType::EDIT:
		m_currentTool = std::make_unique<ToolEdit>(this, m_kernels);
		break;
	case ToolType::GRAPH:
		m_currentTool = std::make_unique<ToolGraph>(this, m_kernels);
		break;
	default: ;
	}
}

QString VectorTerrainRaytracingWidget::getRecordName(const QString& suffix)
{
	if (m_logFolder.isEmpty())
		return "";

	std::time_t now = std::time(nullptr);
	std::tm localTime;
	#ifdef _WIN32
    localtime_s(&localTime, &now);
	#else
    localtime_r(&now, &localTime);
	#endif

	char buffer[80];
	std::strftime(buffer, sizeof(buffer), "%H_%M_%S", &localTime);
	std::string file = buffer;

	const QString filename = m_logFolder + "/" + QString::fromStdString(file) + (suffix.isEmpty() ? "" : "_" + suffix) + ".png";
	std::cout << filename.toStdString() << std::endl;

	return filename;
}

/*
* \brief load the shader
*/
void VectorTerrainRaytracingWidget::loadShader()
{
	QString fullPath = QString::fromStdString(
		std::string(SOLUTION_DIR) + "/shaders/primitives_raytrace.glsl");
	QByteArray ba = fullPath.toLocal8Bit();
	shaderProgram = read_program(ba.data());

	fullPath = QString::fromStdString(std::string(SOLUTION_DIR) + "/shaders/primitives_rasterizer.glsl");
	ba = fullPath.toLocal8Bit();
	m_rasterizerShader = read_program(ba.data());

	fullPath = QString::fromStdString(
		std::string(SOLUTION_DIR) + "/shaders/primitives_acceleration_grid.glsl");
	ba = fullPath.toLocal8Bit();
	m_accelerationGridShader = read_program(ba.data());
}

/*
* \brief Compute the number of primitives to show
*/
void VectorTerrainRaytracingWidget::nbPrimitivesChanged()
{
	if (!m_kernels.empty())
	{
		if (m_influenceRenderer)
		{
			m_influenceRenderer->DeleteBuffers();
			m_influenceRenderer = nullptr;
		}

		m_influenceRenderer = std::make_unique<MayaSimpleRendererColors>(m_influenceCircleSegment * m_kernels.size() * 2, m_influenceCircleSegment * m_kernels.size() * 2, GL_LINES);
		m_influenceRenderer->setDepthTest(false);
		m_influenceRenderer->setLineWidth(1.f);
		updateInfluenceRenderers();

		auto kernelsArray = m_kernels.getArray();
		m_ssbo_primitives.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * kernelsArray.size(), kernelsArray.data(),
		                         GL_STATIC_READ);
		computeAccelerationGrid();
	}
}

void VectorTerrainRaytracingWidget::updatePrimitivesBuffer()
{
	if (!m_kernels.empty())
	{
		auto kernelsArray = m_kernels.getArray();
		m_ssbo_primitives.SetSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(float) * kernelsArray.size(),
		                            kernelsArray.data());
		computeAccelerationGrid();
		updateInfluenceRenderers();
	}
}

void VectorTerrainRaytracingWidget::exportToRes(int res)
{
	const QString filename = QFileDialog::getOpenFileName(this, tr("Open Ground Truth"), QString(),
		tr("Image files(*.jpg, *.png)"));

	if (filename.isEmpty())
	{
		QMessageBox::information(
			this,
			tr("Vector terrain"),
			tr("The export will begin, please do not interact with the ui during the export."));

		delete hf;
		hf = nullptr;
		initHF(res);

		Kernels currentKernels(m_kernels);
		m_kernels = m_originalKernels;

		emit nbPrimitivesChanged(m_kernels.size());
		updatePrimitivesBuffer();
		rasterizePrimitives();

		m_kernels = currentKernels;
		
		emit nbPrimitivesChanged(m_kernels.size());
		updatePrimitivesBuffer();
		rasterizePrimitives();
		recordHF(QString::fromStdString(std::to_string(res)));

		delete hf;
		hf = nullptr;
		initHF();
		emit nbPrimitivesChanged(m_kernels.size());
		updatePrimitivesBuffer();
		rasterizePrimitives();

		QMessageBox::information(
			this,
			tr("Vector terrain"),
			tr("Export finished!"));
		return;
	}
	auto image = QImage(filename);
	if (image.width() != res && image.height() != res)
	{
		switch (QMessageBox::question(
			this,
			tr("Vector terrain"),
			tr("The ground truth you selected is not the same size as the export resolution. Would you like to scale the ground truth?"),

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

	QMessageBox::information(
		this,
		tr("Vector terrain"),
		tr("The export will begin, please do not interact with the ui during the export."));

	auto gt = ScalarField2(Box2(800.0), image);
	//gt.Scale(Vector(1.0, 1.0, 0.001));
	gt.SetRange(0, 1);

	delete hf;
	hf = nullptr;
	initHF(res);
	
	Kernels currentKernels(m_kernels);
	m_kernels = m_originalKernels;
	auto oldDetails = m_details;

	m_details = ScalarField2(m_details, 0.);

	emit nbPrimitivesChanged(m_kernels.size());
	updatePrimitivesBuffer();
	rasterizePrimitives();
	
	setDetails(gt);
	m_kernels = currentKernels;

	emit nbPrimitivesChanged(m_kernels.size());
	updatePrimitivesBuffer();
	rasterizePrimitives();
	recordHF(QString::fromStdString(std::to_string(res)));

	// Reset
	delete hf;
	hf = nullptr;
	initHF(m_hfSize);
	m_details = oldDetails;
	std::vector<float> tmpData(m_details.VertexSize());

	for (int i = 0; i < m_details.VertexSize(); i++)
		tmpData[i] = m_details.at(i);
	m_detailsBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * m_details.VertexSize(), &tmpData.front(), GL_STREAM_READ);
	emit nbPrimitivesChanged(m_kernels.size());
	updatePrimitivesBuffer();
	rasterizePrimitives();

	QMessageBox::information(
		this,
		tr("Vector terrain"),
		tr("Export finished!"));
}

void VectorTerrainRaytracingWidget::computeAccelerationGrid()
{
	if (!m_kernels.empty())
	{
		GLubyte val = 0;
		glClearNamedBufferData(m_gridCellCountsBuffer.GetBuffer(), GL_R32UI, GL_RED_INTEGER, GL_UNSIGNED_BYTE, &val);

		glUseProgram(m_accelerationGridShader);

		glUniform2i(glGetUniformLocation(m_accelerationGridShader, "gridResolution"), m_gridSize, m_gridSize);
		glUniform1i(glGetUniformLocation(m_accelerationGridShader, "maxPerCell"), m_maxPerCell);
		glUniform1i(glGetUniformLocation(m_accelerationGridShader, "primitivesOffset"), utils::sizeKernel());
		m_ssbo_primitives.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
		m_gridCellCountsBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 1);
		m_gridCellMappingsBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 2);

		constexpr int workGroupSize = 64;
		int numWorkGroups = (m_kernels.size() + workGroupSize - 1) / workGroupSize;

		glDispatchCompute(numWorkGroups, 1, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
	}
}

void VectorTerrainRaytracingWidget::rasterizePrimitives()
{
	if (m_kernels.size() > 0)
	{
		std::vector<float> data;
		const int totalBufferSize = hf->GetSizeX() * hf->GetSizeY();
		data.resize(totalBufferSize);

		float zmin = 0., zmax = m_maxRange;
		glUseProgram(m_rasterizerShader);

		hfBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
		m_ssbo_primitives.BindAt(GL_SHADER_STORAGE_BUFFER, 1);

		m_gridCellCountsBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 2);
		m_gridCellMappingsBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 3);
		m_detailsBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 4);

		glUniform2i(glGetUniformLocation(m_rasterizerShader, "gridResolution"), m_gridSize, m_gridSize);
		glUniform1i(glGetUniformLocation(m_rasterizerShader, "maxPerCell"), m_maxPerCell);
		glUniform2i(glGetUniformLocation(m_rasterizerShader, "nxy"), hf->GetSizeX(), hf->GetSizeY());
		glUniform2f(glGetUniformLocation(m_rasterizerShader, "zRange"), zmin, zmax);
		glUniform1f(glGetUniformLocation(m_rasterizerShader, "noiseLevel"), m_noiseLevel);
		glUniform1i(glGetUniformLocation(m_rasterizerShader, "primitiveOffset"), utils::sizeKernel());
		glUniform1i(glGetUniformLocation(m_rasterizerShader, "showNbPrimitives"), m_nbPrimitivesToShow);
		glUniform1i(glGetUniformLocation(m_rasterizerShader, "gaussianID"), static_cast<int>(KernelType::GAUSSIAN));
		glUniform1i(glGetUniformLocation(m_rasterizerShader, "detailsID"), static_cast<int>(KernelType::DETAILS));
		glUniform1i(glGetUniformLocation(m_rasterizerShader, "detailsSize"), m_details.GetSizeX());

		int dispatchSize = (max(hf->GetSizeX(), hf->GetSizeY()) / 8) + 1;

		glDispatchCompute(dispatchSize, dispatchSize, 1);
		glMemoryBarrier(GL_ALL_BARRIER_BITS);

		glUseProgram(0);

		hfBuffer.GetData(data);

		for (int i = 0; i < totalBufferSize; i++)
		{
			(*hf)[i] = static_cast<double>(data[i]);
		}

		K = hf->K();
	}
}

const ScalarField2* VectorTerrainRaytracingWidget::getHF() const
{
	if (hf != nullptr)
	{
		std::vector<float> data;
		const int totalBufferSize = hf->GetSizeX() * hf->GetSizeY();
		data.resize(totalBufferSize);
		hfBuffer.GetData(data);

		for (int i = 0; i < totalBufferSize; i++)
		{
			(*hf)[i] = static_cast<double>(data[i]);
		}
	}

	return hf;
}

void VectorTerrainRaytracingWidget::mouseMoveEvent(QMouseEvent* e)
{
	if (m_currentTool && !isShiftAltPressed(e))
		m_currentTool->mouseMoveEvent(e);

	TerrainRaytracingWidget::mouseMoveEvent(e);
}

void VectorTerrainRaytracingWidget::mousePressEvent(QMouseEvent* e)
{
	if (m_currentTool && !isControlShiftAltPressed(e))
		m_currentTool->mousePressEvent(e);

	TerrainRaytracingWidget::mousePressEvent(e);
}

void VectorTerrainRaytracingWidget::mouseReleaseEvent(QMouseEvent* e)
{
	if (m_currentTool && !isControlShiftAltPressed(e))
		m_currentTool->mouseReleaseEvent(e);

	TerrainRaytracingWidget::mouseReleaseEvent(e);
}

void VectorTerrainRaytracingWidget::wheelEvent(QWheelEvent* e)
{
	if (m_currentTool && (e->modifiers() & Qt::ControlModifier) || (e->modifiers() & Qt::AltModifier))
		m_currentTool->mouseWheelEvent(e);
	else
		TerrainRaytracingWidget::wheelEvent(e);
}

void VectorTerrainRaytracingWidget::keyPressEvent(QKeyEvent* e)
{
	if (m_currentTool)
		m_currentTool->keyPressedEvent(e);

	TerrainRaytracingWidget::keyPressEvent(e);
}

void VectorTerrainRaytracingWidget::setDetails(const ScalarField2& gt)
{
	getHF();
	hf->SetRange(0, 1);

	m_details = gt - *hf;

	recordHF("hf");
	recordHF("GT", &gt);
	recordHF("details", &m_details);

	std::vector<float> tmpData(m_details.VertexSize());

	for (int i = 0; i < m_details.VertexSize(); i++)
		tmpData[i] = m_details.at(i);
	m_detailsBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * m_details.VertexSize(), &tmpData.front(), GL_STREAM_READ);
}

void VectorTerrainRaytracingWidget::saveBrush()
{
	const QString dir("Data/brushes/templates/");
	QDir().mkpath(dir);

	const QString filename = QFileDialog::getSaveFileName(this, tr("Save selected brush"), dir, tr("CSV files(*.csv)"));
	if (m_currentTool && m_currentTool->type() == ToolType::MOVE)
	{
		const auto tool = dynamic_cast<ToolBrush*>(m_currentTool.get());
		tool->saveBrush(filename);
	}
}

void VectorTerrainRaytracingWidget::openBrush(QString filename)
{
	if (filename.isEmpty())
		filename = QFileDialog::getOpenFileName(this, tr("Save selected brush"), QString("Data/brushes/"),
			tr("CSV files(*.csv)"));
	if (!filename.isEmpty() && m_currentTool && m_currentTool->type() == ToolType::MOVE)
	{
		const auto tool = dynamic_cast<ToolBrush*>(m_currentTool.get());
		tool->loadBrush(filename);
	}
}

void VectorTerrainRaytracingWidget::updateBrushThreshold(const int val)
{
	m_brushThreshold = static_cast<float>(val) / 100.f;
	if (m_currentTool)
		m_currentTool->update();
}


void VectorTerrainRaytracingWidget::updateDepthGraphTool(const int val) const
{
	if (m_currentTool && m_currentTool->type() == ToolType::GRAPH)
	{
		const auto tool = dynamic_cast<ToolGraph*>(m_currentTool.get());
		tool->setDepth(val);
	}
}

void VectorTerrainRaytracingWidget::updateStiffnessGraphTool(const int val) const
{
	if (m_currentTool && m_currentTool->type() == ToolType::GRAPH)
	{
		const auto tool = dynamic_cast<ToolGraph*>(m_currentTool.get());
		tool->setStiffness(val / 100.);
	}
}

void VectorTerrainRaytracingWidget::updateBlendGraphTool(int val) const
{
	if (m_currentTool && m_currentTool->type() == ToolType::GRAPH)
	{
		const auto tool = dynamic_cast<ToolGraph*>(m_currentTool.get());
		tool->setBlendThreshold(val / 100.);
	}
}

void VectorTerrainRaytracingWidget::setScaleGraphTool(bool scale) const
{
	if (m_currentTool && m_currentTool->type() == ToolType::GRAPH)
	{
		auto tool = dynamic_cast<ToolGraph*>(m_currentTool.get());
		tool->setScaleKernel(scale);
	}
}

void VectorTerrainRaytracingWidget::setInfluenceRegionGraphTool(bool enable) const
{
	if (m_currentTool && m_currentTool->type() == ToolType::GRAPH)
	{
		auto tool = dynamic_cast<ToolGraph*>(m_currentTool.get());
		tool->setInfluenceRegionEnabled(enable);
	}
}

void VectorTerrainRaytracingWidget::setEditMode(const ToolEdit::Mode& mode) const
{
	if (m_currentTool && m_currentTool->type() == ToolType::EDIT)
	{
		auto tool = dynamic_cast<ToolEdit*>(m_currentTool.get());
		tool->setMode(mode);
	}
}

void VectorTerrainRaytracingWidget::clear()
{
	if (hf == nullptr)
	{
		initHF();
		UseGreenBrownYellowShading(true);
		resetCam();
	}

	m_kernels.clear();
	// Dummy kernel with 0 amplitude to render an empty terrain
	m_kernels.add<GaussianKernel>({ 0.05f, 0.1f, 1.f, 2.5f, 0.f, 0.f, 0.f, 1.f });
	updatePrimitivesBuffer();
	nbPrimitivesChanged(0);
	updateInfluenceRenderers();
	computeAccelerationGrid();
	rasterizePrimitives();
	m_kernels.clear();
}

void VectorTerrainRaytracingWidget::setNoiseLevel(const int val)
{
	m_noiseLevel = static_cast<float>(val) / 500.f;
	rasterizePrimitives();
}

void VectorTerrainRaytracingWidget::updateShowGraphTool(const bool show) const
{
	if (m_currentTool && m_currentTool->type() == ToolType::GRAPH)
	{
		const auto tool = dynamic_cast<ToolGraph*>(m_currentTool.get());
		tool->setShowGraph(show);
	}
}

void VectorTerrainRaytracingWidget::updateTranslateOnlyGraphTool(const bool translate) const
{
	if (m_currentTool && m_currentTool->type() == ToolType::GRAPH)
	{
		const auto tool = dynamic_cast<ToolGraph*>(m_currentTool.get());
		tool->setTranslateOnly(translate);
	}
}
