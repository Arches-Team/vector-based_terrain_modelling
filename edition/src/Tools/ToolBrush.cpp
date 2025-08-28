#include "Tools/ToolBrush.h"

#include "libs/heightfield.h"

#include "VectorTerrainRaytracingWidget.h"
#include "utils.h"
#include "Kernels/GaussianKernel.h"
#include "Kernels/FactoryKernel.h"

void ToolBrush::moveGaussians(const Vector2& dest)
{
	const auto offset = utils::toVec2(QPointF(dest[0], dest[1]) - m_oldPos);
	m_lasso.Translate(offset);
	for (auto& k : m_selectedKernels)
	{
		k->translate(offset);
	}
	updateShowKernels();
	m_parent->updatePrimitivesBuffer();
	m_parent->rasterizePrimitives();
}

void ToolBrush::addToSelection()
{
	if (!m_selectedKernels.empty())
	{
		m_kernels.clear(m_selectedKernels);
		m_selectedKernels.clear();
	}

	const int currentSize = m_kernels.size();
	for (int i = 0; i < currentSize; ++i)
	{
		Kernel& k = m_kernels[i];

		if (m_lasso.Inside(Vector2(k.posX(), k.posY())))
		{
			m_selectedKernels.emplace_back(m_kernels.add(k.clone()));
		}
	}
	if (m_kernels.size() != currentSize)
	{
		emit m_parent->nbPrimitivesChanged(m_parent->getNbPrimitives());
		updateShowKernels();
		m_parent->rasterizePrimitives();
	}
}

void ToolBrush::updateRenderer()
{
	deleteRenderer();

	auto boxSize = m_parent->getHF()->GetBox().Size();

	QVector<Vector> points{};
	for (int i = 0; i < m_lasso.Size(); ++i)
	{
		const auto& p = m_lasso.Vertex(i);
		auto pNorm = Vector2(p[0] * boxSize[0], p[1] * boxSize[1]) / 2.;
		points.append(Vector(pNorm[0], pNorm[1], m_parent->getHF()->Value(pNorm)));
	}

	m_renderer = std::make_unique<MayaSimpleRenderer>(points, Color(0.6, 0.2, 0.2), GL_LINE_LOOP);
	m_renderer->setDepthTest(false);
	m_renderer->setLineWidth(5.);
}

void ToolBrush::updateShowKernels() const
{
	for (auto& k : m_kernels)
	{
		k->setModAmplitude(1.);
		const auto dist = m_lasso.Signed(Vector2(k->posX(), k->posY()));
		if (dist < 0)
		{
			k->setModAmplitude(static_cast<float>(1. - Linear::Step(-dist, 0., m_parent->getBrushThreshold())));
		}
	}

	for (auto& k : m_selectedKernels)
	{
		k->setModAmplitude(1.);
	}
}

void ToolBrush::saveBrush(const QString& filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
	{
		std::cout << file.errorString().toStdString() << std::endl;
		return;
	}
	QTextStream out(&file);

	bool isBeginning = true;
	for (auto& k : m_selectedKernels)
	{
		auto data = k->get();
		for (const float& d : data)
		{
			if (!isBeginning)
				out << ",";
			isBeginning = false;
			out << QString::number(d);
		}
	}
	out << "\n";

	for (int i = 0; i < m_lasso.Size(); ++i)
	{
		if (i != 0)
			out << ",";
		out << QString::number(m_lasso.Vertex(i)[0]) << "," << QString::number(m_lasso.Vertex(i)[1]);
	}

	file.close();

}

void ToolBrush::loadBrush(const QString& filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly))
	{
		std::cout << file.errorString().toStdString() << std::endl;
		return;
	}

	QString data = file.readAll();
	auto lines = data.split('\n');

	auto gaussiansStr = lines[0].split(',');
	auto lassoStr = lines[1].split(',');

	file.close();

	apply();
	m_state = State::MOVE;

	int i = 0;
	int id = -1;
	std::vector<float> params;

	for (const auto& gStr : gaussiansStr)
	{
		float current = std::stof(gStr.toStdString());
		if (i == 0)
			id = static_cast<int>(current);
		else
			params.push_back(current);
		i++;

		if (i == utils::sizeKernel())
		{
			m_selectedKernels.emplace_back(m_kernels.add(kernel::create(static_cast<KernelType>(id), params)));
			i = 0;
			params.clear();
		}
	}

	i = 0;
	Vector2 point(0);
	for (const auto& lStr : lassoStr)
	{
		point[i] = std::stof(lStr.toStdString());
		i++;
		if (i == 2)
		{
			i = 0;
			m_lasso.Append(point);
		}
	}

	update();
}

void ToolBrush::reset()
{
	deleteRenderer();
	m_selectedKernels.clear();
	m_lasso = Polygon2();
	m_state = State::START;
}

void ToolBrush::apply()
{
	reset();

	auto sizeBefore = m_kernels.size();

	m_kernels.removeIf(
		[](const std::unique_ptr<Kernel>& k) { return !k->isShown(); }
	);

	if (sizeBefore > m_kernels.size())
	{
		emit m_parent->nbPrimitivesChanged(m_parent->getNbPrimitives());
		m_parent->rasterizePrimitives();
	}
}

void ToolBrush::update()
{
	updateShowKernels();
	updateRenderer();
	emit m_parent->nbPrimitivesChanged(m_parent->getNbPrimitives());
	m_parent->rasterizePrimitives();
}

void ToolBrush::keyPressedEvent(QKeyEvent* e)
{
	switch (e->key())
	{
	case Qt::Key_Delete:
		if (m_state == State::MOVE)
		{
			const auto sizeBefore = m_kernels.size();
			m_kernels.clear(m_selectedKernels);
			if (sizeBefore > m_kernels.size())
			{
				emit m_parent->nbPrimitivesChanged(m_parent->getNbPrimitives());
				m_parent->rasterizePrimitives();
			}

			reset();
		}
		break;
	default: ;
	}
}

void ToolBrush::mouseMoveEvent(QMouseEvent* e)
{
	if (m_parent->getHF() != nullptr)
	{
		const HeightField* h = (HeightField*)m_parent->getHF();
		Vector intersectionPoint;
		double t;
		const Ray ray = m_parent->ConvertPixelToRay(e->pos());

		if (h->Intersect(ray, t, intersectionPoint, h->GetBox(), h->K()))
		{
			auto boxSize = h->GetBox().Size();
			const auto intersectionNormalized = 2 * QPointF(intersectionPoint[0] / boxSize[0],
			                                                intersectionPoint[1] / boxSize[1]);
			if (m_state == State::LASSO)
				m_lasso.Append(utils::toVec2(intersectionNormalized));

			if (m_state == State::MOVE && m_canMove)
				moveGaussians(Vector2(intersectionNormalized.x(), intersectionNormalized.y()));

			updateRenderer();

			m_oldPos = intersectionNormalized;
		}
	}
}

void ToolBrush::mousePressEvent(QMouseEvent* e)
{
	if (m_parent->getHF() != nullptr)
	{
		const HeightField* h = (HeightField*)m_parent->getHF();
		Vector intersectionPoint;
		double t;
		const Ray ray = m_parent->ConvertPixelToRay(e->pos());

		if (h->Intersect(ray, t, intersectionPoint, h->GetBox(), h->K()))
		{
			if (m_state == State::START)
			{
				m_state = State::LASSO;
				mouseMoveEvent(e);
			}
			if (m_state == State::MOVE)
			{
				m_canMove = true;
			}
		}
	}
}

void ToolBrush::mouseReleaseEvent(QMouseEvent* e)
{
	if (m_state == State::LASSO)
	{
		m_state = State::MOVE;
		addToSelection();
	}
	if (m_state == State::MOVE)
	{
		m_canMove = false;
		m_parent->recordHF("brush_move");
	}
		
}

void ToolBrush::mouseWheelEvent(QWheelEvent* e)
{
	if (m_state == State::MOVE)
	{
		const bool isRotation = e->modifiers() & Qt::AltModifier;


		double val = -e->angleDelta().y() * 0.001;

		if (!isRotation)
			val += 1.;

		const auto center = m_lasso.Center();

		m_lasso.Translate(-center);
		if (isRotation)
			m_lasso.Rotate(Matrix2::Rotation(val));
		else
			m_lasso.Scale(val);
		m_lasso.Translate(center);

		for (auto& k : m_selectedKernels)
		{
			if (isRotation)
				k->rotate(static_cast<float>(val), center);
			else
				k->scale(static_cast<float>(val), center);
		}

		updateShowKernels();
		updateRenderer();
		m_parent->updatePrimitivesBuffer();
		m_parent->rasterizePrimitives();

		m_parent->recordHF("brush_" + isRotation ? "rotation" : "scale");
	}
}
