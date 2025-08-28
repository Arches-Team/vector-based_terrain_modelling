#include "Tools/ToolEdit.h"

#include "libs/evector.h"

#include "VectorTerrainRaytracingWidget.h"


/*
 * Erase primitives that are in the circle of radius m_radius around the intersection point
 */
void ToolEdit::mouseMoveEvent(QMouseEvent* e)
{
	if (m_parent->getHF() != nullptr)
	{
		const HeightField* h = (HeightField*)m_parent->getHF();
		double t;
		const Ray ray = m_parent->ConvertPixelToRay(e->pos());

		// Check if the ray intersects the heightfield
		if (h->Intersect(ray, t, m_intersectionPoint, h->GetBox(), h->K()))
		{
			// If the tool is active
			if (m_canEdit)
			{
				Vector2 boxSize = h->Array2::GetBox().Size();
				// [-1 ; 1]
				auto p = 2 * Vector2(m_intersectionPoint[0] / boxSize[0], m_intersectionPoint[1] / boxSize[1]);

				switch (m_mode)
				{
				case Mode::ERASE:
					erasePrimitives(m_intersectionPoint, boxSize);
					defineRenderer(m_parent->getHF(), m_intersectionPoint);
					break;
				case Mode::AMPLITUDE:
					changeAmplitude(utils::toVec2(e->pos()));
					m_oldMouseScreenPos = utils::toVec2(e->pos());
					break;
				case Mode::AMPLITUDEHR:
					changeAmplitude(utils::toVec2(e->pos()),0.,0.03);
					m_oldMouseScreenPos = utils::toVec2(e->pos());
					break;
				case Mode::AMPLITUDELR:
					changeAmplitude(utils::toVec2(e->pos()), 0.03, 1.);
					m_oldMouseScreenPos = utils::toVec2(e->pos());
					break;
				case Mode::WARP:
					warp(utils::toVec2(e->pos()));
					m_oldMouseScreenPos = utils::toVec2(e->pos());
					break;
				case Mode::MOVE:
					move(p);
					m_oldPos = p;
					break;
				default: ;
				}
			}
			else
			{
				defineRenderer(m_parent->getHF(), m_intersectionPoint);
			}
		}
	}
}

/*
 * Start the editing
 */
void ToolEdit::mousePressEvent(QMouseEvent* e)
{
	if (e->button() == Qt::LeftButton)
	{
		if (m_mode != Mode::ERASE)
		{
			const HeightField* h = (HeightField*)m_parent->getHF();
			double t;
			const Ray ray = m_parent->ConvertPixelToRay(e->pos());

			// Check if the ray intersects the heightfield
			if (h->Intersect(ray, t, m_clickedIntersectionPoint, h->GetBox(), h->K()))
			{
				Vector2 boxSize = h->Array2::GetBox().Size();
				const double radius = 2. * m_radius / boxSize[0];

				// [-1 ; 1]
				const auto p = 2 * Vector2(m_clickedIntersectionPoint[0] / boxSize[0],
				                           m_clickedIntersectionPoint[1] / boxSize[1]);

				for (auto& kernel : m_kernels.getKernels())
				{
					if (kernel->distCenter(p) < radius)
					{
						m_selectedPrimitives.emplace_back(kernel.get());
					}
				}

				m_oldPos = p;
			}
			m_oldMouseScreenPos = utils::toVec2(e->pos());
		}

		m_canEdit = true;

		// Remove the primitives
		mouseMoveEvent(e);
	}
}

/*
 * Stop the editing
 */
void ToolEdit::mouseReleaseEvent(QMouseEvent* e)
{
	if (e->button() == Qt::LeftButton)
	{
		QString mode = "";
		switch (m_mode)
		{
		case Mode::ERASE:
			mode = "edit_erase";
			break;
		case Mode::AMPLITUDE:
			m_selectedPrimitives.clear();
			mode = "edit_amplitude";
			break;
		case Mode::WARP:
			m_selectedPrimitives.clear();
			mode = "edit_warp";
			break;
		case Mode::MOVE:
			m_selectedPrimitives.clear();
			mode = "edit_move";
			break;
		default:;
		}
			
		m_canEdit = false;
        
		m_parent->recordHF(mode);
	}
}

/*!
\brief Delete the primitives that are in the circle of radius m_radius around the intersection point
 */
void ToolEdit::mouseWheelEvent(QWheelEvent* e)
{
	const double scaleRadius = e->angleDelta().y() / 5.;

	m_radius = max(m_radius + scaleRadius, 1.);

	// Delete the circle renderer
	deleteRenderer();

	// Update the radius
	if (m_parent->getHF() != nullptr)
	{
		defineRenderer(m_parent->getHF(), m_intersectionPoint);
	}
}

/*!
\brief Erase the primitives that are in the circle of radius m_radius around the intersection point
*/
void ToolEdit::erasePrimitives(const Vector& intersectionPoint, const Vector2& boxSize) const
{
	const double radius = 2. * m_radius / boxSize[0];

	// [-1 ; 1]
	auto p = 2 * Vector2(intersectionPoint[0] / boxSize[0], intersectionPoint[1] / boxSize[1]);

	const auto sizeBefore = m_kernels.size();
	m_kernels.removeIf([&](const std::unique_ptr<Kernel>& k) { return k->distCenter(p) < radius; });

	if (sizeBefore > m_kernels.size())
	{
		emit m_parent->nbPrimitivesChanged(m_parent->getNbPrimitives());
		m_parent->rasterizePrimitives();
	}
}

void ToolEdit::changeAmplitude(const Vector2& mouseScreenPos, double size_min, double size_max)
{
	double variation =  (m_oldMouseScreenPos[1] - mouseScreenPos[1]) / 100.;

	Vector2 boxSize = (m_parent->getHF())->GetBox().Size();
	// [-1 ; 1]
	Vector2 intersectionVector = 2. * Vector2(m_clickedIntersectionPoint[0] / boxSize[0],
	                                            m_clickedIntersectionPoint[1] / boxSize[1]);
	const double radius = 2. * m_radius / boxSize[0];

	for (Kernel* k : m_selectedPrimitives)
	{
		double size = fabs(0.5 * (k->scaleX() + k->scaleY()));
		if (size > size_min && size < size_max) {
			double dist = Norm(intersectionVector - k->pos()) / radius;
			double falloff = Cubic::SmoothStep(dist, 0., 1.);
			k->amplitude() *= 1. + variation * (1. - falloff);
		}
	}
	m_parent->updatePrimitivesBuffer();
	m_parent->rasterizePrimitives();
}

void ToolEdit::warp(const Vector2& mouseScreenPos)
{
	const auto variation = - (m_oldMouseScreenPos[0] - mouseScreenPos[0]) / 100.;

	const HeightField* h = (HeightField*)m_parent->getHF();
	Vector boxSize = h->GetBox().Size();
	// [-1 ; 1]
	const auto p = 2 * Vector2(m_clickedIntersectionPoint[0] / boxSize[0],
		m_clickedIntersectionPoint[1] / boxSize[1]);
	const double radius = 2. * m_radius / boxSize[0];

	for (auto& k : m_selectedPrimitives)
	{
		Vector2 point(k->posX(), k->posY());
		const double dist = k->distCenter(p);
		const auto currentVariation = (1 - (dist / radius)) * variation;

		k->theta() -= static_cast<float>(currentVariation);

		const double cs = cos(currentVariation);
		const double sn = sin(currentVariation);
		Vector2 pointRotate;
		point -= p;
		pointRotate[0] = point[0] * cs - point[1] * sn;
		pointRotate[1] = point[0] * sn + point[1] * cs;
		pointRotate += p;

		k->posX() = static_cast<float>(pointRotate[0]);
		k->posY() = static_cast<float>(pointRotate[1]);
	}

	m_parent->updatePrimitivesBuffer();
	m_parent->rasterizePrimitives();
}

void ToolEdit::move(const Vector2& intersection)
{
	auto variation = intersection - m_oldPos;

	const HeightField* h = (HeightField*)m_parent->getHF();
	Vector boxSize = h->GetBox().Size();
	const auto p = 2 * Vector2(m_clickedIntersectionPoint[0] / boxSize[0],
		m_clickedIntersectionPoint[1] / boxSize[1]);
	const double radius = 2. * m_radius / boxSize[0];

	for (const auto k : m_selectedPrimitives)
	{
		const auto dist = Norm(p - k->pos()) / radius;
		const auto falloff = 1 - dist;

		k->translate(variation * falloff);

		if(Abs(variation) - Vector(utils::eps) > Vector(0))
		{
			const auto newDist = Norm(p - k->pos()) / radius;
			const auto factor = (dist - newDist + 1.);
			k->scale(variation, factor);
		}
	}
	m_parent->updatePrimitivesBuffer();
	m_parent->rasterizePrimitives();
}

/*!
\brief Define the renderer for the circle around the intersection point
 */
void ToolEdit::defineRenderer(const ScalarField2* hf, const Vector& position)
{
	// Delete the circle renderer
	if (m_renderer)
	{
		m_renderer->DeleteBuffers();
		m_renderer = nullptr;
	}

	// Draw a circle around the intersection point
	const Circle2 circle(position, m_radius);
	QVector<Vector> circlePoints{};
	for (int i = 0; i < 100; i++)
	{
		Vector2 pos = circle.Vertex(2. * Math::Pi * static_cast<double>(i) / 100.);
		circlePoints.append(pos.ToVector(hf->Value(pos)));
	}
	m_renderer = std::make_unique<MayaSimpleRenderer>(circlePoints, Color(0.2, 0.6, 0.2), GL_LINE_LOOP);
	m_renderer->setDepthTest(false);
	m_renderer->setLineWidth(5.);

	m_parent->update();
}
