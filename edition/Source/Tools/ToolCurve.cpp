#include "Tools/ToolCurve.h"

#include "GaussianTerrainRaytracingWidget.h"
#include "libs/heightfield.h"
#include "libs/evector.h"


void ToolCurve::moveGaussians(double distMax) const
{
	for (int i = 0; i < m_kernels.size(); i++)
	{
		Kernel& k = m_kernels[i];
		const float posX = k.posX();
		const float posY = k.posY();
		double u, v;
		Vector2 p(posX, posY);
		m_curveFrom.UV(p, u, v);

		if (m_lasso.containsPoint(QPointF(posX, posY), Qt::OddEvenFill))
		{
			if (u == 0 || u == 1)
			{
				auto normal = Normalized(m_curveFrom.Tangent(u).Orthogonal());
				// TODO: make rotation as utils function
				const double angle = normal.Angle(Normalized(p - m_curveFrom.Eval(u)));
				const double cs = cos(angle);
				const double sn = sin(angle);
				auto normalTo = Normalized(m_curveTo.Tangent(u).Orthogonal());
				Vector2 nRotate;
				nRotate[0] = normalTo[0] * cs - normalTo[1] * sn;
				nRotate[1] = normalTo[0] * sn + normalTo[1] * cs;

				if (i == 0)
				{
					std::cout << "Angle " << angle << std::endl << "Normal " << normalTo << " nRotate " << nRotate <<
						std::endl;
				}
			}
			else
			{
				Vector2 newPos;
				auto normalTo = Normalized(m_curveTo.Tangent(u).Orthogonal());
				auto normalFrom = Normalized(m_curveFrom.Tangent(u).Orthogonal());
				const auto angle = normalTo.Angle(normalFrom);
				newPos = m_curveTo.Eval(u) + normalTo * v;

				k.theta() += static_cast<float>(angle);
				k.posX() = static_cast<float>(newPos[0]);
				k.posY() = static_cast<float>(newPos[1]);
			}
		}
	}
}

void ToolCurve::deleteRenderer(std::unique_ptr<MayaSimpleRenderer>& renderer)
{
	if (renderer)
	{
		renderer->DeleteBuffers();
		renderer = nullptr;
	}
}

void ToolCurve::mouseMoveEvent(QMouseEvent* e)
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

			// TODO: Refactor renderer creation: lot (a lot) of copy paste
			switch (m_state)
			{
			case State::FIRST_CURVE:
				if (!m_curvePointsFrom.empty())
				{
					auto back = m_curvePointsFrom.back();
					const auto p1 = Vector2(back[0] * boxSize[0], back[1] * boxSize[1]) / 2.;
					const Vector2 p2(intersectionPoint);

					deleteRenderer(m_renderer);

					QVector<Vector> points{};
					constexpr int nbPoints = 50;

					if (m_curvePointsFrom.size() == 1)
					{
						for (int i = 0; i < nbPoints; i++)
						{
							double l = static_cast<double>(i) / nbPoints;
							points.append(h->Vertex(p1 * (1 - l) + p2 * l));
						}
					}

					if (m_curvePointsFrom.size() == 2)
					{
						const auto pTemp = 2 * Vector2(intersectionPoint[0] / boxSize[0],
						                               intersectionPoint[1] / boxSize[1]);
						const QuadricCurve2 c =
							QuadricCurve2::Bezier(m_curvePointsFrom[0], m_curvePointsFrom[1], pTemp);

						for (int i = 0; i < nbPoints; i++)
						{
							auto tmp = c.Eval(1. - (static_cast<double>(i) / nbPoints));
							tmp = Vector2(tmp[0] * boxSize[0], tmp[1] * boxSize[1]) / 2.;
							points.append(h->Vertex(tmp));
						}
					}

					m_renderer = std::make_unique<MayaSimpleRenderer>(points, Color(0.2, 0.6, 0.2), GL_LINE_STRIP);
					m_renderer->setDepthTest(false);
					m_renderer->setLineWidth(5.);
				}
				break;

			case State::LASSO:
				{
					if (m_editLasso)
					{
						deleteRenderer(m_rendererLasso);
						m_lasso.append(
							2 * QPointF(intersectionPoint[0] / boxSize[0], intersectionPoint[1] / boxSize[1]));
						QVector<Vector> points{};
						for (const auto& p : m_lasso)
						{
							points.append(h->Vertex(Vector2(p.x() * boxSize[0], p.y() * boxSize[0]) / 2.));
						}

						m_rendererLasso = std::make_unique<MayaSimpleRenderer>(
							points, Color(0.6, 0.2, 0.2), GL_LINE_LOOP);
						m_rendererLasso->setDepthTest(false);
						m_rendererLasso->setLineWidth(5.);
					}
				}
				break;
			case State::SECOND_CURVE:
				if (!m_curvePointsTo.empty())
				{
					auto back = m_curvePointsTo.back();
					const auto p1 = Vector2(back[0] * boxSize[0], back[1] * boxSize[1]) / 2.;
					const Vector2 p2(intersectionPoint);

					deleteRenderer(m_rendererCurveTo);

					QVector<Vector> points{};
					constexpr int nbPoints = 50;

					if (m_curvePointsTo.size() == 1)
					{
						for (int i = 0; i < nbPoints; i++)
						{
							double l = static_cast<double>(i) / nbPoints;
							points.append(h->Vertex(p1 * (1 - l) + p2 * l));
						}
					}

					if (m_curvePointsTo.size() == 2)
					{
						const auto pTemp = 2 * Vector2(intersectionPoint[0] / boxSize[0],
						                               intersectionPoint[1] / boxSize[1]);
						const QuadricCurve2 c = QuadricCurve2::Bezier(m_curvePointsTo[0], m_curvePointsTo[1], pTemp);

						for (int i = 0; i < nbPoints; i++)
						{
							auto tmp = c.Eval(1. - (static_cast<double>(i) / nbPoints));
							tmp = Vector2(tmp[0] * boxSize[0], tmp[1] * boxSize[1]) / 2.;
							points.append(h->Vertex(tmp));
						}
					}

					m_rendererCurveTo = std::make_unique<MayaSimpleRenderer>(
						points, Color(0.2, 0.2, 0.6), GL_LINE_STRIP);
					m_rendererCurveTo->setDepthTest(false);
					m_rendererCurveTo->setLineWidth(5.);
				}
				break;
			}
		}
	}
}

void ToolCurve::mousePressEvent(QMouseEvent* e)
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
			// [-1 ; 1]
			const auto p = 2 * Vector2(intersectionPoint[0] / boxSize[0], intersectionPoint[1] / boxSize[1]);

			switch (m_state)
			{
			case State::FIRST_CURVE:
				m_curvePointsFrom.push_back(p);

				if (m_curvePointsFrom.size() == 3)
				{
					m_curveFrom = QuadricCurve2::Bezier(m_curvePointsFrom[0], m_curvePointsFrom[1],
					                                    m_curvePointsFrom[2]);
					m_state = State::LASSO;
				}
				break;

			case State::LASSO:
				{
					if (!m_editLasso)
						m_editLasso = true;
					else
					{
						m_editLasso = false;
						m_state = State::SECOND_CURVE;
					}
				}
				break;

			case State::SECOND_CURVE:
				m_curvePointsTo.push_back(p);

				if (m_curvePointsTo.size() == 3)
				{
					m_curveTo = QuadricCurve2::Bezier(m_curvePointsTo[0], m_curvePointsTo[1], m_curvePointsTo[2]);
					qDebug() << "Curves finished";
					moveGaussians(0.2);
					m_parent->UpdateGaussiansBuffer();
					m_parent->rasterizeGaussians();
					m_curvePointsFrom.clear();
					m_curvePointsTo.clear();
					m_lasso.clear();
					m_state = State::FIRST_CURVE;
				}
				break;
			}
		}
	}
}
