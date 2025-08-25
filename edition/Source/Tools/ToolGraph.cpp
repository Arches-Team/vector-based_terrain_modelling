#include "Tools/ToolGraph.h"

#include <queue>
#include <ranges>

#include "libs/convex.h"
#include "GaussianTerrainRaytracingWidget.h"
#include "libs/heightfield.h"

ToolGraph::ToolGraph(GaussianTerrainRaytracingWidget* parent, Kernels& kernels) : Tool(parent, kernels)
{
	initGraph();
	m_hfBuffer.Generate();
	m_differenceBuffer.Generate();
	m_hfBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * m_parent->getHF()->VertexSize(), nullptr,
	                   GL_STREAM_READ);
	m_differenceBuffer.SetData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * m_parent->getHF()->VertexSize(), nullptr,
	                           GL_STREAM_READ);

	const QString fullPath = QString::fromStdString(
		std::string(SOLUTION_DIR) + "/Shaders/gaussians_holes_detection.glsl");
	std::cout << fullPath.toStdString() << std::endl;
	QByteArray ba = fullPath.toLocal8Bit();
	m_differenceComputeShader = read_program(ba.data());

}

void ToolGraph::associateKernelsToEdges()
{
	m_edgeLinksKernel.clear();
	m_kernelLinksEdge.clear();

	for (auto& kernel : m_kernels)
	{
		Vector2 kernelPos = kernel->pos();
		const Edge* closestEdge = nullptr;
		double minDistance = std::numeric_limits<double>::max();

		auto checkAndUpdateClosestNode = [&](Graph& graph)
		{
			for (Node& node : graph)
			{
				for (const Edge& edge : node.connectedTo())
				{
					Segment2 segment(edge.orig->pos(), edge.dest->pos());
					const double distance = segment.R(kernelPos);
					if (distance < minDistance)
					{
						minDistance = distance;
						closestEdge = &edge;
					}
				}
			}
		};

		checkAndUpdateClosestNode(m_graphCrest);

		if (closestEdge)
		{
			const Vector2& origPos = closestEdge->orig->pos();
			const Vector2& destPos = closestEdge->dest->pos();
			m_edgeLinksKernel[closestEdge].origOldPos = origPos;
			m_edgeLinksKernel[closestEdge].destOldPos = destPos;

			auto vec = destPos - origPos;
			auto vecKernel = kernel->pos() - origPos;
			auto projection = vec * ((vec * vecKernel) / (vec * vec));

			m_edgeLinksKernel[closestEdge].kernels.push_back(KernelData{
				kernel.get(), Norm(projection) / Norm(vec), minDistance
			});

			m_kernelLinksEdge[kernel.get()] = closestEdge;
		}
	}
	createSprings();
}

void ToolGraph::createSprings()
{
	m_springs.clear();

	// Delaunay spring system
	// auto points = m_graphCrest.getNodesPosition();
	// auto mesh = Mesh2::Delaunay(points);
	// std::unordered_set<std::pair<int, int>, pair_hash> nodeIndexes;
	// const auto& indexes = mesh.Indexes();
	//
	// for (int i = 0; i < indexes.size(); i += 3)
	// {
	// 	nodeIndexes.insert(std::pair(indexes[i], indexes[i + 1]));
	// 	nodeIndexes.insert(std::pair(indexes[i + 1], indexes[i + 2]));
	// 	nodeIndexes.insert(std::pair(indexes[i + 2], indexes[i]));
	// }
	//
	// for (const auto& index : nodeIndexes)
	// {
	// 	auto& orig = m_graphCrest[index.first];
	// 	auto& dest = m_graphCrest[index.second];
	//
	// 	m_springs.emplace_back(Spring(&orig, &dest, Norm(dest.pos() - orig.pos()), m_k));
	// }

	for (const auto& node : m_graphCrest)
	{
		for (const auto& edge : node.connectedTo())
		{
			m_springs.emplace_back(Spring(edge.orig, edge.dest, Norm(edge.dest->pos() - edge.orig->pos()), m_k));
		}
	}
}

void ToolGraph::updateSpringsForce(Node* fixedNode) const
{
	constexpr int nbIter = 100;
	constexpr double tolerance = 1e-6;
	for (int iter = 0; iter < nbIter; ++iter)
	{
		bool stable = true;

		for (const auto& spring : m_springs)
		{
			if (spring.dest == fixedNode)
				continue;
			Vector2 totalForce{0.0, 0.0};
			int i = 0;
			Vector2 diff = spring.dest->pos() - spring.orig->pos();
			const double length = Norm(diff);

			if (length > 0)
			{
				const Vector2 force = -m_k * (length - spring.l0) * (diff / length);
				totalForce += force;
			}


			spring.dest->translate(totalForce * 1 / nbIter);
			if (Norm(totalForce) > tolerance)
			{
				stable = false;
			}
		}

		if (stable)
		{
			break;
		}
	}
}

void ToolGraph::moveNode(Node* node, const Vector2& offset)
{
	if (!node)
		return;

	struct NodeMovingData
	{
		Node* node;
		int depth;
	};

	m_movedNodes.clear();

	std::vector<NodeMovingData> nodesToMove;
	std::queue<NodeMovingData> treated;
	std::unordered_set<Node*> visited;

	nodesToMove.emplace_back(NodeMovingData{node, m_depth});
	treated.emplace(NodeMovingData{node, m_depth});
	visited.insert(node);

	while (!treated.empty())
	{
		auto& [currentNode, currentDepth] = treated.front();
		treated.pop();

		if (currentDepth > 0)
		{
			for (const auto& edge : currentNode->connectedTo())
			{
				Node* neighbor = edge.dest;
				if (!visited.contains(neighbor))
				{
					nodesToMove.emplace_back(NodeMovingData{neighbor, currentDepth - 1});
					treated.emplace(NodeMovingData{neighbor, currentDepth - 1});
					visited.insert(neighbor);
				}
			}
		}
	}

	nodesToMove.erase(nodesToMove.begin());
	node->translate(offset);
	m_movedNodes[node] = offset;

	if (m_translateOnly)
	{
		for (const auto& current : nodesToMove)
		{
			auto& [currentNode, currentDepth] = current;
			auto currentOffset = offset * (static_cast<double>((currentDepth + 1)) / (m_depth + 1));
			currentNode->translate(currentOffset);
			m_movedNodes[currentNode] = currentOffset;
		}
	}
	else
	{
		updateSpringsForce(node);
	}

	updateKernels();
	updateRenderer();
	m_parent->UpdateGaussiansBuffer();
	m_parent->rasterizeGaussians();
}

void ToolGraph::updateKernels()
{
	//std::cout << m_kernels[5] << std::endl;
	std::vector<Kernel*> hasMoved;
	QVector<Vector2> hasMovedPointsSample;

	std::vector<Kernel*> kernelsNotMoved;

    for (const auto& edge : m_edgeLinksKernel | std::views::keys)
	{
		auto [origOldPos, destOldPos, kernels] = m_edgeLinksKernel[edge];
		for (auto& [kernel, curvilignAbscissa, _] : kernels)
		{
			auto oldVec = destOldPos - origOldPos;
			auto newVec = edge->dest->pos() - edge->orig->pos();

			auto oldPos = oldVec * curvilignAbscissa;
			auto newPos = newVec * curvilignAbscissa;
			auto previousPos = kernel->pos();
			auto offset = (newPos - oldPos) + (edge->orig->pos() - origOldPos);
			kernel->translate(offset);
			if (previousPos != kernel->pos())
			{
                constexpr int nbSample = 10;
                hasMoved.emplace_back(kernel);
				auto ellipse = kernel->getEllipse(1.);
				for (int j = 1; j < nbSample; j++)
				{
					hasMovedPointsSample.emplace_back(ellipse.Vertex(2. * Math::Pi * static_cast<double>(j - 1) / nbSample));
				}
			}
			else
			{
				kernelsNotMoved.emplace_back(kernel);
			}

			const double angle = oldVec.Angle(newVec);
			kernel->rotate(static_cast<float>(angle), kernel->pos());

			const auto oldLength = Norm(oldVec);
			const auto newLength = Norm(newVec);
			const auto oldAmplitude = kernel->amplitude();

			float factor = (newLength / oldLength);
			constexpr float percentageFactorMax = 0.5f;
			if (abs(factor - 1) - 1e-6 > 0 && m_scaleKernel)
			{
				//factor = Clamp(factor, 1.f - percentageFactorMax, 1.f + percentageFactorMax);
				kernel->scale(newVec, factor);
			}
				
			//kernel->scale(static_cast<float>(newLength) / oldLength, kernel->pos());
			//kernel->setAmplitude(oldAmplitude);
		}
		m_edgeLinksKernel[edge].origOldPos = edge->orig->pos();
		m_edgeLinksKernel[edge].destOldPos = edge->dest->pos();
	}


	if (m_influenceEnabled)
	{
		std::unordered_map<Node*, std::vector<Vector2>> nodesToMove;

		for (auto kernel : kernelsNotMoved)
		{
			for (const auto& node : m_movedNodes | std::views::keys)
			{
				const auto dist = Norm(node->pos() - kernel->pos()) / m_influenceRadius;
				const auto offset = m_movedNodes[node];
				if (dist < 1.)
				{
					const auto falloff = 1. - dist;

					kernel->translate(offset * falloff);
					auto nodeOrig = m_kernelLinksEdge[kernel]->orig;
					auto nodeDest = m_kernelLinksEdge[kernel]->dest;

					nodesToMove[nodeOrig].emplace_back(offset * falloff);
					nodesToMove[nodeDest].emplace_back(offset * falloff);

					if (Abs(offset) - Vector(utils::eps) > Vector(0))
					{
						const auto newDist = Norm(node->pos() - kernel->pos()) / m_influenceRadius;
						const auto factor = (-(newDist - dist) + 1.);
						kernel->scale(offset, factor);
					}

					break;
				}
			}
		}

		for (auto& node : nodesToMove | std::views::keys)
		{
			Vector2 offset(0);
			for (const auto& v : nodesToMove[node])
				offset += v;
			offset /= nodesToMove[node].size();

			node->translate(offset);
		}
		for (const auto& edge : m_edgeLinksKernel | std::views::keys)
		{
			m_edgeLinksKernel[edge].origOldPos = edge->orig->pos();
			m_edgeLinksKernel[edge].destOldPos = edge->dest->pos();
		}
	}
	else
	{
		const auto hull = Convex2::Hull(hasMovedPointsSample);

		// Hide covered kernels
		if (hull.Size() > 0)
		{
			for (auto& k : m_kernels)
			{
				k->setModAmplitude(1.);
				const auto dist = hull.Signed(Vector2(k->posX(), k->posY()));
				if (dist < 0)
				{
					k->setModAmplitude(static_cast<float>(1. - Linear::Step(-dist, 0., m_blendThreshold)));
				}
			}
			for (const auto k : hasMoved)
			{
				k->setModAmplitude(1.);
			}
		}
	}

    
}

void ToolGraph::fillHoles()
{
	//TODO : fix with the new Kernels architecture
	/*glUseProgram(m_differenceComputeShader);
	m_hfBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 0);
	m_parent->getHFBuffer().BindAt(GL_SHADER_STORAGE_BUFFER, 1);
	m_differenceBuffer.BindAt(GL_SHADER_STORAGE_BUFFER, 2);

	int nx = m_parent->getHF()->GetSizeX(), ny = m_parent->getHF()->GetSizeY();
	glUniform1i(glGetUniformLocation(m_differenceComputeShader, "nx"), nx);
	glUniform1i(glGetUniformLocation(m_differenceComputeShader, "ny"), ny);

	int dispatchSize = (max(nx, ny) / 8) + 1;
	glDispatchCompute(dispatchSize, dispatchSize, 1);
	glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
	ScalarField2 difference(m_parent->getHF()->GetBox(), nx, ny, 0.);
	std::vector<float> data;
	const int totalBufferSize = nx * ny;
	data.resize(totalBufferSize);
	m_differenceBuffer.GetData(data);
	glUseProgram(0);

	for (int j = 0; j < totalBufferSize; j++)
	{
		difference[j] = data[j];
	}

	auto box = difference.GetBox();
	auto points = box.Poisson(10, 10000);
	Kernel mean = m_kernels.getMeanParameters() * 2.;
	for (auto& point : points)
	{
		point -= box[0];
		point /= (box[1] - box[0]);
		auto pointNorm = (point * 2.);
		pointNorm[0] -= 1.;
		pointNorm[1] -= 1.;

		point *= Vector2(nx, ny);
		if (difference.at(QPoint(static_cast<int>(point[0]), static_cast<int>(point[1]))) == 1.)
		{
			Kernel kernel{
				mean.scaleX(), mean.scaleY(), mean.modScale(), mean.theta(), mean.amplitude(),
				static_cast<float>(pointNorm[0]), static_cast<float>(pointNorm[1]), mean.beta()
			};
			m_kernels.add(kernel);
		}
	}
	emit m_parent->nbGaussiansChanged(m_parent->getNbGaussians());
	m_parent->rasterizeGaussians();*/
}

/*
\brief Initialize the graph
*/
void ToolGraph::initGraph()
{
	m_graphCrest.clear();
	m_graphRiver.clear();
	if (m_parent->getHF() != nullptr)
	{
		constexpr double distance = 0.1;
		constexpr double threshold = 70.;
		GraphControl::createCrestRiverGraph(m_graphCrest, m_graphRiver, *((HeightField *)(m_parent->getHF())), threshold);
		//GraphControl::createCrestRiverGraph(m_graphRiver, m_graphCrest, HeightField(*sf));
		m_graphCrest.reduceGraph(distance);
		m_graphRiver.reduceGraph(distance);
		associateKernelsToEdges();
		//updateRenderer();

		for (const auto& node : m_graphCrest)
		{
			m_nbLinesRenderer += m_nodeCircleSegment * 2;
			m_nbLinesRenderer += static_cast<int>(node.connectedTo().size()) * 2;
		}

		/*m_graphRenderer = std::make_unique<MayaSimpleRendererColors>(m_nbLinesRenderer, m_nbLinesRenderer, GL_LINES);
		m_graphRenderer->setDepthTest(false);
		m_graphRenderer->setLineWidth(1.);*/

		updateRenderer();
	}

	recordGraph();
}

void ToolGraph::fillSelectionNodes(Node* firstNode)
{
	if (m_nodeToMove)
		firstNode = m_nodeToMove;

	struct NodeMovingData
	{
		Node* node;
		int depth;
	};

	std::queue<NodeMovingData> treated;

	treated.emplace(NodeMovingData{ firstNode, m_depth });
	m_nodesSelection.clear();
	m_nodesSelection.insert(firstNode);

	while (!treated.empty())
	{
		auto& [currentNode, currentDepth] = treated.front();
		treated.pop();

		if (currentDepth > 0)
		{
			for (const auto& edge : currentNode->connectedTo())
			{
				Node* neighbor = edge.dest;
				if (!m_nodesSelection.contains(neighbor))
				{
					treated.emplace(NodeMovingData{ neighbor, currentDepth - 1 });
					m_nodesSelection.insert(neighbor);
				}
			}
		}
	}
}

void ToolGraph::recordGraph()
{
	const auto graphFilename = m_parent->getRecordName("graph_nodes");
	m_graphCrest.print(graphFilename);
}

Node* ToolGraph::getClosestNode(const Vector2& pos)
{
	Node* closest{nullptr};
	double minDistance = std::numeric_limits<double>::max();
	for (Node& node : m_graphCrest)
	{
		const double distance = Norm(pos - node.pos());
		if (distance < minDistance)
		{
			minDistance = distance;
			closest = &node;
		}
	}
	return closest;
}

void ToolGraph::mousePressEvent(QMouseEvent* e)
{
	if (m_parent->getHF() != nullptr)
	{
		auto x = e->globalPosition().x();
		auto y = e->globalPosition().y();

		const HeightField* h = (HeightField*)m_parent->getHF();
		Vector intersectionPoint;
		double t;
		const Ray ray = m_parent->ConvertPixelToRay(e->pos());

		if (h->Intersect(ray, t, intersectionPoint, h->GetBox(), h->K()))
		{
			m_canMove = true;

			auto boxSize = h->GetBox().Size();
			const auto intersectionNormalized = 2 * Vector2(intersectionPoint[0] / boxSize[0],
			                                                intersectionPoint[1] / boxSize[1]);
			m_nodeToMove = getClosestNode(intersectionNormalized);
			createSprings();
			glBindBuffer(GL_COPY_READ_BUFFER, m_parent->getHFBuffer().GetBuffer());
			glBindBuffer(GL_COPY_WRITE_BUFFER, m_hfBuffer.GetBuffer());
			glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0,
			                    sizeof(float) * m_parent->getHF()->VertexSize());
			glBindBuffer(GL_COPY_READ_BUFFER, 0);
			glBindBuffer(GL_COPY_WRITE_BUFFER, 0);

			updateRenderer();
		}
	}
}

void ToolGraph::mouseMoveEvent(QMouseEvent* e)
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
			const auto intersectionNormalized = 2 * Vector2(intersectionPoint[0] / boxSize[0],
			                                                intersectionPoint[1] / boxSize[1]);
			m_lastIntersection = intersectionNormalized;
			if (m_canMove)
			{
				const auto offset = intersectionNormalized - m_oldPos;
				moveNode(m_nodeToMove, offset);
			}

			m_oldPos = intersectionNormalized;

			fillSelectionNodes(getClosestNode(intersectionNormalized));
			updateRenderer();
		}
	}
}

void ToolGraph::mouseReleaseEvent(QMouseEvent* e)
{
	m_canMove = false;
	m_nodesSelection.clear();
	m_nodeToMove = nullptr;
	if(m_needFillHoles)
	    fillHoles();
	updateRenderer();

	m_parent->recordHF("graph");
	recordGraph();


	// Remove 0. amplitude kernels
	// std::vector<Kernel>& kernels = m_kernels.getKernels();
	// const auto newEnd = std::remove_if(kernels.begin(), kernels.end(),
	//                                    [](Kernel& k) { return !k.isShown(); });
	//
	// const auto sizeBefore = kernels.size();
	// kernels.erase(newEnd, kernels.end());
	//
	// if (sizeBefore > kernels.size())
	// {
	// 	emit m_parent->nbGaussiansChanged(m_parent->getNbGaussians());
	// 	m_parent->rasterizeGaussians();
	// }
	// initGraph();
}

void ToolGraph::mouseWheelEvent(QWheelEvent* e)
{
	const double angle = e->angleDelta().y();

	// Influence radius
	if (e->modifiers() & Qt::ControlModifier && e->modifiers() & Qt::ShiftModifier)
	{
		if (m_influenceEnabled)
		{
			const double scaleRadius = e->angleDelta().y() * 0.0001;
			m_influenceRadius = max(m_influenceRadius + scaleRadius, 0.);

			updateRenderer();
		}
	}
	// Amplitude
	else if (e->modifiers() & Qt::ControlModifier)
	{
		auto closestNode = getClosestNode(m_lastIntersection);
		std::queue<std::pair<Node*, int>> treated;
		std::unordered_set<Node*> visited;

		treated.emplace(closestNode, m_depth);
		visited.insert(closestNode);

		
		while (!treated.empty())
		{
			auto& [currentNode, currentDepth] = treated.front();
			treated.pop();

			if (currentDepth > 0)
			{
				for (const auto& edge : currentNode->connectedTo())
				{
					Node* neighbor = edge.dest;
					if (!visited.contains(neighbor))
					{
						treated.emplace(neighbor, currentDepth - 1);
						visited.insert(neighbor);
					}
				}
			}
		}

		// Get the max distance from the closest node
		double maxDist = 0.;
		for (const auto node : visited)
		{
			for (const auto& edge : m_edgeLinksKernel | std::views::keys)
			{
				if (edge->orig != node && edge->dest != node)
					continue;

				auto [origOldPos, destOldPos, kernels] = m_edgeLinksKernel[edge];
				for (const auto& [kernel, curvilignAbscissa, _] : kernels)
				{
					maxDist = std::max(maxDist, Norm(closestNode->pos() - kernel->pos()));
				}
			}
		}

		// Update the amplitude
		const double variation = 1. + angle / 100.;
		//std::cout << variation << std::endl;
		for (const auto node : visited)
		{
			for (const auto& edge : m_edgeLinksKernel | std::views::keys)
			{
				if (edge->orig != node && edge->dest != node)
					continue;

				auto [origOldPos, destOldPos, kernels] = m_edgeLinksKernel[edge];
				for (const auto& [kernel, curvilignAbscissa, _] : kernels)
				{
					const auto dist = Norm(closestNode->pos() - kernel->pos()) / maxDist;
					kernel->amplitude() *= static_cast<float>(variation * (1 - dist) + dist);
				}
			}
		}
		m_parent->UpdateGaussiansBuffer();
		m_parent->rasterizeGaussians();

		m_parent->recordHF("graph_amplitude");

	}
	// Depth change
	else if(e->modifiers() & Qt::ShiftModifier)
	{
		const double variation = angle/100.;

		m_depth += static_cast<int>(variation);
		m_depth = std::max(0, m_depth);

		fillSelectionNodes(getClosestNode(m_lastIntersection));
		updateRenderer();


		emit m_parent->updateDepthGraph(m_depth);
	}
	
}

void ToolGraph::keyPressedEvent(QKeyEvent* e)
{
	if(e->key() == Qt::Key_Delete)
	{
		std::cout << "Delete" << std::endl;
		fillSelectionNodes(getClosestNode(m_lastIntersection));

		auto& kernels = m_kernels.getKernels();
		std::unordered_set<Kernel*> kernelsToRemove;

		for (const Node* node : m_nodesSelection)
		{
			// Iterate over edges where the node is the origin
			for (const Edge& edge : node->connectedTo())
			{
				auto it = m_edgeLinksKernel.find(&edge);
				if (it != m_edgeLinksKernel.end())
				{
					for (KernelData& kd : it->second.kernels)
					{
						kernelsToRemove.insert(kd.kernel);
					}
				}
			}

			// Iterate over edges where the node is the destination
			for (const Edge& edge : node->isConnected())
			{
				auto it = m_edgeLinksKernel.find(&edge);
				if (it != m_edgeLinksKernel.end())
				{
					for (KernelData& kd : it->second.kernels)
					{
						kernelsToRemove.insert(kd.kernel);
					}
					it->second.kernels.clear(); // Clear the kernels vector
				}
			}
		}

		// Remove kernels from the kernels vector
		m_kernels.removeIf([&kernelsToRemove](const std::unique_ptr<Kernel>& kernel)
			{
				return kernelsToRemove.contains(kernel.get());
			});
	    
		m_parent->UpdateGaussiansBuffer();
		m_parent->rasterizeGaussians();

		initGraph();
		updateRenderer();

		m_parent->recordHF("graph_delete");
	}
}


void ToolGraph::setShowGraph(const bool show)
{
	m_showGraph = show;
	updateRenderer();
}

void ToolGraph::updateRenderer()
{
	deleteRenderer();

	if (!m_showGraph)
		return;

	QVector<VectorFloat> points(m_nbLinesRenderer);
	QVector<ColorFloat> colors(m_nbLinesRenderer);

	//const ColorFloat colorEdge(0.2f, 0.2f, 0.2f);
	//const ColorFloat colorNode(0.8f, 0.2f, 0.2f);
	//const ColorFloat colorNodeSelected(0.2f, 0.8f, 0.2f);


	const ColorFloat colorEdge(0.2f, 0.2f, 0.2f);
	const ColorFloat colorNode(0.7f, 0.7f, 0.7f);
	const ColorFloat colorNodeSelected(0.7f, 0.1f, 0.7f);

	const HeightField* hf = (HeightField*)m_parent->getHF();

	// TODO: change 1250. for a dynamic value>
	auto normPos = [&](const Vector2& pos) { return pos * 1250.; };

	int cpt = 0;
	for (const auto& node : m_graphCrest)
	{
		const auto& color = m_nodesSelection.contains(&node) ? colorNodeSelected : colorNode;
		const Circle2 circle(normPos(node.pos()), 10.);
		for (int i = 1; i < m_nodeCircleSegment; i++)
		{
			Vector2 previousPos = circle.Vertex(2. * Math::Pi * static_cast<double>(i - 1) / m_nodeCircleSegment);
			Vector2 pos = circle.Vertex(2. * Math::Pi * static_cast<double>(i) / m_nodeCircleSegment);
			points[cpt] = VectorFloat(hf->Vertex(previousPos));
			colors[cpt] = color;
			cpt++;
			points[cpt] = VectorFloat(hf->Vertex(pos));
			colors[cpt] = color;
			cpt++;
		}
		points[cpt] = points[cpt - m_nodeCircleSegment * 2 + 2];
		colors[cpt] = color;
		cpt++;
		points[cpt] = VectorFloat(hf->Vertex(circle.Vertex(0.)));
		colors[cpt] = color;
		cpt++;

		// Edge
		for (const auto& edge : node.connectedTo())
		{
			points[cpt] = VectorFloat(hf->Vertex(normPos(edge.orig->pos())));
			colors[cpt] = colorEdge;
			cpt++;
			points[cpt] = VectorFloat(hf->Vertex(normPos(edge.dest->pos())));
			colors[cpt] = colorEdge;
			cpt++;
		}

		//std::cout << cpt << "\n";
	    //std::cout << normPos(node.pos()) << "\n";
	}
	
	if (m_influenceEnabled)
	{
		QVector<Vector> influencePoints{};
		for (auto& node : m_nodesSelection)
		{
			const Circle2 circle(normPos(node->pos()), 1250. * m_influenceRadius);
			int nb = m_nodeCircleSegment * static_cast<int>(m_nodesSelection.size());
			for (int i = 0; i < nb; i++)
			{
				Vector2 pos = circle.Vertex(2. * Math::Pi * static_cast<double>(i) / nb);
				Vector2 nextPos = circle.Vertex(2. * Math::Pi * static_cast<double>(i + 1) / nb);
				bool isDrawn = true;
				for (const auto& n : m_nodesSelection)
				{
					if (n == node)
						continue;
					const Circle2 circle(normPos(n->pos()), 1250. * m_influenceRadius);
					if (circle.Inside(pos) || circle.Inside(nextPos))
					{
						isDrawn = false;
						break;
					}
				}
				if (isDrawn)
				{
					influencePoints.append(hf->Vertex(pos));
					influencePoints.append(hf->Vertex(nextPos));
				}
			}
		}

		if (!influencePoints.isEmpty())
		{
			m_renderer = std::make_unique<MayaSimpleRenderer>(influencePoints, Color(0.2, 0.6, 0.2), GL_LINES);
			m_renderer->setDepthTest(false);
			m_renderer->setLineWidth(2.);
		}
		
	}
	
	m_graphRenderer = std::make_unique<MayaSimpleRendererColors>(points, colors, GL_LINES);
	m_graphRenderer->setDepthTest(false);
	m_graphRenderer->setLineWidth(2.0);

	//m_parent->update();
}
