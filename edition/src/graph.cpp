#include "graph.h"

#include <ranges>
#include <QtGui/qpainter.h>


void Node::removeAllConnections()
{
	for (const auto& connected : m_isConnected)
	{
		connected.orig->removeEdge(*this);
	}

	for (const auto& connectedTo : m_connectedTo)
	{
		connectedTo.dest->removeEdge(*this);
	}
}

void Node::removeEdge(const Node& dest)
{
	m_connectedTo.erase(
		std::ranges::remove_if(m_connectedTo, [&dest](const auto& edge)
		{
			return edge.dest == std::addressof(dest);
		}).begin(),
		m_connectedTo.end()
	);

	m_isConnected.erase(
		std::ranges::remove_if(m_isConnected, [&dest](const auto& edge)
		{
			return edge.orig == std::addressof(dest);
		}).begin(),
		m_isConnected.end()
	);
}

void Graph::getThreeNeighborsNeighbors(const Node& node, std::unordered_set<Node*>& neighbors,
                                       const Vector2& origNodePos, double distanceThresh)
{
	for (const auto& edge : node.connectedTo())
	{
		auto& neighbor = edge.dest;
		const double dist = Norm(neighbor->pos() - origNodePos);
		if (neighbor->nbNeighbors() >= 3 && dist < distanceThresh)
		{
			if (!neighbors.contains(neighbor))
			{
				neighbors.insert(neighbor);
				getThreeNeighborsNeighbors(*neighbor, neighbors, origNodePos, distanceThresh);
			}
		}
	}
}

int Graph::addNode(Vector2 pos)
{
	m_nodes.emplace_back(pos);
	return static_cast<int>(m_nodes.size()) - 1;
}

void Graph::removeNode(const int i)
{
	at(i).removeAllConnections();
	auto it = m_nodes.begin();
	std::advance(it, i);
	m_nodes.erase(it);
}

std::list<Node>::iterator Graph::removeNode(const std::list<Node>::iterator& it)
{
	it->removeAllConnections();
	return m_nodes.erase(it);
}

Node& Graph::at(const int i)
{
	auto it = m_nodes.begin();
	std::advance(it, i);
	return *it;
}

int Graph::getId(const Node& node) const
{
	int i = 0;
	for (const auto& n : m_nodes)
	{
		if (std::addressof(node) == std::addressof(n))
			break;
		i++;
	}
	return i;
}

QVector<Vector2> Graph::getNodesPosition() const
{
	QVector<Vector2> points;
	for (const auto& node : m_nodes)
		points.append(node.pos());
	return points;
}

void Graph::addEdge(Node& a, Node& b, const double weight) const
{
	a.addEdge(b, weight);
	if (!m_isOriented)
		b.addEdge(a, weight);
}

// TODO: fix for oriented graph
void Graph::reduceGraph()
{
	// Primitives are evaluated between [-1 ; 1] which results in 4. for the longest edge
	reduceGraph(4.);
}


// TODO fix for oriented graph
void Graph::reduceGraph(const double distance)
{
	auto it = m_nodes.begin();

	while (it != m_nodes.end())
	{
		auto& neighbors = it->connectedTo();
		if (neighbors.size() == 2)
			++it;
		else if (neighbors.size() >= 3)
		{
			std::unordered_set<Node*> connectedNodes;
			connectedNodes.insert(&*it);

			// Compute all neihborgs with 3+ neighbors and merge them
			getThreeNeighborsNeighbors(*it, connectedNodes, it->pos(), distance);
			if (connectedNodes.size() > 1)
			{
				const int newNodeId = addNode(Vector2(0));
				Node& newNode = at(newNodeId);
				Vector2 avg(0.);
				int i = 0;
				for (const auto node : connectedNodes)
				{
					avg += node->pos();
					i++;

					for (const auto& edge : node->connectedTo())
					{
						auto* neighbor = edge.dest;
						if (!connectedNodes.contains(neighbor))
							addEdge(newNode, *neighbor);
					}

					if (std::addressof(*it) != std::addressof(*node))
					{
						removeNode(getId(*node));
					}
				}

				avg /= i;
				newNode.moveTo(avg);

				it = removeNode(it);
			}
			else
				++it;
		}
		else
		{
			++it;
		}
	}

	it = m_nodes.begin();
	while (it != m_nodes.end())
	{
		auto& edges = it->connectedTo();
		if (edges.size() == 2)
		{
			if (Norm(edges[0].dest->pos() - edges[1].dest->pos()) < distance*2.)
			{
				addEdge(*edges[0].dest, *edges[1].dest);
				it = removeNode(it);
			}
			else
			{
				++it;
			}
			
		}
		else
		{
			++it;
		}
	}

	std::vector<int> idOneNeighborToRemove;
	for (int i = 0; i < m_nodes.size(); ++i)
	{
		auto& node = at(i);
		auto& edges = node.connectedTo();
		if ((edges.size() == 0) || 
			(edges.size() == 1 && (Norm(edges[0].dest->pos() - edges[0].orig->pos()) < distance || edges[0].dest->connectedTo().size() == 1)))
		{
			idOneNeighborToRemove.emplace_back(i);
		}
	}
	std::sort(idOneNeighborToRemove.begin(), idOneNeighborToRemove.end(), std::greater<int>());
	for (auto i : idOneNeighborToRemove)
	{
		removeNode(i);
	}
}

void Graph::print(const QString& name) const
{
	// Draw the graph in a QImage with a qgraphicscene
	QGraphicsScene scene;
	draw(scene, QColor(0, 0, 0));

	const auto image = Draw::CreateImage(scene, 2048, QRectF(0, 0,0, 0));
	image.scaled(QSize(2048, 2048)).save(name);
}

void Graph::draw(QGraphicsScene& scene, const QColor& color) const
{
	for (const auto& node : m_nodes)
	{
		Vector2 pos = node.pos();
		pos[0] += 1.;
		pos[1] += 1.;
		pos = (pos / 2.) * 255.;
		for (const auto& edge : node.connectedTo())
		{
			Vector2 destPost = edge.dest->pos();
			destPost[0] += 1.;
			destPost[1] += 1.;
			destPost = (destPost / 2.) * 255.;
			Segment2(pos, destPost).Draw(scene, QPen(QColor(255, 0, 0), 0.75));
		}
	}

	for (const auto& node : m_nodes)
	{
		Vector2 pos = node.pos();
		pos[0] += 1.;
		pos[1] += 1.;
		pos = (pos / 2.) * 255.;
		Circle2(pos, 2.0).Draw(scene, QPen(color, 0.75));
	}
}

void GraphControl::createCrestRiverGraph(Graph& crest, Graph& river, const HeightField& _hf, double threshold)
{
	double riverThreshold = threshold;
	double crestThreshold = threshold;

	// Clean HeightField
	HeightField hf = _hf;
	if (hf.GetSizeX() > 256)
	{
		int diff = std::log2(hf.GetSizeX()) - 8;
		hf = hf.DownSample(std::pow(2, diff));
	}
		
	hf.CompleteBreach();

	constexpr double eps = 1e-2;

	// Compute the skeleton of the river
	ScalarField2 sf2River = hf.StreamAreaSteepest();

	// Get rivers
	for (int i = 0; i < hf.GetSizeX(); i++)
	{
		for (int j = 0; j < hf.GetSizeY(); j++)
		{
			double streamValue = sf2River.at(i, j);
			sf2River(i, j) = 0;

			const auto gradient = hf.Gradient(i, j);

			if (streamValue > riverThreshold && !(abs(gradient[0]) - eps < 0. && abs(gradient[1]) - eps < 0.))
				sf2River(i, j) = 1;
		}
	}

	ScalarField2 skeletonRiver = sf2River.MorphSkeletonConnected(1);

	// Compute the skeleton of the crest
	HeightField hfC = hf;
	hfC.Invert();
	hfC.CompleteBreach();

	// Compute the skeleton of river and crest
	ScalarField2 sf2Crest = hfC.StreamAreaSteepest();
	// Get rivers
	for (int i = 0; i < hf.GetSizeX(); i++)
	{
		for (int j = 0; j < hf.GetSizeY(); j++)
		{
			double streamValue = sf2Crest.at(i, j);
			sf2Crest(i, j) = 0;

			const auto gradient = hf.Gradient(i, j);
			if (streamValue > crestThreshold && !(abs(gradient[0]) - eps < 0. && abs(gradient[1]) - eps < 0.))
				sf2Crest(i, j) = 1;
		}
	}

	ScalarField2 skeletonCrest = sf2Crest.MorphSkeletonConnected(3);

	// Merge the two skeletons in a QImage with 2 different colors
	QImage imgSkeleton(hf.GetSizeX(), hf.GetSizeY(), QImage::Format_RGB32);
	for (int i = 0; i < hf.GetSizeX(); i++)
	{
		for (int j = 0; j < hf.GetSizeY(); j++)
		{
			imgSkeleton.setPixel(i, j, qRgb(0, 0, 0));
			if (skeletonCrest(i, j) == 1) imgSkeleton.setPixel(i, j, qRgb(255, 0, 0));
			if (skeletonRiver(i, j) == 1) imgSkeleton.setPixel(i, j, qRgb(0, 0, 255));
		}
	}

	// Create the graph
	Array2I indexesArrayCrest(Box2(hf.GetBox()), hf.GetSizeX(), hf.GetSizeY(), -1);
	Array2I indexesArrayRiver(Box2(hf.GetBox()), hf.GetSizeX(), hf.GetSizeY(), -1);

	for (int i = 0; i < hf.GetSizeX(); i++)
	{
		for (int j = 0; j < hf.GetSizeY(); j++)
		{
			double x = (2. * (static_cast<float>(i) / hf.GetSizeX())) - 1.;
			double y = (2. * (static_cast<float>(j) / hf.GetSizeY())) - 1.;
			if (skeletonCrest(i, j) == 1)
			{
				indexesArrayCrest(i, j) = crest.addNode(Vector2(x, y));
			}
			if (skeletonRiver(i, j) == 1)
			{
				indexesArrayRiver(i, j) = river.addNode(Vector2(x, y));
			}
		}
	}

	for (int i = 0; i < hf.GetSizeX(); i++)
	{
		for (int j = 0; j < hf.GetSizeY(); j++)
		{
			for (int s = 0; s < 2; ++s)
			{
				auto& skeleton = s == 0 ? skeletonCrest : skeletonRiver;
				auto& graph = s == 0 ? crest : river;
				auto& indexesArray = s == 0 ? indexesArrayCrest : indexesArrayRiver;

				if (skeleton(i, j) == 1)
				{
					int indexNode = indexesArray(i, j);

					if (i > 0 && skeleton(i - 1, j) == 1)
					{
						int indexNodeCurrent = indexesArray(i - 1, j);
						graph.addEdge(indexNode, indexNodeCurrent);
					}
					if (j > 0 && skeleton(i, j - 1) == 1)
					{
						int indexNodeCurrent = indexesArray(i, j - 1);
						graph.addEdge(indexNode, indexNodeCurrent);
					}
					if (i > 0 && j > 0 && skeleton(i - 1, j - 1) == 1)
					{
						if (skeleton(i - 1, j) == 0 && skeleton(i, j - 1) == 0)
						{
							int indexNodeCurrent = indexesArray(i - 1, j - 1);
							graph.addEdge(indexNode, indexNodeCurrent);
						}
					}
					if (i < hf.GetSizeX() - 1 && j > 0 && skeleton(i + 1, j - 1) == 1)
					{
						if (skeleton(i + 1, j) == 0 && skeleton(i, j - 1) == 0)
						{
							int indexNodeCurrent = indexesArray(i + 1, j - 1);
							graph.addEdge(indexNode, indexNodeCurrent);
						}
					}
				}
			}
		}
	}
}

std::ostream& operator<<(std::ostream& s, const Node& n)
{
	s << "\tNode " << n.pos() << std::endl;
	for (const auto& edge : n.connectedTo())
	{
		s << "\t\t Connected to " << edge.dest->pos() << "\n";
	}
	for (const auto& edge : n.isConnected())
	{
		s << "\t\t Is connected " << edge.orig->pos() << "\n";
	}

	return s;
}

std::ostream& operator<<(std::ostream& s, const Graph& g)
{
	s << "Graph\n";
	for (const auto& node : g.m_nodes)
	{
		s << node;
	}
	s << "\n";
	return s;
}
