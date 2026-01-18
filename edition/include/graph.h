#pragma once

#include <list>
#include <unordered_set>

#include "libs/evector.h"
#include "libs/heightfield.h"
#include "libs/draw.h"


class Node;

struct Edge
{
	double weight;
	Node* orig;
	Node* dest;

	Edge(const double weight, Node* orig, Node* dest): weight(weight), orig(orig), dest(dest)
	{
	}
};

class Node
{
private:
	Vector2 m_pos;

	QVector<Edge> m_connectedTo;
	QVector<Edge> m_isConnected;

public:
	explicit Node(const Vector2& pos) : m_pos(pos)
	{
	}

	const Vector2& pos() const { return m_pos; }
	void moveTo(const Vector2& pos) { m_pos = pos; }
	void translate(const Vector2& offset) { m_pos += offset; }
	int nbNeighbors() const { return static_cast<int>(m_connectedTo.size()); }
	void addEdge(Node& dest) { addEdge(dest, Norm(dest.pos() - this->pos())); }

	void addEdge(Node& dest, const double weight)
	{
		m_connectedTo.append(Edge(weight, this, &dest));
		dest.m_isConnected.append(Edge(weight, this, &dest));
	}

	void removeAllConnections();
	void removeEdge(const Node& dest);
	const QVector<Edge>& connectedTo() const { return m_connectedTo; }
	const QVector<Edge>& isConnected() const { return m_isConnected; }

	friend std::ostream& operator<<(std::ostream& s, const Node& n);
};

class Graph
{
private:
	// Need a std::list instead of std::vector because list is not contiguous. 
	// Removing a node within graph with std::vector would invalidate Node& in node.
	std::list<Node> m_nodes;
	bool m_isOriented;

	static void getThreeNeighborsNeighbors(const Node& node, std::unordered_set<Node*>& neighbors,
	                                       const Vector2& origNodePos, double distanceThresh);

public:
	Graph() : m_isOriented(false)
	{
	}

	explicit Graph(const bool isOriented) : m_isOriented(isOriented)
	{
	}

	int addNode(Vector2 pos);
	void removeNode(int i);
	std::list<Node>::iterator removeNode(const std::list<Node>::iterator& it);
	Node& operator[](const int i) { return at(i); }
	Node& at(int i);
	int getId(const Node& node) const;
	int size() const { return static_cast<int>(m_nodes.size()); }

	QVector<Vector2> getNodesPosition() const;

	/*!
	\brief Add an edge to the graph.

	Weight is computed from the points as the Euclidean distance
	\param a,b Indexes.
	*/
	void addEdge(Node& a, Node& b) { addEdge(a, b, Norm(a.pos() - b.pos())); }
	void addEdge(Node& a, Node& b, double weight) const;
	void addEdge(int i, int j) { addEdge(i, j, Norm(at(i).pos() - at(j).pos())); }
	void addEdge(int i, int j, double weight) { addEdge(at(i), at(j), weight); }

	void clear() { m_nodes.clear(); }

	void reduceGraph();
	void reduceGraph(double distance);

	void print(const QString& name) const;
	void draw(QGraphicsScene& scene, const QColor& color = QColor(50, 50, 50)) const;

	std::list<Node>::iterator begin() { return m_nodes.begin(); }
	std::list<Node>::iterator end() { return m_nodes.end(); }
	std::list<Node>::const_iterator begin() const { return m_nodes.begin(); }
	std::list<Node>::const_iterator end() const { return m_nodes.end(); }
	std::list<Node>::const_iterator cbegin() const { return m_nodes.cbegin(); }
	std::list<Node>::const_iterator cend() const { return m_nodes.cend(); }

	friend std::ostream& operator<<(std::ostream& s, const Graph& g);
};


class GraphControl
{
public:
	static void createCrestRiverGraph(Graph& crest, Graph& river, const HeightField& _hf, double threshold = 70.);
};
