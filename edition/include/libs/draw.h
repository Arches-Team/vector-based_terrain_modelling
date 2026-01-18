// Drawing
#pragma once

#include <QtWidgets/QGraphicsScene>

class Draw {
public:
  static void CreateVector(QGraphicsScene&, const QString&);
  static void CreateImage(QGraphicsScene&, const QString&, int, const QRectF& s = QRectF());
  static QImage CreateImage(QGraphicsScene&, int, const QRectF& s = QRectF());
};

