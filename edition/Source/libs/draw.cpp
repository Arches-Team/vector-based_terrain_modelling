// Draw

#include <QtPrintSupport/QPrinter>
#include <QtGui/QPainter>

#include "libs/draw.h"
#include "libs/mathematics.h"

/*!
\brief Draws the graphics scene as a vector image.
\param scene The scene.
\param name Name of the file.
\param size Size in centimeters
*/
void Draw::CreateVector(QGraphicsScene& scene, const QString& name)
{
  // High resolution printer
  QPrinter printer(QPrinter::HighResolution);
  printer.setOutputFileName(name);
  printer.setFullPage(true);

  // Size to A4
  printer.setPageSize(QPageSize::A4);

  // Painter 
  QPainter painter(&printer);
  painter.setRenderHint(QPainter::Antialiasing);

  // Render
  scene.render(&painter);
}

/*!
\brief Draws the graphics scene as an image.
\param scene The scene.
\param name Name of the file.
*/
void Draw::CreateImage(QGraphicsScene& scene, const QString& name, int width, const QRectF& rect)
{
  // Image
  QImage image = Draw::CreateImage(scene, width, rect);

  image.save(name);
}

/*!
\brief Draws the graphics scene as an image.
\param scene The scene.
\param name Name of the file.
*/
QImage Draw::CreateImage(QGraphicsScene& scene, int width, const QRectF& rect)
{
  QSizeF size;
  if (rect.isNull())
    size = scene.sceneRect().size();
  else
    size = rect.size();
  size = width * size / (Math::Max(size.width(), size.height()));

  // Set to integer
  QSize s = size.toSize();

  // Printer
  QImage image(s, QImage::Format_ARGB32);

  // Paint
  QPainter painter(&image);

  painter.setRenderHint(QPainter::Antialiasing);

  scene.render(&painter, QRectF(QPointF(0., 0.), size), rect);

  return image;
}
