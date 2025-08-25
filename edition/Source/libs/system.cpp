// System

#include <QtWidgets/QGraphicsScene>

#include <QtWidgets/QGraphicsPixmapItem>
#include <QtGui/QImage>
#include <QtGui/QPainter>
#include <QtPrintSupport/QPrinter>

#include "libs/mathematics.h"
#include "libs/cpu.h"

/*!
\brief Save a scene as a PDF file.
\param scene The scene.
\param name Filename.
*/
void System::SavePdf(QGraphicsScene& scene, const QString& name)
{
  QSizeF sizebb = scene.sceneRect().size();
  sizebb = (sizebb / (max(sizebb.width(), sizebb.height()))) * 100;

  QPrinter printer(QPrinter::HighResolution);

  printer.setOutputFileName(name);
  printer.setFullPage(true);
  printer.setPageSize(QPageSize::A4);

  QPainter painter2(&printer);
  painter2.setRenderHint(QPainter::Antialiasing);
  painter2.scale(1.0, 1.0);
  scene.render(&painter2, QRect(), QRect());

  // Use the following code to save as SVG
  /*
  //#include <QtSvg/QSvgGenerator>

  QSvgGenerator generator;
  generator.setFileName(name);
  generator.setSize(QSize(scene.width(), scene.height()));

  QPainter painter(&generator);
  scene.render(&painter);
  */
}

/*!
\brief Save the image as two .png and .jpg files.
\param image The image box.
\param name Filename.
*/
void System::SaveImage(const QImage& image, const QString& name)
{
  // Set name according to date and time
  QString fullname = name + '-' + System::DateTime();

  // Save
  image.save(fullname + ".png", "PNG");
  image.save(fullname + ".jpg", "JPG", 96);
}

/*!
\brief Save a scene as a PDF file.
\param box The cutting box.
\param scene The scene.
\param name Filename.
*/
void System::SavePdf(QGraphicsScene& scene, const QString& name, const Box2& box)
{
  scene.setSceneRect(box[0][0], box[0][1], box.Width(), box.Height());

  QSizeF sizebb = scene.sceneRect().size();
  sizebb = (sizebb / (max(sizebb.width(), sizebb.height()))) * 100;

  QPrinter printer(QPrinter::HighResolution);

  printer.setOutputFileName(name);
  printer.setFullPage(true);
  printer.setPageSize(QPageSize::A4);

  QPainter painter2(&printer);
  painter2.setRenderHint(QPainter::Antialiasing);
  painter2.scale(1.0, 1.0);
  scene.render(&painter2, QRect(), QRect());
}

/*!
\brief Flip all elements in the scene horizontally.
\param scene The scene.
\param box Bouding box for performing the flip.
*/
void System::FlipVertical(QGraphicsScene& scene, const Box2& box)
{
  QList<QGraphicsItem*> l = scene.items();
  if (l.size() == 0)
    return;

  for (int i = 0; i < l.size(); i++)
  {
    QTransform trans = l[i]->transform();
    trans.scale(1, -1);
    l[i]->setTransform(trans);
    QPointF p = l[i]->pos();

    QTransform tempT = l[i]->transform();
    tempT.translate(0, 2 * p.ry() - (box[1][1] + box[0][1]));
    l[i]->setTransform(tempT);

  }
}

/*!
\brief Flip all elements in the scene vertically.
\param scene The scene.
\param box Bouding box for performing the flip.
*/
void System::FlipHorizontal(QGraphicsScene& scene, const Box2& box)
{
  QList<QGraphicsItem*> l = scene.items();
  if (l.size() == 0)
    return;

  for (int i = 0; i < l.size(); i++)
  {
    QTransform trans = l[i]->transform();
    trans.scale(-1, 1);
    l[i]->setTransform(trans);
    QPointF p = l[i]->pos();

    QTransform tempT = l[i]->transform();
    tempT.translate(2 * p.rx() - (box[0][0] + box[1][0]), 0);
    l[i]->setTransform(tempT);
  }
}

/*!
\brief Create an image from a scene.
\param scene The scene.
\param box The box which will be rasterized.
\param size The maximum size (in pixels) of the image.
*/
QImage System::Rasterize(QGraphicsScene& scene, const Box2& box, int size)
{
  scene.setSceneRect(box[0][0], box[0][1], box.Width(), box.Height());

  QSizeF sizebb = scene.sceneRect().size();
  sizebb = (sizebb / (max(sizebb.width(), sizebb.height()))) * size;

  QImage image(sizebb.toSize(), QImage::Format_RGB32);
  image.fill(QColor(Qt::white).rgb());

  QPainter painter(&image);
  scene.render(&painter, QRect(), QRect());
  QImage mirrored = image.mirrored(false, true);

  return mirrored;
}


#include <fstream>

/*!
\brief Compute the number of lines of a file.
\param name Filename.
*/
int System::NumberOfLines(const QString& name)
{
  //  std::cout << "NumberOfLines = " << name.toLocal8Bit().constData() << std::endl;


  std::ifstream inputFile(name.toLocal8Bit().constData());  // Opening the file named "test.txt" for reading

  if (inputFile.is_open()) {  // Checking if the file was successfully opened
    std::string line;  // Declaring a string variable to store each line of text
    int lineCount = 0;  // Initializing a variable to count lines

    while (std::getline(inputFile, line)) {  // Loop through each line in the file
      lineCount++;  // Incrementing line count for each line read
    }

    inputFile.close();  // Closing the file after counting lines

    //  std::cout << "Number of lines in the file: " << lineCount << std::endl;  // Outputting the total line count
    return lineCount;
  }
  return 0;
}

#include <QtCore/QFileInfo>
#include <QtCore/QDir>

/*!
\brief Recurse inside directory hierachy
\param name Directory name.
*/
int System::RecurseDirectory(const QString& name, bool code, bool glsl)
{
   std::cout << "RecurseDirectory = " << name.toLocal8Bit().constData()<<std::endl;

  int total = 0;

  QDir dir(name);
  QFileInfoList list = dir.entryInfoList();
  for (int i = 0; i < list.count(); i++)
  {
    QFileInfo info = list[i];

    QString sFilePath = info.filePath();
    if (info.isDir())
    {
      // recursive
      if (info.fileName() != ".." && info.fileName() != ".")
      {
        total += RecurseDirectory(sFilePath, code, glsl);
      }
    }
    else
    {
        // Do something with the file here
      if (code == true)
      {
        if ((info.completeSuffix() == QString("cpp")) || (info.completeSuffix() == QString("h")))
        {
          total += NumberOfLines(info.absoluteFilePath());
        }
      }
      if (glsl == true)
      {
        if ((info.completeSuffix() == QString("glsl")))
        {
          total += NumberOfLines(info.absoluteFilePath());
        }
      }
    }
  }
  return total;
}

/*!
\brief Compute the number of lines of code in the argument directory.
\param name Directory name.
*/
int System::CodeLines(const QString& name, bool code, bool glsl)
{
  return System::RecurseDirectory(name, code, glsl);
}

