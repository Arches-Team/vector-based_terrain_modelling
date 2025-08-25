#include "qte.h"
#include "libs/triangle.h"

int main(int argc, char *argv[])
{
  QApplication app(argc, argv);

  QteWindow w;
  w.show();

  return app.exec();
}


