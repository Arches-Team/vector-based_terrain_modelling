// Fundamentals

#pragma once

#include "libs/random.h"

#include <QtCore/QVector>

class Shuffle
{
protected:
public:
  //! Empty.
  Shuffle() {}
  //! Empty
  ~Shuffle() {}

  static QVector<int> Table(int);

private:
  static RandomFast r;
};
