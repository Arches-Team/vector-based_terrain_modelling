// System

#include <QtCore/QProcessEnvironment>
#include <QtCore/QFileInfo>

#include "libs/cpu.h"


#ifdef _WIN32

// Architecture
#include <windows.h>

/*!
\defgroup SystemCore Generic core classes.

\brief Generic core classes include System, Hashing and other classes.
*/

/*!
\class System cpu.h
\brief A core class providing core system functions.
\ingroup SystemCore
*/

/*!
\brief Get the environment variable.
\param name Name.
*/
const char* System::WinGetEnv(const char* name)
{
  const unsigned long size = 65535;
  static char buffer[size];

  if (GetEnvironmentVariableA(name, buffer, size))
  {
    return buffer;
  }
  else
  {
    return 0;
  }
}

#endif

bool System::avx = true;

/*!
Set the Advanced Vector Extensions flag.
\param a Boolean, set to true to run AVX accelerated algorithms.
*/
void System::SetAvx(bool a)
{
  avx = a;
}

/*!
Get the Advanced Vector Extensions flag.
\param a Boolean, set to true to run AVX accelerated algorithms.
*/
bool System::Avx()
{
  return avx;
}

/*!
\brief Get the environment variable with Qt.
\param name Name.
*/
QString System::GetEnv(const QString& name)
{
  QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
  return env.value(name);
}

/*!
\brief Get the Arches library environment variable.

This is the same as:

\code
QString s=System::GetEnv("ARCHESLIBDIR");
\endcode
*/
QString System::GetArchesLib()
{
  QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
  return env.value("ARCHESLIBDIR");
}

/*!
\brief Get the repository of heightfields.

This is the same as:

\code
QString s=System::GetArchesLib() + QString("/LibTerra/HeightFields/");
\endcode
*/
QString System::GetHeightFieldDir()
{
  return System::GetArchesLib() + QString("/LibTerra/HeightFields/");
}

/*!
\brief Get the desktop environment.
*/
QString System::GetDesktop()
{
  return QString("C:/Users/") + System::GetEnv("USERNAME") + QString("/Desktop/");
}

/*!
\brief Get the resource file name given an ENV variable and the subfolder/name of the file.

Returns an empty string if the resource has not been found. Searches also in the current folder if the env variable is not defined.
\param env Name of the environment variable.
\param name Name of the file.
*/

QString System::GetResource(const QString& env, const QString& name)
{
  QString envvalue = System::GetEnv(env);
  if (envvalue.isEmpty())
  {
    // The env variable is not defined, we will try to find the file in the current directory
    envvalue = ".";
    if (name.isEmpty())
      return envvalue; // a precise filename has not been asked, we just fall back to the current directory

    QString fullPath = envvalue + name;
    QFileInfo check_file(fullPath);
    if (check_file.exists())
    {
      return fullPath;
    }
  }
  else
  {
    return envvalue + name;
  }
  return QString();
}

/*
\brief Create a string with the formatted elapsed time.
\param time Timer.
*/
QString System::Elapsed(const QElapsedTimer& time)
{
  // Elapsed time is in milliseconds
  int e = time.elapsed();

  // Milliseconds
  int ms = e % 1000;

  // Seconds
  e /= 1000;
  int s = e % 60;

  // Minutes
  e /= 60;
  int m = e;
  return QString("%1:%2:%3").arg(m).arg(s, 2, 10, QChar('0')).arg(ms, 3, 10, QChar('0'));
}

/*
\brief Show elapsed time on console.
\param time Timer.
*/
void System::ShowElapsed(const QElapsedTimer& time)
{
  QString qs = Elapsed(time);
  std::string s = qs.toLocal8Bit().constData();
  std::cout << s << std::endl;
}

/*
\brief Format an integer with spaces.
\param n Integer.
*/
QString System::LongInteger(int n)
{
  if (n < 1000)
  {
    return QString("%1").arg(n);
  }
  else if (n < 1000000)
  {
    return QString("%1.%2").arg(n / 1000).arg(n % 1000, 3, 10, QChar('0'));
  }
  else if (n < 1000000000)
  {
    return QString("%1.%2.%3").arg(n / 1000000).arg((n / 1000) % 1000, 3, 10, QChar('0')).arg(n % 1000, 3, 10, QChar('0'));
  }
  else
  {
    return QString("%1.%2.%3.%4").arg(n / 1000000000).arg((n / 1000000) % 1000, 3, 10, QChar('0')).arg((n / 1000) % 1000, 3, 10, QChar('0')).arg(n % 1000, 3, 10, QChar('0'));
  }
}

/*
\brief Format an 64 bit integer with spaces.
\param n Integer.
*/
QString System::VeryLongInteger(long long n)
{
  if (n < 1000)
  {
    return QString("%1").arg(n);
  }
  else if (n < 1000000)
  {
    return QString("%1.%2").arg(n / 1000).arg(n % 1000, 3, 10, QChar('0'));
  }
  else if (n < 1000000000)
  {
    return QString("%1.%2.%3").arg(n / 1000000).arg((n / 1000) % 1000, 3, 10, QChar('0')).arg(n % 1000, 3, 10, QChar('0'));
  }
  else
  {
    return QString("%1.%2.%3.%4").arg(n / 1000000000).arg((n / 1000000) % 1000, 3, 10, QChar('0')).arg((n / 1000) % 1000, 3, 10, QChar('0')).arg(n % 1000, 3, 10, QChar('0'));
  }
}

/*
\brief Set s string according to date and time.

This function is useful for labeling screenshots.
\param n Integer.
*/
QString System::DateTime()
{
  // Date and time
  QDate date = QDate::currentDate();
  QTime time = QTime::currentTime();

  QString s = QString("%1%2%3-%4%5%6")
    .arg(date.year(), 4)
    .arg(date.month(), 2, 10, QChar('0'))
    .arg(date.day(), 2, 10, QChar('0'))
    .arg(time.hour(), 2, 10, QChar('0'))
    .arg(time.minute(), 2, 10, QChar('0'))
    .arg(time.second(), 2, 10, QChar('0'));
  return s;
}

/*!
\brief Compose two images.
\param image Image that will be modified.
\param over Overlay image.
\parame CompositionMode
*/
void System::Compose(QImage& image, const QImage& over, QPainter::CompositionMode mode)
{
  // Compositing
  QPainter p(&image);
  p.setCompositionMode(mode);
  p.drawImage(0, 0, over);
}

/*!
\brief Compose two images.
\param image Image that will be modified.
\param over Overlay image.
\parame CompositionMode
*/
QImage System::Compose(const QImage& image, const QImage& over, QPainter::CompositionMode mode)
{
  QImage out = image;
  // Compositing
  QPainter p(&out);
  p.setCompositionMode(mode);
  p.drawImage(0, 0, over);
  return out;
}

/*!
\brief Check available memory.
*/
QString System::Memory()
{

#ifdef _WIN32

    const int k = 1024;

    MEMORYSTATUSEX statex;
    statex.dwLength = sizeof(statex);
    GlobalMemoryStatusEx(&statex);

    QString s =
        "Memory   " + LongInteger(statex.dwMemoryLoad / k) +
        "\nPhysical " + LongInteger(statex.ullTotalPhys / k) +
        "\n    Free " + LongInteger(statex.ullAvailPhys / k) +
        "\nVirtual  " + LongInteger(statex.ullTotalVirtual / k) +
        "\n    Free " + LongInteger(statex.ullAvailVirtual / k) +
        "\nExtended " + LongInteger(statex.ullAvailExtendedVirtual / k);

    return s;

#else //make a Linux/Mac version

    return QString("Memory info not available on Linux/MacOS");

#endif
}
