// Camera

#pragma once

#include "libs/box.h"

// Core camera class
class Camera
{
protected:
  Vector eye = -Vector::X;       //!< Eye.
  Vector at = Vector::Null;        //!< Look at point.
  Vector up = Vector::Z;        //!< Up vector.
  Vector view = Vector::X;      //!< View unit vector.

  double width = 1.0;     //!< Screen width.
  double height = 1.0;    //!< Screen height.

  double cah;       //!< Camera aperture horizontal. 
  double cav;       //!< Camera aperture vertical.
  double fl;        //!< Focal length.

  double nearplane = 1.0; //!< Near plane.
  double farplane = 100000.0;  //!< Far plane.

  // Vector horizontal = Normalized(view / Up());
   //Vector vertical = Normalized(horizontal / view);

public:
  Camera();
  explicit Camera(const Vector&, const Vector & = Vector::Null, const Vector & = Vector::Z, const double& = 1.0, const double& = 1.0, const double& = 1.0, const double& = 100000.0);
  explicit Camera(const Vector&, const Vector&, const Vector&, const double&, const double& = 1.0, const double& = 100000.0);

  void Translate(const Vector&);
  void Rotate(const Matrix&);

  Vector At() const;
  Vector Eye() const;
  Vector Up() const;
  Vector View() const;

  Frame GetFrame() const;
  double GetNear() const;
  double GetFar() const;
  double GetCameraApertureH() const;
  double GetCameraApertureV() const;
  double GetFocalLength() const;
  double GetAngleOfViewH() const;
  double GetAngleOfViewV(double, double) const;

  void Step(const double&);

  bool InsideFrustum(const Vector&) const;
  bool InsideFrustum(const Box&) const;

  void LeftRight(const double&);
  void Vertical();
  void UpDown(const double&);
  void SideWay(const double&);
  void BackForth(const double&, bool = false);
  void LeftRightRound(const double&);
  void UpDownRound(const double&);

  void LeftRightFPS(const double&);
  void UpDownRoundFPS(const double&);

  void SlideHorizontal(const double&);
  void SetAt(const Vector&);
  void SetEye(const Vector&);
  void SetPlanes(const double&, const double&);

  // Move camera in a plane
  void UpDownVertical(const double&);
  void LeftRightHorizontal(const double&);

  QString ToString(int = 6) const;

  friend std::ostream& operator<<(std::ostream&, const Camera&);

  // Pixel and sub-pixel sampling
  Ray PixelToRay(int, int, int, int) const;
  Ray PixelToRay(int, int, int, int, int, int, int) const;
  bool VectorToPixel(const Vector&, double&, double&, int, int) const;
public:
  static Camera View(const Box&);
};

//! Returns the look-at point.
inline Vector Camera::At() const
{
  return at;
}

//! Returns the eye point.
inline Vector Camera::Eye() const
{
  return eye;
}

//! Returns the up point.
inline Vector Camera::Up() const
{
  return up;
}

/*!
\brief Returns the viewing direction (not normalized).
*/
inline Vector Camera::View() const
{
  return at - eye;
}
