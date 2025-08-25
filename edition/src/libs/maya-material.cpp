// MayaMaterial
#include "libs/mayageometry.h"

/*!
\class MayaMaterial maya.h
\brief Set of materials handled by the Maya framework.

\ingroup MayaCore
*/

MayaMaterial MayaMaterial::None = MayaMaterial(ShaderPhong, Color(0.1, 0.5, 0.1, 1.0), Color(0.5, 0.2, 0.2, 1.0), Color(0.1, 0.1, 0.1, 1.0), 50.0);

MayaMaterial MayaMaterial::Normal = MayaMaterial(ShaderNormal, Color(0.3, 0.3, 0.3, 1.0), Color(0.5, 0.4, 0.2, 1.0), Color(0.1, 0.1, 0.1, 1.0), 50.0);

MayaMaterial MayaMaterial::NormalWire = MayaMaterial(ShaderNormalWireframe, Color(0.3, 0.3, 0.3, 1.0), Color(0.5, 0.4, 0.2, 1.0), Color(0.1, 0.1, 0.1, 1.0), 50.0);

MayaMaterial MayaMaterial::MayaTransparent = MayaMaterial(ShaderGreyCoolWarm, Color(0.05, 0.05, 0.05, 0.25), Color(0.5, 0.5, 0.5, 0.25), Color(0.1, 0.1, 0.1, 0.25), 50.0);

MayaMaterial MayaMaterial::MayaTransparentWire = MayaMaterial(ShaderGreyCoolWarmWireframe, Color(0.05, 0.05, 0.05, 0.25), Color(0.5, 0.5, 0.5, 0.25), Color(0.1, 0.1, 0.1, 0.25), 50.0);

MayaMaterial MayaMaterial::WireframeAspectRatio = MayaMaterial(ShaderWireframeAspectRatio, Color(0.05, 0.05, 0.05, 0.25), Color(0.5, 0.5, 0.5, 0.25), Color(0.1, 0.1, 0.1, 0.25), 50.0);

MayaMaterial MayaMaterial::GreenColor = MayaMaterial(ShaderGreyCoolWarm, Color(0.5, 0.8, 0.6, 1.0), Color(0.1, 0.1, 0.1, 1.0), Color(0.0, 0.0, 0.0, 1.0), 50.0);

MayaMaterial MayaMaterial::YellowColor = MayaMaterial(ShaderGreyCoolWarm, Color(0.8, 0.7, 0.5, 1.0), Color(0.1, 0.1, 0.1, 1.0), Color(0.0, 0.0, 0.0, 1.0), 50.0);

MayaMaterial MayaMaterial::RedColor = MayaMaterial(ShaderGreyCoolWarm, Color(0.85, 0.5, 0.4, 1.0), Color(0.1, 0.1, 0.1, 1.0), Color(0.0, 0.0, 0.0, 1.0), 50.0);

MayaMaterial MayaMaterial::BlueColor = MayaMaterial(ShaderGreyCoolWarm, Color(0.3, 0.25, 0.7, 1.0), Color(0.1, 0.1, 0.1, 1.0), Color(0.0, 0.0, 0.0, 1.0), 50.0);

MayaMaterial MayaMaterial::ColdWarm = MayaMaterial(ShaderGooch, Color(0.8, 0.7, 0.5, 1.0), Color(0.1, 0.1, 0.1, 1.0), Color(0.0, 0.0, 0.0, 1.0), 50.0);

MayaMaterial MayaMaterial::TriPlanarRock = MayaMaterial(ShaderTerrainTexture3D, Color(0.05, 0.05, 0.05, 0.25), Color(0.5, 0.5, 0.5, 0.25), Color(0.1, 0.1, 0.1, 0.25), 50.0);

MayaMaterial MayaMaterial::TriPlanar = MayaMaterial(ShaderTriPlanar, Color(0.3, 0.3, 0.3, 1.0), Color(0.5, 0.4, 0.2, 1.0), Color(0.1, 0.1, 0.1, 1.0), 50.0);

/*!
\brief Create a material.

\param shaderNode Shader node type.
\param ambient Ambient color.
\param diffuse Diffuse color.
\param specular Specular color.
\param shininess Shininess, i.e. exponent.
\param textIm Texture image.
\param textImAlpha Transparency texture image.
*/
MayaMaterial::MayaMaterial(ShaderNode shaderNode, const Color& ambient, const Color& diffuse, const Color& specular, const double& shininess, const QImage textIm, const QImage& textImAlpha)
{
  MayaMaterial::shaderNode = shaderNode;
  MayaMaterial::materialNode = MaterialNode::None;

  MayaMaterial::ambient = ambient;
  MayaMaterial::diffuse = diffuse;
  MayaMaterial::specular = specular;
  MayaMaterial::shininess = shininess;

  MayaMaterial::texture.albedo = textIm;
  MayaMaterial::texture.opacity = textImAlpha;

  MayaMaterial::name = QString("undefined");
}

/*!
\brief Create a simple color material.

Diffuse and specular colors are set to 0.1.

\param c Ambient color.
*/
MayaMaterial MayaMaterial::SimpleColor(const Color& c)
{
  return MayaMaterial(ShaderPhong, c, Color(0.1, 0.1, 0.1, 1.0), Color(0.1, 0.1, 0.1, 1.0), 24.0);
}
