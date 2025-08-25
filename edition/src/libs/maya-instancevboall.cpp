// Maya

#include "libs/maya.h"

/*!
\class MayaGpuAll maya.h
\brief This class implements a set of mesh instances optimized for the GPU.

\ingroup MayaGpuGroup
*/

/*!
\brief Create a set of instances on the GPU from a set of instances.
\param mia Set of instances.
*/
MayaGpuAll::MayaGpuAll(const MayaGeometryAll& mia)
{
  for (QMap<QString, MayaGeometrySet>::const_iterator i = mia.instances.begin(); i != mia.instances.end(); i++)
  {
    instances[i.key()] = MayaGpuSet(*i);
  }
}

/*!
\brief Render the set of instances.
\param gpu Set of parameters.
*/
void MayaGpuAll::Render(GpuParameter& gpu)
{
  // Activation du shader
  glUseProgram(gpu.shader_program->GetProgram());

  for (QMap<QString, MayaGpuSet>::iterator i = instances.begin(); i != instances.end(); i++)
  {
    i->Render(gpu);
  }

  // Desactivation Shader (pour rendre Background, ...)
  glUseProgram(0);
}

/*!
\brief Render the set of instances.
*/
void MayaGpuAll::RenderBBox()
{
  for (QMap<QString, MayaGpuSet>::iterator i = instances.begin(); i != instances.end(); i++)
  {
    i->RenderBBox();
  }
}

/*!
\brief Clear all instances.
*/
void MayaGpuAll::Clear()
{
  for (QMap<QString, MayaGpuSet>::iterator i = instances.begin(); i != instances.end(); i++)
  {
    i->DeleteBuffers();
  }
  instances.clear();
}

/*!
\brief clear all frames, preserving the instances.
*/
void MayaGpuAll::ClearFrames()
{
  for (QMap<QString, MayaGpuSet>::iterator i = instances.begin(); i != instances.end(); i++)
  {
    i->ClearFrames();
  }
}

/*!
\brief Add an instance set to the scene.

Checks if the object already exists in the collection of models: if so,
it simply creates a new reference to it, otherwise it creates a new entry.
\param ins The instance.
*/
void MayaGpuAll::Append(const MayaGeometrySet& ins)
{
  // Check if static shape already exists
  if (instances.find(ins.GetName()) != instances.end())
  {
    for (int j = 0; j < ins.count(); j++)
    {
      instances[ins.GetName()].Append(ins.GetFrameScaled(j));
    }
    instances[ins.GetName()].RefreshFrames();
  }
  else
  {
    // Add static shape to the set of static objects and creates an instance
    instances[ins.GetName()] = MayaGpuSet(ins);
  }
}

/*!
\brief Append a set of instances with the existing one.
\param scene The set of instances.
*/
void MayaGpuAll::Append(const MayaGeometryAll& scene)
{
  for (QMap<QString, MayaGeometrySet>::const_iterator i = scene.instances.begin(); i != scene.instances.end(); i++)
  {
    Append(*i);
  }
}

/*!
\brief Replace a set of instances with the existing one.
\param scene The set of instances.
*/
void MayaGpuAll::Replace(const MayaGeometryAll& scene)
{
  for (QMap<QString, MayaGeometrySet>::const_iterator i = scene.instances.begin(); i != scene.instances.end(); i++)
  {
    Remove(i.key());
    Append(*i);
  }
}

/*!
\brief Clears one instance by name.
\param n Name of the instance.
*/
void MayaGpuAll::Remove(const QString& n)
{
  if (instances.find(n) != instances.end())
  {
    instances[n].DeleteBuffers();
    instances.remove(n);
  }
}

/*!
\brief Update the material.
*/
void MayaGpuAll::UpdateMaterial(MayaGeometryAll& mga)
{
  for (QMap<QString, MayaGeometrySet>::const_iterator i = mga.instances.begin(); i != mga.instances.end(); i++)
  {
    if (instances.find(i.key()) != instances.end())
    {
      instances[i.key()].SetMaterial(i->GetMaterial());
    }
  }
}

/*!
\brief Compute some statistics.

This function gathers some statistics and pack them into the structure.
*/
MayaStatistics MayaGpuAll::GetStatistics() const
{
  MayaStatistics s;
  for (QMap<QString, MayaGpuSet>::const_iterator i = instances.begin(); i != instances.end(); i++)
  {
    s += i->GetStatistics();
  }
  return s;
}

/*!
\brief Retourne la boite englobante de la scène
*/
Box MayaGpuAll::getBox()
{
  if (instances.size() > 0)
  {
    Box b = instances.begin()->getBox();
    QMap<QString, MayaGpuSet>::iterator i = instances.begin();
    i++;
    for (; i != instances.end(); i++)
    {
      b = Box(b, i->getBox());
    }
    return b;
  }
  else
  {
    return Box::Null;
  }
}
