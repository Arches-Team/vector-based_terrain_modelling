// Maya

#include "libs/maya.h"

/*!
\class MayaGpuSet maya.h
\brief This class implements a geometric instance linked to a set of frames.

\ingroup MayaGpuGroup
*/

/*!
\brief Create a set of geometric instances on the GPU.
\param mis Set of geometric objects.
*/
MayaGpuSet::MayaGpuSet(const MayaGeometrySet& mis) : MayaGpu(mis)
{
  // Init
  frames = mis.GetFrames();

  // Init Buffers  
  buffer_instances_position = 0;
  buffer_instances_scale = 0;
  buffer_instances_rotation = 0;

  buffer_Vao = 0;

  // Statistics from argument geometry set (which include instance stats)
  statistics = mis.GetStatistics();

  InitInstances();
}

/*!
\brief Delete buffers on the GPU.
*/
void MayaGpuSet::DeleteBuffers()
{
  if (buffer_Vao != 0) glDeleteVertexArrays(1, &buffer_Vao);

  // Free the buffer memory
  if (buffer_instances_position != 0) glDeleteBuffers(1, &buffer_instances_position);
  if (buffer_instances_scale != 0)    glDeleteBuffers(1, &buffer_instances_scale);
  if (buffer_instances_rotation != 0) glDeleteBuffers(1, &buffer_instances_rotation);

  MayaGpu::DeleteBuffers();
}

/*!
\brief Get pointer.
*/
void glVertexAttribPointerMatrix(int loc, int line, int col, int stride, int offset)
{
  if (stride == 0)
  {
    stride = line * col * sizeof(GLfloat);
  }
  for (int l = 0; l < line; l++)
  {
    glVertexAttribPointer(loc + l, col, GL_FLOAT, GL_FALSE, stride, (GLvoid*)(offset + (l * col * sizeof(GLfloat)))); // size, type, stride, pointer/offset
  }
}

/*!
\brief Modify the rate at which vertex attributes advance during instanced rendering.
\param loc The index of the generic vertex attribute.
\param line Number of lines.
\param divisor Specify the number of instances that will pass between updates of the generic attribute at slot index.
*/
void glVertexAttribDivisorMatrix(int loc, int line, int divisor)
{
  for (int l = 0; l < line; l++)
  {
    glVertexAttribDivisor(loc + l, divisor);
  }
}

/*!
\brief Initialize the instances.
*/
void MayaGpuSet::InitInstances()
{
  // Generation des buffers
  glGenBuffers(1, &buffer_instances_position);
  glGenBuffers(1, &buffer_instances_scale);
  glGenBuffers(1, &buffer_instances_rotation);

  if (!(mat.materialNode == VertexColor))
    glGenBuffers(1, &buffer_color);

  // Couleur des deux instances.
  GLfloat* inst_pos = new GLfloat[frames.size() * 3];
  GLfloat* inst_scl = new GLfloat[frames.size() * 3];
  GLfloat* inst_rot = new GLfloat[frames.size() * 3 * 3];
  GLfloat* colo = nullptr;
  if (!(mat.materialNode == VertexColor)) colo = new GLfloat[frames.size() * 3];

  for (int i = 0; i < frames.size(); i++)
  {
    Vector translate = frames[i].T();
    inst_pos[3 * i + 0] = GLfloat(translate[0]);
    inst_pos[3 * i + 1] = GLfloat(translate[1]);
    inst_pos[3 * i + 2] = GLfloat(translate[2]);

    Vector scale = frames[i].S();
    inst_scl[3 * i + 0] = GLfloat(scale[0]);
    inst_scl[3 * i + 1] = GLfloat(scale[1]);
    inst_scl[3 * i + 2] = GLfloat(scale[2]);

    Matrix rotation = frames[i].R();
    inst_rot[9 * i + 0] = rotation(0, 0);
    inst_rot[9 * i + 1] = rotation(0, 1);
    inst_rot[9 * i + 2] = rotation(0, 2);
    inst_rot[9 * i + 3] = rotation(1, 0);
    inst_rot[9 * i + 4] = rotation(1, 1);
    inst_rot[9 * i + 5] = rotation(1, 2);
    inst_rot[9 * i + 6] = rotation(2, 0);
    inst_rot[9 * i + 7] = rotation(2, 1);
    inst_rot[9 * i + 8] = rotation(2, 2);

    if (!(mat.materialNode == VertexColor))
    {
      // Declare one color/object if we don't use the vertexcolor
      colo[3 * i + 0] = GLfloat(mat.diffuse[0]);
      colo[3 * i + 1] = GLfloat(mat.diffuse[1]);
      colo[3 * i + 2] = GLfloat(mat.diffuse[2]);
    }
  }

  // Instances Position
  glBindBuffer(GL_ARRAY_BUFFER, buffer_instances_position);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * frames.size() * 3, inst_pos, GL_DYNAMIC_DRAW);
  // Instances Scale
  glBindBuffer(GL_ARRAY_BUFFER, buffer_instances_scale);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * frames.size() * 3, inst_scl, GL_DYNAMIC_DRAW);
  // Instances Rotation
  glBindBuffer(GL_ARRAY_BUFFER, buffer_instances_rotation);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * frames.size() * 3 * 3, inst_rot, GL_DYNAMIC_DRAW);
  // Instances Color
  if (!(mat.materialNode == VertexColor))
  {
    glBindBuffer(GL_ARRAY_BUFFER, buffer_color);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * frames.size() * 3, colo, GL_DYNAMIC_DRAW);
  }

  delete[]inst_pos;
  delete[]inst_scl;
  delete[]inst_rot;
  if (!(mat.materialNode == VertexColor)) delete[]colo;
}



/*!
\brief Initialize vertex array object.
\param GP Parameters.
*/
void MayaGpuSet::InitVAO(GpuParameter& GP)
{
  glGenVertexArrays(1, &buffer_Vao);
  glBindVertexArray(buffer_Vao);

  // Ouverture des positions des attributs
  glEnableVertexAttribArray(GP.uniforms.loc_position);
  glEnableVertexAttribArray(GP.uniforms.loc_normal);
  glEnableVertexAttribArray(GP.uniforms.loc_color);
  if (mat.materialNode == UVMapping && (mat.shaderNode == ShaderTextureUV || mat.shaderNode == ShaderTextureUVWireframe)) glEnableVertexAttribArray(GP.uniforms.loc_UV);
  glEnableVertexAttribArray(GP.uniforms.loc_inst_pos);
  glEnableVertexAttribArray(GP.uniforms.loc_inst_scl);
  glEnableVertexAttribArray(GP.uniforms.loc_inst_rot + 0);
  glEnableVertexAttribArray(GP.uniforms.loc_inst_rot + 1);
  glEnableVertexAttribArray(GP.uniforms.loc_inst_rot + 2);

  // Bind et config du buffer des positions, normals et couleurs
  glBindBuffer(GL_ARRAY_BUFFER, this->buffer_vertex);
  glVertexAttribPointer(GP.uniforms.loc_position, 3, GL_FLOAT, GL_FALSE, 0, nullptr); // size, type, stride, pointer/offset
  glBindBuffer(GL_ARRAY_BUFFER, this->buffer_normal);
  glVertexAttribPointer(GP.uniforms.loc_normal, 3, GL_FLOAT, GL_FALSE, 0, nullptr); // size, type, stride, pointer/offset
  glBindBuffer(GL_ARRAY_BUFFER, this->buffer_color);
  glVertexAttribPointer(GP.uniforms.loc_color, 3, GL_FLOAT, GL_FALSE, 0, nullptr); // size, type, stride, pointer/offset

  if (mat.materialNode == UVMapping && (mat.shaderNode == ShaderTextureUV || mat.shaderNode == ShaderTextureUVWireframe))
  {
    glBindBuffer(GL_ARRAY_BUFFER, this->buffer_UVmap);
    glVertexAttribPointer(GP.uniforms.loc_UV, 2, GL_FLOAT, GL_FALSE, 0, nullptr); // size, type, stride, pointer/offset
  }

  // Bind des Instances 
  glBindBuffer(GL_ARRAY_BUFFER, buffer_instances_position);
  glVertexAttribPointer(GP.uniforms.loc_inst_pos, 3, GL_FLOAT, GL_FALSE, 0, nullptr); // size, type, stride, pointer/offset
  glBindBuffer(GL_ARRAY_BUFFER, buffer_instances_scale);
  glVertexAttribPointer(GP.uniforms.loc_inst_scl, 3, GL_FLOAT, GL_FALSE, 0, nullptr); // size, type, stride, pointer/offset
  glBindBuffer(GL_ARRAY_BUFFER, buffer_instances_rotation);
  glVertexAttribPointerMatrix(GP.uniforms.loc_inst_rot, 3, 3, 0, 0);

  // Ici on met 1 pour que l'accés au cases des tableaux se fassent 
  // par cran d'instance et non d'un cran par vertex (defaut)
  if (mat.materialNode != VertexColor) { glVertexAttribDivisor(GP.uniforms.loc_color, 1); }
  else { glVertexAttribDivisor(GP.uniforms.loc_color, 0); }
  glVertexAttribDivisor(GP.uniforms.loc_inst_pos, 1);
  glVertexAttribDivisor(GP.uniforms.loc_inst_scl, 1);
  glVertexAttribDivisorMatrix(GP.uniforms.loc_inst_rot, 3, 1);

  // Pour les indices, pas besoin de gl****Pointer car c'est bindé sur GL_ELEMENT_ARRAY_BUFFER
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_index);

  glBindVertexArray(0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

/*!
\brief Refresh frames Instances.
*/
void MayaGpuSet::RefreshFrames()
{
  // Couleur des deux instances.
  GLfloat* inst_pos = new GLfloat[frames.size() * 3];
  GLfloat* inst_scl = new GLfloat[frames.size() * 3];
  GLfloat* inst_rot = new GLfloat[frames.size() * 3 * 3];
  GLfloat* colo = nullptr;

  if (!(mat.materialNode == VertexColor))
  {
    colo = new GLfloat[frames.size() * 3];
  }
  // Static analysis raises warning because allocation is far away from usage
  for (int i = 0; i < frames.size(); i++)
  {
    Vector translate = frames[i].T();
    inst_pos[3 * i + 0] = translate[0];
    inst_pos[3 * i + 1] = translate[1];
    inst_pos[3 * i + 2] = translate[2];

    Vector scale = frames[i].S();
    inst_scl[3 * i + 0] = scale[0];
    inst_scl[3 * i + 1] = scale[1];
    inst_scl[3 * i + 2] = scale[2];

    Matrix rotation = frames[i].R();
    inst_rot[9 * i + 0] = rotation(0, 0);
    inst_rot[9 * i + 1] = rotation(0, 1);
    inst_rot[9 * i + 2] = rotation(0, 2);
    inst_rot[9 * i + 3] = rotation(1, 0);
    inst_rot[9 * i + 4] = rotation(1, 1);
    inst_rot[9 * i + 5] = rotation(1, 2);
    inst_rot[9 * i + 6] = rotation(2, 0);
    inst_rot[9 * i + 7] = rotation(2, 1);
    inst_rot[9 * i + 8] = rotation(2, 2);

    if (!(mat.materialNode == VertexColor))
    {
      colo[3 * i + 0] = 0.8f;
      colo[3 * i + 1] = 0.5f;
      colo[3 * i + 2] = 0.8f;
    }
  }

  // Color
  if (!(mat.materialNode == VertexColor)) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer_color);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * frames.size() * 3, colo, GL_DYNAMIC_DRAW);
  }

  // Instances
  glBindBuffer(GL_ARRAY_BUFFER, buffer_instances_position);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * frames.size() * 3, inst_pos, GL_DYNAMIC_DRAW);
  // Instances
  glBindBuffer(GL_ARRAY_BUFFER, buffer_instances_scale);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * frames.size() * 3, inst_scl, GL_DYNAMIC_DRAW);
  // Instances
  glBindBuffer(GL_ARRAY_BUFFER, buffer_instances_rotation);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * frames.size() * 3 * 3, inst_rot, GL_DYNAMIC_DRAW);

  delete[]inst_pos;
  delete[]inst_scl;
  delete[]inst_rot;
  if (!(mat.materialNode == VertexColor))
  {
    delete[]colo;
  }
}

/*!
\brief Render all instances of an object (Shaders are used).
\param GP Parameters for the GPU.
*/
void MayaGpuSet::Render(GpuParameter& GP)
{
  glUseProgram(GP.shader_program->GetProgram());

  if (frames.size() == 0) return;

  if (buffer_Vao == 0) InitVAO(GP);

  glBindVertexArray(buffer_Vao);

  // Lien vers le matériau
  glBindBuffer(GL_UNIFORM_BUFFER, buffer_material);
  glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(material), &material);
  glBindBufferBase(GL_UNIFORM_BUFFER, 0, buffer_material);

  if (mat.materialNode == UVMapping && (mat.shaderNode == ShaderTextureUV || mat.shaderNode == ShaderTextureUVWireframe))
  {
    glAlphaFunc(GL_GREATER, 0.5);
    glEnable(GL_ALPHA_TEST);

    glActiveTexture(GL_TEXTURE11);// ENV
    glBindTexture(GL_TEXTURE_2D, buffer_texture);
    glProgramUniform1i(GP.shader_program->GetProgram(), GP.uniforms.loc_texture, 11);  // ENV
  }
  if (mat.shaderNode == ShaderWireframe)
  {
    glAlphaFunc(GL_GREATER, 0.5);
    glEnable(GL_ALPHA_TEST);
  }

  // Select shader function
  SetSubroutine(GP);

  // Affichage des triangles 
  glDrawArraysInstanced(GL_TRIANGLES, 0, nt, frames.size()); // type, start, count, prims (ie, deux primitives de trois points)

  if (mat.materialNode == UVMapping && (mat.shaderNode == ShaderTextureUV || mat.shaderNode == ShaderTextureUVWireframe) || mat.shaderNode == ShaderWireframe)
  {
    glDisable(GL_ALPHA_TEST);
    glActiveTexture(0);
  }

  glBindVertexArray(0);
}

/*!
\brief Define the shading program.
\param GP Parameters for the GPU.
*/
void MayaGpuSet::SetSubroutine(GpuParameter& GP) const
{
  // Gestion Wireframe
  if (isWireframe())
  {
    glUniform1i(GP.uniforms.loc_wireframe, 1);
  }
  else
  {
    glUniform1i(GP.uniforms.loc_wireframe, 0);
  }

  switch (mat.shaderNode) {
  case ShaderPhong:
  case ShaderPhongWireframe:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.phongIndex);
    break;
  case ShaderPhongVertexColor:
  case ShaderPhongVertexColorWireframe:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.phongVCIndex);
    break;
  case ShaderTextureUV:
  case ShaderTextureUVWireframe:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.UVTextureModel);
    break;
  case ShaderTerrainTexture3D:
  case ShaderTerrainTexture3DWireframe:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.terrain3DIndex);
    break;
  case ShaderTerrainIGN:
  case ShaderTerrainIGNWireframe:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.ignIndex);
    break;
  case ShaderGrid:
  case ShaderGridWireframe:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.gridIndex);
    break;
  case ShaderNormal:
  case ShaderNormalWireframe:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.normalIndex);
    break;
  case ShaderGooch:
  case ShaderGoochWireframe:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.goochIndex);
    break;
  case ShaderBumpNormal:
  case ShaderBumpNormalWireframe:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.bumpNormalIndex);
    break;
  case ShaderThinSnow:
  case ShaderThinSnowWireframe:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.ThinSnowIndex);
    break;
  case ShaderWireframe:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.WireframeIndex);
    break;
  case ShaderWireframeAspectRatio:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.WireframeAspectIndex);
    break;
  case ShaderGreyCoolWarm:
  case ShaderGreyCoolWarmWireframe:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.GreyCoolWarmIndex);
    break;
  case ShaderAxelTerrain:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.AxelTerrainIndex);
    //glUniform1i(GP.uniforms.loc_applyNormalMap, mat.useNormalMapping);
    //glUniform1f(GP.uniforms.loc_normalMapStrength, mat.normalMappingStrength);
    //glUniform1f(GP.uniforms.loc_textureScaling, mat.textureScaling);
    break;
  case ShaderTriPlanar:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.TriPlanarIndex);
    break;
  default:
    glUniformSubroutinesuiv(GL_FRAGMENT_SHADER, 1, &GP.uniforms.phongIndex);
  }
}
/*!
\brief Render the bounding box.
*/
void MayaGpuSet::RenderBBox()
{
  if (frames.size() == 0) return;

  // Render bounding Box
  for (int i = 0; i < frames.size(); i++)
  {
    glPushMatrix();

    Vector t = frames[i].T();
    glTranslated(t[0], t[1], t[2]);

    Vector r = frames[i].R().GetRotationAngles();
    glRotated(Math::RadianToDegree(r[0]), 1.0, 0.0, 0.0);
    glRotated(Math::RadianToDegree(r[1]), 0.0, 1.0, 0.0);
    glRotated(Math::RadianToDegree(r[2]), 0.0, 0.0, 1.0);

    Vector s = frames[i].S();
    glScaled(s[0], s[1], s[2]);

    MayaGpu::RenderBBox();

    glPopMatrix();
  }
}

/*
void MayaGpuSet::Clear()
{
// statistics=MayaStatistics();

frames.clear();

// Free the memory
/*glDeleteBuffers( 1, &buffer_vertex);
glDeleteBuffers( 1, &buffer_normal);
glDeleteBuffers( 1, &buffer_color);
glDeleteBuffers( 1, &buffer_UVmap);
glDeleteBuffers( 1, &buffer_texture);
glDeleteBuffers( 1, &buffer_index);
glDeleteBuffers( 1, &buffer_instances_position);
glDeleteBuffers( 1, &buffer_instances_scale);
glDeleteBuffers( 1, &buffer_instances_rotation);
}*/

/*!
\brief Clear all frames.
*/
void MayaGpuSet::ClearFrames()
{
  // Statistics
  statistics.ClearInstances();

  frames.clear();
}

/*!
\brief Compute statistics.
*/
MayaStatistics MayaGpuSet::GetStatistics() const
{
  // return MayaStatistics(1,MayaGpuSet::nt,MayaGpuSet::nv,frames.size(),MayaGpuSet::nt*frames.size(),MayaGpuSet::nv*frames.size());
  return statistics;
}

/*!
\brief Check if the set should be rendered as wireframe.
*/
bool MayaGpuSet::isWireframe() const
{
  return mat.shaderNode == ShaderPhongWireframe ||
    mat.shaderNode == ShaderPhongVertexColorWireframe ||
    mat.shaderNode == ShaderTextureUVWireframe ||
    mat.shaderNode == ShaderTerrainTexture3DWireframe ||
    mat.shaderNode == ShaderTerrainIGNWireframe ||
    mat.shaderNode == ShaderToonWireframe ||
    mat.shaderNode == ShaderGoochWireframe ||
    mat.shaderNode == ShaderGridWireframe ||
    mat.shaderNode == ShaderNormalWireframe ||
    mat.shaderNode == ShaderBumpNormalWireframe ||
    mat.shaderNode == ShaderGreyCoolWarmWireframe ||
    mat.shaderNode == ShaderThinSnowWireframe ||
    mat.shaderNode == ShaderWireframeAspectRatio;
}