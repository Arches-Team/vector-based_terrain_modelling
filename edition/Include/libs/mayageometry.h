// Maya  

#pragma once

#include <QtCore/QTextStream> 
#include <QtCore/QVector>
#include <QtCore/QMap>


#include "libs/plane.h"
#include "libs/color.h"
#include "libs/scalarfield.h"
#include "libs/vectorfield.h"

#include "libs/mesh.h"
#include "libs/meshcolor.h"

class Voxel;

class MayaStatistics {
protected:
  int s[6]; //!< Array of integers storing statistics.
public:
  explicit MayaStatistics(int = 0, int = 0, int = 0, int = 0, int = 0, int = 0);
  //! Empty
  ~MayaStatistics() {}

  int GetObjects() const;
  int GetObjectsVertices() const;
  int GetObjectsTriangles() const;
  int GetInstances() const;
  int GetInstancesVertices() const;
  int GetInstancesTriangles() const;

  // Assignment operators
  MayaStatistics& operator+= (const MayaStatistics&);

  // Binary operators for combining statistics
  friend MayaStatistics operator+ (const MayaStatistics&, const MayaStatistics&);

  QString GetText() const;

  // Utilities
  void ClearInstances();
  void AddInstances(int = 1);
};

//! Destructive addition.
inline MayaStatistics& MayaStatistics::operator+= (const MayaStatistics& statistics)
{
  for (int i = 0; i < 6; i++)
  {
    s[i] += statistics.s[i];
  }

  return *this;
}

/*!
\brief Return the number of unique objects.
*/
inline int MayaStatistics::GetObjects() const
{
  return s[0];
}
/*!
\brief Return the number of vertices of unique objects.
*/
inline int MayaStatistics::GetObjectsVertices() const
{
  return s[1];
}
/*!
\brief Return the number of triangles of unique objects.
*/
inline int MayaStatistics::GetObjectsTriangles() const
{
  return s[2];
}
/*!
\brief Return the number of instances.
*/
inline int MayaStatistics::GetInstances() const
{
  return s[3];
}
/*!
\brief Return the number of vertices.
*/
inline int MayaStatistics::GetInstancesVertices() const
{
  return s[4];
}
/*!
\brief Return the number of triangless.
*/
inline int MayaStatistics::GetInstancesTriangles() const
{
  return s[5];
}

class MayaIndexVertexData
{
protected:
  int iVertex;//!< Index of Vertex.
  int iNormal;//!< Index of Normal
  int iColor; //!< Index of Color
  int iUVCoord;//!< Index of UVCoord

public:
  MayaIndexVertexData(int = -1, int = -1, int = -1, int = -1);
  //! Empty
  ~MayaIndexVertexData() {}

  // Setters
  void setIV(int i) { iVertex = i; }
  void setIN(int i) { iNormal = i; }
  void setIC(int i) { iColor = i; }
  void setIUV(int i) { iUVCoord = i; }

  // Getters
  int getIV() const { return iVertex; }
  int getIN() const { return iNormal; }
  int getIC() const { return iColor; }
  int getIUV() const { return iUVCoord; }
};


/*!
\brief Create a index structure.
\param v, c, n, uv Vertex, normal, color and texture indexes.
*/
inline MayaIndexVertexData::MayaIndexVertexData(int v, int n, int c, int uv)
{
  iVertex = v;
  iNormal = n;
  iColor = c;
  iUVCoord = uv;
}

enum ShaderNode
{
  ShaderPhong = 0,
  ShaderPhongVertexColor,
  ShaderTextureUV,
  ShaderTerrainTexture3D,
  ShaderTerrainIGN,
  ShaderToon,
  ShaderGooch,
  ShaderGrid,
  ShaderWireframe,
  ShaderNormal,
  ShaderBumpNormal,
  ShaderThinSnow,
  ShaderGreyCoolWarm,
  ShaderGreyCoolWarmWireframe,
  ShaderWireframeAspectRatio,
  ShaderAxelTerrain,

  // Wireframe
  ShaderPhongWireframe,
  ShaderPhongVertexColorWireframe,
  ShaderTextureUVWireframe,
  ShaderTerrainTexture3DWireframe,
  ShaderTerrainIGNWireframe,
  ShaderToonWireframe,
  ShaderGoochWireframe,
  ShaderGridWireframe,
  ShaderNormalWireframe,
  ShaderBumpNormalWireframe,
  ShaderThinSnowWireframe,

  // Triplanar
  ShaderTriPlanar
};

enum MaterialNode
{
  None = 0,
  VertexColor,
  UVMapping
};


// PBR Texture https://www.a23d.co/blog/different-maps-in-pbr-textures/
class MayaTexturePBR {
public:
  QImage albedo;
  QImage normal;
  QImage opacity;

  QImage roughness;
  QImage metalness;
  QImage specular;
  QImage height;
  QImage ambientOcclusion;
  QImage refraction;
  QImage emissive;

public:
  void clear()
  {
    albedo = QImage();
    normal = QImage();
    opacity = QImage();

    roughness = QImage();
    metalness = QImage();
    specular = QImage();
    height = QImage();
    ambientOcclusion = QImage();
    refraction = QImage();
    emissive = QImage();
  }
};

class MayaMaterial
{
public:
  ShaderNode shaderNode;

  QString name;         //!< name of the material (for export purpose)

  Color ambient;				//!< Ambient color.
  Color diffuse;				//!< Diffuse color.
  Color specular;				//!< Specular color.
  double shininess;			//!< Specular exponent.

  MayaTexturePBR texture;

  MaterialNode materialNode;
public:
  //! Empty
  MayaMaterial() {}
  MayaMaterial(ShaderNode, const Color&, const Color&, const Color & = Color(0.1, 0.1, 0.1, 1.0), const double& = 20.0, const QImage = QImage(), const QImage & = QImage());
  void SetName(const QString& materialName) { name = materialName; }

  //! Empty
  ~MayaMaterial() {}

  static MayaMaterial SimpleColor(const Color&);
public:
  static MayaMaterial None;						//!< Default material.
  static MayaMaterial Normal;					//!< Normal shading based material.
  static MayaMaterial NormalWire;				//!< Normal shading based material.
  static MayaMaterial MayaTransparent;			//!< Else.
  static MayaMaterial MayaTransparentWire;		//!< Else.
  static MayaMaterial WireframeAspectRatio;		//!< A Wireframe shading with an emphasis on the degenerate triangles (in red).
  static MayaMaterial GreenColor;				//!< Normal shading based material.
  static MayaMaterial YellowColor;				//!< Normal shading based material.
  static MayaMaterial RedColor;				//!< Normal shading based material.
  static MayaMaterial BlueColor;				//!< Normal shading based material.
  static MayaMaterial ColdWarm;				//!< Gooch shading based material.
  static MayaMaterial TriPlanarRock;					//!< Rock.
  static MayaMaterial TriPlanar;					//!< Tri-planar projection.
};

// Forward declaration
class MayaGeometryAll;

// Base geometry class
class MayaGeometry
{
public:
  QString name;            //!< Signature of the object.
  QVector<Vector> vertex;  //!< Array of vertices.
  QVector<Vector> normal;  //!< Array of normals.
  QVector<Vector> color;   //!< Array of colors.
  QVector<Vector2> UVmap;   //!< Array of UVs.
  QVector<MayaIndexVertexData> indexes; //!< Array of integers referencing vertices, normals, colors and UV. 
  MayaMaterial mat;        //!< Material of the object.

public:
  // Constructors
  explicit MayaGeometry(const QString& = QString("None"), const MayaMaterial & = MayaMaterial::None);
  MayaGeometry(const QString&, const QVector<Vector>&, const QVector<Vector>&, const QVector<MayaIndexVertexData>&, const MayaMaterial & = MayaMaterial::None);
  MayaGeometry(const QString&, const QVector<Vector>&, const QVector<Vector>&, const QVector<Vector>&, const QVector<MayaIndexVertexData>&, const MayaMaterial & = MayaMaterial::None);
  MayaGeometry(const QString&, const QVector<Vector>&, const QVector<Vector>&, const QVector<Vector2>&, const QVector<MayaIndexVertexData>&, const MayaMaterial & = MayaMaterial::None);
  MayaGeometry(const QString&, const Mesh&, const MayaMaterial & = MayaMaterial::None);
  MayaGeometry(const QString&, const MeshColor&, const MayaMaterial & = MayaMaterial::None);
  MayaGeometry(const QString&, const Mesh2&, const MayaMaterial & = MayaMaterial::None);
  MayaGeometry(const QString&, const QVector<Triangle>&, const MayaMaterial & = MayaMaterial::None);
  ~MayaGeometry();

  MayaGeometry(const QString&, const Vector*, int, const Vector*, int, const int*, int, const MayaMaterial & = MayaMaterial::None);

  // Incremental geometry creation
  void AddTriangle(const Vector&, const Vector&, const Vector&);
  void AddSmoothTriangle(const Vector&, const Vector&, const Vector&, const Vector&, const Vector&, const Vector&);

  void SetName(const QString&);
  void SetMaterial(const MayaMaterial&);

  void Merge(const MayaGeometry&);

  // Getters
  Triangle GetTriangle(int) const;
  QVector<Triangle> GetTriangles() const;

  QString GetName() const;
  MayaMaterial GetMaterial() const;
  Box GetBox() const;

  MayaGeometry& Transform(const FrameScaled&);
  MayaGeometry& InverseTransform(const FrameScaled&);

  // Static objects

  static MayaGeometryAll CreateVectorField(const VectorField&, const double& = 1.0);

  void WriteCpp(const QString & = "Class");

  // Load and save functions
  bool Load_OBJ(const QString&);
  bool find_MTL_Address(const QString& url, const QString& name, QString& nameR, QString& nameA);
  bool Save_OBJ(const QString&) const;
  bool Save_OBJ_forVUE(const QString& url) const;
  bool Save_OBJ(QTextStream&, int = 1, int = 1, int = 1) const;
  bool Save_PLY(const QString&) const;
  bool Save_PBRT(QTextStream&) const;

  MayaStatistics GetStatistics() const;
  QString GetText(int = 0, bool = false) const;
  const Vector& GetVertex(int, int) const;
  const Vector& GetVertex(int) const;
  const Vector& GetNormal(int, int) const;
  const Vector& GetNormal(int) const;
  const Vector& GetColor(int, int) const;
  const Vector& GetColor(int) const;

  const Vector2& GetUV(int, int) const;
  const Vector2& GetUV(int) const;

  int size_vertex() const { return vertex.size(); }
  int size_normal() const { return normal.size(); }
  int size_UV() const { return UVmap.size(); }
  int size_triangles() const { return indexes.size() / 3; }

  void Clear();

  // Generate Mapping
  void generatePlannarZ_Mapping(Box2);
  void generatePlannarMapping(const Vector& dir, const double& step);
private:
  void ConvertToVertexEdgeFaceMaya(QVector<int>&, QVector<int>&, QVector<int>&, QVector<Vector>&) const;

protected:
  void AddTriangle(int, int, int, int);
  void AddSmoothTriangle(int, int, int, int, int, int);
  void AddQuadrangle(int, int, int, int, int);
};

/*!
\brief Get the k-th vertex of the i-th triangle.
\param i Triangle number.
\param k Vertex.
*/
inline const Vector& MayaGeometry::GetVertex(int i, int k) const
{
  return vertex.at(indexes.at(i * 3 + k).getIV());
}

/*!
\brief Get the normal of the k-th vertexof the i-th triangle.
\param i Triangle number.
\param k Vertex.
*/
inline const Vector& MayaGeometry::GetNormal(int i, int k) const
{
  return normal.at(indexes.at(i * 3 + k).getIN());
}

/*!
\brief Get the color of the k-th vertex of the i-th triangle.
\param i Triangle number.
\param k Vertex.
*/
inline const Vector& MayaGeometry::GetColor(int i, int k) const
{
  return color.at(indexes.at(i * 3 + k).getIC());
}

/*!
\brief Get the i-th vertec in the array.
\param i Index.
*/
inline const Vector& MayaGeometry::GetVertex(int i) const
{
  return vertex.at(i);
}

/*!
\brief Get the i-th normal in the array.
\param i Index.
*/
inline const Vector& MayaGeometry::GetNormal(int i) const
{
  return normal.at(i);
}

/*!
\brief Get the i-th color in the array.
\param i Index.
*/
inline const Vector& MayaGeometry::GetColor(int i) const
{
  return color.at(i);
}


/*!
\brief Get the i-th UV in the array.
\param i Index.
*/
inline const Vector2& MayaGeometry::GetUV(int i) const
{
  return UVmap.at(i);
}


/*!
\brief Get the i-th UV in the array.
\param i,k Indexes.
*/
inline const Vector2& MayaGeometry::GetUV(int i, int k) const
{
  return UVmap.at(indexes.at(i * 3 + k).getIUV());
}

/*!
\brief Set the name.
\param n Name.
*/
inline void MayaGeometry::SetName(const QString& n)
{
  name = n;
}

/*!
\brief Set the MayaMaterial.
\param mo MayaMaterial.
*/
inline void MayaGeometry::SetMaterial(const MayaMaterial& mo)
{
  mat = mo;
}
/*!
\brief Get the name of the geometry.
*/
inline QString MayaGeometry::GetName() const
{
  return name;
}

/*!
\brief Get the material of the geometry.
*/
inline MayaMaterial MayaGeometry::GetMaterial() const
{
  return mat;
}

// Set of Maya Geometry Instances
class MayaGeometrySet : public MayaGeometry
{
protected:
  QVector<FrameScaled> frames; //!< Array of frames
public:
  //! Empty
  MayaGeometrySet() {}
  MayaGeometrySet(const QString&, const MayaMaterial & = MayaMaterial::None);
  MayaGeometrySet(const QString&, const FrameScaled&, const MayaMaterial & = MayaMaterial::None);

  MayaGeometrySet(const MayaGeometry&);
  MayaGeometrySet(const MayaGeometry&, const FrameScaled&);
  MayaGeometrySet(const MayaGeometry&, const QVector<FrameScaled>&);
  MayaGeometrySet(const MayaGeometry&, const QVector<Vector>&);

  //! Empty
  ~MayaGeometrySet() {}

  // Setters
  void Append(const FrameScaled&);
  void SetName(const QString&);
  void AppendName(const QString& n);
  void SetFrames(QVector<FrameScaled>);
  void ApplyFrames(QVector<FrameScaled>);

  void Rotate(const Vector&);
  void Translate(const Vector&);
  void Scale(const Vector&);

  void Clear();
  void Remove(int);

  // Getters
  QString GetName() const;
  Box GetBox() const;
  FrameScaled GetFrameScaled(int) const;
  QVector<FrameScaled> GetFrames() const;
  int count() const;

  MayaGeometry Collapse() const;
  bool Save_VueScript(QTextStream& out) const;
  bool Save_VUE(QTextStream& out) const;
  QJsonObject Save_JSON(const QString& url) const;
  bool Save_JSON_SubObject(QTextStream& out, const QString& name, const QString& url, const QString& urlRelTexture = "", const bool collapse = true) const;
  void Save_JSON_ObjectInstances(QTextStream& out) const;


  MayaStatistics GetStatistics() const;
  QString GetText(int = 0, bool = false) const;

  static MayaGeometrySet CreateVoxel(const Voxel&, bool = false);
};

/*!
\brief Create a geometry set.
\param name Name.
\param mo Material.
*/
inline MayaGeometrySet::MayaGeometrySet(const QString& name, const MayaMaterial& mo) :MayaGeometry(name, mo)
{
}

/*!
\brief Get the i-th frame.
\param i Index.
*/
inline FrameScaled MayaGeometrySet::GetFrameScaled(int i) const
{
  return frames[i];
}

/*!
\brief Get the set of frames.
*/
inline QVector<FrameScaled> MayaGeometrySet::GetFrames() const
{
  return frames;
}

/*!
\brief Create a geometry set and one instance.
\param name Name.
\param frame Frame of the instance.
\param mo Material.
*/
inline MayaGeometrySet::MayaGeometrySet(const QString& name, const FrameScaled& frame, const MayaMaterial& mo) :MayaGeometry(name, mo)
{
  frames.append(frame);
}


/*!
\brief Get the number of frames in the instance set.
*/
inline int MayaGeometrySet::count() const
{
  return frames.count();
}
/*!
\brief Remove a given instance from the set.
\param i Instance index.
*/
inline void MayaGeometrySet::Remove(int i)
{
  frames.remove(i);
}

/*!
\brief Set the name.
\param n Name.
*/
inline void MayaGeometrySet::SetName(const QString& n)
{
  name = n;
}

inline void MayaGeometrySet::AppendName(const QString& n)
{
  name += n;
}

/*!
\brief Get the name of the geometry set.
*/
inline QString MayaGeometrySet::GetName() const
{
  return name;
}

/*!
\brief Set the set of frames that define the locations of the instances.
\param frames Set of frames.
*/
inline void MayaGeometrySet::SetFrames(QVector<FrameScaled> frames)
{
  MayaGeometrySet::frames = frames;
}

class MayaResources;
class MayaInstanceSet;
class MayaInstanceAll;

class Camera;

// Set of geometry instances along with their frames
class MayaGeometryAll {
protected:
  QMap<QString, MayaGeometrySet> instances; //!< Map between the name of the geometric instances and the MayaGeometrySet structure. A map is used instead of an array to speed up instance search given an input name.
public:
  //! Empty 
  MayaGeometryAll() {}
  MayaGeometryAll(const MayaGeometry&);
  MayaGeometryAll(const MayaGeometrySet&);
  MayaGeometryAll(const QVector<MayaGeometrySet>&);
  MayaGeometryAll(const MayaResources&, const MayaInstanceSet&);
  MayaGeometryAll(const MayaResources&, const MayaInstanceAll&);

  //! Empty
  ~MayaGeometryAll() {}

  // Setters
  void Append(const MayaGeometry&);
  void Append(const MayaGeometrySet&);
  void Append(const QString&, const FrameScaled&);

  void SetFramesToAll(QVector<FrameScaled>);
  void SetMaterial(const MayaMaterial&);
  void SetMaterial_EqualName(const QString&, const MayaMaterial&);
  void SetMaterial_ContainsName(const QString& n, const MayaMaterial& mo);
  void AppendName(const QString&);

  MayaGeometry GetGeometryByName(const QString& n);


  void Rotate(const Vector&);
  void Translate(const Vector&);
  void Scale(const Vector&);

  void Remove(const QString&);

  void Clear();

  void Append(const MayaGeometryAll&);

  // Getters
  int count() const;
  QList<QString> GetNames() const;
  MayaGeometrySet& operator[] (const QString&);
  int GetSetSize() const;
  Box GetBox() const;

  MayaGeometry Collapse() const;

  bool GetMayaGeometrySet(const QString&, MayaGeometrySet&);
  bool GetMayaGeometrySet(const QString&, QVector<MayaGeometrySet>&);

  // Load and save functions
  bool Load_OBJ(const QString&);
  //bool Load_MAYA_Frames(const QString&);
  bool find_MTL_Address(const QString&, const QString&, QString&, QString&);

  bool Save_OBJ(const QString&, const bool& = true) const;
  bool Save_OBJ_Collapse(const QString&) const;

  bool Save_PBRT(const QString&) const;
  bool Save_VUE(const QString&) const;

  bool Save_XML_Instances(const QString&) const;
  bool Save_XML_Object(const QString&, const QString& id) const;
  bool Save_XML_Object(const QString&, const QString& id, const QString& objbasename, const QString& texturepath, const QString& relPath) const;
  bool Save_XML(const QString&, const Camera&, const QString& instancefilename = QString(""), bool = false, const QString & = QString(""), bool objectonly = false) const;

  bool Save_JSON_Scene(const QString&);
  bool Save_JSON_Object(const QString&);
  bool Save_JSON_ObjectsInstances(const QString&);

  bool exportSceneJson(QTextStream& out, const QString& name, const QString& url) const;


  bool Save_VueScript(const QString&) const;

  MayaStatistics GetStatistics() const;
  QString GetText(int = 0, bool = false) const;

  static MayaGeometryAll CreateVoxelSurface(const Voxel&);

  // Friends
  friend class MayaGpuAll;
};

/*!
\brief Get the number of instances.
*/
inline int MayaGeometryAll::count() const
{
  return instances.count();
}

/*!
\brief Get the list of names of the instances in the geometry.
*/
inline QList<QString> MayaGeometryAll::GetNames() const
{
  return instances.keys();
}
/*!
\brief Get the MayaGeometrySet of the geometry.
\param n Name.
*/
inline MayaGeometrySet& MayaGeometryAll::operator[](const QString& n)
{
  return instances[n];
}

/*!
\brief Get the size of instances.
*/
inline int MayaGeometryAll::GetSetSize() const
{
  return instances.size();
}

// Set of instances along with their frames
class MayaResources {
protected:
  QMap<QString, MayaGeometryAll> atlas; //!< Map of ressources identified by their names
public:
  //! Empty 
  MayaResources() {}
  //! Empty 
  ~MayaResources() {}

  // Setters
  void Append(const QString&, const QString&);
  void Append(const QString&, const MayaGeometry&);
  void Append(const QString&, const MayaGeometrySet&);
  void Append(const QString&, const MayaGeometryAll&);

  void Remove(const QString&);

  // Getters
  bool Exist(const QString&) const;
  bool GetRessource(const QString&, MayaGeometryAll&) const;

  void exportSceneJson(QTextStream& out, const QString& url) const;

  void Clear();
  QString GetText(bool = false) const;
};


