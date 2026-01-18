#pragma once

#include "libs/gpu-shader.h"
#include "libs/evector.h"

class GLBuffer
{
protected:
  GLuint buffer;
public:
  inline GLBuffer() : buffer(0) { }
  inline GLBuffer(GLuint b) : buffer(b) { }
  inline ~GLBuffer() { }

  void Generate();
  void Destroy();

  void SetData(int bufferType, size_t size, void* data, int hint = GL_STATIC_DRAW);
  void SetSubData(int bufferType, size_t offset, size_t size, void* data);
  void Bind(int bufferType) const;
  void BindAt(int bufferType, int bindIndex) const;
  void Unbind(int bufferType) const;

  GLuint GetBuffer() const;
  template<typename T> void GetData(std::vector<T>& data) const
  {
    glGetNamedBufferSubData(buffer, 0, sizeof(T) * data.size(), data.data());
  }
};
inline void GLBuffer::Generate() {
  if (buffer == 0)
    glGenBuffers(1, &buffer);
}

inline void GLBuffer::Destroy() {
  glDeleteBuffers(1, &buffer);
}

inline void GLBuffer::SetData(int bufferType, size_t size, void* data, int hint) {
  Bind(bufferType);
  glBufferData(bufferType, size, data, hint);
}

inline void GLBuffer::Unbind(int bufferType) const
{
  glBindBuffer(bufferType, 0);
}

inline void GLBuffer::SetSubData(int bufferType, size_t offset, size_t size, void* data) {
  Bind(bufferType);
  glBufferSubData(bufferType, offset, size, data);
}

inline void GLBuffer::Bind(int bufferType) const {
  glBindBuffer(bufferType, buffer);
}

inline void GLBuffer::BindAt(int bufferType, int bindIndex) const
{
  glBindBufferBase(bufferType, bindIndex, buffer);
}

inline GLuint GLBuffer::GetBuffer() const {
  return buffer;
}


class GLShader
{
protected:
  GLuint program = 0;
public:
  inline GLShader() { }
  inline ~GLShader() { }

  void Initialize(const char*);
  void Initialize(const char*, const char*);
  void Destroy();
  void Bind() const;
  void Unbind() const;

  void SetUniform(const char*, float) const;
  void SetUniform(const char*, float, float) const;
  void SetUniform(const char*, float, float, float) const;
  void SetUniform(const char*, int) const;
  void SetUniform(const char*, const Vector&) const;
  void SetUniform(const char*, const float* mat4) const;
  void SetUniform(const char* name, double v0, double v1) const;
  void SetUniform(const char* name, int i0, int i1) const;

  void SetUniform(int index, float) const;
  void SetUniform(int index, float, float) const;
  void SetUniform(int index, float, float, float) const;
  void SetUniform(int index, int) const;
  void SetUniform(int index, const Vector&) const;
  void SetUniform(int index, const float* mat4) const;
  void SetUniform(int index, double v0, double v1) const;
};

inline void GLShader::Initialize(const char* path) {
  program = read_program(path);
}

inline void GLShader::Initialize(const char* path, const char* definitions) {
  program = read_program(path, definitions);
}

inline void GLShader::Destroy() {
  release_program(program);
}

inline void GLShader::Bind() const {
  glUseProgram(program);
}

inline void GLShader::Unbind() const {
  glUseProgram(0);
}

inline void GLShader::SetUniform(const char* name, float v) const {
  glUniform1f(glGetUniformLocation(program, name), v);
}

inline void GLShader::SetUniform(const char* name, float v0, float v1) const {
  glUniform2f(glGetUniformLocation(program, name), v0, v1);
}

inline void GLShader::SetUniform(const char* name, double v0, double v1) const {
  glUniform2d(glGetUniformLocation(program, name), v0, v1);
}

inline void GLShader::SetUniform(const char* name, int i0, int i1) const {
  glUniform2i(glGetUniformLocation(program, name), i0, i1);
}

inline void GLShader::SetUniform(const char* name, float v0, float v1, float v2) const {
  glUniform3f(glGetUniformLocation(program, name), v0, v1, v2);
}

inline void GLShader::SetUniform(const char* name, int v) const {
  glUniform1i(glGetUniformLocation(program, name), v);
}

inline void GLShader::SetUniform(const char* name, const Vector& v) const {
  glUniform3f(glGetUniformLocation(program, name), float(v[0]), float(v[1]), float(v[2]));
}

inline void GLShader::SetUniform(const char* name, const float* mat4) const {
  glUniformMatrix4fv(glGetUniformLocation(program, name), 1, GL_FALSE, mat4);
}

inline void GLShader::SetUniform(int index, float v) const {
  glUniform1f(index, v);
}

inline void GLShader::SetUniform(int index, float v0, float v1) const {
  glUniform2f(index, v0, v1);
}

inline void GLShader::SetUniform(int index, float v0, float v1, float v2) const {
  glUniform3f(index, v0, v1, v2);
}

inline void GLShader::SetUniform(int index, int v) const {
  glUniform1i(index, v);
}

inline void GLShader::SetUniform(int index, const Vector& v) const {
  glUniform3f(index, float(v[0]), float(v[1]), float(v[2]));
}

inline void GLShader::SetUniform(int index, const float* mat4) const {
  glUniformMatrix4fv(index, 1, GL_FALSE, mat4);
}

inline void GLShader::SetUniform(int index, double v0, double v1) const {
  glUniform2d(index, v0, v1);
}
