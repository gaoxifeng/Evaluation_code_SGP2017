#ifndef RENDERER_H
#define RENDERER_H

#include "../Config.h"
#include "glsl.h"
#include <vector>
#include <boost/shared_ptr.hpp>

PRJ_BEGIN

class Renderable
{
public:
  virtual void drawObject();
};
class DefaultRenderer
{
public:
  virtual ~DefaultRenderer();
  virtual void init();
  virtual void render(Renderable& rend);
  virtual void reshape(int w,int h);
  virtual void disableShadowReceive() {}
  virtual void enableShadowReceive() {}
  virtual void disableShadowCast() {}
  virtual void enableShadowCast() {}
  void drawLights();
  void setShowLights(bool show,GLfloat rad);
  cwc::glShaderManager& getManager();
  //utility
  static GLuint getMaxNumLight();
  static void checkGLError();
  static void drawQuad();
  static void checkFBO();
  static void genTexture(GLuint& tid,GLsizei w,GLsizei h,GLint internalFormat,GLenum format,GLenum type,GLint filter=GL_NEAREST,GLint wrap=GL_CLAMP);
  static void writeTexture(GLuint tid,GLenum format,std::string path);
  static void writeColorTexture(GLuint tid,std::string path);
protected:
  //data
  static const int MAX_LIGHTS=8;
  GLuint _w,_h;
  GLfloat _rad;
  bool _showLight;
  boost::shared_ptr<cwc::glShaderManager> _SM;
};
class PerPixelRenderer : public DefaultRenderer
{
public:
  virtual ~PerPixelRenderer();
  virtual void init();
  virtual void render(Renderable& rend);
  virtual void reshape(int w,int h);
  virtual void disableShadowReceive();
  virtual void enableShadowReceive();
  virtual void disableShadowCast();
  virtual void enableShadowCast();
protected:
  virtual void release();
  void generateShadowFBO(int w,int h,GLuint& tid,GLuint& fboId) const;
  void setupMatrices(int L) const;
  void setTextureMatrix(int L) const;
  GLuint shadowSz() const {
    return (_w+_h)*_shadowRatio/2;
  }
  virtual void beforeShadowDraw(int L);
  virtual void afterShadowDraw(int L);
  virtual void beforeShadowRender(cwc::glShader* shader);
  //data
  cwc::glShader *_shaderDepth;
  cwc::glShader *_shaderDraw;
  GLuint _tid[MAX_LIGHTS],_ctid[MAX_LIGHTS],_fboId[MAX_LIGHTS],_shadowRatio;
  bool _released;
};
class VSMPerPixelRenderer : public PerPixelRenderer
{
public:
  virtual void init();
  virtual void reshape(int w,int h);
  virtual void enableShadowReceive();
protected:
  virtual void release();
  virtual void afterShadowDraw(int L);
  virtual void beforeShadowRender(cwc::glShader* shader);
  void blurPass(GLfloat dx,GLfloat dy,GLuint fbo,GLuint color) const;
  void generateBlurKernel(GLint nr);
  cwc::glShader *_shaderBlur;
  std::vector<GLfloat> _off,_coef;
  GLuint _ctTmp,_fboTmp;
};

PRJ_END

#endif
