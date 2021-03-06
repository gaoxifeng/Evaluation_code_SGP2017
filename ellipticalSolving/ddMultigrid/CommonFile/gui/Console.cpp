#include "Console.h"
#include <GL/freeglut.h>
#include <time.h>

USE_PRJ_NAMESPACE

DefaultConsole::DefaultConsole()
{
  setColor(1.0f,0.0f,0.0f,1.0f,60);
}
void DefaultConsole::setColor(float r,float g,float b,float a,int sz)
{
  _r=r;
  _g=g;
  _b=b;
  _a=a;
  _sz=sz;
}
void DefaultConsole::addMsg(const std::string& msg,sizeType dur)
{
  if(dur > 0)
    _msgs.push(make_pair(msg,(sizeType)clock()+dur));
  else if(dur < 0)
    _cmsgs[dur]=msg;
}
void DefaultConsole::rmvMsg(sizeType dur)
{
  if(dur < 0)
    _cmsgs.erase(dur);
}
void DefaultConsole::reshape(int w,int h)
{
  _w=w;
  _h=h;
}
void DefaultConsole::render()
{
  initDrawMsg();
  sizeType curr=(sizeType)clock();
  for(std::map<sizeType,std::string>::const_iterator
      beg=_cmsgs.begin(),end=_cmsgs.end(); beg!=end; beg++)
    drawMsg(beg->second);

  std::queue<std::pair<std::string,sizeType> > bk;
  while(!_msgs.empty()) {
    std::pair<std::string,sizeType> top=_msgs.front();
    _msgs.pop();
    if(top.second > curr) {
      drawMsg(top.first);
      bk.push(top);
    }
  }
  std::swap(bk,_msgs);
  finishDrawMsg();
}
void DefaultConsole::initDrawMsg()
{
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0.0f,1.0f,0.0f,1.0f);

  glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
  glColor4f(_r,_g,_b,_a);
  glDisable(GL_LIGHTING);

  GLfloat delta=0.01f;
  _posx=delta;
  _posy=1.0f-(GLfloat)_sz/(GLfloat)_h-delta;
}
void DefaultConsole::finishDrawMsg()
{
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}
void DefaultConsole::drawMsg(const std::string& str)
{
  glRasterPos2f(_posx,_posy);
  glRasterPos2f(_posx,_posy);
  for(sizeType i=0; i<(sizeType)str.size(); i++)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,str[i]);
  _posy-=(GLfloat)_sz/(GLfloat)_h;
}
