#ifndef _GLCanvas_GLCanvas_h_
#define _GLCanvas_GLCanvas_h_

#include <GLCtrl/GLCtrl.h>
#include <Surface/Surface.h>
#include "trackball.h"

class GLCanvas : public GLCtrl {
public:
	typedef GLCanvas CLASSNAME;

	GLCanvas();

private:
	TrackBall trackBall;

	void SetUpLighting();
	void SetCamera();
	
protected:
	virtual void Layout();
	virtual Image MouseEvent(int event, Point p, int zdelta, dword keyflags);
	
public:
	void PaintLine(double x0, double y0, double z0, double x1, double y1, double z1, const Color &color);
	void PaintLine(const Point3D &p0, const Point3D &p1, const Color &color);
	void PaintLine(const Segment3D &p, const Color &color);
	void PaintQuad(Point3D &p0, Point3D &p1, Point3D &p2, Point3D &p3, const Color &color, double multx, double multy);
	void PaintAxis(double x, double y, double z);	
	void PaintSurface(Surface &surf, const Color &linCol = LtGreen());

	void SetZoomFactor(double factor)	{trackBall.SetZoomFactor(factor);}
		
	Function <void()> WhenPaint;	
	
	void OnPaint();	
	
	virtual void GLResize(int w, int h) {
		glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	}
	
private:
	void PaintSurface0(Surface &surf, const Color &linCol, bool simX, bool simY);
};

#endif
