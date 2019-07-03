#include <GLCanvas/GLCanvas.h>

using namespace Upp;


GLCanvas::GLCanvas() {
	WhenGLPaint = THISBACK(OnPaint);
	
	WantFocus();
	
	trackBall.Init(this);
	
	SetCamera();
}

Image GLCanvas::MouseEvent(int event, Point p, int zdelta, dword keyflags) {
	Image img = trackBall.MouseEvent(event, p, zdelta, keyflags);
	Refresh();
	return img;
}

void GLCanvas::SetUpLighting() {
	float light1_ambient[4]  = { 1.0f, 1.0f, 1.0f, 1.0f };
	float light1_diffuse[4]  = { 1.0f, 0.9f, 0.9f, 1.0f };
	float light1_specular[4] = { 1.0f, 0.7f, 0.7f, 1.0f };
	float light1_position[4] = { -1.0, 1.0, 1.0, 0.0f };
	glLightfv(GL_LIGHT1, GL_AMBIENT,  light1_ambient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE,  light1_diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular);
	glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
	glEnable(GL_LIGHT1);
	
	float light2_ambient[4]  = { 0.2f, 0.2f, 0.2f, 1.0f };
	float light2_diffuse[4]  = { 0.9f, 0.9f, 0.9f, 1.0f };
	float light2_specular[4] = { 0.7f, 0.7f, 0.7f, 1.0f };
	float light2_position[4] = { 1.0, -1.0, -1.0, 0.0f };
	glLightfv(GL_LIGHT2, GL_AMBIENT,  light2_ambient);
	glLightfv(GL_LIGHT2, GL_DIFFUSE,  light2_diffuse);
	glLightfv(GL_LIGHT2, GL_SPECULAR, light2_specular);
	glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
	//glEnable(GL_LIGHT2);
	
	float front_emission[4] = { 0.3f, 0.2f, 0.1f, 0.0f };
	float front_ambient[4]  = { 0.2f, 0.2f, 0.2f, 0.0f };
	float front_diffuse[4]  = { 0.95f, 0.95f, 0.8f, 0.0f };
	float front_specular[4] = { 0.6f, 0.6f, 0.6f, 0.0f };
	glMaterialfv(GL_FRONT, GL_EMISSION, front_emission);
	glMaterialfv(GL_FRONT, GL_AMBIENT, front_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, front_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, front_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, 16.0);
	glColor4fv(front_diffuse);
	
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
	glEnable(GL_CULL_FACE);
	glColorMaterial(GL_FRONT, GL_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	
}

void GLCanvas::SetCamera(void) {
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	gluPerspective(trackBall.GetZoomFactor(), (float)GetSize().cx / (float)GetSize().cy, 1, 10000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0, 0.0, (float)GetSize().cy/2., 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
}

void GLCanvas::Layout() {
	GLCtrl::Layout();
	glMatrixMode (GL_MODELVIEW);
	glViewport (0, 0, GetSize().cx, GetSize().cy);
	glLoadIdentity();
	SetCamera();
	trackBall.Reshape(GetSize().cx, GetSize().cy);
}

void GLCanvas::OnPaint() {	
	MemoryIgnoreLeaksBlock __; 
	
	glClearColor(1, 1, 1, 0);
	glClearDepth(0); 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	SetCamera();
	glTranslatef(0, 0, 0);
	trackBall.Matrix();
	
	glPushMatrix();
	glDisable(GL_BLEND);
	SetUpLighting();
	WhenPaint();
	glPopMatrix();
}

void GLCanvas::PaintLine(double x0, double y0, double z0, double x1, double y1, double z1, const Color &color) {
	glBegin(GL_LINES);
		glColor4d(color.GetR()/255., color.GetG()/255., color.GetB()/255., 1);
		glVertex3d(x0, y0, z0);
		glVertex3d(x1, y1, z1); // 
	glEnd();
}

void GLCanvas::PaintLine(const Point3D &p0, const Point3D &p1, const Color &color) {
	glBegin(GL_LINES);
		glColor4d(color.GetR()/255., color.GetG()/255., color.GetB()/255., 1);
		glVertex3d(p0.x, p0.y, p0.z);
		glVertex3d(p1.x, p1.y, p1.z); // 
	glEnd();
}

void GLCanvas::PaintLine(const Line3D &p, const Color &color) {
	glBegin(GL_LINES);
		glColor4d(color.GetR()/255., color.GetG()/255., color.GetB()/255., 1);
		glVertex3d(p.from.x, p.from.y, p.from.z);
		glVertex3d(p.to.x, p.to.y, p.to.z); // 
	glEnd();
}

void GLCanvas::PaintQuad(Point3D &p0, Point3D &p1, Point3D &p2, Point3D &p3, const Color &color, double multx, double multy) {
	glBegin(GL_QUADS);
		glColor4d(color.GetR()/255., color.GetG()/255., color.GetB()/255., 1);
		glVertex3d(p0.x*multx, p0.y*multy, p0.z);
		glVertex3d(p1.x*multx, p1.y*multy, p1.z);
		glVertex3d(p2.x*multx, p2.y*multy, p2.z);
		glVertex3d(p3.x*multx, p3.y*multy, p3.z);
	glEnd();
}

void GLCanvas::PaintAxis(double x, double y, double z) {
	PaintLine(0, 0, 0, x, 0, 0, LtRed());
	PaintLine(0, 0, 0, 0, y, 0, LtGreen());
	PaintLine(0, 0, 0, 0, 0, z, LtBlue());
}

void GLCanvas::PaintSurface(Surface &surf, const Color &linCol) {
	PaintSurface0(surf, linCol, false, false);
	if (surf.x0z)
		PaintSurface0(surf, linCol, false, true);
	if (surf.y0z)
		PaintSurface0(surf, linCol, true, false);
}

template<class T>
inline T Avg(T a, T b) 			{return T(a+b)/2;}
template<class T>
inline T Avg(T a, T b, T c)		{return T(a+b+c)/3;}

void GLCanvas::PaintSurface0(Surface &surf, const Color &linCol, bool simX, bool simY) {
	double xsig = simX ? -1 : 1;
	double ysig = simY ? -1 : 1;
	
	for (int ip = 0; ip < surf.panels.GetCount(); ++ip) {
		Panel &panel = surf.panels[ip];
		Point3D p0 = surf.nodes[panel.id[0]];
		Point3D p1 = surf.nodes[panel.id[1]];
		Point3D p2 = surf.nodes[panel.id[2]];
		Point3D p3 = surf.nodes[panel.id[3]];
		p0.x *= xsig;
		p0.y *= ysig;
		p1.x *= xsig;
		p1.y *= ysig;
		p2.x *= xsig;
		p2.y *= ysig;
		p3.x *= xsig;
		p3.y *= ysig;
	
		//PaintQuad(p0, p1, p2, p3, linCol, xsig, ysig);
		
		/*glBegin(GL_QUADS);
			glColor4d(1.0, 0.5, 0.0, 1);	
			if (!sim) {
				glVertex3d(p0.x, p0.y*ysig, p0.z);
				glVertex3d(p1.x, p1.y*ysig, p1.z);
				glVertex3d(p2.x, p2.y*ysig, p2.z);
				glVertex3d(p3.x, p3.y*ysig, p3.z);
			} else {
				glVertex3d(p3.x, p3.y*ysig, p3.z);
				glVertex3d(p2.x, p2.y*ysig, p2.z);
				glVertex3d(p1.x, p1.y*ysig, p1.z);
				glVertex3d(p0.x, p0.y*ysig, p0.z);
			}
		glEnd();*/	
		
		PaintLine(p0, p1, linCol);
		PaintLine(p1, p2, linCol);
		PaintLine(p2, p3, linCol);
		PaintLine(p3, p0, linCol);
		
		Point3D from = GetCentroid(p0, p2);		
		Point3D pnormal = GetNormal(p0, p1, p2);
		Line3D normal(from, pnormal, 1);
		PaintLine(normal, Blue());
	}
}

Point3D GetCentroid(Point3D &a, Point3D &b) {
	return Point3D(Avg(a.x, b.x), Avg(a.y, b.y), Avg(a.z, b.z));	
}

Point3D GetCentroid(Point3D &a, Point3D &b, Point3D &c) {
	return Point3D(Avg(a.x, b.x,c.x), Avg(a.y, b.y, c.y), Avg(a.z, b.z, c.z));	
}

Vector3D GetNormal(Point3D &a, Point3D &b, Point3D &c) {
	return Vector3D((a - b) % (b - c)).Normalize();
}
