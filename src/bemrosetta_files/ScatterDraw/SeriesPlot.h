#ifndef _ScatterDraw_SeriesPlot_h_
#define _ScatterDraw_SeriesPlot_h_

class SeriesPlot {
public:
	virtual ~SeriesPlot() 	{};	
	virtual void Paint(Draw& w, Vector<Pointf> &p, const double &scale, const double opacity, 
				double thick, const Color &color, String pattern, const Color &background, 
				const Color &fillColor, double fx, double fy, double y0, double width, 
				bool isClosed) const = 0;
	virtual void Paint(Painter& w, Vector<Pointf> &p, const double &scale, const double opacity, 
				double thick, const Color &color, String pattern, const Color &background, 
				const Color &fillColor, double fx, double fy, double y0, double width, 
				bool isClosed) const = 0;
	template<class T>
	static void Register(const String& name)
	{
		classMap().FindAdd(name, __Create<T>);
		typeMap().FindAdd(typeid(T).name(), name);
	}
	static void Unregister(const String& name)
	{
		int i = TypeIndex(name);
		ASSERT(i >= 0);
		classMap().Remove(i);
		typeMap().Remove(i);
	}
	static String         TypeName(int i)               {return typeMap()[i];}
	static int            TypeIndex(const String& name) {return typeMap().Find(name);}
	static int            GetCount()                    {return classMap().GetCount();}
	static SeriesPlot*	  Create(int i)                 {return classMap()[i]();}
	static Vector<String> GetTypes()                    {return clone(typeMap()).PickValues();}
	int GetType()const {return typeMap().Find(typeid(*this).name());}
	SeriesPlot* Copy()const;
	       	
protected:
	typedef SeriesPlot* (*CreateFunc)();
	template<class T>	
	static SeriesPlot*                    __Create() {return new T;}
	static VectorMap<String, CreateFunc>& classMap() {static VectorMap<String, CreateFunc> cMap; return cMap;}
	static VectorMap<String, String>&     typeMap()  {static VectorMap<String, String> tMap; 	 return tMap;}
};

class LineSeriesPlot : public SeriesPlot {
private:
	template <class T>
	void DoPaint(T& w, Vector<Pointf> &p, const double &scale, const double opacity, 
				double thick, const Color &color, String pattern, const Color &background, 
				const Color &fillColor, double y0, bool isClosed) const 
	{
		if (!IsNull(fillColor)) {
			if(isClosed) 
				FillPolylineOpa(w, p, scale, opacity, background, fillColor);
			else {
				double x = p[0].x;
				p.Insert(0, Pointf(x, y0));
				x = p[p.GetCount() - 1].x;
				p.Add(Pointf(x, y0));
				FillPolylineOpa(w, p, scale, opacity, background, fillColor);
				p.Remove(0);
				p.Remove(p.GetCount() - 1);
			}
		}
		DrawPolylineOpa(w, p, scale, 1, thick, color, pattern, background);		
	}
	
public:
	void Paint(Draw& w, Vector<Pointf> &p, const double &scale, const double opacity, 
				double thick, const Color &color, String pattern, const Color &background, 
				const Color &fillColor, double , double , double y0, double , bool isClosed) const 
	{
		DoPaint(w, p, scale, opacity, thick, color, pattern, background, fillColor, y0, isClosed);		
	}
	void Paint(Painter& w, Vector<Pointf> &p, const double &scale, const double opacity, 
				double thick, const Color &color, String pattern, const Color &background, 
				const Color &fillColor, double , double , double y0, double , bool isClosed) const 
	{
		DoPaint(w, p, scale, opacity, thick, color, pattern, background, fillColor, y0, isClosed);		
	}
};


class StaggeredSeriesPlot : public SeriesPlot {
private:
	template <class T>
	void DoPaint(T& w, Vector<Pointf> &p, const double &scale, const double opacity, 
				double thick, const Color &color, String pattern, const Color &background,
				const Color &fillColor, double y0) const 
	{
		ASSERT(p.GetCount() > 1);
		Vector<Pointf> ps;
		ps << Pointf(p[0].x - (p[1].x - p[0].x)/2., p[0].y);
		for (int i = 1; i < p.GetCount(); ++i) {
			double x = (p[i].x + p[i-1].x)/2.;
			ps << Pointf(x, p[i-1].y);
			ps << Pointf(x, p[i].y);
		}
		ps << Pointf(p[p.GetCount() - 1].x + (p[p.GetCount() - 1].x - p[p.GetCount() - 2].x)/2., p[p.GetCount() - 1].y);
		if (!IsNull(fillColor)) {
			double x = ps[0].x;
			ps.Insert(0, Pointf(x, y0));
			x = ps[ps.GetCount() - 1].x;
			ps.Add(Pointf(x, y0));
			FillPolylineOpa(w, ps, scale, opacity, background, fillColor);
			ps.Remove(0);
			ps.Remove(ps.GetCount() - 1);
		}
		DrawPolylineOpa(w, ps, scale, 1, fround(thick), color, pattern, background);
	}
	
public:
	void Paint(Draw& w, Vector<Pointf> &p, const double &scale, const double opacity, 
				double thick, const Color &color, String pattern, const Color &background, 
				const Color &fillColor, double , double , double y0, double , bool ) const 
	{
		DoPaint(w, p, scale, opacity, thick, color, pattern, background, fillColor, y0);
	}
	void Paint(Painter& w, Vector<Pointf> &p, const double &scale, const double opacity, 
				double thick, const Color &color, String pattern, const Color &background, 
				const Color &fillColor, double , double , double y0, double , bool ) const 
	{
		DoPaint(w, p, scale, opacity, thick, color, pattern, background, fillColor, y0);
	}
};

class BarSeriesPlot : public SeriesPlot {
private:
	template <class T>
	void DoPaint(T& w, Vector<Pointf> &p, const double &scale, const double opacity, 
				double thick, const Color &color, String pattern, const Color &background, 
				const Color &fillColor, double fx, double y0, double width) const 
	{
		for (int i = 0; i < p.GetCount(); ++i) {
			FillRectangleOpa(w, p[i].x - width*fx, y0, p[i].x + width*fx, p[i].y, opacity, background, fillColor);
			Vector<Pointf> ps;
			ps << Pointf(fround(p[i].x - width*fx), y0) << Pointf(fround(p[i].x - width*fx), p[i].y) 
			   << Pointf(fround(p[i].x + width*fx), p[i].y) << Pointf(fround(p[i].x + width*fx), y0);
			DrawPolylineOpa(w, ps, scale, 1, fround(thick), color, pattern, background);
		}
	}
	
public:
	void Paint(Draw& w, Vector<Pointf> &p, const double &scale, const double opacity, 
				double thick, const Color &color, String pattern, const Color &background, 
				const Color &fillColor, double fx, double , double y0, double width, bool ) const 
	{
		DoPaint(w, p, scale, opacity, thick, color, pattern, background, fillColor, fx, y0, width);
	}
	void Paint(Painter& w, Vector<Pointf> &p, const double &scale, const double opacity, 
				double thick, const Color &color, String pattern, const Color &background, 
				const Color &fillColor, double fx, double , double y0, double width, bool ) const 
	{
		DoPaint(w, p, scale, opacity, thick, color, pattern, background, fillColor, fx, y0, width);
	}
};	


#endif
