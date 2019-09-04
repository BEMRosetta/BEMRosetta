#include "PieDraw.h"
#include "DrawingFunctions.h"

static void DrawPie(Draw& w, double c_x, double c_y, double r, int start, int alpha, int width = 0, 
					Color fill = Null, Color outline = Black, uint64 pattern = 0, Color background = White) {
	const int dalpha = 1;
	int n = alpha/dalpha;
	
	Vector <Point> vP;
	Point centre = Point(int(c_x), int(c_y));
	vP << centre;
	int ix;
	int iy;
	for (int i = 0; i <= n; i++) {
		double x = c_x + r*cos((start+i*dalpha)*M_PI/1800);
		ix = fround(x);
		double y = c_y + r*sin((start+i*dalpha)*M_PI/1800);
		iy = fround(y);
		double dxy = (x-ix)*(x-ix) + (y-iy)*(y-iy);
		if(dxy < 0.1 || i == 0 || i == n)
			vP << Point(ix,iy);
		if(w.Pixels())
			w.DrawRect(ix, iy, 1, 1, Blend(fill, background, 150));
	}
	vP << centre;
	w.DrawPolygon(vP, fill, width, outline, pattern, Null);
}

void PieDraw::AddCategory(const String& name, const double& value, const Color& catcolor) {
	vNames.Add(name);
	vValues.Add(value);
	vColors.Add(catcolor);
	Refresh();
}

void PieDraw::RemoveCategory(const int& index) {
	ASSERT(!vValues.IsEmpty() && vValues.GetCount() > index);
	vNames.Remove(index);
	vValues.Remove(index);
	vColors.Remove(index);
	Refresh();
}

PieDraw& PieDraw::SetCatValue(const int& index, const double& value) {
	ASSERT(!vValues.IsEmpty() && vValues.GetCount() > index);
	vValues[index] = value;
	Refresh();
	return *this;
}

PieDraw& PieDraw::SetCatName(const int& index, const String& name) {
	ASSERT(!vNames.IsEmpty() && vNames.GetCount() > index);
	vNames[index] = name;
	Refresh();
	return *this;
}

PieDraw& PieDraw::SetCatColor(const int& index, const Color& catcolor) {
	ASSERT(!vColors.IsEmpty() && vColors.GetCount() > index);
	vColors[index] = catcolor;
	Refresh();
	return *this;
}

double PieDraw::GetCatValue(const int& index) const {
	ASSERT(!vValues.IsEmpty() && vValues.GetCount() > index);
	return vValues[index];
}

String PieDraw::GetCatName(const int& index) const {
	ASSERT(!vNames.IsEmpty() && vNames.GetCount() > index);
	return vNames[index];
}

Color PieDraw::GetCatColor(const int& index) const {
	ASSERT(!vColors.IsEmpty() && vColors.GetCount() > index);
	return vColors[index];
}

void PieDraw::PaintPie(Draw& w, int scale) {
	Size sz = scale*GetSize();	
	
	if (!IsNull(backColor))
		w.DrawRect(sz, backColor);	
	
	Size textsize;
	textsize.cx = 0;
	textsize.cy = 0;    
	if(!title.IsEmpty()) {
		Upp::Font FontTitle6;
		FontTitle6 = titleFont;
		FontTitle6.Height(scale*titleFont.GetHeight());
		FontTitle6.Width(scale*titleFont.GetWidth());
		textsize = GetTextSizeSpace(title, FontTitle6);
		if(titlePos == TOP) 
			w.DrawText((scale*GetSize().cx - textsize.cx)/2, scale*titleGap, title, FontTitle6, titleColor);
		else  
			w.DrawText((scale*GetSize().cx - textsize.cx)/2, scale*(GetSize().cy - titleFont.GetHeight() - titleGap), 
						title, FontTitle6,titleColor);
	}
	
	if(vValues.IsEmpty())
		return;
	
	int alfa0 = -900 + static_cast<int>(pieAngle);
	int a0 = 0;
	double sum = 0;
	for(int i = 0; i < vValues.GetCount(); i++)
		sum += vValues[i];
	
	double circWidth = sz.cx - pieMarginLeft - pieMarginRight;
	if (circWidth < 0)
		circWidth = 0;
	double circHeight = sz.cy - pieMarginTop - textsize.cy - pieMarginBottom;
	if (circHeight < 0)
		circHeight = 0;
	double circ_r;
	if (circWidth > circHeight)
		circ_r = circHeight/2.;
	else
		circ_r = circWidth/2.;
	double circ_x = pieMarginLeft + circWidth/2.;
	
	double circ_y;
	if(titlePos == TOP)
		circ_y = pieMarginTop + titleGap + textsize.cy + circHeight/2.;
	else
		circ_y = pieMarginTop + circHeight/2.;

	circ_x *= scale;
	circ_y *= scale;
	circ_r *= scale;
	
	for(int i = 0; i < vValues.GetCount(); i++) {
		DrawPie(w, circ_x*scale, circ_y*scale, circ_r*scale, alfa0, fround(3600.0*vValues.At(i)/sum),
					  1*scale, vColors.At(i), vColors.At(i), 0, backColor);
		alfa0 += fround(3600.0*vValues[i]/sum);
	}
	if(showPercent) {
		alfa0 = -900 + static_cast<int>(pieAngle);
		for(int i = 0; i < vValues.GetCount(); i++) {
			a0 = alfa0;                            		              
			alfa0 += fround(3600.0*vValues[i]/sum);
			String percent = GetPercent(vValues[i],sum);
			Upp::Font scaledFont;
			scaledFont.Height(scale*StdFont().GetHeight());
			scaledFont.Width(scale*StdFont().GetWidth());
			Size szz = GetTextSizeSpace(percent, scaledFont);
		
			int px = int(circ_x + scale*circ_r*cos(M_PI*(alfa0+a0)/3600)/1.3 - szz.cx/2.);
			int py = int(circ_y + scale*circ_r*sin(M_PI*(alfa0+a0)/3600)/1.3 - szz.cy/2.);
			w.DrawRect(px,   py,   			  szz.cx, 	  1, 		percentBack);
			w.DrawRect(px-1, py + 1, 		  szz.cx + 2, szz.cy-2, percentBack);
			w.DrawRect(px,   py + szz.cy - 1, szz.cx, 	  1, 		percentBack);
			w.DrawText(px, py, percent, scaledFont);
		}
	}
	if(showLegend) {
		Upp::Font scaledFont;
		scaledFont.Height(scale*legendFont.GetHeight());
		scaledFont.Width(scale*legendFont.GetWidth());
		int fh = scale*legendFont.GetHeight();
		int nr = GetCatCount();
		int legendWidth = 0;
		int legendHeight = (1 + nr)*scaledFont.GetHeight();
		for(int i = 0; i < nr; i++) 
			legendWidth = max<int>(legendWidth, GetTextSizeSpace(vNames[i], scaledFont).cx);
		legendWidth += fround(2.2*fh);
		
		double leg_x = -legendLeft + sz.cx - legendWidth;
		double leg_y;
		if (IsNull(legendTop))
			leg_y = int(circ_y - scale*legendHeight/2.);
		else
			leg_y = legendTop*scale;
		w.DrawRect(int(leg_x), int(leg_y), scale*legendWidth, scale*legendHeight, legendBackColor);
		
		int dly = scale*legendHeight/nr;
		for(int i = 0; i < nr; i++) {
			w.DrawRect(int(leg_x) + 2*scale, int(leg_y + i*dly + dly/2. - fh/2.), fh, fh, vColors[i]);
			w.DrawText(int(leg_x) + fround(1.8*fh), int(leg_y) + i*dly + int((dly - scaledFont.GetLineHeight())/2),
						vNames[i], scaledFont, legendTextColor);
		}
	}
}		                  

String PieDraw::GetPercent(double a, double total) {
	double p = a*100/total;
	return FormatDouble(p, 1) + '%';
}

Drawing PieDraw::GetDrawing(int scale) {
	DrawingDraw ddw(scale*GetSize());
	PaintPie(ddw, scale);
	return ddw;
}

Image PieDraw::GetImage(int scale) {
	int mode = MODE_ANTIALIASED;
	
	ImageBuffer ib(scale*GetSize());	
	BufferPainter bp(ib, mode);	
	
	bp.LineCap(LINECAP_BUTT);
	bp.LineJoin(LINEJOIN_MITER);
	if(IsNull(backColor))
		bp.Clear(RGBAZero());
	PaintPie(bp, scale);

	return ib;
}

PieDraw::PieDraw(): backColor(White), titleFont(StdFont(16)), titleColor(Black), titlePos(TOP), 
					titleGap(2), showPercent(true), percentBack(Null),
					legendFont(StdFont()), legendTextColor(Black), legendBackColor(Null),
					showLegend(true), legendLeft(10), legendTop(Null), pieAngle(0), 
					pieMarginLeft(40), pieMarginTop(40), 
					pieMarginRight(40), pieMarginBottom(40)
{}

PieDraw::~PieDraw(){}


