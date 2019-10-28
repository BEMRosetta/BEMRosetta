#ifndef _PopUpText_PopUpText_h
#define _PopUpText_PopUpText_h

class PopUpText {
public:
	Point pos;
	Size sz;
	String text;
	Font font = StdFont();
	Color color = SColorText;
	Color background = SColorFace;
	bool show = false;
	int xmargin = StdFont().GetWidth('I');
	
	PopUpText &SetText(String _text);
	PopUpText &SetColor(Color _color)		{color = _color;	return *this;}
	PopUpText &SetBackground(Color _color)	{background = _color;return *this;}
	PopUpText &Show(bool _show = true)		{show = _show;		return *this;}
	PopUpText &Hide()						{
	show = false;		
	return *this;}
	PopUpText &Move(Point& p)				{pos = p;			return *this;}	
	Size GetSize()							{return sz;}
	
	template <class T>
	void DoPaint(T& w) {
		if (!show)
			return;
		Rect rect(pos.x, pos.y, pos.x + sz.cx + 3*xmargin, pos.y + sz.cy);
		FillRectangle(w, rect, background);
		DrawRectangle(w, rect, 1, 1, Black());
		DrawText(w, rect.left + xmargin, rect.top, 0, text, font, color);
	}
};

#endif

