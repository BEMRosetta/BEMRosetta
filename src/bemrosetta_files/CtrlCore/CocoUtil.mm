#include "CocoMM.h"

#ifdef PLATFORM_COCOA

namespace Upp {

WString ToWString(CFStringRef s)
{
	if(!s) return Null;
	CFIndex l = CFStringGetLength(s);
	if(!l) return Null;
	WStringBuffer b(l);
    CFStringGetCharacters(s, CFRangeMake(0, l), (UniChar *)~b);
    return b;
}

String ToString(CFStringRef s)
{
	return ToWString(s).ToString();
}

void DrawDragLine(SystemDraw& w, bool horz, int x, int y, int len, int n, const int *pattern, Color color, int animation)
{
	if(len <= 0)
		return;
	if(horz)
		w.Clip(x, y, len, n);
	else
		w.Clip(x, y, n, len);
	
	(horz ? x : y) -= animation;
	len += animation;
	bool ch = false;
	while(len > 0) {
		int segment = pattern[ch];
		int d = segment + pattern[2];
		if(horz) {
			w.DrawRect(x, y, segment, n, color);
			x += d;
		}
		else {
			w.DrawRect(x, y, n, segment, color);
			y += d;
		}
		len -= d;
		ch = !ch;
	}
	w.End();
}

void DrawDragFrame(SystemDraw& w, const Rect& r, int n, const int *pattern, Color color, int animation)
{
	DrawDragLine(w, true, r.left, r.top, r.GetWidth(), n, pattern, color, animation);
	DrawDragLine(w, false, r.left, r.top + n, r.GetHeight() - 2 * n, n, pattern, color, animation);
	DrawDragLine(w, false, r.right - n, r.top + n, r.GetHeight() - 2 * n, n, pattern, color, animation);
	DrawDragLine(w, true, r.left, r.bottom - n, r.GetWidth(), n, pattern, color, animation);
}

void DrawDragRect(Ctrl& q, const DrawDragRectInfo& f)
{
	ViewDraw w(&q);
	w.Clip(f.clip);
	static int dashes[3][3] = {
		{ 32, 32, 0 },
		{ 1, 1, 1 },
		{ 5, 1, 2 },
	};
	const int *dash = dashes[minmax(f.type, 0, 2)];
	Color color = InvertColor;
	DrawDragFrame(w, f.rect1, f.n, dash, color, f.animation);
	DrawDragFrame(w, f.rect2, f.n, dash, color, f.animation);
	w.End();
}

void DrawDragRect(Ctrl& q, const Rect& rect1, const Rect& rect2, const Rect& clip, int n,
                  Color color, int type, int animation)
{
	Ctrl *top = q.GetTopCtrl();
	if(top && top->top->coco) {
		Rect sv = q.GetScreenView();
		Rect tv = top->GetScreenView();
		Point off = sv.TopLeft() - tv.TopLeft();
		DrawDragRectInfo& f = top->top->dr.Create();
		f.rect1 = rect1 + off;
		f.rect2 = rect2 + off;
		f.clip = (clip & q.GetSize()) + off;
		f.n = n;
		f.color = color;
		f.type = type;
		f.animation = animation;
		DrawDragRect(*top, f);
	}
}

void FinishDragRect(Ctrl& q)
{
	Ctrl *top = q.GetTopCtrl();
	if(top && top->top && top->top->dr) {
		DrawDragRect(*top, *top->top->dr);
		top->top->dr.Clear();
	}
}

/*
Size GetScreenSize()
{
	return ScreenInfo().GetPageSize();
}
*/

};

#endif
