#include <CtrlLib/CtrlLib.h>

using namespace Upp;

#include "CtrlScroll.h"

CtrlScroll::CtrlScroll() {
	AddFrame(scroll);
	scroll.AutoHide();
	scroll.WhenScroll = THISBACK(OnScroll);
	hsizepos = vsizepos = false;
}

CtrlScroll &CtrlScroll::AddPane(Ctrl& c, bool scrollH, bool scrollV) { 
	pane = &c; 
	if (scrollH) 
		scroll.x.Enable();
	else {
		c.HSizePos();
		hsizepos = true;
	}
	if (scrollV) 
		scroll.y.Enable();	
	else {
		c.VSizePos();
		vsizepos = true;
	}
	Add(c); 
	return *this;
}
	
void CtrlScroll::Scroll(const Point& p) {
	if(!HasPane()) 
		return;
	Rect _r = pane->GetRect();
	Rect r(-p, _r.GetSize());
	pane->SetRect(r);
	if (hsizepos)
		pane->HSizePos();
	if (vsizepos)
		pane->VSizePos();
	WhenScrolled();
}

void CtrlScroll::OnScroll() {
	Scroll(scroll.Get());
}

void CtrlScroll::Layout() {
	if(!HasPane()) 
		return;
	Size psz = GetSize();
	Size tsz = pane->GetSize();
	scroll.Set(Point(0, 0), psz, tsz);
}
