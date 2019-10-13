#ifndef _CtrlScroll_CtrlScroll_h_
#define _CtrlScroll_CtrlScroll_h_

using namespace Upp;

class CtrlScroll : public StaticRect {
public:
	typedef CtrlScroll CLASSNAME;
	CtrlScroll();

	CtrlScroll &AddPane(Ctrl& c, bool scrollH = true, bool scrollV = true);
	CtrlScroll &AddPaneH(Ctrl& c) 	{return AddPane(c, true, false);}
	CtrlScroll &AddPaneV(Ctrl& c)	{return AddPane(c, false, true);}
	inline bool HasPane() const 	{return (~pane != NULL);}

	Event<> WhenScrolled;
	
public:
	virtual void Layout();
	
	void Scroll(const Point& p);
	void OnScroll();

	ScrollBars scroll;

	Ptr<Ctrl> pane;
	bool hsizepos, vsizepos;
};

#endif
