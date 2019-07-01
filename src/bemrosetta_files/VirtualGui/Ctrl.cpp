#include "Local.h"

#ifdef VIRTUALGUI

#define LLOG(x) // DLOG(x)

NAMESPACE_UPP

void Ctrl::GuiPlatformConstruct()
{
}

void Ctrl::GuiPlatformRemove()
{
}

void Ctrl::GuiPlatformGetTopRect(Rect& r) const
{
}

bool Ctrl::GuiPlatformRefreshFrameSpecial(const Rect& r)
{
	return false;
}

bool Ctrl::GuiPlatformSetFullRefreshSpecial()
{
	return false;
}

void Ctrl::PaintCaret(SystemDraw& w)
{
}

String GuiPlatformGetKeyDesc(dword key)
{
	return Null;
}

void Ctrl::GuiPlatformSelection(PasteClip&)
{
}

void GuiPlatformAdjustDragImage(ImageBuffer&)
{
}

bool GuiPlatformHasSizeGrip()
{
	return true;
}

void GuiPlatformGripResize(TopWindow *q)
{
	q->GripResize();
}

Color GuiPlatformGetScreenPixel(int x, int y)
{
	return Null;
}

void GuiPlatformAfterMenuPopUp()
{
}

String Ctrl::Name() const {
	GuiLock __;
#ifdef CPU_64
	String s = String(typeid(*this).name()) + " : 0x" + FormatIntHex(this);
#else
	String s = String(typeid(*this).name()) + " : " + Format("0x%x", (int) this);
#endif
	if(IsChild())
		s << "(parent " << String(typeid(*parent).name()) << ")";
	return s;
}

void  Ctrl::SetMouseCursor(const Image& image)
{
	GuiLock __;
	fbCursorImage = image;
	if(SystemCursor)
		VirtualGuiPtr->SetMouseCursor(image);
}

END_UPP_NAMESPACE

#endif
