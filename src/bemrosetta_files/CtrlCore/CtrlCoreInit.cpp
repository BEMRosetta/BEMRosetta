#include "CtrlCore.h"

namespace Upp {

#define  TFILE <CtrlCore/CtrlCore.t>
#include <Core/t.h>

static Image sRenderGlyph(int cx, int x, Font font, int chr, int py, int pcy)
{
	ImageDraw iw(cx, pcy);
	iw.DrawRect(0, 0, cx, pcy, White);
	iw.DrawText(x, -py, WString(chr, 1), font, Black);
	return iw;
}

void SetRenderGlyph(Image (*f)(int cx, int x, Font font, int chr, int py, int pcy));

INITIALIZER(CtrlCore) {
	SetRenderGlyph(sRenderGlyph);
}

}
