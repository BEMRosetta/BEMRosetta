#include <CtrlLib/CtrlLib.h>

using namespace Upp;

#include "PopUpText.h"

static Size GetEditSize(const String &_str, const Font &font) {
	WString str(_str);
	Size ret(0, 0);
	int retx = 0, nlines = 1;
	for (int i = 0; i < str.GetCount(); ++i) {
		int c = str[i];
		if (c == '\n') {
			nlines++;
			ret.cx = max(ret.cx, retx);
			retx = 0;
		} else
			retx += font.GetWidth(c);
	}
	ret.cx = max(ret.cx, retx);
	ret.cy = nlines*font.GetHeight() + font.GetDescent();
	return ret;
}

PopUpText &PopUpText::SetText(String _text) {
	text = _text;		
	sz = GetEditSize(text, font);
	return *this;
}