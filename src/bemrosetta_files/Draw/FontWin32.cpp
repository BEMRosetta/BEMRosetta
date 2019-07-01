#include "Draw.h"

namespace Upp {

#ifndef CUSTOM_FONTSYS

#ifdef PLATFORM_WIN32

#define LLOG(x)     //  LOG(x)
#define LTIMING(x)  //  TIMING(x)

#define FONTCACHE 96

struct HFontEntry {
	Font    font;
	HFONT   hfont;
	int     angle;
};

HFONT GetWin32Font(Font fnt, int angle)
{
	LTIMING("GetWin32Font");
	static HFontEntry cache[FONTCACHE];
	ONCELOCK {
		for(int i = 0; i < FONTCACHE; i++)
			cache[i].font.Height(-30000);
	}
	HFontEntry be;
	be = cache[0];
	for(int i = 0; i < FONTCACHE; i++) {
		HFontEntry e = cache[i];
		if(i)
			cache[i] = be;
		if(e.font == fnt && e.angle == angle) {
			if(i)
				cache[0] = e;
			return e.hfont;
		}
		be = e;
	}
	LTIMING("GetWin32Font2");
	if(be.hfont)
		DeleteObject(be.hfont);

	be.font = fnt;
	be.angle = angle;
	be.hfont = CreateFont(
		fnt.GetHeight() ? -abs(fnt.GetHeight()) : -12,
		fnt.GetWidth(), angle, angle, fnt.IsBold() ? FW_BOLD : FW_NORMAL,
		fnt.IsItalic(), fnt.IsUnderline(), fnt.IsStrikeout(),
		fnt.GetFace() == Font::SYMBOL ? SYMBOL_CHARSET : DEFAULT_CHARSET,
		fnt.IsTrueTypeOnly() ? OUT_TT_ONLY_PRECIS : OUT_DEFAULT_PRECIS,
		CLIP_DEFAULT_PRECIS,
		fnt.IsNonAntiAliased() ? NONANTIALIASED_QUALITY : DEFAULT_QUALITY,
		DEFAULT_PITCH|FF_DONTCARE,
		fnt.GetFaceName()
	);

	cache[0] = be;
	return be.hfont;
}

static void Win32_GetGlyphIndices(HDC hdc, LPCWSTR s, int n, LPWORD r, DWORD flag)
{
	typedef DWORD (WINAPI *GGIW)(HDC, LPCWSTR, int, LPWORD, DWORD);
	static GGIW fn;
	ONCELOCK
		if(HMODULE hDLL = LoadLibrary("gdi32"))
			fn = (GGIW) GetProcAddress(hDLL, "GetGlyphIndicesW");
	if(fn)
		fn(hdc, s, n, r, flag);
	else
		memset(r, 0, n * sizeof(WORD));
}

CommonFontInfo GetFontInfoSys(Font font)
{
	CommonFontInfo fi;
	memset(&fi, 0, sizeof(fi));
	HFONT hfont = GetWin32Font(font, 0);
	if(hfont) {
		HDC hdc = ::CreateIC("DISPLAY", NULL, NULL, NULL);
		HFONT ohfont = (HFONT) ::SelectObject(hdc, hfont);
		TEXTMETRIC tm;
		::GetTextMetrics(hdc, &tm);
		fi.ascent = tm.tmAscent;
		fi.descent = tm.tmDescent;
		fi.external = tm.tmExternalLeading;
		fi.internal = tm.tmInternalLeading;
		fi.overhang = tm.tmOverhang;
		fi.avewidth = tm.tmAveCharWidth;
		fi.maxwidth = tm.tmMaxCharWidth;
		fi.firstchar = tm.tmFirstChar;
		fi.charcount = tm.tmLastChar - tm.tmFirstChar + 1;
		fi.default_char = tm.tmDefaultChar;
		fi.fixedpitch = (tm.tmPitchAndFamily & TMPF_FIXED_PITCH) == 0;
		fi.ttf = fi.scaleable = tm.tmPitchAndFamily & TMPF_TRUETYPE;
		if(fi.scaleable) {
			ABC abc;
			GetCharABCWidths(hdc, 'o', 'o', &abc);
			fi.spacebefore = abc.abcA;
			fi.spaceafter = abc.abcC;
		}
		else
			fi.spacebefore = fi.spaceafter = 0;
		::SelectObject(hdc, ohfont);
		::DeleteDC(hdc);
	}
	return fi;
}

static VectorMap<String, FaceInfo> *sList;

static int CALLBACK Win32_AddFace(const LOGFONT *logfont, const TEXTMETRIC *, dword type, LPARAM param)
{
	const char *facename = (const char *)param;
	if(facename && stricmp(logfont->lfFaceName, facename))
		return 1;

	if(logfont->lfFaceName[0] == '@')
		return 1;
	
	String name = FromSystemCharset(logfont->lfFaceName);

	if(FindIndex(Split("Courier New CE;Courier New CYR;Courier New Greek;"
	                   "Courier New TUR;Courier New Baltic;Arial CE;Arial CYR;"
	                   "Arial Greek;Arial TUR;Arial Baltic;Arial CE;Times New Roman CE;"
	                   "Times New Roman CYR;Times New Roman Greek;Times New Roman TUR;"
	                   "Times New Roman Baltic;Times New Roman CE", ';'), name) >= 0)
		return 1;

	int q = sList->Find(name);
	FaceInfo& f = q < 0 ? sList->Add(logfont->lfFaceName) : (*sList)[q];
	f.name = FromSystemCharset(logfont->lfFaceName);

	if(q < 0)
		f.info = Font::SCALEABLE;
	if((logfont->lfPitchAndFamily & 3) == FIXED_PITCH)
		f.info |= Font::FIXEDPITCH;
	if(logfont->lfPitchAndFamily & FF_ROMAN)
		f.info |= Font::SERIFSTYLE;
	if(logfont->lfPitchAndFamily & FF_SCRIPT)
		f.info |= Font::SCRIPTSTYLE;
	if(!(type & TRUETYPE_FONTTYPE))
		f.info &= ~Font::SCALEABLE;
	if(logfont->lfCharSet == SYMBOL_CHARSET || logfont->lfCharSet == OEM_CHARSET)
		f.info |= Font::SPECIAL;
	return facename ? 0 : 1;
}

static int Win32_EnumFace(HDC hdc, const char *face)
{
	return EnumFontFamilies(hdc, face, Win32_AddFace, (LPARAM)face);
}

static void Win32_ForceFace(HDC hdc, const char *face, const char *aface)
{
	if(!aface)
		aface = "Arial";
	if(Win32_EnumFace(hdc, face) && Win32_EnumFace(hdc, aface))
		Panic("Missing font " + String(face));
}

Vector<FaceInfo> GetAllFacesSys()
{
	VectorMap<String, FaceInfo> list;
	sList = &list;
	list.Add("STDFONT").name = "STDFONT";
	list.Top().info = 0;
	HDC hdc = ::CreateIC("DISPLAY", NULL, NULL, NULL);
	Win32_ForceFace(hdc, "Times New Roman", "Arial");
	Win32_ForceFace(hdc, "Arial", "Arial");
	Win32_ForceFace(hdc, "Courier New", "Arial");
	Win32_ForceFace(hdc, "Symbol", "Arial");
	Win32_ForceFace(hdc, "WingDings", "Arial");
	Win32_ForceFace(hdc, "Tahoma", "Arial");
	Win32_EnumFace(hdc, NULL);
	DeleteDC(hdc);
	return list.PickValues();
}

#define GLYPHINFOCACHE 31

#ifdef flagWINGL
#include <CoreGl/FontGl.h>
#include <CoreGl/ResGl.h>
GlyphInfo  GetGlyphInfoSys(Font font, int chr)
{
	static GlyphInfo gi;
	const OpenGLFont& fi = resources.GetFont(font);
	gi.width = chr < fi.chars.GetCount() ? int(fi.chars[chr].xadvance * fi.scale + 0.5f) : 0;
	gi.lspc = 0;
	gi.rspc = 0;
	return gi;
}
#else
GlyphInfo  GetGlyphInfoSys(Font font, int chr)
{
	static Font      fnt[GLYPHINFOCACHE];
	static int       pg[GLYPHINFOCACHE];
	static GlyphInfo li[GLYPHINFOCACHE][256];
	
	ONCELOCK {
		for(int i = 0; i < GLYPHINFOCACHE; i++)
			pg[i] = -1;
	}
	
	int page = chr >> 8;
	int q = CombineHash(font, page) % GLYPHINFOCACHE;
	
	if(fnt[q] != font || pg[q] != page) {
		fnt[q] = font;
		pg[q] = page;
		HFONT hfont = GetWin32Font(font, 0);
		if(!hfont) {
			GlyphInfo n;
			memset(&n, 0, sizeof(GlyphInfo));
			return n;
		}
		HDC hdc = CreateIC("DISPLAY", NULL, NULL, NULL);
		HFONT ohfont = (HFONT) ::SelectObject(hdc, hfont);
		int from = page << 8;
		GlyphInfo *t = li[q];

		bool abca = false, abcw = false;
		Buffer<ABC> abc(256);
		abcw = ::GetCharABCWidthsW(hdc, from, from + 256 - 1, abc);
		if(!abcw)
			abca = ::GetCharABCWidthsA(hdc, from, from + 256 - 1, abc);
		if(abcw) {
			for(ABC *s = abc, *lim = abc + 256; s < lim; s++, t++) {
				t->width = s->abcA + s->abcB + s->abcC;
				t->lspc = s->abcA;
				t->rspc = s->abcC;
			}
		}
		else {
			Buffer<int> wb(256);
			::GetCharWidthW(hdc, from, from + 256 - 1, wb);
			for(int *s = wb, *lim = wb + 256; s < lim; s++, t++) {
				t->width = *s - GetFontInfoSys(font).overhang;
				if(abca) {
					ABC aa = abc[(byte)ToAscii(from++)];
					t->lspc = aa.abcA;
					t->rspc = aa.abcC;
				}
				else
					t->lspc = t->rspc = 0;
			}
		}

		WORD pos[256];
		WCHAR wch[256];
		for(int i = 0; i < 256; i++)
			wch[i] = from + i;
		Win32_GetGlyphIndices(hdc, wch, 256, pos, 1);
		for(int i = 0; i < 256; i++)
			if(pos[i] == 0xffff) {
				li[q][i].width = (int16)0x8000;
				li[q][i].lspc = li[q][i].rspc = 0;
			}

		::SelectObject(hdc, ohfont);
		::DeleteDC(hdc);
	}
	return li[q][chr & 255];
}
#endif

String GetFontDataSys(Font font)
{
	String r;
	HFONT hfont = GetWin32Font(font, 0);
	if(hfont) {
		HDC hdc = CreateIC("DISPLAY", NULL, NULL, NULL);
		HFONT ohfont = (HFONT) ::SelectObject(hdc, hfont);
		DWORD size = GetFontData(hdc, 0, 0, NULL, 0);
		if(size == GDI_ERROR) {
			LLOG("PdfDraw::Finish: GDI_ERROR on font " << pdffont.GetKey(i));
		}
		else {
			StringBuffer b(size);
			GetFontData(hdc, 0, 0, b, size);
			r = b;
		}
		::SelectObject(hdc, ohfont);
		::DeleteDC(hdc);
	}
	return r;
}

double fx_to_dbl(const FIXED& p) {
	return double(p.value) + double(p.fract) * (1.0 / 65536.0);
}

Pointf fx_to_dbl(const Pointf& pp, const POINTFX& p) {
	return Pointf(pp.x + fx_to_dbl(p.x), pp.y - fx_to_dbl(p.y));
}

void RenderCharPath(const char* gbuf, unsigned total_size, FontGlyphConsumer& sw, double xx, double yy)
{
	const char* cur_glyph = gbuf;
	const char* end_glyph = gbuf + total_size;
	Pointf pp(xx, yy);
	while(cur_glyph < end_glyph) {
		const TTPOLYGONHEADER* th = (TTPOLYGONHEADER*)cur_glyph;
		const char* end_poly = cur_glyph + th->cb;
		const char* cur_poly = cur_glyph + sizeof(TTPOLYGONHEADER);
		sw.Move(fx_to_dbl(pp, th->pfxStart));
		while(cur_poly < end_poly) {
			const TTPOLYCURVE* pc = (const TTPOLYCURVE*)cur_poly;
			if (pc->wType == TT_PRIM_LINE)
				for(int i = 0; i < pc->cpfx; i++)
					sw.Line(fx_to_dbl(pp, pc->apfx[i]));
			if (pc->wType == TT_PRIM_QSPLINE)
				for(int u = 0; u < pc->cpfx - 1; u++) {
					Pointf b = fx_to_dbl(pp, pc->apfx[u]);
					Pointf c = fx_to_dbl(pp, pc->apfx[u + 1]);
					if (u < pc->cpfx - 2)
						c = Mid(b, c);
					sw.Quadratic(b, c);
				}
			cur_poly += sizeof(WORD) * 2 + sizeof(POINTFX) * pc->cpfx;
		}
		sw.Close();
		cur_glyph += th->cb;
    }
}

void RenderCharacterSys(FontGlyphConsumer& sw, double x, double y, int ch, Font fnt)
{
	HFONT hfont = GetWin32Font(fnt, 0);
	VERIFY(hfont);
	if(hfont) {
		HDC hdc = CreateIC("DISPLAY", NULL, NULL, NULL);
		HFONT ohfont = (HFONT) ::SelectObject(hdc, hfont);
		GLYPHMETRICS gm;
		MAT2 m_matrix;
		memset(&m_matrix, 0, sizeof(m_matrix));
		m_matrix.eM11.value = 1;
		m_matrix.eM22.value = 1;
		int gsz = GetGlyphOutlineW(hdc, ch, GGO_NATIVE|GGO_UNHINTED, &gm, 0, NULL, &m_matrix);
		if(gsz >= 0) {
			StringBuffer gb(gsz);
			gsz = GetGlyphOutlineW(hdc, ch, GGO_NATIVE|GGO_UNHINTED, &gm, gsz, ~gb, &m_matrix);
			if(gsz >= 0)
				RenderCharPath(~gb, gsz, sw, x, y + fnt.GetAscent());
		}
		::SelectObject(hdc, ohfont);
		::DeleteDC(hdc);
	}
}

#endif

#endif

}
