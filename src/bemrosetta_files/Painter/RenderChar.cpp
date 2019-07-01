#include "Painter.h"

namespace Upp {

struct GlyphPainter : NilPainter, LinearPathConsumer {
	Vector<float>       glyph;
	double              tolerance;
	Pointf              pos, move;
	
	virtual void LineOp(const Pointf& p, bool);
	virtual void MoveOp(const Pointf& p, bool);
	virtual void QuadraticOp(const Pointf& p1, const Pointf& p, bool);
	virtual void CubicOp(const Pointf& p1, const Pointf& p2, const Pointf& p, bool);
	virtual void CloseOp();

	virtual void Line(const Pointf& p);
	virtual void Move(const Pointf& p);
};

void GlyphPainter::Move(const Pointf& p)
{
	glyph.Add((float)1e31);
	Line(p);
	move = pos;
}

void GlyphPainter::Line(const Pointf& p)
{
	glyph.Add((float)p.x);
	glyph.Add((float)p.y);
	pos = p;
}

void GlyphPainter::MoveOp(const Pointf& p, bool)
{
	Move(p);
}

void GlyphPainter::LineOp(const Pointf& p, bool)
{
	Line(p);
}

void GlyphPainter::QuadraticOp(const Pointf& p1, const Pointf& p, bool)
{
	ApproximateQuadratic(*this, pos, p1, p, tolerance);
	pos = p;
}

void GlyphPainter::CubicOp(const Pointf& p1, const Pointf& p2, const Pointf& p, bool)
{
	ApproximateCubic(*this, pos, p1, p2, p, tolerance);
	pos = p;
}

void GlyphPainter::CloseOp()
{
	if(move != pos && !IsNull(move))
		Line(move);
}

struct GlyphKey {
	Font   fnt;
	int    chr;
	double tolerance;
	
	bool operator==(const GlyphKey& b) const {
		return fnt == b.fnt && chr == b.chr && tolerance == b.tolerance;
	}
	unsigned GetHashValue() const {
		return CombineHash(fnt, chr, tolerance);
	}
};

struct sMakeGlyph : LRUCache<Value, GlyphKey>::Maker {
	GlyphKey gk;

	GlyphKey Key() const     { return gk; }
	int      Make(Value& v) const {
		GlyphPainter gp;
		gp.move = gp.pos = Null;
		gp.tolerance = gk.tolerance;
		PaintCharacter(gp, Pointf(0, 0), gk.chr, gk.fnt);
		int sz = gp.glyph.GetCount() * 4;
		v = RawPickToValue(pick(gp.glyph));
		return sz;
	}
};

void ApproximateChar(LinearPathConsumer& t, Pointf at, int ch, Font fnt, double tolerance)
{
	PAINTER_TIMING("ApproximateChar");
	Value v;
	INTERLOCKED {
		PAINTER_TIMING("ApproximateChar::Fetch");
		static LRUCache<Value, GlyphKey> cache;
		cache.Shrink(500000);
		sMakeGlyph h;
		h.gk.fnt = fnt;
		h.gk.chr = ch;
		h.gk.tolerance = tolerance;
		v = cache.Get(h);
#ifdef _DEBUG0
		DLOG("==== ApproximateChar " << ch << " " << (char)ch << " " << fnt << ", tolerance: " << tolerance);
		DDUMP(ValueTo< Vector<float> >(v));
		GlyphPainter chp;
		chp.move = chp.pos = Null;
		chp.tolerance = tolerance;
		PaintCharacter(chp, Pointf(0, 0), ch, fnt);
		DDUMP(chp.glyph);
		ASSERT(ValueTo< Vector<float> >(v) == chp.glyph);
#endif
	}
	const Vector<float>& g = ValueTo< Vector<float> >(v);
	int i = 0;
	while(i < g.GetCount()) {
		Pointf p;
		p.x = g[i++];
		if(p.x > 1e30) {
			p.x = g[i++];
			p.y = g[i++];
			t.Move(p + at);
		}
		else {
			PAINTER_TIMING("ApproximateChar::Line");
			p.y = g[i++];
			t.Line(p + at);
		}
	}
}

}
