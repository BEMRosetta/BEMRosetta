#include "Painter.h"

namespace Upp {

bool Painter::ReadBool(CParser& p)
{
	while(p.Char(','));
	if(p.Char('1')) return true;
	p.Char('0');
	return false;
}

double Painter::ReadDouble(CParser& p)
{
	while(p.Char(','));
	return p.IsDouble2() ? p.ReadDouble() : 0;
}

Pointf Painter::ReadPoint(CParser& p, Pointf current, bool rel)
{
	Pointf t;
	t.x = ReadDouble(p);
	t.y = ReadDouble(p);
	if(rel)
		t += current;
	return t;
}

Painter& Painter::Path(CParser& p)
{
	Pointf current(0, 0);
	bool done = false;
	while(!p.IsEof()) {
		int c = p.GetChar();
		p.Spaces();
		bool rel = IsLower(c);
		Pointf t, t1, t2;
		switch(ToUpper(c)) {
		case 'M':
			current = ReadPoint(p, current, rel);
			Move(current);
		case 'L':
			while(p.IsDouble2()) {
				current = ReadPoint(p, current, rel);
				Line(current);
			}
			done = true;
			break;
		case 'Z':
			Close();
			done = true;
			break;
		case 'H':
			while(p.IsDouble2()) {
				current.x = p.ReadDouble() + rel * current.x;
				Line(current);
				done = true;
			}
			break;
		case 'V':
			while(p.IsDouble2()) {
				current.y = p.ReadDouble() + rel * current.y;
				Line(current);
				done = true;
			}
			break;
		case 'C':
			while(p.IsDouble2()) {
				t1 = ReadPoint(p, current, rel);
				t2 = ReadPoint(p, current, rel);
				current = ReadPoint(p, current, rel);
				Cubic(t1, t2, current);
				done = true;
			}
			break;
		case 'S':
			while(p.IsDouble2()) {
				t2 = ReadPoint(p, current, rel);
				current = ReadPoint(p, current, rel);
				Cubic(t2, current);
				done = true;
			}
			break;
		case 'Q':
			while(p.IsDouble2()) {
				t1 = ReadPoint(p, current, rel);
				current = ReadPoint(p, current, rel);
				Quadratic(t1, current);
				done = true;
			}
			break;
		case 'T':
			while(p.IsDouble2()) {
				current = ReadPoint(p, current, rel);
				Quadratic(current);
				done = true;
			}
			break;
		case 'A':
			while(p.IsDouble2()) {
				t1 = ReadPoint(p, Pointf(0, 0), false);
				double xangle = ReadDouble(p);
				bool large = ReadBool(p);
				bool sweep = ReadBool(p);
				current = ReadPoint(p, current, rel);
				SvgArc(t1, xangle * M_PI / 180.0, large, sweep, current);
				done = true;
			}
			break;
		default:
			if(!done)
				Move(0, 0); // to clear previous path
			return *this;
		}
	}
	if(!done)
		Move(0, 0); // to clear previous path
	return *this;
}

Painter& Painter::Path(const char *path)
{
	try {
		CParser p(path);
		Path(p);
	}
	catch(CParser::Error) {}
	return *this;
}

}
