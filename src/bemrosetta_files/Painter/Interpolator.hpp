inline
void LinearInterpolator::Dda2::Set(int p1, int p2, int len)
{
	count = len <= 0 ? 1 : len;
	lift = (p2 - p1) / count;
	rem = (p2 - p1) % count;
	mod = rem;
	p = p1;
	if(mod <= 0) {
		mod += count;
		rem += count;
		lift--;
    }
	mod -= count;
}

inline
int LinearInterpolator::Dda2::Get()
{
	int pp = p;
	mod += rem;
	p += lift;
	if(mod > 0) {
		mod -= count;
		p++;
	}
	return pp;
}

inline
void LinearInterpolator::Begin(int x, int y, int len)
{
	Pointf p1 = xform.Transform(Pointf(x, y));
	Pointf p2 = xform.Transform(Pointf(x + len, y));
	ddax.Set(Q8(p1.x), Q8(p2.x), len);
	dday.Set(Q8(p1.y), Q8(p2.y), len);
}

inline
Point LinearInterpolator::Get()
{
	return Point(ddax.Get(), dday.Get());
}
