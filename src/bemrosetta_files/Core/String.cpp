#include "Core.h"

namespace Upp {

#ifdef _DEBUG
void String0::Dsyn()
{
	String *d_str = static_cast<String *>(this);
	d_str->s = Begin();
	d_str->len = GetCount();
}
#endif

String0::Rc String0::voidptr[2];

void String0::LSet(const String0& s)
{
	w[2] = s.w[2];
	w[3] = s.w[3];
	if(s.IsRef()) {
		ptr = s.ptr;
		if(ptr != (char *)(voidptr + 1))
			AtomicInc(s.Ref()->refcount);
	}
	else {
		ptr = (char *)MemoryAlloc32();
//		qword *d = qptr;
//		const qword *q = s.qptr;
//		d[0] = q[0]; d[1] = q[1]; d[2] = q[2]; d[3] = q[3];
		fast_copy128(ptr, s.ptr);
		fast_copy128(ptr + 16, s.ptr + 16);
	}
}

void String0::LFree()
{
	if(IsRef()) {
		if(ptr != (char *)(voidptr + 1)) {
			Rc *rc = Ref();
			ASSERT(rc->refcount > 0);
			if(AtomicDec(rc->refcount) == 0) MemoryFree(rc);
		}
	}
	else
		MemoryFree32(ptr);
}

dword String0::LEqual(const String0& s) const
{
	int l = GetCount();
	if(s.GetCount() != l) return 1;
#ifdef CPU_64
	const qword *qa = (const qword *)Begin();
	const qword *qb = (const qword *)s.Begin();
	while(l >= 8) {
		if(*qa++ != *qb++) return 1;
		l -= 8;
	}
	const dword *da = (const dword *)qa;
	const dword *db = (const dword *)qb;
	if((l & 4) && *da++ != *db++) return 1;
#else
	const dword *da = (const dword *)Begin();
	const dword *db = (const dword *)s.Begin();
	while(l >= 4) {
		if(*da++ != *db++) return 1;
		l -= 4;
	}
#endif
	const word *wa = (const word *)da;
	const word *wb = (const word *)db;
	if((l & 2) && *wa++ != *wb++) return 1;
	return (l & 1) ? *(const char *)wa != *(const char *)wb : 0;
}

unsigned String0::LHashValue() const
{
	int l = LLen();
	if(l < 15) {
		dword w[4];
		w[0] = w[1] = w[2] = w[3] = 0;
		SVO_MEMCPY((char *)w, ptr, l);
		((byte *)w)[SLEN] = l;
		return CombineHash(w[0], w[1], w[2], w[3]);
	}
	return memhash(ptr, l);
}

int String0::LCompare(const String0& s) const
{
	const char *a = Begin();
	int la = GetLength();
	const char *b = s.Begin();
	int lb = s.GetLength();
	int q = fast_memcmp(a, b, min(la, lb));
	return q ? q : SgnCompare(la, lb);
}

char *String0::Alloc(int count, char& kind)
{
	if(count < 32) {
		kind = MEDIUM;
		return (char *)MemoryAlloc32();
	}
	size_t sz = sizeof(Rc) + count + 1;
	Rc *rc = (Rc *)MemoryAllocSz(sz);
	rc->alloc = count == INT_MAX ? INT_MAX : (int)sz - sizeof(Rc) - 1;
	rc->refcount = 1;
	kind = min(rc->alloc, 255);
	return rc->GetPtr();
}

char *String0::Insert(int pos, int count, const char *s)
{
	ASSERT(pos >= 0 && count >= 0 && pos <= GetCount());
	int len = GetCount();
	int newlen = len + count;
	if(newlen < len) // overflow, string >2GB
		Panic("String is too big!");
	char *str = (char *)Begin();
	if(newlen < GetAlloc() && !IsSharedRef() && (!s || s < str || s > str + len)) {
		if(pos < len)
			memmove(str + pos + count, str + pos, len - pos);
		if(IsSmall())
			SLen() = newlen;
		else
			LLen() = newlen;
		str[newlen] = 0;
		if(s)
			SVO_MEMCPY(str + pos, s, count);
		Dsyn();
		return str + pos;
	}
	char kind;
	char *p = Alloc(max(len >= int((int64)2 * INT_MAX / 3) ? INT_MAX : len + (len >> 1), newlen),
	                kind);
	if(pos > 0)
		SVO_MEMCPY(p, str, pos);
	if(pos < len)
		SVO_MEMCPY(p + pos + count, str + pos, len - pos);
	if(s)
		SVO_MEMCPY(p + pos, s, count);
	p[newlen] = 0;
	Free();
	ptr = p;
	LLen() = newlen;
	SLen() = 15;
	chr[KIND] = kind;
	Dsyn();
	return ptr + pos;
}

void String0::UnShare()
{
	if(IsSharedRef()) {
		int len = LLen();
		char kind;
		char *p = Alloc(len, kind);
		SVO_MEMCPY(p, ptr, len + 1);
		Free();
		chr[KIND] = kind;
		ptr = p;
	}
}

void String0::SetSLen(int l)
{
	SLen() = l;
	memset(chr + l, 0, 15 - l);
}

void String0::Remove(int pos, int count)
{
	ASSERT(pos >= 0 && count >= 0 && pos + count <= GetCount());
	UnShare();
	char *s = (char *)Begin();
	memmove(s + pos, s + pos + count, GetCount() - pos - count + 1);
	if(IsSmall())
		SetSLen(SLen() - count);
	else
		LLen() -= count;
	Dsyn();
}

void String0::Set(int pos, int chr)
{
	ASSERT(pos >= 0 && pos < GetCount());
	UnShare();
	Ptr()[pos] = chr;
}

void String0::Trim(int pos)
{
	ASSERT(pos >= 0 && pos <= GetCount());
	if(IsSmall()) {
		chr[pos] = 0;
		SetSLen(pos);
	}
	else {
		UnShare();
		ptr[pos] = 0;
		LLen() = pos;
	}
	Dsyn();
}

void String0::LCat(int c)
{
	if(IsSmall()) {
		qword *x = (qword *)MemoryAlloc32();
		x[0] = q[0];
		x[1] = q[1];
		LLen() = SLen();
		SLen() = 15;
		chr[KIND] = MEDIUM;
		qptr = x;
	}
	int l = LLen();
	if(IsRef() ? !IsShared() && l < (int)Ref()->alloc : l < 31) {
		ptr[l] = c;
		ptr[LLen() = l + 1] = 0;
	}
	else {
		char *s = Insert(l, 1, NULL);
		s[0] = c;
		s[1] = 0;
	}
}

void String0::Cat(const char *s, int len)
{
	if(IsSmall()) {
		if(SLen() + len < 14) {
			SVO_MEMCPY(chr + SLen(), s, len);
			SLen() += len;
			chr[(int)SLen()] = 0;
			Dsyn();
			return;
		}
	}
	else
		if((int)LLen() + len < LAlloc() && !IsSharedRef()) {
			SVO_MEMCPY(ptr + LLen(), s, len);
			LLen() += len;
			ptr[LLen()] = 0;
			Dsyn();
			return;
		}
	Insert(GetCount(), len, s);
}

void String0::Reserve(int r)
{
	int l = GetCount();
	Insert(GetCount(), r, NULL);
	Trim(l);
}

void String0::SetL(const char *s, int len)
{
	char *p = Alloc(len, chr[KIND]);
	memcpy(p, s, len);
	p[len] = 0;
	ptr = p;
	LLen() = len;
	SLen() = 15;
}

void String0::Set0(const char *s, int len)
{
	Zero();
	switch(len) {
	#define MOV(x) case x: chr[x - 1] = s[x - 1];
		MOV(14) MOV(13) MOV(12) MOV(11) MOV(10) MOV(9) MOV(8)
		MOV(7) MOV(6) MOV(5) MOV(4) MOV(3) MOV(2) MOV(1)
	case 0:
		SLen() = len;
		break;
	default:
		SetL(s, len);
	};
	Dsyn();
}

void String::AssignLen(const char *s, int slen)
{
	int  len = GetCount();
	char *str = (char *)Begin();
	if(s >= str && s <= str + len)
		*this = String(s, slen);
	else {
		String0::Free();
		String0::Set0(s, slen);
	}
}

String String::GetVoid()
{
	String s;
	s.ptr = (char *)(voidptr + 1);
	s.LLen() = 0;
	s.SLen() = 15;
	s.chr[KIND] = 50;
	return s;
}

bool String::IsVoid() const
{
	return IsRef() && ptr == (char *)(voidptr + 1);
}

WString String::ToWString() const
{
	return WString(Begin(), GetCount());
}

int String::GetCharCount() const
{
	return GetDefaultCharset() == CHARSET_UTF8 ?  utf8len(Begin(), GetCount()) : GetCount();
}

String::String(StringBuffer& b)
{
	if(b.pbegin == b.buffer) {
		String0::Set0(b.pbegin, (int)(uintptr_t)(b.pend - b.pbegin));
		return;
	}
	int l = b.GetLength();
	if(l <= 14) {
		Zero();
		memcpy(chr, b.pbegin, l);
		SLen() = l;
		b.Free();
	}
	else {
		ptr = b.pbegin;
		ptr[l] = 0;
		SLen() = 15;
		LLen() = l;
		chr[KIND] = min(b.GetAlloc(), 255);
		if(GetAlloc() > 4 * GetLength() / 3)
			Shrink();
	}
	b.Zero();

//	char h[100];
//	DLOG(sprintf(h, "String(StringBuffer) end %p (%p)", ptr, this));
	Dsyn();
//	DLOG(sprintf(h, "String(StringBuffer) end2 %p (%p)", ptr, this));
}

char *StringBuffer::Alloc(int count, int& alloc)
{
	if(count <= 31) {
		char *s = (char *)MemoryAlloc32();
		alloc = 31;
		return s;
	}
	else {
		size_t sz = sizeof(Rc) + count + 1;
		Rc *rc = (Rc *)MemoryAlloc(sz);
		alloc = rc->alloc = (int)min((size_t)INT_MAX, sz - sizeof(Rc) - 1);
		rc->refcount = 1;
		return (char *)(rc + 1);
	}
}

void StringBuffer::Free()
{
	if(pbegin == buffer)
		return;
	int all = (int)(limit - pbegin);
	if(all == 31)
		MemoryFree32(pbegin);
	if(all > 31)
		MemoryFree((Rc *)pbegin - 1);
}

void StringBuffer::Realloc(dword n, const char *cat, int l)
{
	int al;
	size_t ep = pend - pbegin;
	if(n > INT_MAX)
		n = INT_MAX;
	bool realloced = false;
	char *p;
	if((int)(limit - pbegin) > 800) {
		size_t sz = sizeof(Rc) + n + 1;
		Rc *rc = (Rc *)pbegin - 1;
		if(MemoryTryRealloc(rc, sz)) {
			realloced = true;
			al = rc->alloc = (int)min((size_t)INT_MAX, sz - sizeof(Rc) - 1);
			p = pbegin;
		}
	}
	if(!realloced) {
		p = Alloc(n, al);
		memcpy(p, pbegin, min((dword)GetLength(), n));
	}
	if(cat) {
		if(ep + l > INT_MAX)
			Panic("StringBuffer is too big (>2GB)!");
		memcpy(p + ep, cat, l);
		ep += l;
	}
	if(!realloced) {
		Free();
		pbegin = p;
	}
	pend = pbegin + ep;
	limit = pbegin + al;
}

void StringBuffer::Expand()
{
	Realloc(GetLength() * 3 / 2);
	if(pend == limit)
		Panic("StringBuffer is too big!");
}

void StringBuffer::SetLength(int l)
{
	if(l > GetAlloc())
		Realloc(l);
	pend = pbegin + l;
}

void StringBuffer::Shrink()
{
	int l = GetLength();
	if(l < GetAlloc() && l > 14)
		Realloc(l);
	pend = pbegin + l;
}

void StringBuffer::Cat(const char *s, int l)
{
	if(pend + l > limit)
		Realloc(max(GetLength(), l) + GetLength(), s, l);
	else {
		SVO_MEMCPY(pend, s, l);
		pend += l;
	}
}

void StringBuffer::Cat(int c, int l)
{
	if(pend + l > limit)
		Realloc(max(GetLength(), l) + GetLength(), NULL, l);
	memset(pend, c, l);
	pend += l;
}

void StringBuffer::Set(String& s)
{
	s.UnShare();
	int l = s.GetLength();
	if(s.GetAlloc() == 14) {
		pbegin = (char *)MemoryAlloc32();
		limit = pbegin + 31;
		memcpy(pbegin, s.Begin(), l);
		pend = pbegin + l;
	}
	else {
		pbegin = s.ptr;
		pend = pbegin + l;
		limit = pbegin + s.GetAlloc();
	}
	s.Zero();
}

String TrimLeft(const String& str)
{
	const char *s = str;
	if(!IsSpace(*s))
		return str;
	while(IsSpace((byte)*s)) s++;
	return String(s, str.End());
}

String TrimRight(const String& str)
{
	if(str.IsEmpty())
		return str;
	const char *s = str.Last();
	if(!IsSpace(*s))
		return str;
	while(s >= ~str && IsSpace((byte)*s)) s--;
	return String(~str, s + 1);
}

String TrimBoth(const String& str)
{
	return TrimLeft(TrimRight(str));
}

String TrimLeft(const char *sw, int len, const String& s)
{
	return s.StartsWith(sw, len) ? s.Mid(len) : s;
}

String TrimRight(const char *sw, int len, const String& s)
{
	return s.EndsWith(sw, len) ? s.Mid(0, s.GetCount() - len) : s;
}

struct StringICompare__
{
	int encoding;
	int operator()(char a, char b) const { return ToUpper(a, encoding) - ToUpper(b, encoding); }

	StringICompare__(int e) : encoding(e) {}
};

int CompareNoCase(const String& a, const String& b, byte encoding)
{
	if(encoding == CHARSET_DEFAULT) encoding = GetDefaultCharset();
	if(encoding == CHARSET_UTF8) return CompareNoCase(FromUtf8(a), FromUtf8(b));
#ifdef DEPRECATED
	return IterCompare(a.Begin(), a.End(), b.Begin(), b.End(), StringICompare__(encoding));
#else
	return CompareRanges(a, b, StringICompare__(encoding));
#endif
}

int CompareNoCase(const String& a, const char *b, byte encoding)
{
	if(encoding == CHARSET_DEFAULT) encoding = GetDefaultCharset();
	if(encoding == CHARSET_UTF8) return CompareNoCase(FromUtf8(a), FromUtf8(b, (int)strlen(b)));
#ifdef DEPRECATED
	return IterCompare(a.Begin(), a.End(), b, b + strlen(b), StringICompare__(encoding));
#else
	return CompareRanges(a, SubRange(b, b + strlen(b)), StringICompare__(encoding));
#endif
}

}
