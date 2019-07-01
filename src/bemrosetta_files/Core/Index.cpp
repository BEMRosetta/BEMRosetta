#include <Core/Core.h>

namespace Upp {

int IndexCommon::empty[1] = { -1 };

IndexCommon::IndexCommon()
{
	hash = NULL;
	map = empty;
	mask = 0;
	unlinked = -1;
}

void IndexCommon::Pick(IndexCommon& b)
{
	Free();

	hash = b.hash;
	map = b.map;
	mask = b.mask;
	unlinked = b.unlinked;
	
	b.hash = NULL;
	b.map = empty;
	b.mask = 0;
	b.unlinked = -1;
}

void IndexCommon::Copy(const IndexCommon& b, int count)
{
	memcpy(hash, b.hash, sizeof(Hash) * count);
	mask = b.mask;
	unlinked = b.unlinked;

	FreeMap();
	map = (int *)MemoryAlloc((mask + 1) * sizeof(int));
	memcpy(map, b.map, (mask + 1) * sizeof(int));
}

void IndexCommon::Swap(IndexCommon& b)
{
	UPP::Swap(hash, b.hash);
	UPP::Swap(map, b.map);
	UPP::Swap(mask, b.mask);
	UPP::Swap(unlinked, b.unlinked);
}

IndexCommon::~IndexCommon()
{
	Free();
}

void IndexCommon::FreeMap()
{
	if(map != empty)
		MemoryFree(map);
}

void IndexCommon::Free()
{
	if(hash)
		MemoryFree(hash);
	FreeMap();
}

void IndexCommon::Remap(int count)
{
	Fill(map, map + mask + 1, -1);
	for(int i = 0; i < count; i++) // todo: unlinked
		if(hash[i].hash)
			Link(i, hash[i].hash);
}

void IndexCommon::Reindex(int count)
{
	FreeMap();
	map = (int *)MemoryAlloc((mask + 1) * sizeof(int));
	Remap(count);
}

void IndexCommon::Clear()
{
	Free();
	hash = NULL;
	map = empty;
	mask = 0;
	unlinked = -1;
}

void IndexCommon::GrowMap(int count)
{
	mask = (mask << 1) | 3;
	Reindex(count);
}

Vector<int> IndexCommon::GetUnlinked() const
{
	Vector<int> r;
	int i = unlinked;
	if(i >= 0) {
		do {
			i = hash[i].prev;
			r.Add(i);
		}
		while(i != unlinked);
	}
	return r;
}

void IndexCommon::AdjustMap(int count, int alloc)
{
	if(alloc == 0) {
		FreeMap();
		map = empty;
		mask = 0;
		return;
	}
	dword msk = 0;
	while(msk < (dword)alloc)
		msk = (msk << 1) | 3;
	if(msk != mask) {
		mask = msk;
		Reindex(count);
	}
}

void IndexCommon::MakeMap(int count)
{
	mask = 0;
	AdjustMap(count, count);
}

void IndexCommon::Trim(int n, int count)
{
	if(n == 0) {
		int n = (int)(mask + 1);
		for(int i = 0; i < n; i++)
			map[i] = -1;
		unlinked = -1;
		return;
	}
	
	for(int i = n; i < count; i++) { // remove items in trimmed area from buckets / unlinked
		Hash& hh = hash[i];
		if(hh.hash)
			Del(map[hh.hash & mask], hh, i);
		else
			Del(unlinked, hh, i);
	}
}

void IndexCommon::Sweep(int n)
{
	int ti = 0;
	for(int i = 0; i < n; i++)
		if(hash[i].hash)
			hash[ti++].hash = hash[i].hash;
	Remap(ti);
	unlinked = -1;
}

#ifdef CPU_UNALIGNED

NOUBSAN // CPU supports unaligned memory access
unsigned memhash(const void *ptr, size_t count)
{
	unsigned hash = 1234567890U;

	const unsigned *ds = (unsigned *)ptr;
	const unsigned *de = ds + (count >> 2);
	while(ds < de)
		hash = ((hash << 5) - hash) ^ *ds++;

	const byte *s = (byte *)ds;
	const byte *e = s + (count & 3);
	while(s < e)
		hash = ((hash << 5) - hash) ^ *s++;

	return hash;
}

#else

unsigned memhash(const void *ptr, size_t count)
{
	unsigned hash = 1234567890U;

	const byte *s = (byte *)ptr;
	const byte *e = s + count;
	while(s < e)
		hash = ((hash << 5) - hash) ^ *s++;

	return hash;
}

#endif

unsigned GetHashValue0(const double& d)
{
	return memhash(&d, sizeof(double));
}

}