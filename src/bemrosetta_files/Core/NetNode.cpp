#include "Core.h"

#ifdef PLATFORM_WIN32

namespace Upp {

NetNode::NetNode()
{
	memset(&net, 0, sizeof(NETRESOURCE));
}

NetNode& NetNode::operator=(const NetNode& s)
{
	net = s.net;
	local = s.local;
	remote = s.remote;
	comment = s.comment;
	provider = s.provider;
	name = s.name;
	path = s.path;
	SetPtrs();
	return *this;
}

String DosInitCaps(const char *name)
{
	for(const char *s = name; *s; s++)
		if(IsLetter(*s) && IsLower(*s))
			return name;
	return InitCaps(name);
}

void NetNode::SetPtr(String& s, char *& ptr)
{
	if(ptr) ptr = (char *)~s;
}

void NetNode::SetPtrs()
{
	SetPtr(local, net.lpLocalName);
	SetPtr(remote, net.lpRemoteName);
	SetPtr(comment, net.lpComment);
	SetPtr(provider, net.lpProvider);
}

void NetNode::Serialize(Stream& s)
{
	s % net.dwScope % net.dwType % net.dwDisplayType % net.dwUsage;
	s % local % remote % comment % provider % name % path;
	SetPtrs();
}

Array<NetNode> NetNode::Enum() const
{
	HANDLE hEnum;
	if(::WNetOpenEnum(RESOURCE_GLOBALNET, RESOURCETYPE_DISK, 0, (NETRESOURCE *)&net, &hEnum))
		return Array<NetNode>();
	return Enum0(hEnum);
}

Array<NetNode> NetNode::EnumRoot()
{
	HANDLE hEnum;
	if(::WNetOpenEnum(RESOURCE_CONTEXT, RESOURCETYPE_DISK, 0, NULL, &hEnum))
		return Array<NetNode>();
	return Enum0(hEnum);
}

Array<NetNode> NetNode::EnumRemembered()
{
	HANDLE hEnum;
	if(::WNetOpenEnum(RESOURCE_REMEMBERED, RESOURCETYPE_DISK, 0, NULL, &hEnum))
		return Array<NetNode>();
	return Enum0(hEnum);
}

Array<NetNode> NetNode::Enum0(HANDLE hEnum)
{
	Array<NetNode> r;
	DWORD cEntries = (DWORD)-1, cbBuffer = 0x4000;
	Buffer<NETRESOURCE> lpnr(cbBuffer);
	while(::WNetEnumResource(hEnum, &cEntries, lpnr, &cbBuffer) == 0) {
		for(int i = 0; i < (int)cEntries; i++) {
			NETRESOURCE& sn = lpnr[i];
			const char *s = sn.lpRemoteName;
			NetNode& nn = r.Add();
			NETRESOURCE& n = nn.net;
			n = sn;
			nn.local = n.lpLocalName;
			nn.remote = n.lpRemoteName;
			nn.comment = n.lpComment;
			nn.provider = n.lpProvider;
			nn.SetPtrs();
			if(s) {
				if(s[0] == '\\' && s[1] == '\\')
					nn.name = FromSystemCharset(DosInitCaps(GetFileName(s)));
				else
					nn.name = FromSystemCharset(s);
			}
			if(n.lpComment && *n.lpComment) {
				if(nn.name.GetCount())
					nn.name = String().Cat() << FromSystemCharset(n.lpComment)
					                         << " (" << nn.name << ")";
				else
					nn.name = FromSystemCharset(n.lpComment);
			}
			if(!(n.dwUsage & RESOURCEUSAGE_CONTAINER))
				nn.path = FromSystemCharset(n.lpRemoteName);
		}
	}
	::WNetCloseEnum(hEnum);
	return r;
}

int NetNode::GetDisplayType() const
{
	switch(net.dwDisplayType) {
	case RESOURCEDISPLAYTYPE_NETWORK:
		return NETWORK;
	case RESOURCEDISPLAYTYPE_DOMAIN:
	case RESOURCEDISPLAYTYPE_GROUP:
	case RESOURCEDISPLAYTYPE_TREE:
	case /*RESOURCEDISPLAYTYPE_NDSCONTAINER*/0x0000000B:
		return GROUP;
	case RESOURCEDISPLAYTYPE_SERVER:
		return SERVER;
	case RESOURCEDISPLAYTYPE_SHARE:
	case RESOURCEDISPLAYTYPE_SHAREADMIN:
	case RESOURCEDISPLAYTYPE_DIRECTORY:
		return SHARE;
	}
	return UNKNOWN;
}

}

#endif
