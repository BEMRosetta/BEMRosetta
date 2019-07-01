//#BLITZ_PROHIBIT
// Need to prohibit BLITZ for this file as Xos.h defines some nasty macros like 'index'

#include <CtrlCore/CtrlCore.h>

#ifdef PLATFORM_POSIX

#include <X11/Xlib.h>
#include <X11/Xos.h>
#include <X11/Xfuncs.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>

using namespace Upp;

int GetKeyCodeX(int key) {
	_XDisplay *dpy = XOpenDisplay(NULL);
	if (!dpy)
		return Null;

	if (key > 0x00ff)
    	key = key | 0x01000000;
 	
	key = XKeysymToKeycode(dpy, key) + K_DELTA;
	
	XFlush(dpy);
 	XCloseDisplay(dpy);
 	return key;
}
	
#endif