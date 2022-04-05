#include <CtrlLib/CtrlLib.h>

using namespace Upp;

GUI_APP_MAIN
{
	/*FileSel fs;
	
	fs.ExecuteOK();
	*/
	bool large = false;
	bool exe = false;
	SHFILEINFO info = {0};
	SHGetFileInfo("w.lnk", FILE_ATTRIBUTE_NORMAL,
		               &info, sizeof(info),
		               SHGFI_ICON|
		               (large ? SHGFI_LARGEICON : SHGFI_SMALLICON)|
		               (exe ? 0 : SHGFI_USEFILEATTRIBUTES));
	
	Exclamation("It works!");
}
