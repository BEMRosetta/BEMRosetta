#ifndef flagGUI

#include "BEMRosetta.h"

CONSOLE_APP_MAIN {
	const Vector<String>& command = CommandLine();
	
	ConsoleMain(command, false);
}

#endif