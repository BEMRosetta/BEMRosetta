// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2021, the BEMRosetta author and contributors
#ifndef _BEM_Rosetta_GUI_arrange_h_
#define _BEM_Rosetta_GUI_arrange_h_

#include <CtrlLib/CtrlLib.h>

using namespace Upp;

#define LAYOUTFILE <BEMRosetta/premain.lay>
#include <CtrlCore/lay.h>

class ArrangeDOF : public WithArrange<StaticRect> {
public:
	void Init(Hydro &hydro);

private:		
	bool DnDInsert(int line, PasteClip& d);
	bool selecting = false;
};



#endif
