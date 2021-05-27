#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>

#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#define IMAGECLASS ImgMoor
#define IMAGEFILE <BEMRosetta/main.iml>
#include <Draw/iml.h>

#include "main.h"


void MainMoor::Init() {
	CtrlLayout(*this);	

	edVessX <<= 0;
	edVessY <<= 0;

	lineTypes.Init(mooring);
	lineProperties.Init(mooring);
	lineConnections.Init(mooring);

	tab.Add(lineTypes.SizePos(), t_("Types"));
	tab.Add(lineConnections.SizePos(), t_("Connections"));
	tab.Add(lineProperties.SizePos(), t_("Properties"));
	tab.WhenSet = [&] {
		lineTypes.Save();
		lineProperties.Save();
		lineConnections.Save();
		if (tab.IsAt(lineProperties)) 
			lineProperties.LoadDrop();	
		
	};
	
	const String moorFiles = ".json";
	String moorFilesAst = clone(moorFiles);
	moorFilesAst.Replace(".", "*.");
	file.Type(Format("All supported mooring files (%s)", moorFiles), moorFilesAst);
	file.AllFilesType();
	file.BrowseRightWidth(40).UseOpenFolder(true).BrowseOpenFolderWidth(10);
	butLoad.Tip(t_("Loads mooring file")).WhenAction = [&] {OnLoad();};
	butSave.Tip(t_("Saves mooring file")).WhenAction = [&] {OnSave();};
	butUpdate.Tip(t_("Update mooring"))  .WhenAction = [&] {OnUpdate();};
}


bool MainMoor::OnLoad() {
	try {
		if (!mooring.Load(~file)) {
			Exclamation(Format("Problem loading %s file", DeQtf(~file)));
			return false;
		}
		lineTypes.Load();
		lineProperties.Load();
		lineConnections.Load();
	} catch (const Exc &e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	
	return true;
}

bool MainMoor::OnSave() {
	try {
		lineTypes.Save();
		lineProperties.Save();
		lineConnections.Save();
		if (!mooring.Save(~file)) {
			Exclamation(Format("Problem loading %s file", DeQtf(~file)));
			return false;
		}
	} catch (const Exc &e) {
		Exclamation(DeQtfLf(e));
		return false;
	}
	
	return true;
}

void MainMoor::OnUpdate() {
	Exclamation("Pendiente");	
}
