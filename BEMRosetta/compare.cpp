#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <GLCanvas/GLCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>

#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"



void CompareParameters::Init(ScatterDraw& scatter) {
	CtrlLayout(*this);
	
	pscatter = &scatter;
	swRelative <<= 0;
	swRelative.WhenAction = [=]() {Load();};
}

void CompareParameters::Init(Hydro::DataToShow dataToShow) {
	list.Reset();
	list.SetLineCy(EditField::GetStdHeight()).MultiSelect();
	list.WhenBar = [&](Bar &menu) {ArrayCtrlWhenBar(menu, list);};
	
	this->dataToShow = dataToShow;
	
	if (dataToShow == Hydro::DATA_A) {
		list.AddColumn(t_("Name"), 40);
		list.AddColumn(t_("RMS"), 10);
		list.AddColumn(t_("A∞"), 10);
	} else if (dataToShow == Hydro::DATA_B) {
		list.AddColumn(t_("Name"), 40);
		list.AddColumn(t_("RMS"), 10);
	} else
		list.AddColumn(t_("No parameters to compare"));
}

void CompareParameters::Load() {
	ScatterDraw &scatter = *pscatter;
	
	list.Clear();
	int col = 1;
	
	Vector<double> rms;
	for(int i = 0; i < scatter.GetCount(); i++) {
		String leg = scatter.GetLegend(i);
		if (leg.StartsWith("A∞") || leg.StartsWith("A0"))
			continue;
		list.Add(leg);
		DataSource &data = scatter.GetDataSource(i);
		rms << data.RMSY();
	}
	if (swRelative == 1) {
		for(int row = 1; row < rms.size(); row++) 
			rms[row] = rms[row]/rms[0];
		rms[0] = 1;
	}
	for(int row = 0; row < rms.size(); row++) 
		list.Set(row, col, FormatDoubleSize(rms[row], 10, false));	
    
    //col++;
    Vector<Vector<double>> data;
    for(int i = 0; i < scatter.GetCount(); i++) {
		String leg = scatter.GetLegend(i);
		if (leg.StartsWith("A∞") || leg.StartsWith("A0"))
			continue;
		Vector<double> &d = data.Add();
		
		Eigen::VectorXd x, y;
		scatter.GetDataSource(i).CopyXY(x, y);
		CleanNANDupXSort(x, y, x, y);
		Resample(x, y, x, y, Null);		// Set right srate
		Copy(y, d);
	}
    
    ///// XCorr
    
    
 
    if (dataToShow == Hydro::DATA_A) {
        col++;
        Vector<double> ainf;
        //int row = 0;
		for(int i = 0; i < scatter.GetCount(); i++) {
			String leg = scatter.GetLegend(i);
			if (leg.StartsWith("A∞")) 
				ainf << scatter.GetDataSource(i).y(int64(0));
		}
		if (swRelative == 1) {
			for(int row = 1; row < ainf.size(); row++) 
				ainf[row] = ainf[row]/ainf[0];
			ainf[0] = 1;
		}
		for(int row = 0; row < ainf.size(); row++) 
			list.Set(row, col, FormatDoubleSize(ainf[row], 10, false));	
    }

}

