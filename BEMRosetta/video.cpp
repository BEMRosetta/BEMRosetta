// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2025, the BEMRosetta author and contributors
#include <CtrlLib/CtrlLib.h>
#include <Controls4U/Controls4U.h>
#include <ScatterCtrl/ScatterCtrl.h>
#include <SurfaceCanvas/SurfaceCanvas.h>
#include <RasterPlayer/RasterPlayer.h>
#include <TabBar/TabBar.h>
#include <DropGrid/DropGrid.h>

#include <BEMRosetta_cl/BEMRosetta.h>

using namespace Upp;

#include "main.h"

void VideoCtrl::Init(Function <int(UVector<int> &ids)> _GetBodyId, Function <void(int id, const UVector<int> &ids, const Point3D &pos, const Point3D &angle, const Point3D &c0, bool full, bool saveBitmap)> _Action) {
	CtrlLayout(*this);
	
	butPlay <<= THISBACK(OnPlay);
	butPlay.Disable();
	butRecord <<= THISBACK(OnRecord);
	butRecordClear <<= THISBACK(OnRecordClear);
	butRecordClear.Disable();
	butLoad <<= THISBACK(OnLoad);
	butSave <<= THISBACK(OnSave);
	butSave.Disable();
	
	opSaveBitmap.Disable();
	
	deltaT.Enable();
	
	if (IsNull(deltaT)) 
		deltaT <<= 0.5;
	
	editTime <<= "0";
	numFrames <<= 0;
	
	progress.SetTotal(1000);
	
	GetBodyId = _GetBodyId;
	Action = _Action;
}

void VideoCtrl::OnPlay() {
	if ((meshId = GetBodyId(ids)) < 0) 
		return;
	
	playing = !playing;
	
	if (!playing) {
		butPlay.SetLabel("Play");	
		progress.Set(1000);
		
		KillTimeCallback(1);
		meshId = -1;
	} else {	
		if (video->size() == 0)
			return;
		if (!video->HasTime())
			video->SetDeltaTime(~deltaT);
		
		butPlay.SetLabel("Stop");	
		dT = ~deltaT;
		
		SetTimeCallback(-200, THISBACK(TimerFun), 1);
		time.Reset();
		TimerFun();
	} 
	
	butRecord.Enable(!playing);
}

void VideoCtrl::OnRecord() {
	recording = !recording;
	
	butPlay.Enable(!recording);
    butRecordClear.Enable();
    
	if (!recording) {
		butRecord.SetLabel(t_("Record"));
		butSave.Enable();
	} else {
		video.Clear();
		video.Create<BasicVideoSequence>();
		butRecord.SetLabel(t_("..."));
	}
}

void VideoCtrl::OnRecordClear() {
	if (!IsBasicOpened())
		return;
	static_cast<BasicVideoSequence&>(*video).Clear();	
	butRecordClear.Disable();
	butPlay.Disable();
	butSave.Disable();
	editTime <<= "0";
	progress.Set(0);
	numFrames <<= 0;
}

void VideoCtrl::TimerFun() {
	if (!playing)
		return;

	VideoSequence::VideoRecord record;
	int res;
	while((res = video->GetRecord(record, time.Seconds())) > 0) {
		Action(meshId, ids, record.pos, record.angle, record.c0, res < 0, ~opSaveBitmap);
		editTime <<= SecondsToString(record.time, 1, false, false, false, false, true); 
		progress.Set(int(1000.*((record.time+dT)/(dT*video->size()))));
		if (res == 2) {
			OnPlay();
			return;
		}
	}		
}

VideoCtrl::VideoType VideoCtrl::GetVideoType(String name) {
	String ext = ToLower(GetFileExt(name));
	
	if (ext == ".json")
		return JSON;
	else if (ext == ".out" || ext == ".outb")
		return OpenFAST;
	else if (ext == ".csv")
		return CSV;
	return UNKNOWN;
}

void VideoCtrl::OnLoad() {
	try {
		int type = GetVideoType(~editFile);
		
		if (type == JSON) {
			video.Clear();
			video.Create<BasicVideoSequence>();
			if (!video->Load(~editFile))
				throw Exc(t_("Problem loading file"));
		} else
			throw Exc(t_("File type is not already handled"));
		
		numFrames <<= video->size();
		butPlay.Enable();
		butSave.Enable();
		butRecordClear.Enable();
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

void VideoCtrl::OnSave() {
	try {
		if (!IsBasicOpened())
			throw Exc(t_("This sequence cannot be saved"));;
		
		if (!video->Save(~editFile))
			throw Exc(t_("Problem saving file"));
		
	} catch (Exc e) {
		BEM::PrintError(DeQtfLf(e));
	}
}

VideoSequence::~VideoSequence() {}
