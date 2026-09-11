// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2026, the BEMRosetta author and contributors
#include <Core/Core.h>
#include <Surface/Surface.h>
#include <STEM4U/Utility.h>
#include <Hdf5/hdf5.h>
#include "BEMRosetta.h"
#include "BEMRosetta_int.h"

#include "diffrac.h"

using namespace Upp;

		
static bool ReadBool(XmlNode &node, const String &tag, bool &ret) {
	const XmlNode &child = node(tag);
	if (child.GetCount() <= 0) 
		return false;
    String v = ToLower(Trim(child.GatherText()));
    ret = v == "1" || v == "true" || v == "yes";
    return true;
}

static bool ReadString(XmlNode &node, const String &tag, String &ret) {
	const XmlNode &child = node(tag);
	if (child.GetCount() <= 0) 
		return false;
    ret = Trim(child.GatherText());
    return true;
}

static bool ReadInt(XmlNode &node, const String &tag, int &ret) {
	const XmlNode &child = node(tag);
	if (child.GetCount() <= 0) 
		return false;
    ret = ScanInt(child.GatherText());
    return !IsNull(ret);
}

static bool ReadDouble(XmlNode &node, const String &tag, double &ret) {
	const XmlNode &child = node(tag);
	if (child.GetCount() <= 0) 
		return false;
    ret = ScanDouble(child.GatherText());
    return !IsNull(ret);
}

// 0,45,90,135,180,210; or 0(45)180,210
static bool ReadVectorDouble(XmlNode &node, const String &tag, UVector<double> &ret) {
	const XmlNode &child = node(tag);
	if (child.GetCount() <= 0) 
		return false;
	ret.Clear();
	String str = Trim(child.GatherText());
	UVector<String> t = Split(str, ',');
	for (const String &s : t) {
		int posleft = s.FindAfter("(");
		if (posleft > 0) {
			int posright = str.FindAfter(")", posleft);	
			if (posright <= 0)
				return false;
			double from = ScanDouble(str);
			if (!IsNum(from))
				return false;
			double delta = ScanDouble(str.Mid(posleft));
			if (!IsNum(delta))
				return false;
			double to = ScanDouble(str.Mid(posright));
			if (!IsNum(to))
				return false;
			if (from > to)
				return false;
			UVector<double> r;
			Arange(r, from, to, delta);
			ret.Append(r);
		} else {
			double d = ScanDouble(s);
			if (!IsNum(d))
				return false;
			ret << d;
		}
	}
    return true;
}

static bool ReadPointf(XmlNode &node, const String &tag, Pointf &ret) {
	const XmlNode &child = node(tag);
	if (child.GetCount() <= 0) 
		return false;
    UVector<String> t = Split(Trim(child.GatherText()), ',');
    if(t.GetCount() >= 2)
        ret = Pointf(ScanDouble(t[0]), ScanDouble(t[1]));
    return !IsNull(ret);
}

static bool ReadPoint3D(XmlNode &node, const String &tag, Point3D &ret) {
	const XmlNode &child = node(tag);
	if (child.GetCount() <= 0) 
		return false;
    UVector<String> t = Split(Trim(child.GatherText()), ',');
    if(t.GetCount() >= 3)
        ret = Point3D(ScanDouble(t[0]), ScanDouble(t[1]), ScanDouble(t[2]));
    return !IsNull(ret);
}

void DiffracData::LoadXML(const String &xml) {
	XmlNode xn = ParseXML(xml);
	XmlNode &sim = xn("sim");
	
	XmlNode &parVTKDIFx = sim("parVTKDIF");
	if (parVTKDIFx.GetCount())
		ReadBool(parVTKDIFx, "runProgram", this->parVTKDIF.runProgram);
		
	XmlNode &parINIDIFx = sim("parINIDIF");
	if (parINIDIFx.GetCount()) {
		ReadBool(parINIDIFx, "runProgram", this->parINIDIF.runProgram);
		ReadString(parINIDIFx, "dataBaseFn", this->parINIDIF.dataBaseFn);
		XmlNode &quay = parINIDIFx("quay");
		if (quay.GetCount()) {
			ReadBool(quay, "apply", this->parINIDIF.quay.apply);
			ReadDouble(quay, "yQuay", this->parINIDIF.quay.yQuay);
		}
		XmlNode &basinWall = parINIDIFx("basinWall");
		if (basinWall.GetCount()) {
			ReadBool(basinWall,   "apply",  this->parINIDIF.basinWall.apply);
			ReadPointf(basinWall, "origin", this->parINIDIF.basinWall.origin);
		}
		ReadDouble(parINIDIFx, "density", this->parINIDIF.density);
		XmlNode &springMatrix = parINIDIFx("springMatrix");
		if (springMatrix.GetCount()) 
			ReadBool(springMatrix, "fromGeometry", this->parINIDIF.springMatrixfromGeometry);
		XmlNode &dampLids = parINIDIFx("dampLids");
		for (int i = 0; i < dampLids.GetCount(); ++i) {
			XmlNode &dampLid = dampLids.At(i);	
			if (dampLid.GetCount()) {
				DiffracData::ParINIDIF::DampingLid &lid = this->parINIDIF.lids.Add();
				ReadPointf(dampLid, "origin", lid.origin);
				ReadDouble(dampLid, "orientation", lid.orientation);
				ReadDouble(dampLid, "dampingValue", lid.dampingValue);
				ReadDouble(dampLid, "length", lid.length);
				ReadDouble(dampLid, "width", lid.width);
				ReadInt(dampLid, "NrPanelsLength", lid.NrPanelsLength);
				ReadInt(dampLid, "NrPanelsWidth", lid.NrPanelsWidth);
				ReadBool(dampLid, "dampIncWave", lid.dampIncWave);
			}
		}
	}
	XmlNode &parDIFFRACx = sim("parDIFFRAC");
	if (parDIFFRACx.GetCount()) {
		double minFrequency, frequencyStep;
		int nFrequencies;
				
		ReadBool(parDIFFRACx, "runProgram", this->parDIFFRAC.runProgram);
		ReadVectorDouble(parDIFFRACx, "waveDir", this->parDIFFRAC.waveDir);
		ReadDouble(parDIFFRACx, "waterDepth", this->parDIFFRAC.waterDepth);
		XmlNode &current = parDIFFRACx("current");
		if (current.GetCount()) {
			ReadDouble(current, "speed", this->parDIFFRAC.current.speed);
			ReadDouble(current, "direction", this->parDIFFRAC.current.direction);
		}
		ReadString(parDIFFRACx, "irregFreqSuppression", this->parDIFFRAC.irregFreqSuppression);
		ReadDouble(parDIFFRACx, "irregFreqDamping", this->parDIFFRAC.irregFreqDamping);
		ReadVectorDouble(parDIFFRACx, "waveFreq", this->parDIFFRAC.waveFreq);
		ReadDouble(parDIFFRACx, "frequencyStep", frequencyStep);
		ReadInt(parDIFFRACx, "nFrequencies", nFrequencies);
		ReadDouble(parDIFFRACx, "minFrequency", minFrequency);
		ReadBool(parDIFFRACx, "exportKinematicsVTK", this->parDIFFRAC.exportKinematicsVTK);
		
		if (this->parDIFFRAC.waveFreq.IsEmpty()) {
			this->parDIFFRAC.waveFreq.SetCount(nFrequencies);
			for (int i = 0; i < nFrequencies; ++i)
				this->parDIFFRAC.waveFreq[i] = minFrequency + frequencyStep*i;
		}
	}
	XmlNode &parDBRESPx = sim("parDBRESP");
	if (parDBRESPx.GetCount()) {		
		ReadBool(parDBRESPx, "runProgram", this->parDBRESP.runProgram);
		XmlNode &BodyInputs = parDBRESPx("BodyInputs");
		for (int i = 0; i < BodyInputs.GetCount(); ++i) {
			XmlNode &BodyInput = BodyInputs.At(i);
			if (BodyInput.GetCount()) {	
				DiffracData::ParDBRESP::BodyInputDamping &b = this->parDBRESP.inputs.Add();
				ReadInt(BodyInput, "index", b.index);
				XmlNode &totalDampings = BodyInput("totalDampings");
				for (int j = 0; j < totalDampings.GetCount(); ++j) {
					XmlNode &totalDamping = totalDampings.At(j);
					if (totalDamping.GetCount()) {	
						DiffracData::ParDBRESP::BodyInputDamping::BodyTotalDamping &total = b.total.Add();
						ReadBool(totalDamping, "allowNegativeAddedDamping", total.allowNegativeAddedDamping);
						ReadString(totalDamping, "mode", total.mode);
						ReadString(totalDamping, "type", total.type);
						ReadDouble(totalDamping, "value", total.value);
					}
				}
			}
		}
		XmlNode &springMatrix = parDBRESPx("springMatrix");
		if (springMatrix.GetCount()) {	
			ReadString(springMatrix, "fileName", this->parDBRESP.springMatrix.fileName);	
			ReadBool(springMatrix, "isEarthFixed", this->parDBRESP.springMatrix.isEarthFixed);	
			ReadString(springMatrix, "unit", this->parDBRESP.springMatrix.unit);	
		}
		XmlNode &dampingMatrix = parDBRESPx("dampingMatrix");
		if (springMatrix.GetCount()) {	
			ReadString(dampingMatrix, "fileName", this->parDBRESP.dampingMatrix.fileName);	
			ReadBool(dampingMatrix, "isEarthFixed", this->parDBRESP.dampingMatrix.isEarthFixed);	
			ReadString(dampingMatrix, "unit", this->parDBRESP.dampingMatrix.unit);	
		}
	}
	XmlNode &parDRIFTPx = sim("parDRIFTP");
	if (parDRIFTPx.GetCount()) {		
		ReadBool(parDRIFTPx, "runProgram", this->parDRIFTP.runProgram);	
		ReadBool(parDRIFTPx, "exportContribution1", this->parDRIFTP.exportContribution[0]);
		ReadBool(parDRIFTPx, "exportContribution2", this->parDRIFTP.exportContribution[1]);
		ReadBool(parDRIFTPx, "exportContribution3", this->parDRIFTP.exportContribution[2]);
		ReadBool(parDRIFTPx, "exportContribution4", this->parDRIFTP.exportContribution[3]);
		ReadBool(parDRIFTPx, "exportContribution5", this->parDRIFTP.exportContribution[4]);
		ReadDouble(parDRIFTPx, "minFrequency", this->parDRIFTP.minFrequency);	
		ReadDouble(parDRIFTPx, "maxFrequency", this->parDRIFTP.maxFrequency);
		ReadInt(parDRIFTPx, "numberOfWavefrequencyDiagonals", this->parDRIFTP.numberOfWavefrequencyDiagonals);
		ReadBool(parDRIFTPx, "waveDirInteraction", this->parDRIFTP.waveDirInteraction);	
	}
	XmlNode &parSUMFREQUENCYWAVEFORCESx = sim("parSUMFREQUENCYWAVEFORCES");
	if (parSUMFREQUENCYWAVEFORCESx.GetCount()) {		
		ReadBool(parSUMFREQUENCYWAVEFORCESx, "runProgram", this->parSUMFREQUENCYWAVEFORCES.runProgram);	
		ReadBool(parSUMFREQUENCYWAVEFORCESx, "exportContribution1", this->parSUMFREQUENCYWAVEFORCES.exportContribution[0]);
		ReadBool(parSUMFREQUENCYWAVEFORCESx, "exportContribution2", this->parSUMFREQUENCYWAVEFORCES.exportContribution[1]);
		ReadBool(parSUMFREQUENCYWAVEFORCESx, "exportContribution3", this->parSUMFREQUENCYWAVEFORCES.exportContribution[2]);
		ReadBool(parSUMFREQUENCYWAVEFORCESx, "exportContribution4", this->parSUMFREQUENCYWAVEFORCES.exportContribution[3]);
		ReadBool(parSUMFREQUENCYWAVEFORCESx, "exportContribution5", this->parSUMFREQUENCYWAVEFORCES.exportContribution[4]);
		ReadDouble(parSUMFREQUENCYWAVEFORCESx, "minFrequency", this->parSUMFREQUENCYWAVEFORCES.minFrequency);	
		ReadDouble(parSUMFREQUENCYWAVEFORCESx, "maxFrequency", this->parSUMFREQUENCYWAVEFORCES.maxFrequency);
		ReadInt(parSUMFREQUENCYWAVEFORCESx, "numberOfWavefrequencyDiagonals", this->parSUMFREQUENCYWAVEFORCES.numberOfWavefrequencyDiagonals);
		ReadBool(parSUMFREQUENCYWAVEFORCESx, "waveDirInteraction", this->parSUMFREQUENCYWAVEFORCES.waveDirInteraction);	
	}
	XmlNode &parEXPORTx = sim("parEXPORT");
	if (parEXPORTx.GetCount()) {		
		ReadBool(parEXPORTx, "runProgram", this->parEXPORT.runProgram);		
		XmlNode &hydFile = parEXPORTx("hydFile");
		if (hydFile.GetCount()) {	
			ReadBool(hydFile, "export", this->parEXPORT.hyd.exportOn);	
			ReadBool(hydFile, "exportQTFContribution1", this->parEXPORT.hyd.exportQTF[0]);
			ReadBool(hydFile, "exportQTFContribution2", this->parEXPORT.hyd.exportQTF[1]);
			ReadBool(hydFile, "exportQTFContribution3", this->parEXPORT.hyd.exportQTF[2]);
			ReadBool(hydFile, "exportQTFContribution4", this->parEXPORT.hyd.exportQTF[3]);
			ReadBool(hydFile, "exportQTFContribution5", this->parEXPORT.hyd.exportQTF[4]);
			ReadInt(hydFile, "numberOfWavefrequencyDiagonals", this->parEXPORT.hyd.numberOfWavefrequencyDiagonals);
		}
		XmlNode &MonitorRelativeWaveHeights = parEXPORTx("MonitorRelativeWaveHeights");
		for (int i = 0; i < MonitorRelativeWaveHeights.GetCount(); ++i) {
			XmlNode &MonitorRelativeWaveHeight = MonitorRelativeWaveHeights.At(i);
			if (MonitorRelativeWaveHeight.GetCount()) {	
				 DiffracData::ParEXPORT::RelWaveHeight &rel = this->parEXPORT.relHeights.Add();
				 ReadString(MonitorRelativeWaveHeight, "name", rel.name);
				 ReadBool(MonitorRelativeWaveHeight, "includeIncidentWave", rel.incIncident);
				 ReadBool(MonitorRelativeWaveHeight, "includeDiffractedWave", rel.incDiffracted);
				 ReadBool(MonitorRelativeWaveHeight, "includeRadiatedWave", rel.incRadiated);
				 ReadBool(MonitorRelativeWaveHeight, "includeMotions", rel.incMotions);
				 XmlNode &relWaveHeightBodyInput = MonitorRelativeWaveHeight("relWaveHeightBodyInput");
				 if (relWaveHeightBodyInput.GetCount()) {
				 	ReadInt(relWaveHeightBodyInput, "index", rel.index);
				 	XmlNode &referencePoint2D = relWaveHeightBodyInput("referencePoint2D");
					if (referencePoint2D.GetCount()) 
						ReadPointf(referencePoint2D, "coordinate", rel.referencePoint);
				 }
			}
		}
		XmlNode &CGNS = parEXPORTx("CGNS");
		if (CGNS.GetCount()) {
			ReadBool(CGNS, "export", this->parEXPORT.cgns.exportOn);
			ReadDouble(CGNS, "waveFreq", this->parEXPORT.cgns.waveFreq);
			XmlNode &MonitorFlowDatas = CGNS("MonitorFlowDatas");
			for (int i = 0; i < MonitorFlowDatas.GetCount(); ++i) {			
				XmlNode &MonitorFlowData = MonitorFlowDatas.At(i);
				if (MonitorFlowData.GetCount()) {
					DiffracData::ParEXPORT::CGNS::MonitorFlowData &data = this->parEXPORT.cgns.monitorFlowData.Add();
					XmlNode &grids = MonitorFlowData("grids");
					for (int j = 0; j < grids.GetCount(); ++j) {
						XmlNode &grid = grids.At(j);
						if (grid.GetCount()) {
							DiffracData::ParEXPORT::CGNS::MonitorFlowData::Grid &gr = data.grids.Add();	
							ReadInt(grid, "NrOfPointsX", gr.Nx);
							ReadInt(grid, "NrOfPointsY", gr.Ny);
							ReadInt(grid, "NrOfPointsZ", gr.Nz);
							ReadPoint3D(grid, "origin", gr.origin);
							ReadDouble(grid, "spacing_x", gr.sx);
							ReadDouble(grid, "spacing_y", gr.sy);
							ReadDouble(grid, "spacing_z", gr.sz);
							ReadDouble(grid, "orientation", gr.orientation);
						}
					}
				}
			}
			ReadBool(CGNS, "movingBodies", this->parEXPORT.cgns.movingBodies);
			ReadBool(CGNS, "movingFreeSurface", this->parEXPORT.cgns.movingFreeSurface);
			ReadInt(CGNS, "numberOfTimeSteps", this->parEXPORT.cgns.numberOfTimeSteps);
			ReadDouble(CGNS, "waveAmplificationFactor", this->parEXPORT.cgns.waveAmplificationFactor);
			ReadDouble(CGNS, "waveDir", this->parEXPORT.cgns.waveDir);
		}
	}
	XmlNode &bodiesx = sim("bodies");
	for (int i = 0; i < bodiesx.GetCount(); ++i) {
		XmlNode &body = bodiesx.At(i);	
		if (body.GetCount()) {
			DiffracData::Body &b = this->bodies.Add();
			XmlNode &hstat = body("hstat");
			if (hstat.GetCount()) {
				ReadDouble(hstat, "lengthBetweenPerp", b.hstat.lengthBetweenPerp);
				ReadDouble(hstat, "draft", b.hstat.draft);	
			}
			ReadInt(body, "index", b.index);	
			XmlNode &mesh = body("mesh");
			if (mesh.GetCount())
				ReadString(mesh, "meshFn", b.meshFn);
			ReadString(body, "name", b.name);
			XmlNode &position = body("position");
			if (position.GetCount()) {
				ReadPoint3D(position, "translation", b.translation);
				ReadDouble(position, "rotation", b.rotation);	
			}			
			XmlNode &massElements = body("massElements");
			for (int im = 0; im < massElements.GetCount(); ++im) {
				XmlNode &massElement = massElements.At(im);	
				if (massElement.GetCount()) {
					DiffracData::Body::MassElement &m = b.massElements.Add();				
					ReadPoint3D(massElement, "COGwrtKeel", m.COGwrtKeel);
					ReadDouble(massElement, "mass", m.mass);		
					ReadDouble(massElement, "rollRadiusGyr", m.rollRadiusGyr);
					ReadDouble(massElement, "pitchRadiusGyr", m.pitchRadiusGyr);
					ReadDouble(massElement, "yawRadiusGyr", m.yawRadiusGyr);
				}
			}
		}
	}
	ReadInt(sim, "projectNumber", this->projectNumber);
	ReadInt(sim, "nProcs", this->nProcs);
}

static String ToText(bool d) {
	return d ? "true" : "false";
}

static String ToText(int d) {
	return FormatInt(d);
}

static String ToText(double d) {
	//if (abs(d) < 0.000000001)
	//	return "0.0";
	String ret = F("%f", d);
	if (ret.Find(".") >= 0) {
		while (*(ret.Last()) == '0')
			ret.Remove(ret.GetCount()-1);
	}
	if (*(ret.Last()) == '.')
		ret.Remove(ret.GetCount()-1);
	if (ret == "-0")
		ret = "0";
	return ret;
}

static String ToText(const Pointf &d) {
	return F("%s,%s", ToText(d.x), ToText(d.y));
}

static String ToText(const Point3D &d) {
	return F("%s,%s,%s", ToText(d.x), ToText(d.y), ToText(d.z));
}

bool AnalyzeStep(const UVector<double>& data, double tolerance, double &minVal, double &step, double &maxVal) {
    int n = data.GetCount();
	
	ASSERT(n > 0);
    if (n < 2) {
    	minVal = maxVal = data[0];
    	step = 0;
        return true;
    }

    minVal = data[0];
    maxVal = data.Top();
    step = (maxVal - minVal)/(n-1);

    for (int i = 0; i < n-1; i++) {
        double currentStep = data[i+1] - data[i];
        if (abs(currentStep - step) > tolerance)
            return false;
    }
    return true;
}

template <class Range>
String ToText(const Range& a) {
	double minVal, step, maxVal;
	if (AnalyzeStep(a, 1e-5, minVal, step, maxVal)) {
		if (step == 0)
			return F("%s", ToText(minVal));
		else
			return F("%s(%s)%s", ToText(minVal), ToText(step), ToText(maxVal));
	} else {
		String ret;
		for (int i = 0; i < a.size(); i++) {
			if (i > 0)
				ret << ",";
			ret << a[i]; 
		}
		return ret;
	}
}

String DiffracData::SaveXML() {
	String out;
	out << "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n"
		   "<sim xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"I:\\Applications\\Diffrac\\releases\\v3.5.6\\XML\\XML_input\\DIFFRAC_input.xsd\">\n";
    
	out << "  <parVTKDIF>\n"
		<< F("    <runProgram>%s</runProgram>\n", ToText(parVTKDIF.runProgram))
		<< "  </parVTKDIF>\n"
	;
	out << "  <parINIDIF>\n"
		<< F("    <runProgram>%s</runProgram>\n", ToText(parINIDIF.runProgram))
		<< F("    <dataBaseFn>%s</dataBaseFn>\n", parINIDIF.dataBaseFn)
		<< "    <quay>\n"
		<< F("      <apply>%s</apply>\n", ToText(parINIDIF.quay.apply))
		<< F("      <yQuay>%s</yQuay>\n", ToText(parINIDIF.quay.yQuay))
		<< "    </quay>\n"
		<< "    <basinWall>\n"
		<< F("      <apply>%s</apply>\n", ToText(parINIDIF.basinWall.apply))
		<< "      <origin>\n"
		<< F("            %s\n", ToText(parINIDIF.basinWall.origin))
		<< "      </origin>\n"
		<< "    </basinWall>\n"
		<< "    <density>\n"
		<< F("          %s\n", ToText(parINIDIF.density))
		<< "    </density>\n"
		<< "    <springMatrix>\n"
		<< F("      <fromGeometry>%s</fromGeometry>\n", ToText(parINIDIF.springMatrixfromGeometry))
		<< "      <springMatrixFn></springMatrixFn>\n"
		<< "    </springMatrix>\n"
    ;
    if (!parINIDIF.lids.IsEmpty())
		out	<< "    <dampLids>\n";
	for (int i = 0; i < parINIDIF.lids.size(); ++i) {
		out << "      <dampLid>\n";
		out << F("        <origin>%s</origin>\n", ToText(parINIDIF.lids[i].origin));
		out << F("        <orientation>%s</orientation>\n", ToText(parINIDIF.lids[i].orientation));
		out << F("        <dampingValue>%s</dampingValue>\n", ToText(parINIDIF.lids[i].dampingValue));
		out << F("        <length>%s</length>\n", ToText(parINIDIF.lids[i].length));
		out << F("        <width>%s</width>\n", ToText(parINIDIF.lids[i].width));
		out << F("        <NrPanelsLength>%s</NrPanelsLength>\n", ToText(parINIDIF.lids[i].NrPanelsLength));
		out << F("        <NrPanelsWidth>%s</NrPanelsWidth>\n", ToText(parINIDIF.lids[i].NrPanelsWidth));
		out << F("        <dampIncWave>%s</dampIncWave>\n", ToText(parINIDIF.lids[i].dampIncWave));
		out << "      </dampLid>\n";
	}
	if (!parINIDIF.lids.IsEmpty())
		out << "    </dampLids>\n";
	out << "  </parINIDIF>\n";	
	
	out << "  <parDIFFRAC>\n"
		<< F("    <runProgram>%s</runProgram>\n", ToText(parDIFFRAC.runProgram))
		<< F("    <waveDir>%s</waveDir>\n", ToText(parDIFFRAC.waveDir))
		<< F("    <waterDepth>%s</waterDepth>\n", ToText(parDIFFRAC.waterDepth))
		<< "    <current>\n"
		<< "      <speed>\n"
		<< F("            %s\n", ToText(parDIFFRAC.current.speed))
		<< "      </speed>\n"
		<< "      <direction>\n"
		<< F("            %s\n", ToText(parDIFFRAC.current.direction))
		<< "      </direction>\n"
		<< "    </current>\n"
		<< F("    <irregFreqSuppression>%s</irregFreqSuppression>\n", parDIFFRAC.irregFreqSuppression)
		<< F("    <irregFreqDamping>%s</irregFreqDamping>\n", ToText(parDIFFRAC.irregFreqDamping))
		<< F("    <waveFreq>%s</waveFreq>\n", ToText(parDIFFRAC.waveFreq))
		<< F("    <exportKinematicsVTK>%s</exportKinematicsVTK>\n", ToText(parDIFFRAC.exportKinematicsVTK))
		<< "  </parDIFFRAC>\n"
	;

	out << "  <parDBRESP>\n"
		<< F("    <runProgram>%s</runProgram>\n", ToText(parDBRESP.runProgram));
	
	if (!parDBRESP.inputs.IsEmpty())
		out << "    <BodyInputs>\n";
	for (int i = 0; i < parDBRESP.inputs.size(); ++i) {
		out << "      <BodyInput>\n"
			<< F("        <index>%s</index>\n", ToText(parDBRESP.inputs[i].index))
			<< "        <totalDampings>\n"
		;
		for (int j = 0; j < parDBRESP.inputs[i].total.size(); ++j) {
			out	<< "          <totalDamping>\n"
				<< F("            <allowNegativeAddedDamping>%s</allowNegativeAddedDamping>\n", ToText(parDBRESP.inputs[i].total[j].allowNegativeAddedDamping))
				<< F("            <mode>%s</mode>\n", parDBRESP.inputs[i].total[j].mode)
				<< F("            <type>%s</type>\n", parDBRESP.inputs[i].total[j].type)
				<< F("            <value>%s</value>\n", ToText(parDBRESP.inputs[i].total[j].value))
				<< "          </totalDamping>\n"
			;
		}
		out	<< "        </totalDampings>\n"
			<< "      </BodyInput>\n"
		;
	}
	if (!parDBRESP.inputs.IsEmpty())
		out	<< "    </BodyInputs>\n";
	if (!parDBRESP.springMatrix.fileName.IsEmpty()) {
		out	<< "    <springMatrix>\n"
			<< F("      <fileName>%s</fileName>\n", parDBRESP.springMatrix.fileName)
			<< F("      <isEarthFixed>%s</isEarthFixed>\n", ToText(parDBRESP.springMatrix.isEarthFixed))
			<< F("      <unit>%s</unit>\n", parDBRESP.springMatrix.unit)
			<< "    </springMatrix>\n";
	}
	if (!parDBRESP.dampingMatrix.fileName.IsEmpty()) {
		out	<< "    <dampingMatrix>\n"
			<< F("      <fileName>%s</fileName>\n", parDBRESP.dampingMatrix.fileName)
			<< F("      <isEarthFixed>%s</isEarthFixed>\n", ToText(parDBRESP.dampingMatrix.isEarthFixed))
			<< F("      <unit>%s</unit>\n", parDBRESP.dampingMatrix.unit)
			<< "    </dampingMatrix>\n";
	}
	out	<< "  </parDBRESP>\n";

	if (parDRIFTP.runProgram) {
		out << "  <parDRIFTP>\n"
			<< F("    <runProgram>%s</runProgram>\n", ToText(parDRIFTP.runProgram))
			<< F("    <exportContribution1>%s</exportContribution1>\n", ToText(parDRIFTP.exportContribution[0]))
			<< F("    <exportContribution2>%s</exportContribution2>\n", ToText(parDRIFTP.exportContribution[1]))
			<< F("    <exportContribution3>%s</exportContribution3>\n", ToText(parDRIFTP.exportContribution[2]))
			<< F("    <exportContribution4>%s</exportContribution4>\n", ToText(parDRIFTP.exportContribution[3]))
			<< F("    <exportContribution5>%s</exportContribution5>\n", ToText(parDRIFTP.exportContribution[4]))
			<< F("    <minFrequency>%s</minFrequency>\n", ToText(parDRIFTP.minFrequency))
			<< F("    <maxFrequency>%s</maxFrequency>\n", ToText(parDRIFTP.maxFrequency))
			<< F("    <numberOfWavefrequencyDiagonals>%s</numberOfWavefrequencyDiagonals>\n", ToText(parDRIFTP.numberOfWavefrequencyDiagonals))
			<< F("    <waveDirInteraction>%s</waveDirInteraction>\n", ToText(parDRIFTP.waveDirInteraction))
			<< "  </parDRIFTP>\n"
		;
	}
	if (parSUMFREQUENCYWAVEFORCES.runProgram) {
		out << "  <parSUMFREQUENCYWAVEFORCES>\n"
			<< F("    <runProgram>%s</runProgram>\n", ToText(parSUMFREQUENCYWAVEFORCES.runProgram))
			<< F("    <exportContribution1>%s</exportContribution1>\n", ToText(parSUMFREQUENCYWAVEFORCES.exportContribution[0]))
			<< F("    <exportContribution2>%s</exportContribution2>\n", ToText(parSUMFREQUENCYWAVEFORCES.exportContribution[1]))
			<< F("    <exportContribution3>%s</exportContribution3>\n", ToText(parSUMFREQUENCYWAVEFORCES.exportContribution[2]))
			<< F("    <exportContribution4>%s</exportContribution4>\n", ToText(parSUMFREQUENCYWAVEFORCES.exportContribution[3]))
			<< F("    <exportContribution5>%s</exportContribution5>\n", ToText(parSUMFREQUENCYWAVEFORCES.exportContribution[4]))
			<< F("    <minFrequency>%s</minFrequency>\n", ToText(parSUMFREQUENCYWAVEFORCES.minFrequency))
			<< F("    <maxFrequency>%s</maxFrequency>\n", ToText(parSUMFREQUENCYWAVEFORCES.maxFrequency))
			<< F("    <numberOfWavefrequencyDiagonals>%s</numberOfWavefrequencyDiagonals>\n", ToText(parSUMFREQUENCYWAVEFORCES.numberOfWavefrequencyDiagonals))
			<< F("    <waveDirInteraction>%s</waveDirInteraction>\n", ToText(parSUMFREQUENCYWAVEFORCES.waveDirInteraction))
			<< "  </parSUMFREQUENCYWAVEFORCES>\n"
		;
	}
	out << "  <parEXPORT>\n"
		<< F("    <runProgram>%s</runProgram>\n", ToText(parEXPORT.runProgram));
	out	<< "    <hydFile>\n"
		<< F("      <export>%s</export>\n", ToText(parEXPORT.hyd.exportOn));
	if (parSUMFREQUENCYWAVEFORCES.runProgram) 
		out	<< F("      <exportQTFContribution1>%s</exportQTFContribution1>\n", ToText(parEXPORT.hyd.exportQTF[0]))
			<< F("      <exportQTFContribution2>%s</exportQTFContribution2>\n", ToText(parEXPORT.hyd.exportQTF[1]))
			<< F("      <exportQTFContribution3>%s</exportQTFContribution3>\n", ToText(parEXPORT.hyd.exportQTF[2]))
			<< F("      <exportQTFContribution4>%s</exportQTFContribution4>\n", ToText(parEXPORT.hyd.exportQTF[3]))
			<< F("      <exportQTFContribution5>%s</exportQTFContribution5>\n", ToText(parEXPORT.hyd.exportQTF[4]));
	out	<< F("      <numberOfWavefrequencyDiagonals>%s</numberOfWavefrequencyDiagonals>\n", ToText(parEXPORT.hyd.numberOfWavefrequencyDiagonals))
		<< "    </hydFile>\n"
	;
	if (!parEXPORT.relHeights.IsEmpty())
		out << "    <MonitorRelativeWaveHeights>\n";
	for (int i = 0; i < parEXPORT.relHeights.size(); ++i) {
		out << "      <MonitorRelativeWaveHeight>\n"
			<< F("        <name>%s</name>\n", parEXPORT.relHeights[i].name)
			<< F("        <includeIncidentWave>%s</includeIncidentWave>\n", ToText(parEXPORT.relHeights[i].incIncident))
			<< F("        <includeDiffractedWave>%s</includeDiffractedWave>\n", ToText(parEXPORT.relHeights[i].incDiffracted))
			<< F("        <includeRadiatedWave>%s</includeRadiatedWave>\n", ToText(parEXPORT.relHeights[i].incRadiated))
			<< F("        <includeMotions>%s</includeMotions>\n", ToText(parEXPORT.relHeights[i].incMotions))
			<< "        <relWaveHeightBodyInput>\n"
			<< F("          <index>%s</index>\n", ToText(parEXPORT.relHeights[i].index))
			<< "          <referencePoint2D>\n"
			<< F("            <coordinate>%s</coordinate>\n", ToText(parEXPORT.relHeights[i].referencePoint))
			<< "          </referencePoint2D>\n"
			<< "        </relWaveHeightBodyInput>\n"
			<< "      </MonitorRelativeWaveHeight>\n"
		;		
	}
	if (!parEXPORT.relHeights.IsEmpty())
		out << "    </MonitorRelativeWaveHeights>\n";
	
	out << "    <CGNS>\n"
		<< F("      <export>false</export>\n", ToText(parEXPORT.cgns.exportOn))
		<< F("      <waveFreq>0.8</waveFreq>\n", ToText(parEXPORT.cgns.exportOn));
	
	if (!parEXPORT.cgns.monitorFlowData.IsEmpty())
		out << "      <MonitorFlowDatas>\n";
	for (int i = 0; i < parEXPORT.cgns.monitorFlowData.size(); ++i) {
		out	<< "        <MonitorFlowData>\n"
			<< "          <grids>\n";
		for (int j = 0; j < parEXPORT.cgns.monitorFlowData[i].grids.size(); ++j) {
			out	<< "            <grid>\n"
				<< F("              <NrOfPointsX>201</NrOfPointsX>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].Nx))
				<< F("              <NrOfPointsY>101</NrOfPointsY>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].Ny))
				<< F("              <origin>0.0,0.0,0.0</origin>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].origin))
				<< F("              <spacing_x>5</spacing_x>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].sx))
				<< F("              <spacing_y>5</spacing_y>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].sy))
				<< F("              <orientation>0.0</orientation>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].orientation))
				<< F("              <NrOfPointsZ>1</NrOfPointsZ>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].Nz))
				<< F("              <spacing_z>1e-6</spacing_z>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].sz))
				<< "            </grid>\n"
			;
		}
		out	<< "          </grids>\n"
		    << "        </MonitorFlowData>\n";
	}
	if (!parEXPORT.cgns.monitorFlowData.IsEmpty())
		out	<< "      </MonitorFlowDatas>\n";
	
	out	<< F("      <movingBodies>%s</movingBodies>\n", ToText(parEXPORT.cgns.movingBodies))
		<< F("      <movingFreeSurface>%s</movingFreeSurface>\n", ToText(parEXPORT.cgns.movingFreeSurface))
		<< F("      <numberOfTimeSteps>%s</numberOfTimeSteps>\n", ToText(parEXPORT.cgns.numberOfTimeSteps))
		<< F("      <waveAmplificationFactor>%s</waveAmplificationFactor>\n", ToText(parEXPORT.cgns.waveAmplificationFactor))
		<< F("      <waveDir>%s</waveDir>\n", ToText(parEXPORT.cgns.waveDir))
		<< "    </CGNS>\n"
    ;
    out	<< "    <H5M>\n"
      	<< "      <export>true</export>\n"
    	<< "    </H5M>\n"
    	<< "    <monitorAmassDampAtCoG>false</monitorAmassDampAtCoG>\n";
	out << "  </parEXPORT>\n";
	
	out	<< "  <bodies>\n";
	for(int i = 0; i < bodies.size(); ++i) {
		out << "    <body>\n"
			<< "      <hstat>\n"
			<< F("        <lengthBetweenPerp>%s</lengthBetweenPerp>\n", ToText(bodies[i].hstat.lengthBetweenPerp))
			<< F("        <draft>%s</draft>\n", ToText(bodies[i].hstat.draft))
			<< "      </hstat>\n"
			<< F("      <index>%s</index>\n", ToText(bodies[i].index))
			<< "      <mesh>\n"
			<< F("        <meshFn>%s</meshFn>\n", bodies[i].meshFn)
			<< "      </mesh>\n"
			<< F("      <name>%s</name>\n", bodies[i].name)
			<< "      <position>\n"
			<< F("        <translation>%s</translation>\n", ToText(bodies[i].translation))
			<< F("        <rotation>0</rotation>\n", ToText(bodies[i].rotation))
			<< "      </position>\n"
			<< "      <massElements>\n"
		;
		for (int j = 0; j < bodies[i].massElements.size(); ++j) {
			out << "        <massElement>\n"
				<< F("          <COGwrtKeel>%s</COGwrtKeel>\n", ToText(bodies[i].massElements[j].COGwrtKeel))
				<< F("          <mass>%s</mass>\n", ToText(bodies[i].massElements[j].mass))
				<< F("          <rollRadiusGyr>%s</rollRadiusGyr>\n", ToText(bodies[i].massElements[j].rollRadiusGyr))
				<< F("          <pitchRadiusGyr>%s</pitchRadiusGyr>\n", ToText(bodies[i].massElements[j].pitchRadiusGyr))
				<< F("          <yawRadiusGyr>%s</yawRadiusGyr>\n", ToText(bodies[i].massElements[j].yawRadiusGyr))
				<< "        </massElement>\n"
				<< "      </massElements>\n"
			;
		}
		out << "    </body>\n";
	}
    out	<< "  </bodies>\n";	
	
	out << F("  <projectNumber>98800</projectNumber>\n", ToText(projectNumber))
		<< F("  <nProcs>4</nProcs>\n", ToText(nProcs)) 
		<< "</sim>"
	;
	return out;
}

String Diffrac::LoadCase(String file) {
	dt.file = file;
	dt.name = GetFileTitle(file);
	dt.dimen = true;
	dt.len = 1;
	dt.solver = Hydro::DIFFRAC;
	dt.Nb = Null;
	dt.x_w = dt.y_w = 0;
	
	String folder = GetFileFolder(file);
	
	try {
		BEM::Print("\n\n" + F(t_("Loading '%s'"), file));
		
	    String xml = LoadFile(file);
	    if(IsNull(xml))
			return F(t_("Cannot read file %s"), file);

		DiffracData data;
		data.LoadXML(xml);
		
		dt.name = data.parINIDIF.dataBaseFn;
		
		dt.rho = data.parINIDIF.density*1000;
		dt.h = data.parDIFFRAC.waterDepth;
		dt.w = pick(data.parDIFFRAC.waveFreq);
		dt.Nf = dt.w.size();
		dt.head = pick(data.parDIFFRAC.waveDir);
		dt.Nh = dt.head.size();
		
		dt.Nb = data.bodies.size(); 	
		dt.msh.SetCount(dt.Nb);
		for (int ib = 0; ib < dt.Nb; ++ib) {
			DiffracData::Body &b = data.bodies[ib];		
			Body &msh = dt.msh[ib];
			
			Body::Load(msh, AFX(folder, b.meshFn), dt.rho, Bem().g, Null, Null, false);
			if (msh.dt.mesh.IsEmpty())
				return F(t_("Impossible to load mesh file %s"), b.meshFn);
			
			msh.dt.name = b.name;
			
			DiffracData::Body::MassElement &m = b.massElements[0];
			msh.dt.M = MatrixXd::Zero(6, 6);
			msh.dt.M(0, 0) = msh.dt.M(1, 1) = msh.dt.M(2, 2) = m.mass;
			msh.dt.M(3, 3) = m.mass*sqr(m.rollRadiusGyr);
			msh.dt.M(4, 4) = m.mass*sqr(m.pitchRadiusGyr);
			msh.dt.M(5, 5) = m.mass*sqr(m.yawRadiusGyr);
			msh.dt.M *= 1000.;
			
			double minZ = msh.dt.mesh.env.minZ;
			msh.dt.cg = pick(m.COGwrtKeel);
			msh.dt.cg.z += minZ;
			msh.dt.c0 = clone(msh.dt.cg);
			
			msh.dt.cg += b.translation;
			msh.dt.c0 += b.translation;
			msh.dt.mesh.Translate(b.translation);
			
			Surface lid;
			if (lid.GetPanels(msh.dt.mesh, false, true, false, Null, Null)) {
				dt.lids.SetCount(dt.Nb);
				dt.lids[ib].dt.mesh = pick(lid);
				dt.lids[ib].AfterLoad(Bem().rho, Bem().g, false, true);
			}
		}
		if (!data.parDBRESP.springMatrix.fileName.IsEmpty()) {
			MatrixXd d(6*dt.Nb, 6*dt.Nb);
			FileInLine in(AFX(folder, data.parDBRESP.springMatrix.fileName));
			if (!in.IsOpen())
				return F(t_("Imposssible to open stiffness matrix file %s"), data.parDBRESP.springMatrix.fileName);
			
			LineParserWamit f(in);
			f.IsSeparator = IsTabSpace;
	
			for (int row = 0; row < 6*dt.Nb; ++row) {
				f.GetLine();
				for (int col = 0; col < 6*dt.Nb; ++col)
					d(row, col) = f.GetDouble(col);
			}
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.Cadd = d.block(6*ib, 6*ib, 6, 6);
		}
		if (!data.parDBRESP.dampingMatrix.fileName.IsEmpty()) {
			MatrixXd d(6*dt.Nb, 6*dt.Nb);
			FileInLine in(AFX(folder, data.parDBRESP.dampingMatrix.fileName));
			if (!in.IsOpen())
				return F(t_("Imposssible to open damping matrix file %s"), data.parDBRESP.dampingMatrix.fileName);
			
			LineParserWamit f(in);
			f.IsSeparator = IsTabSpace;
	
			for (int row = 0; row < 6*dt.Nb; ++row) {
				f.GetLine();
				for (int col = 0; col < 6*dt.Nb; ++col)
					d(row, col) = f.GetDouble(col);
			}
			for (int ib = 0; ib < dt.Nb; ++ib)
				dt.msh[ib].dt.Dlin = d.block(6*ib, 6*ib, 6, 6);
		}
	} catch (Exc e) {
		return e;
	}
	
	return String();
}

void Diffrac::SaveCase(String folder, int numThreads, bool withPotentials, bool x0z, bool y0z, 
					   const UVector<bool> &listDOF, bool irregular, int qtfType) const {
	bool onlyMeanDrift = qtfType > 10;
	qtfType %= 10;
		       
	if (!DirectoryCreateX(folder))
		throw Exc(Format(t_("Problem creating '%s' folder"), folder));

	DiffracData data;
	
	data.nProcs = numThreads;
	
	data.parINIDIF.dataBaseFn = dt.name;
	
	data.parINIDIF.density = dt.rho/1000;
	data.parDIFFRAC.waterDepth = dt.h;
	data.parDIFFRAC.waveFreq = clone(dt.w);
	data.parDIFFRAC.waveDir = clone(dt.head);
	if (irregular) {
		data.parDIFFRAC.irregFreqSuppression = "RIGID LID";
		data.parDIFFRAC.irregFreqDamping = 0.03;
	} else
		data.parDIFFRAC.irregFreqSuppression = "NONE";
	
	bool thereIsStiffness = false;
	for (int ib = 0; ib < dt.Nb; ++ib) {
		if (dt.msh[ib].dt.Cadd.size() == 36 || dt.msh[ib].dt.Cmoor.size() == 36) {
			thereIsStiffness = true;
			break;
		}
	}
	if (thereIsStiffness) {
		FileOut fs(AFX(folder, "springMatrix.txt"));
		if (!fs.IsOpen())
			throw Exc("Cannot save springMatrix.txt");
		for (int row = 0; row < 6*dt.Nb; ++row) {
			int ib = row/6;
			for (int col = 0; col < 6*dt.Nb; ++col) {
				if (row >= 6*ib && row < 6*(ib+1) && col >= 6*ib && col < 6*(ib+1)) {
					double val = 0;
					if (dt.msh[ib].dt.Cadd.size() == 36)
						val = dt.msh[ib].dt.Cadd(row - ib*6, col - ib*6);
					if (dt.msh[ib].dt.Cmoor.size() == 36)
						val +=dt.msh[ib].dt.Cmoor(row - ib*6, col - ib*6);
					fs << val << " ";
				} else
					fs << "0 ";
			}
			fs << "\n";
		}
		data.parDBRESP.springMatrix.fileName = "springMatrix.txt";
	}
	bool thereIsDamping = false;
	for (int ib = 0; ib < dt.Nb; ++ib) {
		if (dt.msh[ib].dt.Dlin.size() == 36) {
			thereIsDamping = true;
			break;
		}
	}
	if (thereIsDamping) {
		FileOut fs(AFX(folder, "dampingMatrix.txt"));
		if (!fs.IsOpen())
			throw Exc("Cannot save dampingMatrix.txt");
		for (int row = 0; row < 6*dt.Nb; ++row) {
			int ib = row/6;
			for (int col = 0; col < 6*dt.Nb; ++col) {
				if (dt.msh[ib].dt.Dlin.size() == 36 && row >= 6*ib && row < 6*(ib+1) && col >= 6*ib && col < 6*(ib+1))
					fs << dt.msh[ib].dt.Dlin(row - ib*6, col - ib*6) << " ";
				else
					fs << "0 ";
			}
			fs << "\n";
		}
		data.parDBRESP.dampingMatrix.fileName = "dampingMatrix.txt";
	}
	
	if (qtfType <= 0) {
		data.parDRIFTP.runProgram = data.parSUMFREQUENCYWAVEFORCES.runProgram = false;
		data.parEXPORT.hyd.exportOn = false;
	} else {
		data.parEXPORT.hyd.exportOn = true;
		
		data.parDRIFTP.runProgram = true;
		for (int i = 0; i < 5; ++i)
			data.parDRIFTP.exportContribution[i] = true;
		
		data.parDRIFTP.minFrequency = First(dt.w);
		data.parDRIFTP.maxFrequency = Last(dt.w);
		data.parDRIFTP.numberOfWavefrequencyDiagonals = dt.Nf;
		data.parDRIFTP.waveDirInteraction = true;
		
		if (!onlyMeanDrift) {
			data.parSUMFREQUENCYWAVEFORCES.runProgram = true;
			for (int i = 0; i < 5; ++i)
				data.parSUMFREQUENCYWAVEFORCES.exportContribution[i] = true;
			
			data.parSUMFREQUENCYWAVEFORCES.minFrequency = First(dt.w);
			data.parSUMFREQUENCYWAVEFORCES.maxFrequency = Last(dt.w);
			data.parSUMFREQUENCYWAVEFORCES.numberOfWavefrequencyDiagonals = dt.Nf;
			data.parSUMFREQUENCYWAVEFORCES.waveDirInteraction = true;
		} else {
			data.parSUMFREQUENCYWAVEFORCES.runProgram = false;
			data.parEXPORT.hyd.numberOfWavefrequencyDiagonals = 0;
		}
		for (int i = 0; i < 5; ++i)
			data.parEXPORT.hyd.exportQTF[i] = !onlyMeanDrift;
	}

	data.bodies.SetCount(dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) {
		DiffracData::Body &b = data.bodies[ib];
		const Body &msh = dt.msh[ib];
		
		b.index = ib+1;
		b.name = msh.dt.name;
		
		DiffracData::Body::MassElement &m = b.massElements.Add();
		MatrixXd M = msh.dt.M/1000.;
		
		m.mass = M(0, 0);		
		Surface::TranslateInertia66(M, msh.dt.cg, msh.dt.c0, msh.dt.cg);
		m.rollRadiusGyr = sqrt(M(3, 3)/m.mass);
		m.pitchRadiusGyr = sqrt(M(4, 4)/m.mass);
		m.yawRadiusGyr = sqrt(M(5, 5)/m.mass);
		
		m.COGwrtKeel = clone(msh.dt.cg);
		m.COGwrtKeel.z -= msh.dt.mesh.env.minZ;
		
		b.hstat.draft = abs(msh.dt.mesh.env.minZ);
		b.hstat.lengthBetweenPerp = max(msh.dt.mesh.env.maxX, msh.dt.mesh.env.maxY) - min(msh.dt.mesh.env.minX, msh.dt.mesh.env.minY);
	}
	
	UArray<Body> bodies(dt.Nb);
	for (int ib = 0; ib < dt.Nb; ++ib) {
		DiffracData::Body &b = data.bodies[ib];
		
		b.meshFn = F("Body_%d.vtk", ib+1);
		String dest = AFX(folder, b.meshFn);
		
		bodies[ib] = clone(dt.msh[ib]);
		bool isLid = irregular && dt.lids.size() > ib && !dt.lids[ib].dt.mesh.panels.IsEmpty();
		if (isLid) {
			Surface lid = clone(dt.lids[ib].dt.mesh);
			if (x0z) {			// Assures that even with symmetry, no triangle is in the lid
				Surface nlid;			
				nlid.CutY(lid);
				nlid.TrianglesToFalseQuads();
				nlid.DeployYSymmetry();
				lid = pick(nlid);				
			} else
				lid.TrianglesToFalseQuads();
			bodies[ib].Append(lid, rho_ndim(), g_ndim());
		}
	}
	for (int ib = 0; ib < dt.Nb; ++ib) {
		DiffracData::Body &b = data.bodies[ib];
		
		bool save = true;
		if (ib > 0) {
			Value3D dist;
			for (int iib = 0; iib < ib; ++iib) {
				if (bodies[ib].dt.mesh.CompareDistance(bodies[iib].dt.mesh, 0.0001, dist)) {
					save = false;
					b.translation = -dist;
					b.massElements[0].COGwrtKeel += dist;
					b.meshFn = data.bodies[iib].meshFn;
					break;
				}
			}
		}
		if (save) {
			String dest = AFX(folder, F("Body_%d.vtk", ib+1));
			Body::SaveAs(bodies[ib], dest, Body::VTK_ASCII_4, Body::ALL, rho_ndim(), g_ndim(), y0z, x0z);
		}
	}
	SaveFile(AFX(folder, "diffrac.xml"), data.SaveXML());				
}

String Diffrac::Load(String file, double) {
	dt.file = file;
	dt.name = GetFileTitle(file);
	dt.dimen = true;
	dt.len = 1;
	dt.solver = Hydro::DIFFRAC_H5;
	dt.Nb = Null;
	dt.x_w = dt.y_w = 0;
	
	try {
		BEM::Print("\n\n" + F(t_("Loading '%s'"), file));

		BEM::Print("\n- " + F(t_("H5m file")));
		
		Load_H5();
		
		if (dt.Nb == 0)
			return t_("No data found");
	
	} catch (Exc e) {
		return e;
	}
	
	return String();
}

void Diffrac::Load_H5() {
	String fileName = ForceExtSafer(dt.file, ".h5m");
	
	Hdf5File hfile;
	hfile.Open(fileName, H5F_ACC_RDONLY);

	for (dt.Nb = 0; hfile.ExistGroup(F("Body_Nr_%d", dt.Nb+1)); dt.Nb++) 
		;	
		
	if (dt.Nb == 0)
		return;
	
	dt.msh.SetCount(dt.Nb);
	
	{	
		hfile.ChangeGroup("Fluid");
		
		dt.h = hfile.GetFloat("water depth");
		if (dt.h == 0)
			dt.h = -1;
		
		if (hfile.ExistDataset("wave direction")) {
			UVector<float> head;
			hfile.GetDouble("wave direction", head);
			if (head.IsEmpty())
				throw Exc("No headings found");
			dt.head.SetCount(head.size());
			for (int i = 0; i < head.size(); ++i)
				dt.head[i] = head[i];
			dt.Nh = head.size();
		}
		if (hfile.ExistDataset("wave frequency")) {
			UVector<float> w;
			hfile.GetDouble("wave frequency", w);
			if (w.IsEmpty())
				throw Exc("No frequencies found");
			dt.w.SetCount(w.size());
			for (int i = 0; i < w.size(); ++i)
				dt.w[i] = w[i];
			dt.Nf = w.size();
		}
		hfile.UpGroup();
	}
	
	for (int ib = 0; ib < dt.Nb; ++ib) {
		hfile.ChangeGroup(F("Body_Nr_%d", ib+1));
		
		MultiDimMatrixRowMajor<float> data;
		for (int idof1 = 0; idof1 < 6; ++idof1)
			for (int idof2 = 0; idof2 < 6*dt.Nb; ++idof2) {
				String svar = F("A_m") + FormatInt(ib*6 + idof1 + 1) + F("m") + FormatInt(idof2 + 1);
				if (hfile.ExistDataset(svar)) {
					if (!IsLoadedA())
						Initialize_AB(dt.A);
					hfile.GetDouble(svar, data);
					for (int ifr = 0; ifr < dt.Nf; ++ifr)
						dt.A[ib*6 + idof1][idof2](ifr) = data.begin()[ifr];
				}
			}
		for (int idof1 = 0; idof1 < 6; ++idof1)
			for (int idof2 = 0; idof2 < 6*dt.Nb; ++idof2) {
				String svar = F("B_m") + FormatInt(ib*6 + idof1 + 1) + F("m") + FormatInt(idof2 + 1);
				if (hfile.ExistDataset(svar)) {
					if (!IsLoadedB())
						Initialize_AB(dt.B);
					hfile.GetDouble(svar, data);
					for (int ifr = 0; ifr < dt.Nf; ++ifr)
						dt.B[ib*6 + idof1][idof2](ifr) = data.begin()[ifr];
				}
			}		
		
		MultiDimMatrixRowMajor<std::complex<float>> f;
		for (int idof = 0; idof < 6; ++idof) {
			String svar = "F_exc_m" + FormatInt(idof+1);
			if (hfile.ExistDataset(svar)) {
				if (!IsLoadedFex())
					Initialize_Forces(dt.ex);
				hfile.GetComplex(svar, f);
				for (int ifr = 0; ifr < dt.Nf; ++ifr)
					for (int ih = 0; ih < dt.Nh; ++ih)
						dt.ex[ib][ih](ifr, idof) = f(0, 0, 0, ih, ifr);
			}
		}
		for (int idof = 0; idof < 6; ++idof) {
			String svar = "F_dif_m" + FormatInt(idof+1);
			if (hfile.ExistDataset(svar)) {
				if (!IsLoadedFsc())
					Initialize_Forces(dt.sc);
				hfile.GetComplex(svar, f);
				for (int ifr = 0; ifr < dt.Nf; ++ifr)
					for (int ih = 0; ih < dt.Nh; ++ih)
						dt.sc[ib][ih](ifr, idof) = f(0, 0, 0, ih, ifr);
			}
		}
		for (int idof = 0; idof < 6; ++idof) {
			String svar = "F_inc_m" + FormatInt(idof+1);
			if (hfile.ExistDataset(svar)) {
				if (!IsLoadedFfk())
					Initialize_Forces(dt.fk);
				hfile.GetComplex(svar, f);
				for (int ifr = 0; ifr < dt.Nf; ++ifr)
					for (int ih = 0; ih < dt.Nh; ++ih)
						dt.fk[ib][ih](ifr, idof) = f(0, 0, 0, ih, ifr);
			}
		}
		for (int idof = 0; idof < 6; ++idof) {
			String svar = F("%s Motion", BEM::strDOFtext[idof]);
			if (hfile.ExistDataset(svar)) {
				if (!IsLoadedRAO())
					Initialize_Forces(dt.rao);
				hfile.GetComplex(svar, f);
				for (int ifr = 0; ifr < dt.Nf; ++ifr)
					for (int ih = 0; ih < dt.Nh; ++ih)
						dt.rao[ib][ih](ifr, idof) = f(0, 0, 0, ih, ifr);
			}
		}
		for (int idof = 0; idof < 6; ++idof) {
			String svar = "QTF_m" + FormatInt(idof+1);
			if (hfile.ExistDataset(svar)) {
				if (!IsLoadedQTF(false))
					Hydro::Initialize_QTF(dt.qtfdif, dt.Nb, dt.Nh, dt.Nf);
				hfile.GetComplex(svar, f);
				for (int ifr1 = 0; ifr1 < dt.Nf; ++ifr1) {
					for (int ifr2 = 0; ifr2 < dt.Nf; ++ifr2)
						for (int ih = 0; ih < dt.Nh; ++ih)
							dt.qtfdif[ib][ih][idof](ifr1, ifr2) = f(0, 0, 0, ih, ifr1, ifr2);			
				}
			}
		}
		if (hfile.ExistGroup("Metadata")) {
		 	hfile.ChangeGroup("Metadata");
			if (hfile.ExistGroup("Particulars")) {
			 	hfile.ChangeGroup("Particulars");	
				if (hfile.ExistDataset("COB")) {
					UVector<float> COB;
					hfile.GetDouble("COB", COB);
					dt.msh[ib].dt.cb.Set(COB);
				}
				if (hfile.ExistDataset("COG")) {
					UVector<float> COG;
					hfile.GetDouble("COG", COG);
					dt.msh[ib].dt.cg.Set(COG);
				}				
				
				hfile.UpGroup();
			}
			hfile.UpGroup();
		}
		hfile.UpGroup();
	}
	if (IsLoadedQTF(false)) {
		dt.qtftype = 9;
		dt.qw = Get_w();
		dt.qhead.resize(dt.Nh);
		for (int ih = 0; ih < dt.Nh; ++ih)
			dt.qhead[ih] = std::complex<double>(dt.head[ih], dt.head[ih]);
	}
}
