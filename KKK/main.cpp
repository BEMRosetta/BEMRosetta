#include <Core/Core.h>
#include <Surface/Surface.h>
#include <STEM4U/Utility.h>

#include "DiffRacModel.h"

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

static bool ReadVectorDouble(XmlNode &node, const String &tag, UVector<double> &ret) {
	const XmlNode &child = node(tag);
	if (child.GetCount() <= 0) 
		return false;
	ret.Clear();
	String str = Trim(child.GatherText());
	int posleft = str.FindAfter("(");
	if (posleft > 0) {
		int posright = str.FindAfter(")", posleft);	
		if (posright <= 0)
			return false;
		double from = ScanDouble(str);
		double delta = ScanDouble(str.Mid(posleft));
		double to = ScanDouble(str.Mid(posright));
		if (from > to)
			return false;
		Arange(ret, from, to, delta);
	} else {
		UVector<String> t = Split(str, ',');
		for (const String s : t) {
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
	
	XmlNode &parVTKDIF = sim("parVTKDIF");
	if (parVTKDIF.GetCount())
		ReadBool(parVTKDIF, "runProgram", this->parVTKDIF.runProgram);
		
	XmlNode &parINIDIF = sim("parINIDIF");
	if (parINIDIF.GetCount()) {
		ReadBool(parINIDIF, "runProgram", this->parINIDIF.runProgram);
		ReadString(parINIDIF, "dataBaseFn", this->parINIDIF.dataBaseFn);
		XmlNode &quay = parINIDIF("quay");
		if (quay.GetCount()) {
			ReadBool(quay, "apply", this->parINIDIF.quay.apply);
			ReadDouble(quay, "yQuay", this->parINIDIF.quay.yQuay);
		}
		XmlNode &basinWall = parINIDIF("basinWall");
		if (basinWall.GetCount()) {
			ReadBool(basinWall,   "apply",  this->parINIDIF.basinWall.apply);
			ReadPointf(basinWall, "origin", this->parINIDIF.basinWall.origin);
		}
		ReadDouble(parINIDIF, "density", this->parINIDIF.density);
		XmlNode &springMatrix = parINIDIF("springMatrix");
		if (springMatrix.GetCount()) 
			ReadBool(springMatrix, "fromGeometry", this->parINIDIF.springMatrixfromGeometry);
		XmlNode &dampLids = parINIDIF("dampLids");
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
	XmlNode &parDIFFRAC = sim("parDIFFRAC");
	if (parDIFFRAC.GetCount()) {		
		ReadBool(parDIFFRAC, "runProgram", this->parDIFFRAC.runProgram);
		ReadVectorDouble(parDIFFRAC, "waveDir", this->parDIFFRAC.waveDir);
		ReadDouble(parDIFFRAC, "waterDepth", this->parDIFFRAC.waterDepth);
		XmlNode &current = parDIFFRAC("current");
		if (current.GetCount()) {
			ReadDouble(current, "speed", this->parDIFFRAC.current.speed);
			ReadDouble(current, "direction", this->parDIFFRAC.current.direction);
		}
		ReadString(parDIFFRAC, "irregFreqSuppression", this->parDIFFRAC.irregFreqSuppression);
		ReadDouble(parDIFFRAC, "irregFreqDamping", this->parDIFFRAC.irregFreqDamping);
		ReadDouble(parDIFFRAC, "frequencyStep", this->parDIFFRAC.frequencyStep);
		ReadInt(parDIFFRAC, "nFrequencies", this->parDIFFRAC.nFrequencies);
		ReadDouble(parDIFFRAC, "minFrequency", this->parDIFFRAC.minFrequency);
		ReadBool(parDIFFRAC, "exportKinematicsVTK", this->parDIFFRAC.exportKinematicsVTK);
	}
	XmlNode &parDBRESP = sim("parDBRESP");
	if (parDBRESP.GetCount()) {		
		ReadBool(parDBRESP, "runProgram", this->parDBRESP.runProgram);
		XmlNode &BodyInputs = parDBRESP("BodyInputs");
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
	}
	XmlNode &parDRIFTP = sim("parDRIFTP");
	if (parDRIFTP.GetCount()) {		
		ReadBool(parDRIFTP, "runProgram", this->parDRIFTP.runProgram);	
		ReadBool(parDRIFTP, "exportContribution1", this->parDRIFTP.exportContribution[0]);
		ReadBool(parDRIFTP, "exportContribution2", this->parDRIFTP.exportContribution[1]);
		ReadBool(parDRIFTP, "exportContribution3", this->parDRIFTP.exportContribution[2]);
		ReadBool(parDRIFTP, "exportContribution4", this->parDRIFTP.exportContribution[3]);
		ReadBool(parDRIFTP, "exportContribution5", this->parDRIFTP.exportContribution[4]);
		ReadDouble(parDRIFTP, "minFrequency", this->parDRIFTP.minFrequency);	
		ReadDouble(parDRIFTP, "maxFrequency", this->parDRIFTP.maxFrequency);
		ReadInt(parDRIFTP, "numberOfWavefrequencyDiagonals", this->parDRIFTP.numberOfWavefrequencyDiagonals);
		ReadBool(parDRIFTP, "waveDirInteraction", this->parDRIFTP.waveDirInteraction);	
	}
	XmlNode &parEXPORT = sim("parEXPORT");
	if (parEXPORT.GetCount()) {		
		ReadBool(parEXPORT, "runProgram", this->parEXPORT.runProgram);		
		XmlNode &hydFile = parEXPORT("hydFile");
		if (hydFile.GetCount()) {	
			ReadBool(hydFile, "export", this->parEXPORT.hyd.exportOn);	
			ReadBool(hydFile, "exportQTFContribution1", this->parEXPORT.hyd.exportQTF[0]);
			ReadBool(hydFile, "exportQTFContribution2", this->parEXPORT.hyd.exportQTF[1]);
			ReadBool(hydFile, "exportQTFContribution3", this->parEXPORT.hyd.exportQTF[2]);
			ReadBool(hydFile, "exportQTFContribution4", this->parEXPORT.hyd.exportQTF[3]);
			ReadBool(hydFile, "exportQTFContribution5", this->parEXPORT.hyd.exportQTF[4]);
			ReadInt(hydFile, "numberOfWavefrequencyDiagonals", this->parEXPORT.hyd.numberOfWavefrequencyDiagonals);
		}
		XmlNode &MonitorRelativeWaveHeights = parEXPORT("MonitorRelativeWaveHeights");
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
		XmlNode &CGNS = parEXPORT("CGNS");
		if (CGNS.GetCount()) {
			ReadBool(CGNS, "export", this->parEXPORT.cgns.exportOn);
			ReadDouble(CGNS, "waveFreq", this->parEXPORT.cgns.waveFreq);
			XmlNode &MonitorFlowDatas = CGNS("MonitorFlowDatas");
			for (int i = 0; i < MonitorFlowDatas.GetCount(); ++i) {			
				XmlNode &MonitorFlowData = MonitorFlowDatas.At(i);
				if (MonitorFlowData.GetCount()) {
					DiffracData::ParEXPORT::CGNS::MonitorFlowData &data = this->parEXPORT.cgns.monitorFlowData.Add();
					XmlNode &grids = MonitorFlowData("grids");
					for (int i = 0; i < grids.GetCount(); ++i) {
						XmlNode &grid = grids.At(i);
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
	XmlNode &bodies = sim("bodies");
	for (int i = 0; i < bodies.GetCount(); ++i) {
		XmlNode &body = bodies.At(i);	
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
	String ret = Format("%f", d);
	if (ret.Find(".") >= 0) {
		while (*(ret.Last()) == '0')
			ret.Remove(ret.GetCount()-1);
	}
	if (*(ret.Last()) == '.')
		ret.Remove(ret.GetCount()-1);
	return ret;
}

static String ToText(const Pointf &d) {
	return Format("%s,%s", ToText(d.x), ToText(d.y));
}

static String ToText(const Point3D &d) {
	return Format("%s,%s,%s", ToText(d.x), ToText(d.y), ToText(d.z));
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
    step = data[1] - data[0];

    for (int i = 1; i < n-1; i++) {
        double currentStep = data[i+1] - data[i];
        if (abs(currentStep - step) > tolerance)
            return false;
    }
    maxVal = data.Top();
    return true;
}

template <class Range>
String ToText(const Range& a) {
	double minVal, step, maxVal;
	if (AnalyzeStep(a, 1e-7, minVal, step, maxVal)) {
		return Format("%s(%s)%s", ToText(minVal), ToText(step), ToText(maxVal));
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
		<< Format("    <runProgram>%s</runProgram>\n", ToText(parVTKDIF.runProgram))
		<< "  </parVTKDIF>\n"
	;
	out << "  <parINIDIF>\n"
		<< Format("    <runProgram>%s</runProgram>\n", ToText(parINIDIF.runProgram))
		<< Format("    <dataBaseFn>%s</dataBaseFn>\n", parINIDIF.dataBaseFn)
		<< "    <quay>\n"
		<< Format("      <apply>%s</apply>\n", ToText(parINIDIF.quay.apply))
		<< Format("      <yQuay>%s</yQuay>\n", ToText(parINIDIF.quay.yQuay))
		<< "    </quay>\n"
		<< "    <basinWall>\n"
		<< Format("      <apply>%s</apply>\n", ToText(parINIDIF.basinWall.apply))
		<< "      <origin>\n"
		<< Format("            %s\n", ToText(parINIDIF.basinWall.origin))
		<< "          </origin>\n"
		<< "    </basinWall>\n"
		<< "    <density>\n"
		<< Format("          %s\n", ToText(parINIDIF.density))
		<< "        </density>\n"
		<< "    <springMatrix>\n"
		<< Format("      <fromGeometry>%s</fromGeometry>\n", ToText(parINIDIF.springMatrixfromGeometry))
		<< "      <springMatrixFn></springMatrixFn>\n"
		<< "    </springMatrix>\n"
    ;
	out	<< "    <dampLids>\n";
	for (int i = 0; i < parINIDIF.lids.size(); ++i) {
		out << "      <dampLid>\n";
		out << Format("        <origin>%s</origin>\n", ToText(parINIDIF.lids[i].origin));
		out << Format("        <orientation>%s</orientation>\n", ToText(parINIDIF.lids[i].orientation));
		out << Format("        <dampingValue>%s</dampingValue>\n", ToText(parINIDIF.lids[i].dampingValue));
		out << Format("        <length>%s</length>\n", ToText(parINIDIF.lids[i].length));
		out << Format("        <width>%s</width>\n", ToText(parINIDIF.lids[i].width));
		out << Format("        <NrPanelsLength>%s</NrPanelsLength>\n", ToText(parINIDIF.lids[i].NrPanelsLength));
		out << Format("        <NrPanelsWidth>%s</NrPanelsWidth>\n", ToText(parINIDIF.lids[i].NrPanelsWidth));
		out << Format("        <dampIncWave>%s</dampIncWave>\n", ToText(parINIDIF.lids[i].dampIncWave));
		out << "      </dampLid>\n";
	}
	out << "    </dampLids>\n";
	out << "  </parINIDIF>\n";	
	
	out << "  <parDIFFRAC>\n"
		<< Format("    <runProgram>%s</runProgram>\n", ToText(parDIFFRAC.runProgram))
		<< Format("    <waveDir>%s</waveDir>\n", ToText(parDIFFRAC.waveDir))
		<< Format("    <waterDepth>%s</waterDepth>\n", ToText(parDIFFRAC.waterDepth))
		<< "    <current>\n"
		<< "      <speed>\n"
		<< Format("            %s\n", ToText(parDIFFRAC.current.speed))
		<< "          </speed>\n"
		<< "      <direction>\n"
		<< Format("            %s\n", ToText(parDIFFRAC.current.direction))
		<< "          </direction>\n"
		<< "    </current>\n"
		<< Format("    <irregFreqSuppression>%s</irregFreqSuppression>\n", parDIFFRAC.irregFreqSuppression)
		<< Format("    <irregFreqDamping>%s</irregFreqDamping>\n", ToText(parDIFFRAC.irregFreqDamping))
		<< Format("    <frequencyStep>%s</frequencyStep>\n", ToText(parDIFFRAC.frequencyStep))
		<< Format("    <nFrequencies>%s</nFrequencies>\n", ToText(parDIFFRAC.nFrequencies))
		<< Format("    <minFrequency>%s</minFrequency>\n", ToText(parDIFFRAC.minFrequency))
		<< Format("    <exportKinematicsVTK>%s</exportKinematicsVTK>\n", ToText(parDIFFRAC.exportKinematicsVTK))
		<< "  </parDIFFRAC>\n"
	;

	out << "  <parDBRESP>\n"
		<< Format("    <runProgram>%s</runProgram>\n", ToText(parDBRESP.runProgram))
		<< "    <BodyInputs>\n";
	for (int i = 0; i < parDBRESP.inputs.size(); ++i) {
		out << "      <BodyInput>\n"
			<< Format("        <index>%s</index>\n", ToText(parDBRESP.inputs[i].index))
			<< "        <totalDampings>\n"
		;
		for (int j = 0; j < parDBRESP.inputs[i].total.size(); ++j) {
			out	<< "          <totalDamping>\n"
				<< Format("            <allowNegativeAddedDamping>%s</allowNegativeAddedDamping>\n", ToText(parDBRESP.inputs[i].total[j].allowNegativeAddedDamping))
				<< Format("            <mode>%s</mode>\n", parDBRESP.inputs[i].total[j].mode)
				<< Format("            <type>%s</type>\n", parDBRESP.inputs[i].total[j].type)
				<< Format("            <value>%s</value>\n", ToText(parDBRESP.inputs[i].total[j].value))
				<< "          </totalDamping>\n"
			;
		}
		out	<< "        </totalDampings>\n"
			<< "      </BodyInput>\n"
		;
	}
	out	<< "    </BodyInputs>\n"
		<< "  </parDBRESP>\n"
	;

	out << "  <parDRIFTP>\n"
		<< Format("    <runProgram>%s</runProgram>\n", ToText(parDRIFTP.runProgram))
		<< Format("    <exportContribution1>%s</exportContribution1>\n", ToText(parDRIFTP.exportContribution[0]))
		<< Format("    <exportContribution2>%s</exportContribution2>\n", ToText(parDRIFTP.exportContribution[1]))
		<< Format("    <exportContribution3>%s</exportContribution3>\n", ToText(parDRIFTP.exportContribution[2]))
		<< Format("    <exportContribution4>%s</exportContribution4>\n", ToText(parDRIFTP.exportContribution[3]))
		<< Format("    <exportContribution5>%s</exportContribution5>\n", ToText(parDRIFTP.exportContribution[4]))
		<< Format("    <minFrequency>%s</minFrequency>\n", ToText(parDRIFTP.minFrequency))
		<< Format("    <maxFrequency>%s</maxFrequency>\n", ToText(parDRIFTP.maxFrequency))
		<< Format("    <numberOfWavefrequencyDiagonals>%s</numberOfWavefrequencyDiagonals>\n", ToText(parDRIFTP.numberOfWavefrequencyDiagonals))
		<< Format("    <waveDirInteraction>%s</waveDirInteraction>\n", ToText(parDRIFTP.waveDirInteraction))
		<< "  </parDRIFTP>\n"
	;
	
	out << "  <parEXPORT>\n"
		<< Format("    <runProgram>%s</runProgram>\n", ToText(parEXPORT.runProgram))
		<< "    <hydFile>\n"
		<< Format("      <export>%s</export>\n", ToText(parEXPORT.hyd.exportOn))
		<< Format("      <exportQTFContribution1>%s</exportQTFContribution1>\n", ToText(parEXPORT.hyd.exportQTF[0]))
		<< Format("      <exportQTFContribution2>%s</exportQTFContribution2>\n", ToText(parEXPORT.hyd.exportQTF[1]))
		<< Format("      <exportQTFContribution3>%s</exportQTFContribution3>\n", ToText(parEXPORT.hyd.exportQTF[2]))
		<< Format("      <exportQTFContribution4>%s</exportQTFContribution4>\n", ToText(parEXPORT.hyd.exportQTF[3]))
		<< Format("      <exportQTFContribution5>%s</exportQTFContribution5>\n", ToText(parEXPORT.hyd.exportQTF[4]))
		<< Format("      <numberOfWavefrequencyDiagonals>%s</numberOfWavefrequencyDiagonals>\n", ToText(parEXPORT.hyd.numberOfWavefrequencyDiagonals))
		<< "    </hydFile>\n"
	;
	out << "    <MonitorRelativeWaveHeights>\n";
	for (int i = 0; i < parEXPORT.relHeights.size(); ++i) {
		out << "      <MonitorRelativeWaveHeight>\n"
			<< Format("        <name>%s</name>\n", parEXPORT.relHeights[i].name)
			<< Format("        <includeIncidentWave>%s</includeIncidentWave>\n", ToText(parEXPORT.relHeights[i].incIncident))
			<< Format("        <includeDiffractedWave>%s</includeDiffractedWave>\n", ToText(parEXPORT.relHeights[i].incDiffracted))
			<< Format("        <includeRadiatedWave>%s</includeRadiatedWave>\n", ToText(parEXPORT.relHeights[i].incRadiated))
			<< Format("        <includeMotions>%s</includeMotions>\n", ToText(parEXPORT.relHeights[i].incMotions))
			<< "        <relWaveHeightBodyInput>\n"
			<< Format("          <index>%s</index>\n", ToText(parEXPORT.relHeights[i].index))
			<< "          <referencePoint2D>\n"
			<< Format("            <coordinate>%s</coordinate>\n", ToText(parEXPORT.relHeights[i].referencePoint))
			<< "          </referencePoint2D>\n"
			<< "        </relWaveHeightBodyInput>\n"
			<< "      </MonitorRelativeWaveHeight>\n"
		;		
	}
	out << "    </MonitorRelativeWaveHeights>\n";
	
	out << "    <CGNS>\n"
		<< Format("      <export>false</export>\n", ToText(parEXPORT.cgns.exportOn))
		<< Format("      <waveFreq>0.8</waveFreq>\n", ToText(parEXPORT.cgns.exportOn))
		<< "      <MonitorFlowDatas>\n"
	;
	for (int i = 0; i < parEXPORT.cgns.monitorFlowData.size(); ++i) {
		out	<< "        <MonitorFlowData>\n"
			<< "          <grids>\n";
		for (int j = 0; j < parEXPORT.cgns.monitorFlowData[i].grids.size(); ++j) {
			out	<< "            <grid>\n"
				<< Format("              <NrOfPointsX>201</NrOfPointsX>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].Nx))
				<< Format("              <NrOfPointsY>101</NrOfPointsY>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].Ny))
				<< Format("              <origin>0.0,0.0,0.0</origin>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].origin))
				<< Format("              <spacing_x>5</spacing_x>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].sx))
				<< Format("              <spacing_y>5</spacing_y>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].sy))
				<< Format("              <orientation>0.0</orientation>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].orientation))
				<< Format("              <NrOfPointsZ>1</NrOfPointsZ>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].Nz))
				<< Format("              <spacing_z>1e-6</spacing_z>\n", ToText(parEXPORT.cgns.monitorFlowData[i].grids[j].sz))
				<< "            </grid>\n"
			;
		}
		out	<< "          </grids>\n"
			<< "        </MonitorFlowData>\n";
	}
	out	<< "      </MonitorFlowDatas>\n"
		<< Format("      <movingBodies>%s</movingBodies>\n", ToText(parEXPORT.cgns.movingBodies))
		<< Format("      <movingFreeSurface>%s</movingFreeSurface>\n", ToText(parEXPORT.cgns.movingFreeSurface))
		<< Format("      <numberOfTimeSteps>%s</numberOfTimeSteps>\n", ToText(parEXPORT.cgns.numberOfTimeSteps))
		<< Format("      <waveAmplificationFactor>%s</waveAmplificationFactor>\n", ToText(parEXPORT.cgns.waveAmplificationFactor))
		<< Format("      <waveDir>%s</waveDir>\n", ToText(parEXPORT.cgns.waveDir))
		<< "    </CGNS>\n"
    ;	
	out << "  </parEXPORT>\n";
	
	out	<< "  <bodies>\n";
	for(int i = 0; i < bodies.size(); ++i) {
		out << "    <body>\n"
			<< "      <hstat>\n"
			<< Format("        <lengthBetweenPerp>%s</lengthBetweenPerp>\n", ToText(bodies[i].hstat.lengthBetweenPerp))
			<< Format("        <draft>%s</draft>\n", ToText(bodies[i].hstat.draft))
			<< "      </hstat>\n"
			<< Format("      <index>%s</index>\n", ToText(bodies[i].index))
			<< "      <mesh>\n"
			<< Format("        <meshFn>%s</meshFn>\n", bodies[i].meshFn)
			<< "      </mesh>\n"
			<< Format("      <name>%s</name>\n", bodies[i].name)
			<< "      <position>\n"
			<< Format("        <translation>%s</translation>\n", ToText(bodies[i].translation))
			<< Format("        <rotation>0</rotation>\n", ToText(bodies[i].rotation))
			<< "      </position>\n"
			<< "      <massElements>\n"
		;
		for (int j = 0; j < bodies[i].massElements.size(); ++j) {
			out << "        <massElement>\n"
				<< Format("          <COGwrtKeel>%s</COGwrtKeel>\n", ToText(bodies[i].massElements[j].COGwrtKeel))
				<< Format("          <mass>%s</mass>\n", ToText(bodies[i].massElements[j].mass))
				<< Format("          <rollRadiusGyr>%s</rollRadiusGyr>\n", ToText(bodies[i].massElements[j].rollRadiusGyr))
				<< Format("          <pitchRadiusGyr>%s</pitchRadiusGyr>\n", ToText(bodies[i].massElements[j].pitchRadiusGyr))
				<< Format("          <yawRadiusGyr>%s</yawRadiusGyr>\n", ToText(bodies[i].massElements[j].yawRadiusGyr))
				<< "        </massElement>\n"
				<< "      </massElements>\n"
			;
		}
		out << "    </body>\n";
	}
    out	<< "  </bodies>\n";	
	
	out << Format("  <projectNumber>98800</projectNumber>\n", ToText(projectNumber))
		<< Format("  <nProcs>4</nProcs>\n", ToText(nProcs)) 
		<< "</sim>"
	;
	return out;
}


CONSOLE_APP_MAIN
{
    StdLogSetup(LOG_COUT|LOG_FILE);
    String path = "C:\\Desarrollo\\BEMRosetta\\KKK\\SBS.input.xml";

    String xml = LoadFile(path);
    if(IsNull(xml)) {
        Cout() << "Cannot read file: " << path << '\n';
        return;
    }
	
	DiffracData data;
	data.LoadXML(xml);
	
	SaveFile("C:\\Desarrollo\\BEMRosetta\\KKK\\SBS.processed.xml", data.SaveXML());
	
	Cout() << "\nEnd:";
	ReadStdIn();
}