#pragma once
#include <Core/Core.h>

using namespace Upp;

struct DiffracData {
	struct ParVTKDIF {
		bool runProgram = false;
	};
	struct ParINIDIF {
		struct Quay {
		    bool apply = false;
		    double yQuay = 0;
		};
		struct BasinWall {
		    bool apply = false;
		    Pointf origin = Pointf(0, 0);
		};
		struct DampingLid {
		    Pointf origin = Pointf(0, 0);
		    double orientation = 0.0;
		    double dampingValue = 0.0;
		    double length = 0.0;
		    double width = 0.0;
		    int NrPanelsLength = 0;
		    int NrPanelsWidth = 0;
		    bool dampIncWave = false;
		};
		bool runProgram = false;
	    String dataBaseFn;
	    Quay quay;
	    BasinWall basinWall;
	    double density = 1025./1000.;
	    bool springMatrixfromGeometry = false;
	    UArray<DampingLid> lids;
	};
	struct ParDIFFRAC {
		struct Current {
		    double speed = 0;
		    double direction = 0;
		};
		bool runProgram = false;
	    UVector<double> waveDir;
	    double waterDepth = 0.0;
	    Current current;
	    String irregFreqSuppression;
	    double irregFreqDamping = 0.0;
	    double frequencyStep = 0.0;
	    int nFrequencies = 0;
	    double minFrequency = 0.0;
	    bool exportKinematicsVTK = false;
	};
	struct ParDBRESP {
		struct BodyInputDamping {
			struct BodyTotalDamping {
			    String mode;
			    String type;
			    bool allowNegativeAddedDamping = false;
			    double value = 0;
			};
		    int index = 0;
		    UArray<BodyTotalDamping> total;
		};
		bool runProgram = false;
	    UArray<BodyInputDamping> inputs;
	};
	struct ParDRIFTP {
		bool runProgram = false;
	    bool exportContribution[5] = {false,false,false,false,false};
	    double minFrequency = 0;
	    double maxFrequency = 0;
	    int numberOfWavefrequencyDiagonals = 0;
	    bool waveDirInteraction = false;
	};
	struct ParEXPORT {
		struct HydFileExport {
		    bool exportOn = true;
		    bool exportQTF[5] = {false,false,false,false,false};
		    int numberOfWavefrequencyDiagonals = 0;
		};
		struct RelWaveHeight {
		    String name;
		    bool incIncident = true;
		    bool incDiffracted = true;
		    bool incRadiated = true;
		    bool incMotions = false;
		    int index = 0;
		    Pointf referencePoint;
		};
		struct CGNS {
			struct MonitorFlowData {
				struct Grid {
					int Nx = 0, Ny = 0, Nz = 0;
				    double sx = 0, sy = 0, sz = 0;
				    Point3D origin = Point3D(0,0,0);
				    double orientation = 0.0;	
				};
				UArray<Grid> grids;	
			};
		    bool exportOn = false;
		    double waveFreq = 0.0;
			UArray<MonitorFlowData> monitorFlowData;		    
		    bool movingBodies = false;
		    bool movingFreeSurface = false;
		    int numberOfTimeSteps = 0;
		    double waveAmplificationFactor = 1.0;
		    double waveDir = 0.0;
		};
	    bool runProgram = true;
	    HydFileExport hyd;
	    UArray<RelWaveHeight> relHeights;
	    CGNS cgns;
	};
	struct Body {
		struct HStat {
		    double lengthBetweenPerp = 0.0;
		    double draft = 0.0;
		};
		struct MassElement {
		    Point3D COGwrtKeel = Point3D(0,0,0);
		    double mass = 0.0;
		    double rollRadiusGyr = 0.0, pitchRadiusGyr = 0.0, yawRadiusGyr = 0.0;
		};
	    HStat hstat;
	    int index = 0;
	    String meshFn;
	    String name;
	    Point3D translation = Point3D(0,0,0);
	    double rotation = 0.0;
	    UArray<MassElement> massElements;
	};
    ParVTKDIF parVTKDIF;
    ParINIDIF parINIDIF;
    ParDIFFRAC parDIFFRAC;
    ParDBRESP parDBRESP;
    ParDRIFTP parDRIFTP;
    ParEXPORT parEXPORT;
    UArray<Body> bodies;
    int projectNumber;
    int nProcs = 1;
    
    void LoadXML(const String &xml);
    String SaveXML();
};

