#include <Core/Core.h>
#include <Functions4U/Functions4U.h>
#include <plugin/Eigen/Eigen.h>
#include <STEM4U/Sundials.h>

using namespace Upp;
using namespace Eigen;


bool FuerzaAmarre(double lmb_seco, double dens, double bl, double LTOT, double XTOT, double ha, double hb, double Fbuque, double &Fx, double &Fya, double &Fyb, bool &roto, double *x = 0, double *y = 0, int numData = 0) {
	double g = 9.81;       	// m/s2. Aceleración de la gravedad.
	double dagua = 1020;	// Kg/m3
	// lmb 		kg/m. Peso sumergido por metro lineal de cadena (=(6.85/7.85)*lmb_seco)
	// dens 	kg/m3. Densidad de la cadena.
	// bl		N. Carga de rotura (Break Load) del amarre. 
	double lmb = lmb_seco*(dens - dagua)/dens;
	// LTOT 	m. Longitud total de cadena.
	// ha 		m. Altura desde el fondo del ancla/muerto.
	// hb 		m. Altura desde el fondo del buque
	// XTOT 	m. Distancia horizontal entre el ancla/muerto y el buque.
	// Reacciones sobre el buque
	// Fx 		N. Fuerza horizontal sobre el buque
	// Fy 		N. Fuerza vertical descendente sobre el buque
	
  	Buffer<double> udata(1);
  	udata[0] = 1;	
	SolveNonLinearEquationsSun(udata, 1, [&](const double y[], double residuals[])->bool {
		double Blimit = y[0];
		residuals[0] =  sqrt(ha*(ha + 2*Blimit)) + sqrt(hb*(hb + 2*Blimit)) - LTOT; 
		return true;
	});
	double Blimit = udata[0];

    double Xcata_limit = Blimit * acosh(ha/Blimit + 1);
    double Xcatb_limit = Blimit * acosh(hb/Blimit + 1);
    double XTOT_limit = Xcata_limit + Xcatb_limit;
	
	for (int i = 0; i < numData; ++i) 
    	x[i] = y[i] = 0;
	
	if (XTOT < LTOT - ha - hb) {				
		printf("\nCaso AA");
		
        Fx = 0;
        Fya = lmb*g*ha;
        Fyb = lmb*g*hb;
        //double LTOT_res = LTOT;   
        //double h_res = hb - ha;
	} else if (XTOT == LTOT - ha - hb) {		
		printf("\nCaso BB");
		
        Fx = 0;
        Fya = lmb*g*ha;
        Fyb = lmb*g*hb;
        //double LTOT_res = LTOT;   
        //double h_res = hb - ha;
    } else if (XTOT < XTOT_limit) {
        printf("\nCaso CC");
		int neq = 1;
	
	  	Buffer<int> consdata(1);
	  	consdata[0] = 2;
	  	
	  	Buffer<double>udata(1);
	  	udata[0] = Blimit/2.;	
		SolveNonLinearEquationsSun(udata, neq, [&](const double y[], double residuals[])->bool {
			double B = y[0];
			residuals[0] = XTOT - B*acosh(ha/B + 1) - B*acosh(hb/B + 1) + sqrt(ha*(ha + 2*B)) + sqrt(hb*(hb + 2*B)) - LTOT; 
			return true;
		}, consdata);		
		double B = udata[0];

        double LTOT_res = XTOT - B*acosh(ha/B + 1) - B*acosh(hb/B + 1) + sqrt(ha*(ha + 2*B)) + sqrt(hb*(hb + 2*B));
        double Xcata = B*acosh(ha/B + 1);
        double Xcatb = B*acosh(hb/B + 1);
        double h_res = B*(cosh(Xcatb/B) - 1) - B*(cosh(Xcata/B) - 1);
		if (abs(LTOT_res - LTOT) > 0.0001 || abs(h_res - (hb - ha)) > 0.0001) {
			printf("\nError en calculo. LTOT_res: %f, LTOT: %f, h_res: %f, ha: %f", LTOT_res, LTOT, h_res, ha);
			return false;
		}
        Fx = lmb*g*B;          	// N. Fuerza horizontal sobre el buque
        Fya = Fx*sinh(Xcata/B); // N. Fuerza vertical descendente sobre el ancla
        Fyb = Fx*sinh(Xcatb/B); // N. Fuerza vertical descendente sobre el buque
        // DIBUJOS Ejes con origen en el ancla/muerto
       	if (numData != 0) {
            double delta = XTOT/(numData - 1);
            for (int i = 0; i < numData; ++i) {
                x[i] = i*delta;
                if (x[i] < Xcata)
	        		y[i] = B*(cosh((Xcata - x[i])/B) - 1) - ha;
                else if (x[i] > (XTOT - Xcatb))
                    y[i] = B*(cosh((max(0., x[i] - (XTOT - Xcatb)))/B) - 1) - ha;
                else
                    y[i] = 0;
	        }
        }
    } else if (XTOT < sqrt(pow2(LTOT) - pow2(hb - ha))) {
        printf("\nCaso DD");
		int neq = 2;
        double h = hb - ha;
        double Lcat = LTOT;
        double Xcat = XTOT;

	  	Buffer<double> udata(neq);
	  	Buffer<int> consdata(2);
	  	consdata[0] = 2;
	  	consdata[1] = 0;
	  	
	  	double x1_0, B1_0, x1, B;
	  	bool done = false;
	  	for (x1_0 = 0; x1_0 < abs(XTOT) && !done; x1_0 += abs(XTOT)/4) {
	  		for (B1_0 = fabs(Blimit); B1_0 < 10*B1_0 && !done; B1_0 += 5*fabs(Blimit)/4) {		
		  		udata[0] = B1_0;	
		  		udata[1] = x1_0;
				SolveNonLinearEquationsSun(udata, neq, [&](const double y[], double residuals[])->bool {
					double B = y[0];
					double x1 = y[1];
					
					residuals[0] = B*(sinh((x1 + Xcat)/B) - sinh(x1/B)) - Lcat;
					residuals[1] = B*(cosh((x1 + Xcat)/B) - cosh(x1/B)) - h;
					
					return true;
				}, consdata);
				B = udata[0];
				x1 = udata[1];		
				if (fabs((x1 + Xcat)/B) < 3)
					done = true;
	  		}
	  	}
	  	if (!done) {
			printf("\nError en resolucion de ecuaciones");
			return false;
		}

     	double Lcat_res = B*(sinh((x1 + Xcat)/B) - sinh(x1/B));
     	double h_res    = B*(cosh((x1 + Xcat)/B) - cosh(x1/B));
		if (abs(Lcat_res - Lcat) > 0.0001 || abs(h_res - h) > 0.0001) {
			printf("\nError en calculo. Lcat_res: %f, Lcat: %f, h_res: %f, h: %f", Lcat_res, Lcat, h_res, h);
			return false; 
		}
        Fx = lmb*g*fabs(B);           			// N. Fuerza horizontal sobre el buque
        Fya = -Fx*sinh(x1/B);  			
        Fyb = Fx*sinh((x1 + Xcat)/B);  	
        //double LTOT_res = Lcat_res;
        // DIBUJOS Ejes con origen en el ancla/muerto
        if (numData != 0) {
            double delta = Xcat/(numData - 1);
            for (int i = 0; i < numData; ++i) {
                x[i] = i*delta;
	        	y[i] = B*(cosh((x[i] + x1)/B) - 1) - B*(cosh(x1/B) - 1);
	        }
        }
    } else if (XTOT == sqrt(pow2(LTOT) - pow2(hb - ha))) {
        printf("\nCaso EE");
           
        double h = hb - ha;
        Fx = Fbuque;        // N. Fuerza horizontal sobre el buque
        Fya = -Fx*h/XTOT;	// N. Fuerza vertical descendente sobre el ancla
        Fyb = -Fya;  		// N. Fuerza vertical descendente sobre el buque

        //double h_res = h;    
        //double LTOT_res = LTOT;
    } else {
		printf("\nCadena rota por exceso de longitud!");
		return false;
    }

	roto = max(sqrt(pow2(Fx) + pow2(Fya)), sqrt(pow2(Fx) + pow2(Fyb))) >= bl ;
	return true;
}

#define NUMDATA 30

void Ensayo_IHC() {
	printf("\nCalculo amarre IHC");
	
	double Fx, Fya, Fyb;
	bool roto;
	double lmb_seco = 0.166;   
	double dens = 7850;		// Kg/m3
	double bl = 3965000;	// N Carga de rotura Grado R3S
	double Fbuque = 0;
	double x[NUMDATA], y[NUMDATA];
	
	double ha = 0;
	double hb = 1.17;
	double LTOT = 4.167;
	double XTOT = 3.617;	// 3.8
	
	if (FuerzaAmarre(lmb_seco, dens, bl, LTOT, XTOT, ha, hb, Fbuque, Fx, Fya, Fyb, roto, x, y, NUMDATA))  {
		double Fa = sqrt(pow2(Fx) + pow2(Fya));
		double Fb = sqrt(pow2(Fx) + pow2(Fyb));
		printf("\tXTOT=%.3f, Fx=%.2f\nFya=%.2f, Fa=%.2f, ang_a=%.1f grad\nFyb=%.2f, Fb=%.2f, ang_b=%.1f grad\n(%.0f%%)", 
				XTOT, Fx, 
				-Fya, Fa, atan2(Fya, Fx)*180./M_PI, 
				-Fyb, Fb, atan2(Fyb, Fx)*180./M_PI,
				100*max(Fa, Fb)/bl);
		printf("\nX\t\t;Y");
		for (int i = 0; i < NUMDATA; ++i) 
			printf("\n%f\t;%f", x[i], y[i]);
	}

	printf("\nPulsa una tecla para finalizar");
	getchar();	
	
}

void Ensayo_Nautilus() {
	printf("\nCalculo amarre Nautilus BE");
	
	double Fx, Fya, Fyb;
	bool roto;
	double lmb_seco = 346.74;   
	double dens = 7850;		// Kg/m3
	double bl = 1.49E9;		// N Carga de rotura Grado R3S
	double Fbuque = 0;
	double x[NUMDATA], y[NUMDATA];
	
	double ha = 0;
	double hb = 45.5+8.4;
	double LTOT = 480;
	double XTOT = sqrt(2*sqr(356.37-28.5)) + 11.25;	
		
	if (FuerzaAmarre(lmb_seco, dens, bl, LTOT, XTOT, ha, hb, Fbuque, Fx, Fya, Fyb, roto, x, y, NUMDATA))  {
		double Fa = sqrt(pow2(Fx) + pow2(Fya));
		double Fb = sqrt(pow2(Fx) + pow2(Fyb));
		printf("\nFx=%.2f\nFya=%.2f, Fa=%.2f, ang_a=%.1f grad\nFyb=%.2f, Fb=%.2f, ang_b=%.1f", 
				Fx, 
				-Fya, Fa, atan2(Fya, Fx)*180./M_PI, 
				-Fyb, Fb, atan2(Fyb, Fx)*180./M_PI);
		printf("\nX\t;Y");
		for (int i = 0; i < NUMDATA; ++i) 
			printf("\n%.2f\t;%.2f", x[i], y[i]);
	}

	printf("\nPulsa una tecla para finalizar");
	getchar();	
	
}

void Amarre() {
	Ensayo_Nautilus();
	//Ensayo_IHC();
	return;
	
	printf("\nCalculo de fuerza de amarre");
	
	double Fx, Fya, Fyb;
	bool roto;
	double lmb_seco = 90;   // 64 mm diam
	double dens = 7850;		// Kg/m3
	double bl = 3965000;	// N Carga de rotura Grado R3S
	double Fbuque = 0;
	double x[NUMDATA], y[NUMDATA];
	
	double ha = 0;
	double hb = 2.245;
	double LTOT = 3;
	//double XTOT = 1.708;
	lmb_seco = 0.16;
	
	/*if (FuerzaAmarre(lmb_seco, dens, bl, LTOT, XTOT, ha, hb, Fbuque, Fx, Fya, Fyb, roto, x, y, NUMDATA))  {
		printf("\tXTOT=%.3f, Fx=%.2f, Fya=%.2f, Fyb=%.2f (%.0f%%)", XTOT, Fx, Fya, Fyb, 100*max(sqrt(pow2(Fx) + pow2(Fya)), sqrt(pow2(Fx) + pow2(Fyb)))/bl);
		printf("\nX\t\t;Y");
		for (int i = 0; i < NUMDATA; ++i) 
			printf("\n%f\t;%f", x[i], y[i]);
	}*/
	for (double XTOT = 1; XTOT < 2; XTOT += 0.01) {
		if (FuerzaAmarre(lmb_seco, dens, bl, LTOT, XTOT, ha, hb, Fbuque, Fx, Fya, Fyb, roto, x, y, NUMDATA))  {
			printf("\tXTOT=%.2f, Fx=%.2f, Fya=%.2f, Fyb=%.2f (%.0f%%)", XTOT, Fx, Fya, Fyb, 100*max(sqrt(pow2(Fx) + pow2(Fya)), sqrt(pow2(Fx) + pow2(Fyb)))/bl);
//			printf("\nX\t\t;Y");
//			for (int i = 0; i < NUMDATA; ++i) 
//				printf("\n%f\t;%f", x[i], y[i]);
		}
	}
	printf("\nPulsa una tecla para finalizar");
	getchar();
}
