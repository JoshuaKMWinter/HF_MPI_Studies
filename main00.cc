#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <string> 
#include <cmath> 

using namespace Pythia8;

double Delta_Phi(double phi1, double phi2) {				
	const double pi = 3.14159265358979323846;
	double dphi;
	if ((phi1>=0 && phi2>=0) || (phi1<0 && phi2<0)) dphi = std::abs(phi1-phi2);
	else dphi = 2*pi - std::abs(phi1-phi2);
	return dphi;
}

int main() {

	std::ofstream datafile; datafile.open("EventData.csv");
	datafile << "D_pT,D_phi,D_eta,Dbar_pT,Dbar_phi,Dbar_eta,c_pT,c_phi,c_eta,cbar_pT,cbar_phi,cbar_eta\n";

	Pythia pythia;
	pythia.readString("Beams:eCM = 13000.");
	pythia.readString("HardQCD:hardccbar = on");
	/*pythia.readString("4:onMode = 0"); 
	pythia.readString("4:onIfAny = 421"); this will not work as the confinement is not registered as a decay channel of the c-quark. Can confinement to a D0 be forced?*/ 

	pythia.init();

	//Hist h_nD("Number of D0",4,-0.5,3.5);

	// Begin event loop. Generate event. Skip if error. List first one.
	for (int iEvent = 0; iEvent < 200000; ++iEvent) { //100k
		if (!pythia.next()) continue;
		
		double c_pt  = 0, cbar_pt  = 0, D_pt  = 0, Dbar_pt  = 0;
		double c_phi = 0, cbar_phi = 0, D_phi = 0, Dbar_phi = 0;
		double c_eta = 0, cbar_eta = 0, D_eta = 0, Dbar_eta = 0;
		int num_D = 0, num_Dbar = 0, num_c = 0, num_cbar = 0;
		int D_ind = -1, Dbar_ind = -1;
		
		for (int i = 0; i < pythia.event.size(); ++i) {
			
			//find status and id of current particle and the indices and ids of mother particles
			int status = pythia.event[i].status(),   id = pythia.event[i].id();

		       	/*int mother1 = pythia.event[i].mother1(), mother2 = pythia.event[i].mother2();
		        int m1_id = pythia.event[mother1].id(),  m2_id = pythia.event[mother2].id();*/

			//kinematics of final charm quarks prior to confinement
			if (id == 4 && status/10 == -7) { //71-79 preparation for hadronisation
				c_pt  = pythia.event[i].pT();
				c_phi = pythia.event[i].phi();
				c_eta = pythia.event[i].eta();
				num_c++;
			} 
			if (id == -4 && status/10 == -7) {
				cbar_pt  = pythia.event[i].pT();
				cbar_phi = pythia.event[i].phi();
				cbar_eta = pythia.event[i].eta();
				num_cbar++;
			}

			if (id ==  421) { num_D++; D_ind = i; }
			if (id == -421) { num_Dbar++; Dbar_ind = i; }
							
		}

		//choose only events with a single D0-D0bar pair
		if (num_D == 1 && num_Dbar == 1 && num_c == 1 && num_cbar == 1) {
			D_pt  = pythia.event[D_ind].pT();  Dbar_pt  = pythia.event[Dbar_ind].pT();
			D_phi = pythia.event[D_ind].phi(); Dbar_phi = pythia.event[Dbar_ind].phi();
			D_eta = pythia.event[D_ind].eta(); Dbar_eta = pythia.event[Dbar_ind].eta();

			//calculate angular separations	
			//dphi_cc = Delta_Phi(c_phi,cbar_phi); dphi_DD = Delta_Phi(D_phi,Dbar_phi); 
			//deta_cc = std::abs(c_eta - cbar_eta); deta_DD = std::abs(D_eta - Dbar_eta);

			//write to datafile
			datafile <<D_pt<<","<<D_phi<<","<<D_eta
				 <<","<<Dbar_pt<<","<<Dbar_phi<<","<<Dbar_eta
				 <<","<<c_pt<<","<<c_phi<<","<<c_eta
				 <<","<<cbar_pt<<","<<cbar_phi<<","<<cbar_eta<<"\n";
		}
	// End of event loop.
	}
	datafile.close();

	pythia.stat();
//	cout << h_nD;
	
	return 0;
}
	

