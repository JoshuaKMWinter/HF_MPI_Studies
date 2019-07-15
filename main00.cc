#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <string> 
#include <cmath> 

using namespace Pythia8;

int main() {

	std::ofstream datafile; datafile.open("EventData.csv");
	datafile << "cLO_pT,cLO_phi,cLO_eta,c_pT,c_phi,c_eta,D_pT,D_phi,D_eta,cbarLO_pT,cbarLO_phi,cbarLO_eta,cbar_pT,cbar_phi,cbar_eta,Dbar_pT,Dbar_phi,Dbar_eta,Kpl_pT,Kpl_phi,Kpl_eta,Kmi_pT,Kmi_phi,Kmi_eta,pipl_pT,pipl_phi,pipl_eta,pimi_pT,pimi_phi,pimi_eta,multiplicity,D_cone_mult,D_ptcone,Dbar_cone_mult,Dbar_ptcone\n";
	
	//std::ofstream multfile; multfile.open("NDs_multData.csv");
	//multfile << "num_D,num_Dbar,num_Ds,num_Dsbar,num_lambdac,num_lambdacbar,multiplicity\n";

	Pythia pythia;
	pythia.readString("Beams:eCM = 13000.");
	pythia.readString("HardQCD:hardccbar = on");
	pythia.readString("421:onMode = off"); 
	pythia.readString("421:onIfMatch = 211 321");  

	pythia.init();

	//bool print = true; //boolean to print first selected event
	
	// Begin event loop. Generate event. Skip if error. List first one.
	for (int iEvent = 0; iEvent < 200000; ++iEvent) { //100k
		if (!pythia.next()) continue;
		
		double cLO_pt = 0,  cLO_phi = 0,  cLO_eta = 0,  cbarLO_pt = 0, cbarLO_phi = 0, cbarLO_eta = 0;
		double c_pt = 0,    c_phi = 0,    c_eta = 0,    cbar_pt = 0,   cbar_phi = 0,   cbar_eta = 0;
		double D_pt = 0,    D_phi = 0,    D_eta = 0,    Dbar_pt = 0,   Dbar_phi = 0,   Dbar_eta = 0;
		double Kpl_pt = 0,  Kpl_phi = 0,  Kpl_eta = 0,  Kmi_pt = 0,    Kmi_phi = 0,    Kmi_eta = 0;
		double pipl_pt = 0, pipl_phi = 0, pipl_eta = 0, pimi_pt = 0,   pimi_phi = 0,   pimi_eta = 0;
		double D_ptcone = 0, Dbar_ptcone = 0;

		int num_D = 0, num_Dbar = 0, num_c = 0, num_cbar = 0, num_pipl = 0, num_pimi = 0, num_Kpl = 0, num_Kmi = 0, num_Ds=0, num_Dsbar=0, num_lambdac=0, num_lambdacbar=0;
		int D_ind = -1, Dbar_ind = -1;
		int multiplicity = 0, D_cone_mult = 0, Dbar_cone_mult = 0;
		
		for (int i = 0; i < pythia.event.size(); ++i) {
			
			//use all final particles as an estimate of multiplicity
			if (pythia.event[i].isFinal()) multiplicity++;

			//find status and id of current particle and the indices and ids of mother particles
			int status = pythia.event[i].status(),   id = pythia.event[i].id();
			int mother1 = pythia.event[i].mother1(), mother2 = pythia.event[i].mother2();


			if (id == 4 && status == -23) { //outgoing from hardest subprocess
				cLO_pt  = pythia.event[i].pT();
				cLO_phi = pythia.event[i].phi();
				cLO_eta = pythia.event[i].eta();
			} 
			if (id == -4 && status == -23) {
				cbarLO_pt  = pythia.event[i].pT();
				cbarLO_phi = pythia.event[i].phi();
				cbarLO_eta = pythia.event[i].eta();
			}

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

			if (id ==   421) { num_D++; D_ind = i; }
			if (id ==  -421) { num_Dbar++; Dbar_ind = i; }
			//if (id ==   431)  num_Ds++;
			//if (id ==  -431)  num_Dsbar++;
			//if (id ==  4122)  num_lambdac++;
			//if (id == -4122)  num_lambdacbar++;

			if (id == 211 && (mother1 == D_ind || mother2 == D_ind)) {
				pipl_pt  = pythia.event[i].pT();
				pipl_phi = pythia.event[i].phi();
				pipl_eta = pythia.event[i].eta();
				num_pipl++;
			}
			if (id == -211 && (mother1 == Dbar_ind || mother2 == Dbar_ind)) {
				pimi_pt  = pythia.event[i].pT();
				pimi_phi = pythia.event[i].phi();
				pimi_eta = pythia.event[i].eta();
				num_pimi++;
			}
			if (id == 321 && (mother1 == Dbar_ind || mother2 == Dbar_ind)) {
				Kpl_pt  = pythia.event[i].pT();
				Kpl_phi = pythia.event[i].phi();
				Kpl_eta = pythia.event[i].eta();
				num_Kpl++;
			}
			if (id == -321 && (mother1 == D_ind || mother2 == D_ind)) {
				Kmi_pt  = pythia.event[i].pT();
				Kmi_phi = pythia.event[i].phi();
				Kmi_eta = pythia.event[i].eta();
				num_Kmi++;
			}
		}
		
		//multfile << num_D << "," <<  num_Dbar << "," <<  num_Ds << "," <<  num_Dsbar << "," <<  num_lambdac << "," <<  num_lambdacbar << "," <<  multiplicity <<"\n";

		//choose only events with 1 of each particle, as generation will ensure correct matching
		if (num_D == 1 && num_Dbar == 1 && num_c == 1 && num_cbar == 1 && num_Kpl == 1 && num_Kmi == 1 && num_pipl == 1 && num_pimi == 1) {
			//if (print) {pythia.event.list(); print = false;} //print first selected event

			D_pt  = pythia.event[D_ind].pT();  Dbar_pt  = pythia.event[Dbar_ind].pT();
			D_phi = pythia.event[D_ind].phi(); Dbar_phi = pythia.event[Dbar_ind].phi();
			D_eta = pythia.event[D_ind].eta(); Dbar_eta = pythia.event[Dbar_ind].eta();
						
			//Find total pT and number of final particles inside cone of radius R around each D meson
			const double R = 0.7;
			for (int i = 0; i < pythia.event.size(); ++i) {
				if (!pythia.event[i].isFinal()) continue; //only look at final particles
				double eta = pythia.event[i].eta(), phi = pythia.event[i].phi(), pt = pythia.event[i].pT();
				double D_deltaR = std::pow( std::pow(phi - D_phi, 2) + std::pow(eta - D_eta, 2) , 0.5); //Distance from D0
				double Dbar_deltaR = std::pow( std::pow(phi - Dbar_phi, 2) + std::pow(eta - Dbar_eta, 2) , 0.5); //Distance from D0bar

				if (D_deltaR < R) {D_cone_mult++; D_ptcone += pt;}
				if (Dbar_deltaR < R) {Dbar_cone_mult++; Dbar_ptcone += pt;}
			}

			//write to datafile
			datafile <<cLO_pt<<","<<cLO_phi<<","<<cLO_eta<<","<<c_pt<<","<<c_phi<<","<<c_eta<<","<<D_pt<<","<<D_phi<<","<<D_eta<<","
				 <<cbarLO_pt<<","<<cbarLO_phi<<","<<cbarLO_eta<<","<<cbar_pt<<","<<cbar_phi<<","<<cbar_eta<<","<<Dbar_pt<<","<<Dbar_phi<<","<<Dbar_eta<<","
				 <<Kpl_pt<<","<<Kpl_phi<<","<<Kpl_eta<<","<<Kmi_pt<<","<<Kmi_phi<<","<<Kmi_eta<<","
				 <<pipl_pt<<","<<pipl_phi<<","<<pipl_eta<<","<<pimi_pt<<","<<pimi_phi<<","<<pimi_eta<<","
				 <<multiplicity<<","<<D_cone_mult<<","<<D_ptcone<<","<<Dbar_cone_mult<<","<<Dbar_ptcone<<"\n";
			
		}
	// End of event loop.
	}
	//multfile.close();

	datafile.close();

	pythia.stat();
		
	return 0;
}
	

