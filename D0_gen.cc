#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <string> 
#include <cmath>
#include <vector>

using namespace Pythia8;

int main() {

	const int chadron_id = 421; 	//PDG ID of charm meson for selected decay
	const int d1_id = 211;		//Decay product ids. pi+
	const int d2_id = 321;		//K+

	//outfile for selected event kinematic analysis
	std::ofstream datafile; datafile.open("CSVs/D0_EventData.csv");
	datafile << "cLO_pT,cLO_phi,cLO_eta,c_pT,c_phi,c_eta,chad_pT,chad_phi,chad_eta,cbarLO_pT,cbarLO_phi,cbarLO_eta,cbar_pT,cbar_phi,cbar_eta,cbarhad_pT,cbarhad_phi,cbarhad_eta,d1_pT,d1_phi,d1_eta,d1bar_pT,d1bar_phi,d1bar_eta,d2_pT,d2_phi,d2_eta,d2bar_pT,d2bar_phi,d2bar_eta,multiplicity,chad_cone_mult,chad_ptcone,cbarhad_cone_mult,cbarhad_ptcone\n";
	
	Pythia pythia;
	pythia.readString("Beams:eCM = 13000.");
	pythia.readString("HardQCD:hardccbar = on");
	pythia.readString("421:onMode = off"); 
	pythia.readString("421:onIfMatch = 211 321"); 
	//pythia.readString("PhaseSpace:pTHatMin = 4.");
	pythia.init();

	//bool print = true; //boolean to print first selected event
	
	// Begin event loop. Generate event. Skip if error. List first one.
	for (int iEvent = 0; iEvent < 200000; ++iEvent) { //100k
		if (!pythia.next()) continue;
		
		double cLO_pt = 0,  cLO_phi = 0,  cLO_eta = 0,  cbarLO_pt = 0, cbarLO_phi = 0, cbarLO_eta = 0; //Leading Order charm kinematics
		double c_pt = 0,    c_phi = 0,    c_eta = 0,    cbar_pt = 0,   cbar_phi = 0,   cbar_eta = 0;   //Charm kin. at hadronisation
		double chad_pt = 0, chad_phi = 0, chad_eta = 0, cbarhad_pt = 0,   cbarhad_phi = 0,   cbarhad_eta = 0;   //Charm kin. in hadron
		double d1_pt = 0,   d1_phi = 0,   d1_eta = 0,   d1bar_pt = 0,  d1bar_phi = 0,  d1bar_eta = 0;  //kinematics of daughter particles
		double d2_pt = 0,   d2_phi = 0,   d2_eta = 0,   d2bar_pt = 0,  d2bar_phi = 0,  d2bar_eta = 0;  
		double chad_ptcone = 0, cbarhad_ptcone = 0;

		int n_chad = 0, n_cbarhad = 0, n_c = 0, n_cbar = 0, n_d1 = 0, n_d1bar = 0, n_d2 = 0, n_d2bar = 0; //counters for number of selected particles
		int chad_ind = -1, cbarhad_ind = -1; //integers to track charm meson index
		int mult = 0, chad_cone_mult = 0, cbarhad_cone_mult = 0; //multiplicity counters.

		for (int i = 0; i < pythia.event.size(); ++i) {
			
			//use all final particles as an estimate of multiplicity
			if (pythia.event[i].isFinal()) mult++;

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
				n_c++;
			} 
			if (id == -4 && status/10 == -7) {
				cbar_pt  = pythia.event[i].pT();
				cbar_phi = pythia.event[i].phi();
				cbar_eta = pythia.event[i].eta();
				n_cbar++;
			}
                        //Identify charmed hadrons. (Kinematics determined after event selection)
			if (id ==  chadron_id) { n_chad++; chad_ind = i;}	
                        if (id == -chadron_id) { n_cbarhad++; cbarhad_ind = i;} 

			//d1 = pi+
			if (id == d1_id && (mother1 == chad_ind || mother2 == chad_ind)) {
				d1_pt  = pythia.event[i].pT();
				d1_phi = pythia.event[i].phi();
				d1_eta = pythia.event[i].eta();
				n_d1++;
			}
			//d1bar = pi-
			if (id == -d1_id && (mother1 == cbarhad_ind || mother2 == cbarhad_ind)) {
				d1bar_pt  = pythia.event[i].pT();
				d1bar_phi = pythia.event[i].phi();
				d1bar_eta = pythia.event[i].eta();
				n_d1bar++;
			}
			//d2 = K+
			if (id == d2_id && (mother1 == cbarhad_ind || mother2 == cbarhad_ind)) {
				d2_pt  = pythia.event[i].pT();
				d2_phi = pythia.event[i].phi();
				d2_eta = pythia.event[i].eta();
				n_d2++;
			}
			//d2bar = K-
			if (id == -d2_id && (mother1 == chad_ind || mother2 == chad_ind)) {
				d2bar_pt  = pythia.event[i].pT();
				d2bar_phi = pythia.event[i].phi();
				d2bar_eta = pythia.event[i].eta();
				n_d2bar++;
			}
		}
		
		//choose only events with 1 of each particle, as generation will ensure correct matching
		if (n_chad == 1 && n_cbarhad == 1 && n_c == 1 && n_cbar == 1 && n_d1 == 1 && n_d1bar == 1 && n_d2 == 1 && n_d2bar == 1) {
			
			chad_pt  = pythia.event[chad_ind].pT();  cbarhad_pt  = pythia.event[cbarhad_ind].pT();
			chad_phi = pythia.event[chad_ind].phi(); cbarhad_phi = pythia.event[cbarhad_ind].phi();
			chad_eta = pythia.event[chad_ind].eta(); cbarhad_eta = pythia.event[cbarhad_ind].eta();
						
			//Find total pT and number of final particles inside cone of radius R around each D meson
			const double R = 0.7;
			for (int i = 0; i < pythia.event.size(); ++i) {
				if (!pythia.event[i].isFinal()) continue; //only look at final particles	
                                double eta = pythia.event[i].eta(), phi = pythia.event[i].phi(), pt = pythia.event[i].pT();
                                double chad_deltaR = std::pow( std::pow(phi - chad_phi, 2) + std::pow(eta - chad_eta, 2) , 0.5); //Distance from c hadron
                                double cbarhad_deltaR = std::pow( std::pow(phi - cbarhad_phi, 2) + std::pow(eta - cbarhad_eta, 2) , 0.5); //Distance from cbar hadron

                                if (chad_deltaR < R) {chad_cone_mult++; chad_ptcone += pt;}
                                if (cbarhad_deltaR < R) {cbarhad_cone_mult++; cbarhad_ptcone += pt;}
                        }
			//write to datafile
			datafile <<cLO_pt<<","<<cLO_phi<<","<<cLO_eta<<","<<c_pt<<","<<c_phi<<","<<c_eta<<","<<chad_pt<<","<<chad_phi<<","<<chad_eta<<","
				 <<cbarLO_pt<<","<<cbarLO_phi<<","<<cbarLO_eta<<","<<cbar_pt<<","<<cbar_phi<<","<<cbar_eta<<","<<cbarhad_pt<<","<<cbarhad_phi<<","<<cbarhad_eta<<","
				 <<d1_pt<<","<<d1_phi<<","<<d1_eta<<","<<d1bar_pt<<","<<d1bar_phi<<","<<d1bar_eta<<","
				 <<d2_pt<<","<<d2_phi<<","<<d2_eta<<","<<d2bar_pt<<","<<d2bar_phi<<","<<d2bar_eta<<","
				 <<mult<<","<<chad_cone_mult<<","<<chad_ptcone<<","<<cbarhad_cone_mult<<","<<cbarhad_ptcone<<"\n";

		}
	// End of event loop.
	}
	datafile.close();

	pythia.stat();
		
	return 0;
}
	

