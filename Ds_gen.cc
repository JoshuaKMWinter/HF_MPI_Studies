//Ds_gen.cc
//Code for generation of ccbar pairs, with selection of events based on the confinement of c-quarks in Ds mesons
//Following Ds+ > pi+ phi, phi > K+ K- and complex conjugate for cbar.

#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <string> 
#include <cmath> 

using namespace Pythia8;

int main(int argc, char* argv[]) {

	const int chadron_id = 431; 	//PDG ID of charm meson for selected decay
	const int d1_id = 211;		//Decay product ids. pi+
	const int d2_id = 333;		//phi0
	const int d3_id = 321;		//K+
	const int d4_id = 321;		//K+

        char* outfile;

	if (argc == 1) outfile = (char*)"CSVs/Ds_EventData.csv";
        else if (argc !=3) {cout<<"Wrong number of arguments. One output file and on pythia tune expected. Program stopped."<<endl; return(1);}
        else {
                outfile = argv[1];
                ifstream is(argv[2]);
                if (!is) {cout<<"Pythia tune not found! Program stopped."; return(1);}
        }

	std::ofstream datafile; datafile.open(outfile);
	datafile <<"cLO_pT,cLO_phi,cLO_eta,cLO_rap,c_pT,c_phi,c_eta,c_rap,chad_pT,chad_phi,chad_eta,chad_rap,cbarLO_pT,cbarLO_phi,cbarLO_eta,cbarLO_rap,cbar_pT,cbar_phi,cbar_eta,cbar_rap,cbarhad_pT,cbarhad_phi,cbarhad_eta,cbarhad_rap,d1_pT,d1_phi,d1_eta,d1_rap,d1bar_pT,d1bar_phi,d1bar_eta,d1bar_rap,d2_pT,d2_phi,d2_eta,d2_rap,d2bar_pT,d2bar_phi,d2bar_eta,d2bar_rap,d3_pT,d3_phi,d3_eta,d3_rap,d3bar_pT,d3bar_phi,d3bar_eta,d3bar_rap,d4_pT,d4_phi,d4_eta,d4_rap,d4bar_pT,d4bar_phi,d4bar_eta,d4bar_rap,multiplicity,chad_cone_mult,chad_ptcone,cbarhad_cone_mult,cbarhad_ptcone\n";

	//Set-up event properties
	Pythia pythia;

	//Tune of Pythia
        if (argc ==3) {
                pythia.readFile(argv[2]);
                cout<<"Reading Pythia Tune: "<<argv[2]<<endl;
        }

	pythia.readString("Beams:eCM = 13000.");      //set collision centre of mass energy
	pythia.readString("HardQCD:hardccbar = on");  //Select ccbar production
	pythia.readString("431:onMode = off");        //Turn off decay modes for selected charmed hadron
	pythia.readString("431:onIfMatch = 211 333"); //Selected desired decays for charmed hadron
	pythia.readString("333:onMode = off");
	pythia.readString("333:onIfMatch = 321 321"); //Force phi > K+ K- decay
	//pythia.readString("PhaseSpace:pTHatMin = 4.");
	pythia.init();
	
	// Begin event loop. Generate event. Skip if error. List first one.
	for (int iEvent = 0; iEvent < 400000; ++iEvent) { 
		if (!pythia.next()) continue;
		//if(iEvent == 419) pythia.event.list();

		double cLO_pt = 0,  cLO_phi = 0,  cLO_eta = 0,  cLO_rap = 0,  cbarLO_pt = 0,  cbarLO_phi = 0,  cbarLO_eta = 0,  cbarLO_rap = 0; //Leading Order charm kinematics
                double c_pt = 0,    c_phi = 0,    c_eta = 0,    c_rap = 0,    cbar_pt = 0,    cbar_phi = 0,    cbar_eta = 0,    cbar_rap = 0;      //Charm kin. at hadronisation
                double chad_pt = 0, chad_phi = 0, chad_eta = 0, chad_rap = 0, cbarhad_pt = 0, cbarhad_phi = 0, cbarhad_eta = 0, cbarhad_rap = 0;//Charm kin. in hadron
                double d1_pt = 0,   d1_phi = 0,   d1_eta = 0,   d1_rap = 0,   d1bar_pt = 0,   d1bar_phi = 0,   d1bar_eta = 0,   d1bar_rap = 0;  //kinematics of daughter particles
                double d2_pt = 0,   d2_phi = 0,   d2_eta = 0,   d2_rap = 0,   d2bar_pt = 0,   d2bar_phi = 0,   d2bar_eta = 0,   d2bar_rap = 0;
		double d3_pt = 0,   d3_phi = 0,   d3_eta = 0,   d3_rap = 0,   d3bar_pt = 0,  d3bar_phi = 0,  d3bar_eta = 0,  d3bar_rap = 0;  
		double d4_pt = 0,   d4_phi = 0,   d4_eta = 0,   d4_rap = 0,   d4bar_pt = 0,  d4bar_phi = 0,  d4bar_eta = 0,  d4bar_rap = 0;  
		double chad_ptcone = 0, cbarhad_ptcone = 0;		

		int n_chad = 0, n_cbarhad = 0, n_c = 0, n_cbar = 0, n_d1 = 0, n_d1bar = 0, n_d2 = 0, n_d2bar = 0; //counters for number of selected particles
		int chad_ind = -1, cbarhad_ind = -1, cphi_ind = -1, cbarphi_ind = -1; //integers to track charm meson index
		int mult = 0, chad_cone_mult = 0, cbarhad_cone_mult = 0; //multiplicity counters.
	
                int cind = -1, cbarind = -1;

		for (int i = 0; i < pythia.event.size(); ++i) {
			
			//use all final particles as an estimate of multiplicity
			if (pythia.event[i].isFinal()) mult++;

			//find status and id of current particle and the indices and ids of mother particles
			int status = pythia.event[i].status(),   id = pythia.event[i].id();
			int mother1 = pythia.event[i].mother1(), mother2 = pythia.event[i].mother2();

			////////////////////////////////////Obtain particle kinematic properties////////////////////////////////////
			//Charm Leading Order kinematics
			if (id == 4 && status == -23) { //outgoing from hardest subprocess
				cLO_pt  = pythia.event[i].pT();
				cLO_phi = pythia.event[i].phi();
				cLO_eta = pythia.event[i].eta();
				cLO_rap = pythia.event[i].y();
				cind = i;
				while (pythia.event[cind].id() == 4) {
                                        if (pythia.event[pythia.event[cind].daughter1()].id() == 4) cind = pythia.event[cind].daughter1();
                                        else break;
                                }
			} 
			if (id == -4 && status == -23) {
				cbarLO_pt  = pythia.event[i].pT();
				cbarLO_phi = pythia.event[i].phi();
				cbarLO_eta = pythia.event[i].eta();
				cbarLO_rap = pythia.event[i].y();
				cbarind = i;
				while (pythia.event[cbarind].id() == -4) {
                                        if (pythia.event[pythia.event[cbarind].daughter1()].id() == -4) cbarind = pythia.event[cbarind].daughter1();
                                        else break;
                                }
			}

			//Charm kinematics prior to hadronisation
			if (id == 4 && status/10 == -7 && cind == i) { //71-79 preparation for hadronisation
				c_pt  = pythia.event[i].pT();
				c_phi = pythia.event[i].phi();
				c_eta = pythia.event[i].eta();
				c_rap = pythia.event[i].y();
				n_c++;
			} 
			if (id == -4 && status/10 == -7 && cbarind == i) {
				cbar_pt  = pythia.event[i].pT();
				cbar_phi = pythia.event[i].phi();
				cbar_eta = pythia.event[i].eta();
				cbar_rap = pythia.event[i].y();
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
				d1_rap = pythia.event[i].y();
				n_d1++;
			}
			//d1bar = pi-
			if (id == -d1_id && (mother1 == cbarhad_ind || mother2 == cbarhad_ind)) {
				d1bar_pt  = pythia.event[i].pT();
				d1bar_phi = pythia.event[i].phi();
				d1bar_eta = pythia.event[i].eta();
				d1bar_rap = pythia.event[i].y();
				n_d1bar++;
			}
			//d2 = phi from Ds+
			if (id == d2_id && (mother1 == chad_ind || mother2 == chad_ind)) { 
				cphi_ind = i;
				d2_pt  = pythia.event[i].pT();
				d2_phi = pythia.event[i].phi();
				d2_eta = pythia.event[i].eta();
				d2_rap = pythia.event[i].y();
				n_d2++;
			}
			//d2bar = phi from Ds-
			if (id == d2_id && (mother1 == cbarhad_ind || mother2 == cbarhad_ind)) { 
				cbarphi_ind = i;
				d2bar_pt  = pythia.event[i].pT();
				d2bar_phi = pythia.event[i].phi();
				d2bar_eta = pythia.event[i].eta();
				d2bar_rap = pythia.event[i].y();
				n_d2bar++;
			}
			//d3 = K+ from d2
			if (id == d3_id && (mother1 == cphi_ind || mother2 == cphi_ind)) {
				d3_pt  = pythia.event[i].pT();
				d3_phi = pythia.event[i].phi();
				d3_eta = pythia.event[i].eta();
				d3_rap = pythia.event[i].y();
			}
			//d3bar = K- from d2
			if (id == -d3_id && (mother1 == cphi_ind || mother2 == cphi_ind)) {
				d3bar_pt  = pythia.event[i].pT();
				d3bar_phi = pythia.event[i].phi();
				d3bar_eta = pythia.event[i].eta();
				d3bar_rap = pythia.event[i].y();
			}
			//d4 = K+ from d2bar
			if (id == d4_id && (mother1 == cbarphi_ind || mother2 == cbarphi_ind)) {
				d4_pt  = pythia.event[i].pT();
				d4_phi = pythia.event[i].phi();
				d4_eta = pythia.event[i].eta();
				d4_rap = pythia.event[i].y();
			}
			//d4bar = K- from d2bar
			if (id == -d4_id && (mother1 == cbarphi_ind || mother2 == cbarphi_ind)) {
				d4bar_pt  = pythia.event[i].pT();
				d4bar_phi = pythia.event[i].phi();
				d4bar_eta = pythia.event[i].eta();
				d4bar_rap = pythia.event[i].y();
			}
		}
	
		//choose only events with 1 of each particle, as generation will ensure correct matching
		if (n_chad == 1 && n_cbarhad == 1 && n_c == 1 && n_cbar == 1 && n_d1 == 1 && n_d1bar == 1 && n_d2 == 1 && n_d2bar == 1) {
			//pythia.event.list();

			chad_pt  = pythia.event[chad_ind].pT();  cbarhad_pt  = pythia.event[cbarhad_ind].pT();
			chad_phi = pythia.event[chad_ind].phi(); cbarhad_phi = pythia.event[cbarhad_ind].phi();
			chad_eta = pythia.event[chad_ind].eta(); cbarhad_eta = pythia.event[cbarhad_ind].eta();
                        chad_rap = pythia.event[chad_ind].y();   cbarhad_rap = pythia.event[cbarhad_ind].y();

 			//Find total pT and number of final particles inside cone of radius R around each charmed hadron
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
                        datafile <<cLO_pt<<","<<cLO_phi<<","<<cLO_eta<<","<<cLO_rap<<","<<c_pt<<","<<c_phi<<","<<c_eta<<","<<c_rap<<","<<chad_pt<<","<<chad_phi<<","<<chad_eta<<","<<chad_rap<<","
                                 <<cbarLO_pt<<","<<cbarLO_phi<<","<<cbarLO_eta<<","<<cbarLO_rap<<","<<cbar_pt<<","<<cbar_phi<<","<<cbar_eta<<","<<cbar_rap<<","<<cbarhad_pt<<","<<cbarhad_phi<<","<<cbarhad_eta<<","<<cbarhad_rap<<","
                                 <<d1_pt<<","<<d1_phi<<","<<d1_eta<<","<<d1_rap<<","<<d1bar_pt<<","<<d1bar_phi<<","<<d1bar_eta<<","<<d1bar_rap<<","
                                 <<d2_pt<<","<<d2_phi<<","<<d2_eta<<","<<d2_rap<<","<<d2bar_pt<<","<<d2bar_phi<<","<<d2bar_eta<<","<<d2bar_rap<<"," 
				 <<d3_pt<<","<<d3_phi<<","<<d3_eta<<","<<d3_rap<<","<<d3bar_pt<<","<<d3bar_phi<<","<<d3bar_eta<<","<<d3bar_rap<<","
				 <<d4_pt<<","<<d4_phi<<","<<d4_eta<<","<<d4_rap<<","<<d4bar_pt<<","<<d4bar_phi<<","<<d4bar_eta<<","<<d4bar_rap<<","
				 <<mult<<","<<chad_cone_mult<<","<<chad_ptcone<<","<<cbarhad_cone_mult<<","<<cbarhad_ptcone<<"\n";

		}

	// End of event loop.
	}
	datafile.close();

	pythia.stat();
		
	return 0;
}
	

