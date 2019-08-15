#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <string> 
#include <cmath>
#include <vector>

using namespace Pythia8;

int main() {

	//outfile for multiplicity analysis
	//std::ofstream multfile; multfile.open("CSVs/NDs_multData.csv");
	//multfile << "num_D,num_Dbar,num_Ds,num_Dsbar,num_lambdac,num_lambdacbar,D_leadpt,Ds_leadpt,Lc_leadpt,multiplicity\n";

	//outfile for pt of all hadrons analysis
	std::ofstream ptfile; ptfile.open("CSVs/pThad_Data_WithMPI_NoCR.csv");
	ptfile << "pt_had,multiplicity,hadron_pdg\n"; 

	Pythia pythia;

	pythia.readFile("../PythiaTunes/WithMPI_NoCR.cmnd");

	pythia.readString("Beams:eCM = 13000.");
	pythia.readString("HardQCD:hardccbar = on");

	pythia.init();
	
	// Begin event loop. Generate event. Skip if error. List first one.
	for (int iEvent = 0; iEvent < 500000; ++iEvent) { //100k
		if (!pythia.next()) continue;
		
		//iint num_D = 0, num_Dbar = 0, num_Ds=0, num_Dsbar=0, num_lambdac=0, num_lambdacbar=0;
		
		int multiplicity = 0;
		
		vector<double> D0_ptvec, Ds_ptvec, Lc_ptvec;

		for (int i = 0; i < pythia.event.size(); ++i) {
			
			//use all final particles as an estimate of multiplicity
			if (pythia.event[i].isFinal()) multiplicity++;

			//find status and id of current particle and the indices and ids of mother particles
			int status = pythia.event[i].status(),   id = pythia.event[i].id();

			if (abs(id) ==   421) { 
				//num_D++; D_ind = i;
				D0_ptvec.push_back(pythia.event[i].pT());
			}
			if (abs(id) ==   431)  {
				//num_Ds++;
				Ds_ptvec.push_back(pythia.event[i].pT());
			}
			if (abs(id) ==  4122)  {
			 	//num_lambdac++;
				Lc_ptvec.push_back(pythia.event[i].pT());
			}	

		}
		for (size_t i=0; i<D0_ptvec.size(); i++) ptfile << D0_ptvec[i] <<","<< multiplicity <<","<< 421 << "\n";
		for (size_t i=0; i<Ds_ptvec.size(); i++) ptfile << Ds_ptvec[i] <<","<< multiplicity <<","<< 431 << "\n";
		for (size_t i=0; i<Lc_ptvec.size(); i++) ptfile << Lc_ptvec[i] <<","<< multiplicity <<","<< 4122 << "\n";
		
		//multfile << num_D << "," <<  num_Dbar << "," <<  num_Ds << "," <<  num_Dsbar << "," <<  num_lambdac << "," <<  num_lambdacbar<< "," <<  multiplicity <<"\n";

	// End of event loop.
	}
	ptfile.close();
	//multfile.close();
	//datafile.close();

	pythia.stat();
		
	return 0;
}
	

