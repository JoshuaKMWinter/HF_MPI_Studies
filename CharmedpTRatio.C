#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLegend.h>

void CharmedpTRatio() {
	TFile *outFile = new TFile("charm_hadptratio.root","RECREATE");
	TTree *tree = new TTree("data","charmed hadron pT");
	tree->ReadFile("pThad_Data.csv","",',');	
	
	TTreeReader reader("data");
	TTreeReaderValue<float> had_pt(reader,"pt_had");
	TTreeReaderValue<float> mult(reader,"multiplicity");
	TTreeReaderValue<float> had_pdg(reader,"hadron_pdg");
	
	float bins[] = {0,1,2,3,5,10,20,30};

	TH1F *h_D0_pt_bin1 = new TH1F("h_D0_pt_bin1","",7,bins);
	TH1F *h_D0_pt_bin2 = new TH1F("h_D0_pt_bin2","",7,bins);
	TH1F *h_D0_pt_bin3 = new TH1F("h_D0_pt_bin3","",7,bins);
	TH1F *h_Ds_pt_bin1 = new TH1F("h_Ds_pt_bin1","",7,bins);
	TH1F *h_Ds_pt_bin2 = new TH1F("h_Ds_pt_bin2","",7,bins);
	TH1F *h_Ds_pt_bin3 = new TH1F("h_Ds_pt_bin3","",7,bins);
	TH1F *h_Lc_pt_bin1 = new TH1F("h_Lc_pt_bin1","",7,bins);
	TH1F *h_Lc_pt_bin2 = new TH1F("h_Lc_pt_bin2","",7,bins);
	TH1F *h_Lc_pt_bin3 = new TH1F("h_Lc_pt_bin3","",7,bins);
	
	h_D0_pt_bin1->Sumw2();
	h_D0_pt_bin2->Sumw2();
	h_D0_pt_bin3->Sumw2();
	h_Ds_pt_bin1->Sumw2();
	h_Ds_pt_bin2->Sumw2();
	h_Ds_pt_bin3->Sumw2();
	h_Lc_pt_bin1->Sumw2();
	h_Lc_pt_bin2->Sumw2();
	h_Lc_pt_bin3->Sumw2();

	while(reader.Next()) {
		if (*had_pdg == 421) {
			if (*mult < 100) h_D0_pt_bin1->Fill(*had_pt);
			else if (*mult < 300) h_D0_pt_bin2->Fill(*had_pt);
			else h_D0_pt_bin3->Fill(*had_pt);
		}
		else if (*had_pdg == 431) {
			if (*mult < 100) h_Ds_pt_bin1->Fill(*had_pt);
			else if (*mult < 300) h_Ds_pt_bin2->Fill(*had_pt);
			else h_Ds_pt_bin3->Fill(*had_pt);
		}
		else if (*had_pdg == 4122){
			if (*mult < 100) h_Lc_pt_bin1->Fill(*had_pt);
			else if (*mult < 300) h_Lc_pt_bin2->Fill(*had_pt);
			else h_Lc_pt_bin3->Fill(*had_pt);
		}
	}
	
	TH1F *h_Dsratio_bin1 = (TH1F*)h_Ds_pt_bin1->Clone("h_Dsratio_bin1");
	TH1F *h_Dsratio_bin2 = (TH1F*)h_Ds_pt_bin2->Clone("h_Dsratio_bin2");
	TH1F *h_Dsratio_bin3 = (TH1F*)h_Ds_pt_bin3->Clone("h_Dsratio_bin3");
	TH1F *h_Lcratio_bin1 = (TH1F*)h_Lc_pt_bin1->Clone("h_Lcratio_bin1");
	TH1F *h_Lcratio_bin2 = (TH1F*)h_Lc_pt_bin2->Clone("h_Lcratio_bin2");
	TH1F *h_Lcratio_bin3 = (TH1F*)h_Lc_pt_bin3->Clone("h_Lcratio_bin3");
	
	h_Dsratio_bin1->Divide(h_D0_pt_bin1);
	h_Dsratio_bin2->Divide(h_D0_pt_bin2);
	h_Dsratio_bin3->Divide(h_D0_pt_bin3);
	h_Lcratio_bin1->Divide(h_D0_pt_bin1);
	h_Lcratio_bin2->Divide(h_D0_pt_bin2);
	h_Lcratio_bin3->Divide(h_D0_pt_bin3);

	h_Dsratio_bin1->SetLineColor(kBlue);
	h_Dsratio_bin2->SetLineColor(kBlack);
	h_Dsratio_bin3->SetLineColor(kRed);
	h_Lcratio_bin1->SetLineColor(kBlue);
	h_Lcratio_bin2->SetLineColor(kBlack);
	h_Lcratio_bin3->SetLineColor(kRed);

	TCanvas *c1 = new TCanvas("c1","Ds to D0 pt ratio");
	TCanvas *c2 = new TCanvas("c2","Lc to D0 pt ratio");

	TLegend *leg = new TLegend();
	leg->SetHeader("Event multiplicity ranges","C");
	leg->AddEntry(h_Dsratio_bin1,"[0,100)");
	leg->AddEntry(h_Dsratio_bin2,"[100,300)");
	leg->AddEntry(h_Dsratio_bin3,"300+");
	
	c1->cd();
	//h_Dsratio_bin1->SetMaximum(0.5);
	h_Dsratio_bin1->Draw("e1");
	h_Dsratio_bin2->Draw("e1 same");
	h_Dsratio_bin3->Draw("e1 same");
	leg->Draw();
	c2->cd();
	h_Lcratio_bin1->SetMaximum(0.15);
	h_Lcratio_bin1->Draw("e1");
	h_Lcratio_bin2->Draw("e1 same");
	h_Lcratio_bin3->Draw("e1 same");
	leg->Draw();

	h_D0_pt_bin1->Write();
	h_D0_pt_bin2->Write();
	h_D0_pt_bin3->Write();
	h_Ds_pt_bin1->Write();
	h_Ds_pt_bin2->Write();
	h_Ds_pt_bin3->Write();
	h_Lc_pt_bin1->Write();
	h_Lc_pt_bin2->Write();
	h_Lc_pt_bin3->Write();
	h_Dsratio_bin1->Write();
	h_Dsratio_bin2->Write();
	h_Dsratio_bin3->Write();
	h_Lcratio_bin1->Write();
	h_Lcratio_bin2->Write();
	h_Lcratio_bin3->Write();
	c1->Write();
	c2->Write();

	outFile->Close();
}
