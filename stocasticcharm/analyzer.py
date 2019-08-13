#############################################################################
##  © Copyright CERN 2018. All rights not expressly granted are reserved.  ##
##                 Author: Gian.Michele.Innocenti@cern.ch                  ##
## This program is free software: you can redistribute it and/or modify it ##
##  under the terms of the GNU General Public License as published by the  ##
## Free Software Foundation, either version 3 of the License, or (at your  ##
## option) any later version. This program is distributed in the hope that ##
##  it will be useful, but WITHOUT ANY WARRANTY; without even the implied  ##
##     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    ##
##           See the GNU General Public License for more details.          ##
##    You should have received a copy of the GNU General Public License    ##
##   along with this program. if not, see <https://www.gnu.org/licenses/>. ##
#############################################################################
import os
# pylint: disable=unused-wildcard-import, wildcard-import
from array import *
# pylint: disable=import-error, no-name-in-module, unused-import
import pandas as pd
import pickle
from root_numpy import fill_hist
from ROOT import TFile, TH1F, TCanvas
from ROOT import gStyle, TLegend
from ROOT import gROOT
from ROOT import TStyle
import numpy as np

def Delta_Phi(phi1, phi2):
    pi = 3.14159265358979323846
    dphi = abs(phi1 - phi2)
    return np.where(dphi > pi, 2*pi - dphi, dphi)
def Delta_R(phi1, eta1, phi2, eta2):
    dphi = Delta_Phi(phi1,phi2)
    deta = abs( eta1 - eta2 )
    return pow( dphi*dphi + deta*deta, 0.5 )

def Find_Hist_Max(histlist):
    hmax = 0
    for i in range(len(histlist)):
        if histlist[i].GetMaximum() > hmax : hmax = histlist[i].GetMaximum()
    return hmax

# pylint: disable=too-few-public-methods, too-many-instance-attributes, too-many-statements
class Analyzer:
    species = "analyzer"
    def __init__(self, datap, case):
        self.case = case
        self.mass = datap["mass"]
        self.pdg = datap["pdg"]
        self.rapiditymax = datap["rapiditymax"]
        self.rapiditymin = datap["rapiditymin"]
        self.var_all = datap["variables"]["var_all"]
        self.inputfile = datap["fileinput"]
        
        self.sel_1d_ptchad = datap["sel_1d_ptchad"]

        self.var_1d_distr   = datap["var_1d_distr"]
        self.nbins_1d_distr = datap["nbins_1d_distr"]
        self.min_1d_distr   = datap["min_1d_distr"]
        self.max_1d_distr   = datap["max_1d_distr"]
        self.title_1d_distr = datap["title_1d_distr"]
        self.leg_1d_distr   = datap["leg_1d_distr"]
        self.logy_1d_distr  = datap["logy_1d_distr"]
        self.n_1d = len(self.max_1d_distr) 
        self.colours = datap["colours"]

        self.dfm = None
        self.dfm = pd.read_csv(self.inputfile)
        self.add_derived()

    def add_derived(self):
        if self.case == "D0" or "Ds" or "Lc":
            self.dfm["dphi_cc"] = Delta_Phi( self.dfm["c_phi"], self.dfm["cbar_phi"] )
            self.dfm["dphi_ccLO"] = Delta_Phi( self.dfm["cLO_phi"], self.dfm["cbarLO_phi"] )
            self.dfm["dphi_cchad"] = Delta_Phi( self.dfm["chad_phi"], self.dfm["cbarhad_phi"] )
            self.dfm["dphi_chad"] = Delta_Phi( self.dfm["c_phi"], self.dfm["chad_phi"] )
            self.dfm["dphi_cbarhad"] = Delta_Phi( self.dfm["cbar_phi"], self.dfm["cbarhad_phi"] )
            self.dfm["deta_cc"] = abs(self.dfm["c_eta"] - self.dfm["cbar_eta"])
            self.dfm["deta_ccLO"] = abs(self.dfm["cLO_eta"] - self.dfm["cbarLO_eta"])
            self.dfm["deta_cchad"] = abs(self.dfm["chad_eta"] - self.dfm["cbarhad_eta"])
            self.dfm["deta_chad"] = abs(self.dfm["c_eta"] - self.dfm["chad_eta"])
            self.dfm["deta_cbarhad"] = abs(self.dfm["cbar_eta"] - self.dfm["cbarhad_eta"])
            self.dfm["FFc"] = self.dfm["chad_pT"] / self.dfm["chad_ptcone"]
            self.dfm["FFcbar"] = self.dfm["cbarhad_pT"] / self.dfm["cbarhad_ptcone"]
            self.dfm["cLO_ptratio"] = self.dfm["cLO_pT"] / self.dfm["cbarLO_pT"]
            self.dfm["c_ptratio"] = self.dfm["c_pT"] / self.dfm["cbar_pT"]
            self.dfm["chad_ptratio"] = self.dfm["chad_pT"] / self.dfm["cbarhad_pT"]

    def plot(self):
        print("Running plotter")
       
        outfile = TFile("data/" + self.case + "_hists.root", "RECREATE")
       
        dfsel = self.dfm.query(self.sel_1d_ptchad)
        
        for index in range(self.n_1d):    
            c = TCanvas('c_'+self.var_1d_distr[index][0])
            legend = TLegend(.77,.6,.98,.75)
            histlist = []
            for i in range( len(self.var_1d_distr[index]) ):
                hist = TH1F("h_"+self.var_1d_distr[index][i],
                            "h_"+self.var_1d_distr[index][i],
                            self.nbins_1d_distr[index],
                            self.min_1d_distr[index],
                            self.max_1d_distr[index])
                hist.Sumw2()
                hist.SetLineColor(self.colours[i])
                fill_hist(hist, dfsel[self.var_1d_distr[index][i]])
                histlist.append(hist)
                hist.Write()
                legend.AddEntry(hist, self.leg_1d_distr[index][i])
            hmax = Find_Hist_Max(histlist) 
            if self.logy_1d_distr[index] :
                c.SetLogy()
                histlist[0].SetMinimum(0.5)
            histlist[0].SetMaximum(hmax*1.2)
            histlist[0].SetTitle(self.case + ": " + self.title_1d_distr[index] )
            histlist[0].Draw()
            for ihist in range( len(histlist)-1 ):
                histlist[ihist+1].Draw("SAME")
            legend.Draw()
            c.SaveAs("plots/%s/c_%s_%s.eps" % (self.case, self.case, self.var_1d_distr[index][0]))
        outfile.Close()

    def hadron_ptratio(self):
        self.df = pd.read_csv("../input/pThad_Data.csv")

        multbins = [self.df['multiplicity']<100, self.df['multiplicity']>=100 and self.df['multiplicity']<300, self.df['multiplicity']>=300]

        dfD0 = self.df.query(hadron_pdg == 421)
        df_chad = self.df.query(hadron_pdg == self.pdg)
        
        D0_pt = [[],[],[]]
        chad_pt = [[],[],[]]

        for i in range(len(multbins)):
            D0_pt[i] = [x for x in dfD0['pt_had'] if multbins[i]]
            chad_pt[i] = [x for x in df_chad['pt_had'] if multbins[i]]
        
        print(D0_pt)
