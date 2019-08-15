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
import copy

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
        self.colours = datap["colours"]
        
        self.sel_1d_ptchad = datap["sel_1d_ptchad"]

        self.var_1d_distr   = datap["var_1d_distr"]
        self.nbins_1d_distr = datap["nbins_1d_distr"]
        self.min_1d_distr   = datap["min_1d_distr"]
        self.max_1d_distr   = datap["max_1d_distr"]
        self.title_1d_distr = datap["title_1d_distr"]
        self.leg_1d_distr   = datap["leg_1d_distr"]
        self.xlab_1d_distr  = datap["xlab_1d_distr"]
        self.ylab_1d_distr  = datap["ylab_1d_distr"]
        self.logy_1d_distr  = datap["logy_1d_distr"]
        self.n_1d = len(self.max_1d_distr) 
        self.dfm = None
        #self.dfm = pd.read_csv(self.inputfile)
        #self.add_derived(self.dfm)

        self.tune_fileinput   = datap["tune_fileinput"]
        self.var_tune_distr   = datap["var_tune_distr"]
        self.nbins_tune_distr = datap["nbins_tune_distr"]
        self.min_tune_distr   = datap["min_tune_distr"]
        self.max_tune_distr   = datap["max_tune_distr"]
        self.title_tune_distr = datap["title_tune_distr"]
        self.leg_tune_distr   = datap["leg_tune_distr"]
        self.xlab_tune_distr  = datap["xlab_tune_distr"]
        self.ylab_tune_distr  = datap["ylab_tune_distr"]

        self.dfpt = None
        self.pthad_fileinput  = datap["pthad_fileinput"]
        self.multi_bins = datap["multi_bins"]

#    def add_derived(self): #change this function so it can be applied to ANY dataframe
#        if self.case == "D0" or "Ds" or "Lc":
#            self.dfm["dphi_cc"] = Delta_Phi( self.dfm["c_phi"], self.dfm["cbar_phi"] )
#            self.dfm["dphi_ccLO"] = Delta_Phi( self.dfm["cLO_phi"], self.dfm["cbarLO_phi"] )
#            self.dfm["dphi_cchad"] = Delta_Phi( self.dfm["chad_phi"], self.dfm["cbarhad_phi"] )
#            self.dfm["dphi_chad"] = Delta_Phi( self.dfm["c_phi"], self.dfm["chad_phi"] )
#            self.dfm["dphi_cbarhad"] = Delta_Phi( self.dfm["cbar_phi"], self.dfm["cbarhad_phi"] )
#            self.dfm["deta_cc"] = abs(self.dfm["c_eta"] - self.dfm["cbar_eta"])
#            self.dfm["deta_ccLO"] = abs(self.dfm["cLO_eta"] - self.dfm["cbarLO_eta"])
#            self.dfm["deta_cchad"] = abs(self.dfm["chad_eta"] - self.dfm["cbarhad_eta"])
#            self.dfm["deta_chad"] = abs(self.dfm["c_eta"] - self.dfm["chad_eta"])
#            self.dfm["deta_cbarhad"] = abs(self.dfm["cbar_eta"] - self.dfm["cbarhad_eta"])
#            self.dfm["FFc"] = self.dfm["chad_pT"] / self.dfm["chad_ptcone"]
#            self.dfm["FFcbar"] = self.dfm["cbarhad_pT"] / self.dfm["cbarhad_ptcone"]
#            self.dfm["cLO_ptratio"] = self.dfm["cLO_pT"] / self.dfm["cbarLO_pT"]
#            self.dfm["c_ptratio"] = self.dfm["c_pT"] / self.dfm["cbar_pT"]
#            self.dfm["chad_ptratio"] = self.dfm["chad_pT"] / self.dfm["cbarhad_pT"]
    def add_derived(self,df): #change this function so it can be applied to ANY dataframe
        if self.case == "D0" or "Ds" or "Lc":
            df["dphi_cc"] = Delta_Phi( df["c_phi"], df["cbar_phi"] )
            df["dphi_ccLO"] = Delta_Phi( df["cLO_phi"], df["cbarLO_phi"] )
            df["dphi_cchad"] = Delta_Phi( df["chad_phi"], df["cbarhad_phi"] )
            df["dphi_chad"] = Delta_Phi( df["c_phi"], df["chad_phi"] )
            df["dphi_cbarhad"] = Delta_Phi( df["cbar_phi"], df["cbarhad_phi"] )
            df["deta_cc"] = abs(df["c_eta"] - df["cbar_eta"])
            df["deta_ccLO"] = abs(df["cLO_eta"] - df["cbarLO_eta"])
            df["deta_cchad"] = abs(df["chad_eta"] - df["cbarhad_eta"])
            df["deta_chad"] = abs(df["c_eta"] - df["chad_eta"])
            df["deta_cbarhad"] = abs(df["cbar_eta"] - df["cbarhad_eta"])
            df["FFc"] = df["chad_pT"] / df["chad_ptcone"]
            df["FFcbar"] = df["cbarhad_pT"] / df["cbarhad_ptcone"]
            df["cLO_ptratio"] = df["cLO_pT"] / df["cbarLO_pT"]
            df["c_ptratio"] = df["c_pT"] / df["cbar_pT"]
            df["chad_ptratio"] = df["chad_pT"] / df["cbarhad_pT"]

    def plot(self):
        print("Running analyzer")
        
        for ifile in range(len(self.tune_fileinput)):
            outfile = TFile("data/" + self.case+"_"+self.leg_tune_distr[ifile] + "_hists.root", "RECREATE")
            
            self.dfm = pd.read_csv(self.tune_fileinput[ifile])
            self.add_derived(self.dfm)
            dfsel = self.dfm.query(self.sel_1d_ptchad)
            
            for index in range(self.n_1d):    
                c = TCanvas('c_'+self.var_1d_distr[index][0],'',600,600)
                legend = TLegend(.7,.75,.9,.9)
                histlist = []
                for i in range( len(self.var_1d_distr[index]) ):
                    hist = TH1F("h_"+self.var_1d_distr[index][i],
                                "h_"+self.var_1d_distr[index][i],
                                self.nbins_1d_distr[index],
                                self.min_1d_distr[index],
                                self.max_1d_distr[index])
                    hist.Sumw2()
                    hist.SetLineColor(self.colours[i])
                    hist.SetTitle(self.case + ": " + self.title_1d_distr[index] )
                    hist.SetXTitle(self.xlab_1d_distr[index])
                    hist.SetYTitle(self.ylab_1d_distr[index])
                    hist.SetStats(0)
                    fill_hist(hist, dfsel[self.var_1d_distr[index][i]])
                    histlist.append(hist)
                    hist.Write()
                    legend.AddEntry(hist, self.leg_1d_distr[index][i])
                hmax = Find_Hist_Max(histlist) 
                   
                for ihist in range( len(histlist) ):
                    if ihist == 0:
                        histlist[ihist].SetMaximum(hmax*1.3)
                        if self.logy_1d_distr[index] :
                            c.SetLogy()
                            histlist[ihist].SetMinimum(0.5)
                            histlist[ihist].SetMaximum(hmax*10)
                        histlist[ihist].Draw()
                    else: histlist[ihist].Draw("SAME")
                legend.Draw()
                c.SaveAs("plots/%s/%s/c_%s_%s.eps" % (self.leg_tune_distr[ifile], self.case, self.case, self.var_1d_distr[index][0]))
            outfile.Close()

    def hadron_ptratio(self):
       
        for ifile in range(len(self.pthad_fileinput)):
            self.dfpt = pd.read_csv(self.pthad_fileinput[ifile])
            
            dfD0 = self.dfpt[self.dfpt.hadron_pdg == 421]
            dfhad = self.dfpt[self.dfpt.hadron_pdg == self.pdg]

            pt_bins = [0,1,2,3,5,10,20,50]
            c = TCanvas('c_'+self.case+'pt_over_D0pt','',600,600)
            legend = TLegend(.7,.75,.9,.9)
            histlist = []
            
            for i in range(len(self.multi_bins)):
                legstr = ''
                if i<len(self.multi_bins)-1: 
                    dfD0sel  =  dfD0[np.logical_and(dfD0.multiplicity >= self.multi_bins[i],   dfD0.multiplicity < self.multi_bins[i+1])]
                    dfhadsel = dfhad[np.logical_and(dfhad.multiplicity >= self.multi_bins[i], dfhad.multiplicity < self.multi_bins[i+1])]
                    legstr = str(self.multi_bins[i]) + "<= multiplicity <" + str(self.multi_bins[i+1])
                else: 
                    dfD0sel  =  dfD0[dfD0.multiplicity >=  self.multi_bins[i]]
                    dfhadsel = dfhad[dfhad.multiplicity >= self.multi_bins[i]]
                    legstr = "multiplicity >=" + str(self.multi_bins[i])
                h_D0 = TH1F("h_D0_pt_multbin"+str(i),"h_D0_pt_multbin"+str(i),7,array('d',pt_bins))
                h_had = TH1F("h_"+self.case+"_pt_multbin"+str(i),"h_"+self.case+"_pt_multbin"+str(i),7,array('d',pt_bins))
                h_D0.Sumw2()
                h_had.Sumw2()
                fill_hist(h_D0, dfD0sel['pt_had'])
                fill_hist(h_had, dfhadsel['pt_had'])
                h_had.Divide(h_D0)
                h_had.SetTitle("Histogram ratio of " + self.case + " pT over D0 pT spectra")
                h_had.SetXTitle("p_{T} (GeV)")
                h_had.SetYTitle("# "+ self.case + "/# D^{0}")
                h_had.SetLineColor(self.colours[i])
                h_had.SetStats(0)
                legend.AddEntry(h_had, legstr)
                histlist.append(h_had)
            hmax = Find_Hist_Max(histlist)
            for ihist in range(len(histlist)):
                if ihist == 0:
                    histlist[ihist].SetMaximum(hmax*1.5)
                    histlist[ihist].Draw()
                else: histlist[ihist].Draw("SAME")
            legend.Draw()
            c.SaveAs("plots/%s/%s/c_%spt_over_D0pt.eps" % (self.leg_tune_distr[ifile], self.case, self.case))

    def hadron_multiratio(self):

        for ifile in range(len(self.pthad_fileinput)):
            self.dfpt = pd.read_csv(self.pthad_fileinput[ifile])

            dfD0 = self.dfpt[self.dfpt.hadron_pdg == 421]
            dfhad = self.dfpt[self.dfpt.hadron_pdg == self.pdg]
            multibins = self.multi_bins
            multibins.append(1000)

            c = TCanvas('c_'+self.case+'_multiratio','',600,600)
            h_D0 = TH1F("h_D0_multiratio","h_D0_multiratio",3,array('d',multibins))
            h_had = TH1F("h_"+self.case+"_multiratio","h_"+self.case+"_multiratio",3,array('d',multibins))
            h_D0.Sumw2()
            h_had.Sumw2()
            fill_hist(h_D0, dfD0['multiplicity'])
            fill_hist(h_had, dfhad['multiplicity'])
            h_had.Divide(h_D0)
            h_had.SetTitle("Histogram ratio of " + self.case + " multiplicity over D0 multiplicity")
            h_had.SetXTitle("# final particles")
            h_had.SetYTitle("# "+ self.case + "/# D^{0}")
            h_had.SetStats(0)
            h_had.Draw()
            c.SaveAs("plots/%s/%s/c_%s_multiratio.eps" % (self.leg_tune_distr[ifile], self.case, self.case))

    def plot_tunes(self):
        print("plotting tunes of Pythia")
        
        for index in range(len(self.var_tune_distr)):
            histlist = [ [] for x in range(len(self.var_tune_distr[index])) ]
            for ifile in range(len(self.tune_fileinput)):
                dftune = pd.read_csv(self.tune_fileinput[ifile])
                self.add_derived(dftune)
                for i in range(len(self.var_tune_distr[index])):
                    hist = TH1F("h_"+self.leg_tune_distr[ifile]+"_"+self.var_tune_distr[index][i],
                            "h_"+self.leg_tune_distr[ifile]+"_"+self.var_tune_distr[index][i],
                            self.nbins_tune_distr[index],
                            self.min_tune_distr[index],
                            self.max_tune_distr[index])
                    hist.Sumw2()
                    hist.SetLineColor(self.colours[ifile])
                    hist.SetTitle(self.case + ": " + self.title_tune_distr[index][i] )
                    hist.SetXTitle(self.xlab_tune_distr[index])
                    hist.SetYTitle(self.ylab_tune_distr[index])
                    hist.SetStats(0)
                    fill_hist(hist, dftune[self.var_tune_distr[index][i]])
                    histlist[i].append(hist)
            c = TCanvas('c_'+self.var_tune_distr[index][0],'',600*len(self.var_tune_distr[index]),600)
            c.Divide(len(histlist))
            for icanv in range(len(histlist)):
                legend = TLegend(.7,.75,.9,.9)
                c.cd(icanv+1)
                hmax = Find_Hist_Max(histlist[icanv])
                for ihist in range(len(histlist[icanv])):
                    legend.AddEntry(histlist[icanv][ihist], self.leg_tune_distr[ihist])
                    if ihist == 0:
                        histlist[icanv][ihist].SetMaximum(hmax*1.3)
                        histlist[icanv][ihist].Draw()
                    else: histlist[icanv][ihist].Draw("SAME")
                legend.Draw()
            c.SaveAs("plots/TuneComparison/%s/c_%s_%s.eps" % (self.case, self.case, self.var_tune_distr[index][0]))
