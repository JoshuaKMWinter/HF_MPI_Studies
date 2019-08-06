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

# pylint: disable=too-few-public-methods, too-many-instance-attributes, too-many-statements
class Analyzer:
    species = "analyzer"
    def __init__(self, datap, case):
        self.case = case
        self.mass = datap["mass"]
        self.rapiditymax = datap["rapiditymax"]
        self.rapiditymin = datap["rapiditymin"]
        self.var_all = datap["variables"]["var_all"]
        self.inputfile = datap["fileinput"]
        self.sel_1d_ptchad = datap["sel_1d_ptchad"]
        self.var_1d_distr = datap["var_1d_distr"]
        self.nbins_1d_distr = datap["nbins_1d_distr"]
        self.min_1d_distr = datap["min_1d_distr"]
        self.max_1d_distr = datap["max_1d_distr"]
        self.n_1d = len(self.max_1d_distr)
        self.dfm = None
        self.dfm = pd.read_csv(self.inputfile)
        self.add_derived()

    def add_derived(self):
        if self.case == "Ds" or "Lc":
            self.dfm["deltaphi"] = self.dfm["c_phi"] - self.dfm["cbar_phi"]

    def plot(self):
        print("I am running plotter")
        dfsel = self.dfm.query(self.sel_1d_ptchad)
        cdistr = TCanvas('cEff', 'The Fit Canvas', 1000,200)
        cdistr.Divide(self.n_1d, 1)
        histlist = []
        for index in range(self.n_1d):
            hist = TH1F("hist" + self.var_1d_distr[index],
                        "hist" + self.var_1d_distr[index],
                        self.nbins_1d_distr[index],
                        self.min_1d_distr[index],
                        self.max_1d_distr[index])
            fill_hist(hist, dfsel[self.var_1d_distr[index]])
            histlist.append(hist)

        for index in range(self.n_1d):
            cdistr.cd(index+1)
            histlist[index].Draw()
        cdistr.SaveAs("cdistr%s.eps" % self.case)
