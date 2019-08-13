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
from ROOT import TFile, TObject, TList, TClass, TKey, TTree, TTreeReader, TH1F, TCanvas
from ROOT import gStyle, TLegend
from ROOT import gROOT
from ROOT import TStyle
import numpy as np
import copy

# pylint: disable=too-few-public-methods, too-many-instance-attributes, too-many-statements

def Find_Hist_Max(histlist, scaled = False):
    hmax = 0
    for i in range(len(histlist)):
        norm = histlist[i].GetEntries()
        normmax = histlist[i].GetMaximum() / norm;
        if normmax > hmax : hmax = normmax
    return hmax

class Plotter:
    species = "analyzer"
    def __init__(self):
        self.cases = ["D0", "Ds", "Lc"]

    def get_num_objects(self, infile, objclass):
        n_obj = 0
        for key in infile.GetListOfKeys():
            tclass = gROOT.GetClass(key.GetClassName())
            if (tclass.InheritsFrom(objclass)): n_obj+=1
        return n_obj

    def get_objects(self, infile, objclass):
        objarray = []
        for key in infile.GetListOfKeys():
            tclass = gROOT.GetClass(key.GetClassName())
            if (tclass.InheritsFrom(objclass)):
                obj = key.ReadObj()
                objarray.append(obj)
        return objarray

    def plotcomparison(self):
        print("Running plotter")
        tempfile = TFile.Open("data/"+self.cases[0]+"_hists.root")
        histlist = [ [] for i in range(len(self.cases))]
        colours = [4,2,1]

        for i in range(len(self.cases)):
            infile = TFile.Open("data/" + self.cases[i] + "_hists.root")
            caselist = self.get_objects(infile, "TH1")
            histlist[i] = copy.deepcopy(caselist)
            infile.Close()
        for icanv in range(len(histlist[0])):
            c = TCanvas("c_"+histlist[0][icanv].GetName())
            leg = TLegend(.87,.6,.98,.75)
            temparr = []
            for i in range(len(histlist)):
                temparr.append(histlist[i][icanv])
            hmax = Find_Hist_Max(temparr)
            for i in range(len(histlist)):
                histlist[i][icanv].SetLineColor(colours[i])
                norm = histlist[i][icanv].GetEntries()
                histlist[i][icanv].Scale(1./norm)
                leg.AddEntry(histlist[i][icanv], self.cases[i])
                if i==0 : 
                    histlist[i][icanv].SetMaximum(hmax*1.2)
                    histlist[i][icanv].Draw()
                else: histlist[i][icanv].Draw("same")
            leg.Draw()
            c.SaveAs("plots/Comparison/%s.eps" % histlist[0][icanv].GetName())
        
