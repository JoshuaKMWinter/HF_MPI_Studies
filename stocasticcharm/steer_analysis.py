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

"""
main script for doing data processing, machine learning and analysis
"""

#import os
import subprocess
import argparse
from os.path import exists
import yaml
from pkg_resources import resource_stream
from analyzer import Analyzer

def steer_analysis(case):

    with open("data/database.yml", 'r') as datafile:
        data_param = yaml.load(datafile, Loader=yaml.FullLoader)

    myan = Analyzer(data_param[case], case)
    do_distributions = True
    if do_distributions is True:
        myan.plot()

steer_analysis("Ds")