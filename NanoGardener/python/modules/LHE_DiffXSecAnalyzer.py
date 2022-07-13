
import ROOT
import math 
import numpy
import ctypes
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.common.collectionMerger import collectionMerger

import os.path


class LHE_DiffXSecAnalyzer(Module):
    def __init__(self, sample, year):
        print '####################', sample
        self.sample = sample
        self.year = year
        self.cmssw_base = os.getenv('CMSSW_BASE')
        self.cmssw_arch = os.getenv('SCRAM_ARCH')
        self.writeHistFile = True

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        
        self.out = wrappedOutputTree

        self.diffXSec = ROOT.TH1F('diffXSec', '$d\sigma / dE$ vs $m_{WW}$$', 5000, -0.5, 4999.5)
        #Adding histogram as an Object to be saved later
        self.addObject(self.diffXSec)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        genPart = Collection(event,"GenPart")

        
        nWs = 0
        w_ids = []

        for i in range(0, len(genPart)):
            part = genPart[i]
            if(abs(part.pdgId) == 24 and abs(genPart[genPart[i].genPartIdxMother].pdgId) == 25):
                nWs += 1
                w_ids.append(i)

        if(nWs != 2):
            return True

        m_ww = genPart[w_ids[0]].p4() + genPart[w_ids[1]].p4()
        genWeight = event.genWeight
        self.diffXSec.Fill(m_ww.M(), genWeight)

        return True
