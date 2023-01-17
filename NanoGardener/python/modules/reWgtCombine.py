from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.common.collectionMerger import collectionMerger
import ROOT
import os.path
ROOT.PyConfig.IgnoreCommandLineOptions = True


class ReWgtCombineSamples(Module):
    def __init__(self, sample, year, sample_type, renormWgt_path, combineWgt_path):
        print("########################", sample)
        self.sample = sample
        self.year = year
        self.cmssw_base = os.getenv('CMSSW_BASE')
        self.cmssw_arch = os.getenv('SCRAM_ARCH')

        self.sample_type = sample_type

        self.sampleMass = str(self.sample)[str(self.sample).find('_M') + 2:]
        print("SAMPLE MASS: " + str(self.sampleMass) +  " | SAMPLE TYPE: " + str(self.sample_type))

        self.renormWgt =  1
        renorm_wgt_file = open(renormWgt_path)
        for line in renorm_wgt_file.readlines():
            line_comps = line.split("\t")
            if(int(line_comps[0]) == int(self.sampleMass)):
                self.renormWgt = float(line_comps[1])
        print("RENORM_WGT: ", self.renormWgt)

        self.combineWgts =  []
        combine_wgt_file = open(combineWgt_path)
        for line in combine_wgt_file.readlines(): 
            line_comps = line.split(":")
            if(int(line_comps[0]) == int(self.sampleMass)):
                self.combineWgts = [float(wgt) for wgt in line_comps[1].split(",")]
        print("COMBINE_WGTS: ", self.combineWgts)

        pass

    def beginJob(self, histFile=None, histDirName=None):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.newbranches =  ['HWWOffshell_combWeight_' + str(self.sample_type)]
        for nameBranches in self.newbranches :
            self.out.branch(nameBranches  ,  "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        genPart = Collection(event, "GenPart")


        W_BOSON_pdgID = 24
        H_BOSON_pdgID = 25
        
        nWs = 0
        w_ids = []

        for i in range(0, len(genPart)):
            part = genPart[i]
            if(abs(part.pdgId) == W_BOSON_pdgID and abs(genPart[genPart[i].genPartIdxMother].pdgId) == H_BOSON_pdgID):
                nWs += 1
                w_ids.append(i)

        if(nWs != 2):
            return True

        mass_window_index = 0
        m_ww = genPart[w_ids[0]].p4() + genPart[w_ids[1]].p4()
        m_ww_mass = m_ww.M()

        mass_window_edges = [0,136.7,148.3,165,175,185,195,205,220,240,260,285,325,275,450,550,650,750,850,950,1250,1750,2250,2750,14000]

        for i in range(0, len(mass_window_edges)-1):
            if(m_ww_mass > mass_window_edges[i] and m_ww_mass < mass_window_edges[i+1]):
                mass_window_index = i
                break

        self.out.fillBranch('HWWOffshell_combWeight_' + str(self.sample_type), self.combineWgts[mass_window_index]*self.renormWgt)
        return True
