from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.common.collectionMerger import collectionMerger
import ROOT
import json
import os.path
ROOT.PyConfig.IgnoreCommandLineOptions = True


class ReWgtCombineSamples(Module):
    def __init__(self, sample, year, hww_wgt_path):
        print("########################", sample)
        self.sample = sample
        self.year = year
        self.cmssw_base = os.getenv('CMSSW_BASE')
        self.cmssw_arch = os.getenv('SCRAM_ARCH')

        
        self.sample_type = "ggH"
        if(str(self.sample).__contains__("VBF")):
            self.sample_type = "VBF"

        self.sampleMass = str(self.sample)[str(self.sample).find('_M') + 2:]

        sample_json_id = "HM_" + str(self.sampleMass)
        print("SAMPLE MASS: " + str(self.sampleMass) +  " | SAMPLE TYPE: " + str(self.sample_type))

        hww_wgt_file = open(self.cmssw_base + '/src/' + hww_wgt_path + "_" + self.sample_type + ".json")
        hww_wgt_contents = hww_wgt_file.read()
        hww_wgts = json.loads(hww_wgt_contents)

        self.renormWgt = hww_wgts[sample_json_id]["RENORM_WGT"]
        self.combineWgt = hww_wgts[sample_json_id]["COMB_WGTS"]

        self.lossCompWgt = hww_wgts[sample_json_id]["LOSS_COMP"]

        self.cutoffSIG = hww_wgts[sample_json_id]["SIG_CUTOFF"]
        self.cutoffCONT = hww_wgts[sample_json_id]["CONT_CUTOFF"]
        self.cutoffSIGplusCONT = hww_wgts[sample_json_id]["SIGplusCONT_CUTOFF"]    
        print("RENORM WGT: " + str(self.renormWgt))
        print("COMB WGT:" + str(self.combineWgt))
        print("LOSS COMP:" + str(self.lossCompWgt))
        pass

    def beginJob(self, histFile=None, histDirName=None):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.newbranches =  ['HWWOffshell_combineWgt']
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

        mass_window_edges = [0,136.7,148.3,165,175,185,195,205,220,240,260,285,325,375,425,475,525,575,650,750,850,950,1250,1750,2250,2750,14000]

        for i in range(0, len(mass_window_edges)-1):
            if(m_ww_mass > mass_window_edges[i] and m_ww_mass < mass_window_edges[i+1]):
                mass_window_index = i
                break

        if(mass_window_index >= len(self.combineWgt)):
            print("EVENT OUT OF BOUNDS >>> MASS: " + str(m_ww_mass) + ", WINDOW: " + str(mass_window_index))
            print("Adding to overflow bin...")
            mass_window_index = len(self.combineWgt)-1

        mc_wgt_SIG = event.p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM
        mc_wgt_CONT = event.p_Gen_GG_BKG_MCFM
        mc_wgt_SIGplusCONT = event.p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM
        if(self.sample_type == "VBF"):
            mc_wgt_SIG = event.p_Gen_JJEW_SIG_ghv1_1_MCFM
            mc_wgt_CONT = event.p_Gen_JJEW_BKG_MCFM
            mc_wgt_SIGplusCONT = event.p_Gen_JJEW_BSI_ghv1_1_MCFM

        SIG_wgt = mc_wgt_SIG*event.p_Gen_CPStoBWPropRewgt*event.XSWeight
        CONT_wgt = mc_wgt_CONT*event.p_Gen_CPStoBWPropRewgt*event.XSWeight 
        SIGplusCONT_wgt = mc_wgt_SIGplusCONT*event.p_Gen_CPStoBWPropRewgt*event.XSWeight
 
        if(SIG_wgt > self.cutoffSIG[mass_window_index] or CONT_wgt > self.cutoffCONT[mass_window_index] or SIGplusCONT_wgt > self.cutoffSIGplusCONT[mass_window_index]):
            self.out.fillBranch('HWWOffshell_combineWgt', 0)
        elif(SIG_wgt == 0 and CONT_wgt == 0 and SIGplusCONT_wgt == 0):
            self.out.fillBranch('HWWOffshell_combineWgt', 0)
        else:
            self.out.fillBranch('HWWOffshell_combineWgt', self.renormWgt * self.combineWgt[mass_window_index] * self.lossCompWgt[mass_window_index])

        return True

