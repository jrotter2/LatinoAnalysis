import ROOT
import math 
import numpy
import ctypes
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.common.collectionMerger import collectionMerger

import os.path

# THIS MODULE NEEDS TO BE MODIFIED IF ZH AND/OR WH ARE CONSIDERED, IN ORDER TO ACCOUNT FOR THE ASSOCIATED W/Z BOSONS

# In that case, one should check if NanoAODv9 keeps intermediate H bosons. Need to determine unambiguously which WW pair is assumed to be the Higgs boson for WH in POWHEG, so that we can reweight s-channel WH to the s+u configuration properly.

class LHE_MEMaker(Module):
    def __init__(self, sample, year, cfg_path, xsecs_path):
        print '####################', sample
        self.sample = sample
        self.year = year
        self.cmssw_base = os.getenv('CMSSW_BASE')
        self.cmssw_arch = os.getenv('SCRAM_ARCH')

        ROOT.gSystem.AddIncludePath(self.cmssw_base+"/lib/")
        ROOT.gSystem.Load("./libJHUGenMELAMELA.so") 
        ROOT.gSystem.Load(self.cmssw_base+"/src/JHUGenMELA/MELA/data/"+self.cmssw_arch+"/libmcfm_707.so")
        ROOT.gSystem.Load("./libIvyFrameworkIvyDataTools.so") 
        ROOT.gSystem.Load("./libIvyFrameworkIvyAutoMELA.so") 
        ROOT.gSystem.Load("./libMelaAnalyticsGenericMEComputer.so") 
        ROOT.gSystem.Load("./libMelaAnalyticsEventContainer.so") 
        ROOT.gSystem.Load("./libMelaAnalyticsCandidateLOCaster.so") 
        ROOT.gInterpreter.Declare('#include "/eos/user/f/fernanpe/CMSSW_10_6_4/src/MelaAnalytics/EventContainer/interface/MELAEvent.h"')
        ROOT.gInterpreter.Declare('#include "/eos/user/f/fernanpe/CMSSW_10_6_4/src/IvyFramework/IvyAutoMELA/interface/IvyMELAHelpers.h"')
        ROOT.gInterpreter.Declare('#include "/eos/user/f/fernanpe/CMSSW_10_6_4/src/JHUGenMELA/MELA/interface/PDGHelpers.h"')
        ROOT.gInterpreter.Declare('#include "/eos/user/f/fernanpe/CMSSW_10_6_4/src/MelaAnalytics/EventContainer/interface/HiggsComparators.h"')

        # LOAD THE theLHEProbabilities COLLECTION

        lheinfo = {}
        execfile(self.cmssw_base+'/src/'+cfg_path, lheinfo) # REF: https://github.com/usarica/CMS3NtupleTools/tree/combined/NtupleMaker/data/LHEProbabilities 

        # REPLACE '<HMASS>' BY THE SPECIFIC MASS VALUE OF THE INPUT SAMPLE

        self.theLHEProbabilities = ROOT.vector('string')() #IvyMELAHelpers::GMECBlock::buildMELABranches(const vector<string>& MElist, bool isGen) CORRECT FORMAT TO MAKE THE CALL
        
        rawLHEProbabilities = lheinfo['theLHEProbabilities']
        self.mass = str(self.sample)[str(self.sample).find('_M') + 2:] # Followed the v7 name convention for the HM samples
        for idx, item in enumerate(rawLHEProbabilities):
            if '<HMASS>' in item:
                rawLHEProbabilities[idx] = item.replace('<HMASS>', self.mass)
            self.theLHEProbabilities.push_back(rawLHEProbabilities[idx])
            
        # Get the Xsec and BR for this particular mass value from the txt file, and Build the binning

        xsecs = open(self.cmssw_base+'/src/'+xsecs_path) # REF: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBSMAt13TeV
        central_values = []
        
        for line in xsecs:
            central_values.append(float(line.split()[0]))
            if int(line.split()[0]) == int(self.mass):
                self.xsec = float(line.split()[1])
                self.br = float(line.split()[2])
                
        self.final_binning = []

        for idx,value in enumerate(central_values):
            if(idx < len(central_values)-1):
                edgevalue = (value+central_values[idx+1])/2 
                if idx == 0:
                    self.final_binning.append(70.) #from Ulascan [70,.....,13000]
                    #self.final_binning.append(value - abs(value-edgevalue)) #another option?
                self.final_binning.append(edgevalue)
            else:
                self.final_binning.append(13000.)
                #self.final_binning.append(abs(value-edgevalue) + value)
                
        #print(self.final_binning) #Is this the binning we want to have?

        # DEFINE THE 'VVMode' AND 'VVDecaymode' FLAGS (https://github.com/usarica/CMS3NtupleTools/blob/combined/NtupleMaker/test/samples_MC_signals_WW_2018.csv)

        self.VVMode = ROOT.MELAEvent.getCandidateVVModeFromString("WW")
        self.isGen = True

        if "GluGluHToWW" in self.sample :
          self.VVDecayMode = 0
        elif "VBF" in self.sample :
          self.VVDecayMode = 0
        elif "WminusHToWW" in self.sample or "WplusHToWW" in self.sample:
          self.VVDecayMode = -1
        elif "HZJ_HToWW" in self.sample or "ZH_HToWW" in self.sample:
          self.VVDecayMode = -1
        else:
          raise NameError(self.sample, "is an unrecognised simulation")
          # ONLY VALID FOR WWT02L2Nu. FULL MAP: https://github.com/MELALabs/MelaAnalytics/blob/master/EventContainer/src/MELAEvent.cc#L177-L187

        ROOT.IvyMELAHelpers.setupMela(self.year, float(self.mass), ROOT.TVar.ERROR)

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        
        self.out = wrappedOutputTree
        
        self.newbranches =  ['LHECandMass']

        for nameBranches in self.newbranches :
          self.out.branch(nameBranches  ,  "F");

        self.lheMEblock = ROOT.IvyMELAHelpers.GMECBlock()
        self.lheMEblock.addRefTree(self.out.tree()) # to be able to push the MEs to the output trees directly 
        self.lheMEblock.buildMELABranches(self.theLHEProbabilities, self.isGen)
        
        # Get new ME branches names
        self.strMEs = []
        mes = ROOT.std.unordered_map('string','float')() 
        self.lheMEblock.getBranchValues(mes)
        for me in mes:
            self.strMEs.append(me.first)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        gen_dme_hsm   = -999

        self.LHE = Collection(event,"LHEPart")
        Gen = Collection(event,"GenPart")

        incoming = ROOT.vector('TLorentzVector')()
        incomingIDs = ROOT.vector('int')()
        outgoing   = ROOT.vector('TLorentzVector')()
        outgoingIDs = ROOT.vector('int')()
        
        # COVERT THE P4 INCOMING & OUTOGOING PARTICLES INTO MELA PARTICLES

        melaParticlesList = ROOT.vector('MELAParticle')()
        
        for partid, part in enumerate(self.LHE):
            temp = ROOT.TLorentzVector()
            if part.status == 1: # LHE particle status; -1:incoming, 1:outgoing
                temp.SetPtEtaPhiM(part.pt, part.eta, part.phi, part.mass)
                outgoing.push_back(temp)
                outgoingIDs.push_back(part.pdgId)
            elif part.status == -1:
                temp.SetPxPyPzE(0.,0., abs(part.incomingpz)*6500, abs(part.incomingpz)*6500)
                incoming.push_back(temp)
                incomingIDs.push_back(part.pdgId)

        daughter_coll = ROOT.SimpleParticleCollection_t() 
        mother_coll = ROOT.SimpleParticleCollection_t()

        for idx, dau in enumerate(outgoing):
         tempPart = ROOT.SimpleParticle_t(outgoingIDs[idx], dau)         
         tempMela = ROOT.MELAParticle(tempPart.first, tempPart.second)
         tempMela.setGenStatus(1) # Status 2 for associated W production in WH? Is this info available in nanoAOD?
         melaParticlesList.push_back(tempMela)
         
        
        for idx, mot in enumerate(incoming):
         tempPart = ROOT.SimpleParticle_t(incomingIDs[idx], mot)         
         tempMela = ROOT.MELAParticle(tempPart.first, tempPart.second)
         tempMela.setGenStatus(-1) #DANGER: GENSTATUS = 1 IF STABLE?
         melaParticlesList.push_back(tempMela)


        # BUILD THE MELA EVENT FROM THE MELA PARTICLES
         
        LHEEvent = ROOT.MELAEvent()
        writtenGenCands = ROOT.vector('MELAParticle')()
        writtenGenTopCands = ROOT.vector('MELAParticle')()
        ThereIsHiggs = False

        for part in melaParticlesList:
            if part.genStatus == -1:
                LHEEvent.addMother(part)
            else:
                if ROOT.PDGHelpers.isALepton(part.id):
                    LHEEvent.addLepton(part)
                elif ROOT.PDGHelpers.isANeutrino(part.id):
                    LHEEvent.addNeutrino(part)
                elif ROOT.PDGHelpers.isAKnownJet(part.id) and ROOT.PDGHelpers.isATopQuark(part.id) == False:
                    LHEEvent.addJet(part)
                elif ROOT.PDGHelpers.isAPhoton(part.id):
                    LHEEvent.addPhoton(part)
                elif ROOT.PDGHelpers.isAHiggs(part.id):
                    writtenGenCands.push_back(part)
                    ThereIsHiggs = True
                    if (self.VVMode==ROOT.MELAEvent.UndecayedMode and (part.genStatus==1 or part.genStatus==2)):
                        LHEEvent.addIntermediate(part);
                elif ROOT.PDGHelpers.isATopQuark(part.id):
                    writtenGenTopCands.push_back(part)
                    if (part.genStatus==1):
                        LHEEvent.addIntermediate(part);
                else:
                    print("FATAL: UNIDENTIFIED PARTICLE IN THE FINAL STATE")
                    print("ID: ", part.id)
                
        
        LHEEvent.constructTopCandidates()
        # FOLLOWING BLOCK DOES NOT SEEM TO BE NEEDED IN OUR CASE: https://github.com/usarica/CommonLHETools/blob/master/LHEHandler/src/LHEHandler.cc#L183-L198
        # matchedTops = ROOT.vector('MELATopCandidate_t')()
        # for writtenGenTopCand in writtenGenTopCands:
        #     tmpCand = ROOT.TopComparators.matchATopToParticle(LHEEvent, writtenGenTopCand)
        #     if tmpCand:
        #         matchedTops.push_back(tmpCand)
        # for tmpCand in LHEEvent.getTopCandidates():
        #     if (ROOT.std.find(matchedTops.begin, matchedTops.end, tmpCand)==matchedTops.end):
        #         tmpCand.setSelected(False)

        LHEEvent.constructVVCandidates(self.VVMode, self.VVDecayMode)
        LHEEvent.addVVCandidateAppendages()
        
        ThereIsCand = False
        if ThereIsHiggs:
            for writtenGenCand in writtenGenCands:
                tmpCand = ROOT.HiggsComparators.matchAHiggsToParticle(LHEEvent, writtenGenCand)
                if tmpCand:
                    if not ThereIsCand:
                        genCand = tmpCand    
                    else: 
                        genCand = ROOT.HiggsComparators.candComparator(genCand, tmpCand, ROOT.HiggsComparators.BestZ1ThenZ2ScSumPt, self.VVMode);
        else:
            genCand = ROOT.HiggsComparators.candidateSelector(LHEEvent, ROOT.HiggsComparators.BestZ1ThenZ2ScSumPt, self.VVMode);

        # MAKE THE IvyAutoMELA CALL 

        ROOT.IvyMELAHelpers.melaHandle.setCurrentCandidate(genCand)
        self.lheMEblock.computeMELABranches()
        

        ################# REWEIGHTING STEP ##########################
        # Needed inputs: binning, xsec, BR, gen MC weight, LHE MEs weights

        # event.genWeight -> nanoAOD genWeight branch: nominal generator-level weight, computed as the product of the nominal weight from LHE and the nominal weight from Pythia
        # self.final_binning 
        # self.xsec
        # self.br
        # L119-123 -> do something similar to get the computed LHE ME weights from me.second



        # FILL THE OUTPUT BRANCHES


        self.out.fillBranch('LHECandMass', genCand.m())
        self.lheMEblock.pushMELABranches()
        
        # CLEAR MELA

        ROOT.IvyMELAHelpers.melaHandle.resetInputEvent()


        return True   
        


