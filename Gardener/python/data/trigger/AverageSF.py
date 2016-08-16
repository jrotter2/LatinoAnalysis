#!/usr/bin/env python

import json
import sys
from sys import exit
import ROOT
import optparse
import os.path
import string
import logging
import traceback
from array import array
from collections import OrderedDict
import math

import os




if __name__ == '__main__':
    print '''
--------------------------------------------------------------------------------------------------

             \                                            \ \     /       |                    
            _ \ \ \   /  _ \   __|  _` |   _` |   _ \      \ \   /  _` |  |  |   |   _ \   __| 
           ___ \ \ \ /   __/  |    (   |  (   |   __/       \ \ /  (   |  |  |   |   __/ \__ \ 
         _/    _\ \_/  \___| _|   \__,_| \__, | \___|        \_/  \__,_| _| \__,_| \___| ____/ 
                                         |___/                                                 
 
--------------------------------------------------------------------------------------------------
'''    

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)

    parser.add_option('--outputFile'     , dest='outputFile'     , help='output file with averaged values'           , default='test.txt')
          
    # read default parsing options as well
    #hwwtools.addOptions(parser)
    #hwwtools.loadOptDefaults(parser)
    (opt, args) = parser.parse_args()

    sys.argv.append( '-b' )
    ROOT.gROOT.SetBatch()


    print " outputFile       = ", opt.outputFile

     
      
    print "list of trigger files"
    
    files_names_trigger = {}
    files_names_trigger["ICHEP2016/HLT_DoubleMuLegLowPt.txt"] = 1.32
    files_names_trigger["ICHEP2016/HLT_DoubleMuLegLowPt_BeforeRun274094.txt"] = 1.60
    #files_names_trigger["ICHEP2016/HLT_EleMuLegLowPt_BeforeRun274094.txt"] = 4.50
    
    print " files_names_trigger = ", files_names_trigger
    
    files_trigger = {}    
    for names in files_names_trigger :
       #print names
       ##print files_names_trigger[names]       
       
       files_trigger [names] = {
         'file' : open (names),
         'lumi' : files_names_trigger[names]
         }

      
   
    list_triggers = {}

    for names in files_trigger :
      list_triggers[names]  = {  
             'values' : [line.rstrip().split() for line in files_trigger[names]['file']        if '#' not in line] ,
             'lumi'   : files_trigger[names]['lumi']
          }
      




    outFile = open(opt.outputFile, 'w')
    outFile.write('# etamin  etamax  ptmin   ptmax   eff     deff_high       deff_low\n')

    first_key = list_triggers.keys()[0] 
    
    for point in list_triggers[first_key]['values']:
      outFile.write( '  {0:.2f}  {1:.2f}  {2:.2f}  {3:.2f} '.format(float(point[0]), float(point[1]), float(point[2]), float(point[3])) )
      
      eff = float(point[4])
      error_eff = float(point[5])
      
      error_eff_up = eff + error_eff
      error_eff_lo = eff - error_eff

      # muons and electrons have provided different formats!
      if len(point) > 6 :
         error_eff_lo = float(point[6])
         error_eff_lo = eff - error_eff_lo                  
      
      luminosity = list_triggers[first_key]['lumi']
     
      all_lumi       =   luminosity
      eff            =   eff          * luminosity
      error_eff_up   =   error_eff_up * luminosity
      error_eff_lo   =   error_eff_lo * luminosity
      
      
      for names in list_triggers :
        if names != first_key:
          for point_other in list_triggers[names]['values']:
            if    float(point_other[0]) == float(point[0]) \
              and float(point_other[1]) == float(point[1]) \
              and float(point_other[2]) == float(point[2]) \
              and float(point_other[3]) == float(point[3]) :

              new_eff = float(point_other[4])
              new_error_eff = float(point_other[5])
              
              new_error_eff_up = new_eff + new_error_eff
              new_error_eff_lo = new_eff - new_error_eff

              # muons and electrons have provided different formats!
              if len(point_other) > 6 :
                 new_error_eff_lo = float(point_other[6])
                 new_error_eff_lo = new_eff - new_error_eff_lo                  
        
        
              luminosity = list_triggers[names]['lumi']

              all_lumi       =   all_lumi     + luminosity
              eff            =   eff          +  new_eff          * luminosity
              error_eff_up   =   error_eff_up + new_error_eff_up * luminosity
              error_eff_lo   =   error_eff_lo + new_error_eff_lo * luminosity
                      
      eff            =   eff /  all_lumi
      error_eff_up   =   error_eff_up / all_lumi
      error_eff_lo   =   error_eff_lo / all_lumi
              
      outFile.write( '        {0:.4f}  {1:.4f}  {2:.4f}  \n'.format(float(eff), float(error_eff_up - eff), float(eff-error_eff_lo)) )
    
          
        
        

    #for point in self.list_triggers[whichTrigger] :
           
           #if ( eta >= float(point[0]) and eta <= float(point[1]) and         # the "=" in both directions is only used by the overflow bin
                #pt  >= float(point[2]) and pt  <= float(point[3]) ) :         # in other cases the set is (min, max]
               
               #eff       = float(point[4])
               #error_eff = float(point[5])
               
               #error_eff_up = eff + error_eff
               #error_eff_lo = eff - error_eff
               
               ## muons and electrons have provided different formats!
               #if len(point) > 6 :
                  #error_eff_lo = float(point[6])
                  #error_eff_lo = eff - error_eff_lo                  
               
               ## correction not to have >1 probability!!!
               #if error_eff_up > 1.0 : 
                 #error_eff_up = 1.0
               #if error_eff_lo < 0.0 :
                 #error_eff_lo = 0.0
                  
               #return eff, error_eff_lo, error_eff_up
        
        ##print " is 0 ... ", pt, "  ", eta, "   " , whichTrigger
        ##if pt > 100:
          ##for point in self.list_triggers[whichTrigger] :
            ##print " ---> ", float(point[0]), " ", float(point[1]), " ", float(point[2]), " ", float(point[3])
        
        #return 0. , 0., 0.
           

    print " outputFile       = ", opt.outputFile


