#
# Example script to load a LZAP ROOT output, loop over the events
# and then loop over pulses from each of the detectors.
#
# Setup:
# $ source /cvmfs/lz.opensciencegrid.org/LzBuild/latest/setup.sh 
#
# Before running for the first time you need to generate some libraries
# to be able to read the LZAP RQ variables. Open the LZAP file in a root 
# session: 
# root /unix/darkmatter/jim/fordaniel/reconstructed/lz_20180622000_lzap.root
# [root] _file0->MakeProject("rqlib", "*", "update++")
# [root] .q 
#
# Then run the script as:
# $ python radon_lzap_loop.py --help 
# $ python radon_lzap_loop.py --nevents 10 --verbose /unix/darkmatter/jim/fordaniel/reconstructed/lz_20180622000_lzap.root 
# $ radon_lzap_loop.py --nevents 10 --verbose /unix/darkmatter/jim/fordaniel/reconstructed/lz_20180622000_lzap.root --plot
# $ radon_lzap_loop.py --nevents 20000 /unix/darkmatter/jim/fordaniel/pb212/LZAP-5.0.1-reconstructed/lz_20180622\*_lzap.root

import sys
import argparse
from typing import List

from pandas.core.algorithms import diff
from pandas.core.base import IndexOpsMixin
from pandas.core.indexing import IndexSlice
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

def main():    
    parser = argparse.ArgumentParser(description='Example loop over events for BACCARAT ROOT output.')
    parser.add_argument('-n', '--nevents', type=int, help='limit number of events to loop over.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show more verbose output to screen.')
    parser.add_argument('-p', '--plot', action='store_true', help='plot energy deposits for each event.')
    parser.add_argument('input', help='path to input file')
    args = parser.parse_args()
    
    print("Loading events from:\n {}".format(args.input))
    if args.nevents:
        print("will limit loop to {} events".format(args.nevents))

    verbose = args.verbose
    plot = args.plot 
 
    # load ROOT after argparse as takes a while
    import ROOT
    
    # Load the LZAP RQ libraries so that can access the TTree branches
    if ROOT.gSystem.Load("rqlib/rqlib.so") != 0:
        print("Unable to load rqlib/rqlib.so. Make sure you have generated the .so file")
        sys.exit(1)

    #S1S2Correction 2DHistograms:

    calib_file = ROOT.TFile("s1_s2_xy_corrections.root","OPEN")
    s1corrections_h2 = calib_file.Get('s1corrections_h2')
    s2corrections_h2 = calib_file.Get('s2corrections_h2')
    assert s2corrections_h2 and s1corrections_h2

    # Multi Files: Add many files to a TChain of the "Events" tree

    events = ROOT.TChain("Events")
    n_files_added = events.Add(args.input)
    print('Found {} inputs file(s) matching {}'.format(n_files_added,args.input))
    if n_files_added <=0:
        print('Cannot find input files, exiting!')
        sys.exit(1)
    if 1:
        # SpeedupExample: disable all branches 
        events.SetBranchStatus("*", 0)
        # and then enable just the ones you need
        events.SetBranchStatus("eventHeader.eventID", 1)
        events.SetBranchStatus("eventHeader.runID", 1)
        events.SetBranchStatus("eventHeader.rawFileName", 1)
        events.SetBranchStatus("pulsesTPC.nPulses", 1)
        events.SetBranchStatus("pulsesSkin.nPulses", 1)
        events.SetBranchStatus("pulsesODLG.nPulses", 1)
        events.SetBranchStatus("pulsesTPC.pulseArea_phd", 1)
        events.SetBranchStatus("pulsesTPC.bottomArea_phd", 1)
        events.SetBranchStatus("pulsesTPC.classification", 1)
        events.SetBranchStatus("pulsesTPC.pulseStartTime_ns", 1)
        events.SetBranchStatus("pulsesTPC.pulseEndTime_ns", 1)
        events.SetBranchStatus("pulsesTPC.s2Xposition_cm", 1)
        events.SetBranchStatus("pulsesTPC.s2Yposition_cm", 1)
        events.SetBranchStatus("pulsesTPC.topBottomAsymmetry", 1)
        events.SetBranchStatus("pulsesTPC.promptFraction50ns", 1)
        events.SetBranchStatus("pulsesTPC.pulseID", 1)
        events.SetBranchStatus("pulsesTPC.peakAmp",1)
        events.SetBranchStatus("pulsesSkin.pulseID",1)
        events.SetBranchStatus("pulsesSkin.pulseArea_phd",1)
        events.SetBranchStatus("pulsesSkin.pulseStartTime_ns",1)
        events.SetBranchStatus("pulsesODLG.pulseID",1)
        events.SetBranchStatus("pulsesODLG.pulseArea_phd",1)
        events.SetBranchStatus("pulsesODLG.pulseStartTime_ns",1)



    #tfile = ROOT.TFile.Open(args.input)
    #if not tfile:
    #    print("Cannot open input file, exiting!")
    #    sys.exit(1)

    # Set ROOT plotting style
    ROOT.gROOT.SetStyle("Plain")

    # If requested cap the number events
    #events = tfile.Get("Events")
    nevents = events.GetEntries()
    print("Total number of events {}".format(nevents))
    if args.nevents and args.nevents < nevents: 
        print("limiting number events from {} to {}".format(nevents, args.nevents))
        nevents = args.nevents
    #Recording general information about event types to print as summary at end of program

    overlap_event_id_array = []


    types_of_scatter_test = {
   
    }
    #defining variables to store data which are used in scatter plots to see classification performance
    spe_pulse_area=[]
    spe_pulse_width=[]
    spe_tba=[]
    spe_pf50ns=[]
    s1_pulse_area=[]
    s1_pulse_width=[]
    s1_tba=[]
    s1_pf50ns=[]
    s2_pulse_area=[]
    s2_pulse_width=[]
    s2_tba=[]
    s2_pf50ns=[]
    se_pulse_area=[]
    se_pulse_width=[]
    se_tba=[]
    se_pf50ns=[]
    mpe_pulse_area=[]
    mpe_pulse_width=[]
    mpe_tba=[]
    mpe_pf50ns=[]
    event_type=''

    #Arrays to store area of first S2 and second S2 in S1=1, S2=2 events:
    
    first_s2_areas=[]
    second_s2_areas=[]
    s1_area_list_corrected=[]
    s1_bottom_corrected_areas_list=[]
    first_s2_bottom_area_list_corrected=[]
    second_s2_bottom_area_list_corrected=[]


 
    






    overlaps_in_TPC_and_outer_detector ={
        'Overlaps in skin' : 0,
        'Overlaps in OD' : 0,
    }
    #Variables to be used

    coincidence_window = 4000000  # used to search for overlaps in S1 in TPC and OD regions
    area_decrease=0.05 # When looking at coincidence window want to see if this is overlap with TPC for similar areas but want to remove background -ie the dark rate.

    threshold = float(input("-> Enter Energy threshold to remove background pulses..... example:  0.10 = 10%, 0.20 = 20% \n-> Suggest around 0.10 \n-> "))

    # Loop over events in file 
    for entry in range(nevents):
        events.GetEntry(entry)

        tpc_pulse_counter = {
            'S1' : 0,
            'S2' : 0
            }

        tpc_pulse_counter_threshold = {
         'S1' : 0,
         'S2' : 0
        }

        s1_area= 0
        s2_area = 0
        s1_area_array =[0]

        #Arrays to gather indicies for each event
        s1_pulse_indices = -1
        s2_first_pulse_indices=-1
        s2_second_pulse_indices=-1

        if entry % 1000 == 0:
            print("-> Processed {} events of {}".format(entry,nevents))
        else:
            pass
        

        # Print event information
        if verbose: 
            print("===========================")
            print("eventId {}, run {} from {}:".\
            format(events.eventHeader.eventID, events.eventHeader.runID, 
                   events.eventHeader.rawFileName)) 

        # For convenience rename some variables
        tpc = events.pulsesTPC
        skin = events.pulsesSkin
        od = events.pulsesODLG
        if verbose: print(" -> event has {} TPC pulses, {} Skin pulses, {} OD pulses".format(
            tpc.nPulses, skin.nPulses, od.nPulses))

        # Loop over TPC detector pulses
        tpc = events.pulsesTPC # for convenience
        if verbose: print(" -> looping over {} TPC pulses:".format(tpc.nPulses))
        for ipulse in range(tpc.nPulses):
            id = tpc.pulseID[ipulse] # just number of pulse
            area = tpc.pulseArea_phd[ipulse] # total pulse area in detected photons (phd)
            classification = tpc.classification[ipulse] # classification: S1, S2, SE, SPE, MPE, Other
            start_time_ns = tpc.pulseStartTime_ns[ipulse] # pulse start
            end_time_ns = tpc.pulseEndTime_ns[ipulse] # pulse end
            width_ns =  end_time_ns-start_time_ns # pulse length in ns
            tba = tpc.topBottomAsymmetry[ipulse] # -1 all light in bottom array, +1 all light in top array 
            pf50ns = tpc.promptFraction50ns[ipulse] # fraction of pulse area that was in first 50 ns of pulse

            s1_area_array.append(area)

            if verbose:
                print("   -> {}: {:0.2f} phd, {}, {} ns, {:0.2f} tba, {:0.2f} pf50ns,  {} ns".format(
                    id, area, classification, width_ns, tba, pf50ns,start_time_ns))

            if classification == 'S1':
                tpc_pulse_counter['S1'] +=1
                if area > s1_area:
                    s1_area = area

            elif classification == 'S2':
                tpc_pulse_counter['S2'] +=1
                if area > s2_area:
                    s2_area = area 
            else:
                pass
        

         #Finding 10% of area for S1 and S2 signals in TPC to use as a threshold

        s1_threshold = s1_area*threshold*2
        s2_threshold = s2_area*threshold

       


            


        #Finding the S1 pulse with the largest area 
        max_s1_area = max(s1_area_array)
        #Additional cut to distinguish between SE and S2
        area_cut_s2= (0.188 * width_ns) + 624

        for ipulse in range(tpc.nPulses):
            classification = tpc.classification[ipulse] # classification: S1, S2, SE, SPE, MPE, Other
            area = tpc.pulseArea_phd[ipulse] # total pulse area in detected photons (phd)
            start_time_ns = tpc.pulseStartTime_ns[ipulse] # pulse start
            end_time_ns = tpc.pulseEndTime_ns[ipulse] # pulse end
            width_ns =  end_time_ns-start_time_ns # pulse length in ns
            tba = tpc.topBottomAsymmetry[ipulse] #-1 all light in bottom array, +1 all light in top array 
            pf50ns = tpc.promptFraction50ns[ipulse] # fraction of pulse area that was in first 50 ns of pulse
            peak_height = events.pulsesTPC.peakAmp[ipulse]
            


            if classification == 'S1' and area >= s1_threshold:
                tpc_pulse_counter_threshold['S1'] +=1
            if classification == 'S2' and area >= s2_threshold and area >= area_cut_s2:
                tpc_pulse_counter_threshold['S2'] +=1
            #Finding S1 pulse start time which has the largest area (Assumed to be the main signal and storing the start time to search for similar event in OS/OD)
            if area == max_s1_area:
                start_time_of_s1=start_time_ns
            
            # 2D plots using Pulse wdith/ Promt Fraction / asymmetry.
            if classification == 'SPE': 
                spe_pulse_area.append(area)
                spe_pulse_width.append(width_ns)
                spe_tba.append(tba)
                spe_pf50ns.append(pf50ns)
            if classification =='S1' and area >= s1_threshold:
                s1_pulse_area.append(area)
                s1_pulse_width.append(width_ns)
                s1_tba.append(tba)
                s1_pf50ns.append(pf50ns)
            if classification =='S2'and area >= s2_threshold and area >= area_cut_s2 :
                s2_pulse_area.append(area)
                s2_pulse_width.append(width_ns)
                s2_tba.append(tba)
                s2_pf50ns.append(pf50ns)
            if classification =='SE':
                se_pulse_area.append(area)
                se_pulse_width.append(width_ns)
                se_tba.append(tba)
                se_pf50ns.append(pf50ns)
            if classification == 'MPE':
                mpe_pulse_area.append(area)
                mpe_pulse_width.append(width_ns)
                mpe_tba.append(tba)
                mpe_pf50ns.append(pf50ns)

                
        # Obtaining S1=1 and S2=2 to plot against each other to see symmetry and potential features

        s2_area_first=[0]
        s2_bottom_area_first =[0]
        s2_area_second=[0]
        s2_bottom_area_second=[0]
        s1_start_time= 0
        s2_start_time_1= 0
        s2_start_time_2= 0
        first_s2_start_time = 0
        veto_time = 0
        area_of_s1 = 0


        # Code to identify first and second S2 pulses and record the area. 

        if tpc_pulse_counter_threshold['S1'] == 1 and tpc_pulse_counter_threshold['S2'] == 2:
            for ipulse in range(tpc.nPulses):
                area = tpc.pulseArea_phd[ipulse] # total pulse area in detected photons (phd)
                classification = tpc.classification[ipulse] # classification: S1, S2, SE, SPE, MPE, Other
                start_time_ns = tpc.pulseStartTime_ns[ipulse] # pulse start
                end_time_ns = tpc.pulseEndTime_ns[ipulse] # pulse end
                bottom_area = tpc.bottomArea_phd[ipulse]

                if classification == 'S1' and area >= s1_threshold:
                    s1_start_time = start_time_ns
                    s1_pulse_indices=ipulse
                    area_of_s1 = area
                    s1_bottom_area = bottom_area
                if classification == 'S2' and area >= s2_threshold and s2_start_time_1 == 0 and area >= area_cut_s2 :
                    s2_start_time_1 = start_time_ns
                    s2_first_pulse_indices=ipulse
                if classification == 'S2' and area >= s2_threshold and s2_start_time_1 != 0 and area>= area_cut_s2:
                    s2_start_time_2 = start_time_ns 
                    s2_second_pulse_indices =ipulse            
                if classification == 'S2' and area >= s2_threshold and area >= area_cut_s2:
                    if s2_area_first[0] == 0:
                        s2_area_first[0] = area
                        s2_bottom_area_first[0] = bottom_area     
                    else:
                        s2_area_second[0] = area
                        s2_bottom_area_second[0] = bottom_area
        # Removing S1 events which occur close to S2 or after and S2 
        if s2_start_time_2 >= s2_start_time_1:
            first_s2_start_time = s2_start_time_1
        else:
            first_s2_start_time = s2_start_time_2
        
        veto_time = first_s2_start_time - s1_start_time
        if veto_time > 50:
            #Correcting areas to append to the area arrays for S2
            #start times
            s1_start_time_us = tpc.pulseStartTime_ns[s1_pulse_indices]/1000.0
            s2_first_start_time_us = tpc.pulseStartTime_ns[s2_first_pulse_indices]/1000.0
            s2_second_start_time_us = tpc.pulseStartTime_ns[s2_second_pulse_indices]/1000.0
            #drift times
            drift_time_s2_first_us = s2_first_start_time_us - s1_start_time_us
            drift_time_s2_second_us = s2_second_start_time_us - s1_start_time_us
            #position in x-y of S2 
            s2_x_cm_first = tpc.s2Xposition_cm[s2_first_pulse_indices]
            s2_y_cm_first = tpc.s2Yposition_cm[s2_first_pulse_indices]
            s2_x_cm_second = tpc.s2Xposition_cm[s2_second_pulse_indices]
            s2_y_cm_second = tpc.s2Yposition_cm[s2_second_pulse_indices]
            s2_r2_cm_first = (s2_x_cm_first*s2_x_cm_first)+(s2_y_cm_first*s2_y_cm_first)
            s2_r2_cm_second = (s2_x_cm_second*s2_x_cm_second)+(s2_y_cm_second*s2_y_cm_second)

            #Correcting S1 area
            # S1 is a mix of light from both S2 as gamma hit PMT instantly
            #Thus to correct for S1 area need to use weighted average of the drift time between
            #Both S2 and S1. Thus going to use a weighted average.
            #i.e. (drift_s2_1*s2_1_area + drift_s2_2*s2_s2_area )/ (s2_1_area+s2_2_area)
            #Using s2 as a relative proxy for S1 as each contribution to S1 is unknown

            s1_drift_proxy = ((drift_time_s2_first_us*s2_area_first[0])*(drift_time_s2_second_us*s2_area_second[0]))/(s2_area_first[0]+s2_area_second[0])
            s2_pos_proxy = ((s2_r2_cm_first*s2_area_first[0])*(s2_r2_cm_second*s2_area_second[0]))/(s2_area_first[0]+s2_area_second[0])

            s1corr_bin = s1corrections_h2.FindBin(s2_pos_proxy,s1_drift_proxy)
            s1corr_factor = s1corrections_h2.GetBinContent(s1corr_bin)
            if s1corr_factor == 0:
                s1corr_factor = 1  # i.e. no corrections
            s1_area_corrected = area_of_s1/s1corr_factor
            s1_bottom_area_corrected = s1_bottom_area/s1corr_factor

        


            # Correcting S2 areas

            s2corr_bin_first = s2corrections_h2.FindBin(s2_x_cm_first,s2_y_cm_first)
            s2corr_factor_first = s2corrections_h2.GetBinContent(s2corr_bin_first)
            # 1026.1 is the attenuation length in us
            s2eff_factor_first = ROOT.TMath.Exp(drift_time_s2_first_us/1026.1)
            if s2corr_factor_first == 0:
                s2corr_factor_first =1
            s2_area_correction_first = (s2eff_factor_first/s2corr_factor_first)

            
            s2corr_bin_second = s2corrections_h2.FindBin(s2_x_cm_second,s2_y_cm_second)
            s2corr_factor_second = s2corrections_h2.GetBinContent(s2corr_bin_second)
            s2_eff_factor_second = ROOT.TMath.Exp(drift_time_s2_second_us/1026.1)
            if s2corr_factor_second == 0:
                s2corr_factor_second =1
            s2_area_correction_second = (s2_eff_factor_second/s2corr_factor_second)





            #Appending to outside main entry loop to store and use to plot graphs.
            if s2_area_first[0] != 0:  
                first_s2_areas.append(s2_area_first[0] * s2_area_correction_first)
                s2_bottom_area_first_corrected = s2_bottom_area_first[0]*s2_area_correction_first
                if math.isnan(s2_bottom_area_first_corrected):
                    s2_bottom_area_first_corrected = 0
                first_s2_bottom_area_list_corrected.append(s2_bottom_area_first_corrected)
                s1_area_list_corrected.append(s1_area_corrected)
                if math.isnan(s1_bottom_area_corrected):
                    s1_bottom_area_corrected = 0
                s1_bottom_corrected_areas_list.append(s1_bottom_area_corrected)


            if s2_area_second[0] != 0:
                second_s2_areas.append(s2_area_second[0] * s2_area_correction_second)
                s2_bottom_area_second_corrected =s2_bottom_area_second[0]*s2_area_correction_second
                if math.isnan(s2_bottom_area_second_corrected):
                    s2_bottom_area_second_corrected = 0
                second_s2_bottom_area_list_corrected.append(s2_bottom_area_second_corrected)
                









            #Now all areas variables to be used in plots have been corrected for 
            #S1_Area,S2_areas for first and second, and the bottoms areas for each of the pulses just listed.
            #The corrected S2_areas were appended back into the inital arrays to can use them to plot

            



        
            




    
        #Creating dict keys for end summary print out.
        dict_key = "{} S1 and {} S2".format(tpc_pulse_counter_threshold['S1'], tpc_pulse_counter_threshold['S2'])
        if dict_key in types_of_scatter_test:
            types_of_scatter_test[dict_key] +=1
        else:
            types_of_scatter_test[dict_key] = 1
          
        #Printing event TPC information and labeling as single scatter / multi-scatter   
        if verbose:
            print("-> Total number of TPC S1/S2 pulses: {}".format(tpc_pulse_counter))
            print('-> The total number of TPC S1/S2 pulses above the 10% thrshold is: {}'.format(tpc_pulse_counter_threshold))
            if tpc_pulse_counter_threshold['S1'] == 1 and tpc_pulse_counter_threshold['S2'] == 1:
                print("-> This event is likely to be a single scatter ")
                event_type = 'TPC Single Scatter'
            elif tpc_pulse_counter_threshold['S1'] ==1 and tpc_pulse_counter_threshold['S2'] > 1:
                print('-> This event is likely to be a multi-scatter')
                event_type = 'TPC Multi-scatter'
            elif tpc_pulse_counter_threshold['S1'] == 0 and tpc_pulse_counter_threshold['S2'] == 0:
                print("-> No pulses were recorded in TPC")
                event_type = "TPC No Pulse" 

            







        #types_of_scatter_test['Single Scatter'] = types_of_scatter_test['S1 and S2']



        # Loop over Skin detector pulses
        skin = events.pulsesSkin
        if verbose: print(" -> looping over {} Skin pulses:".format(skin.nPulses))
        for ipulse in range(skin.nPulses):
            id = skin.pulseID[ipulse] # just number of pulse
            area = skin.pulseArea_phd[ipulse] # total pulse area in detected photons (phd)
            start_time_ns = skin.pulseStartTime_ns[ipulse] # pulse start
            # only need the area and start time of pulses to decide if should veto event
            #Checking coincidence window
            if event_type == 'TPC Single Scatter' or event_type == 'TPC Multi-scatter':
                if (start_time_of_s1-coincidence_window) <=  start_time_ns <=  (start_time_of_s1+coincidence_window) and area >= max_s1_area*area_decrease:  
                    overlap ='Overlap of S1 signal in TPC and Skin'
                    overlaps_in_TPC_and_outer_detector['Overlaps in skin'] +=1
                    overlap_event_id_array.append(id)

                else:
                    overlap ='No overlap in S1 or Skin'
            if verbose:
                print("   -> {}: {:0.2f} phd, start {} ns".format(
                    id, area, start_time_ns))


        # Loop over OD detector pulse
        od = events.pulsesODLG
        if verbose: print(" -> looping over {} OS pulses:".format(od.nPulses))
        for ipulse in range(od.nPulses):
            id = od.pulseID[ipulse] # just number of pulse
            area = od.pulseArea_phd[ipulse] # total pulse area in detected photons (phd)
            start_time_ns = od.pulseStartTime_ns[ipulse] # pulse start
            # only need the area and start time of pulses to decide if should veto event
            if event_type == 'TPC Single Scatter' or event_type == 'TPC Multi-scatter':
                if (start_time_of_s1-coincidence_window) <=  start_time_ns <=  (start_time_of_s1+coincidence_window) and area >= max_s1_area*area_decrease:
                    overlap ='Overlap of S1 signal in TPC and Outer Detector'
                    overlaps_in_TPC_and_outer_detector['Overlaps in OD'] +=1
                    overlap_event_id_array.append(id)
                else:
                    overlap ='No overlap in S1 and outer detector'
            if verbose:
                print("   -> {}: {:0.2f} phd, start {} ns".format(
                    id, area, start_time_ns))



    #Dataframe for df1. Breaks down pulse catagories and frequency 
    df = pd.DataFrame(types_of_scatter_test.items(),columns=['Pulses', 'Frequency'])

    #Generating lower dimension df2 from df1 for easy reading ( I.e- collecting multiscatters and single scatters and showing only S1 and only S2 events)
    pdToListKeys = list(df['Pulses'])
    pdToListValues=list(df['Frequency'])

    dict_for_df2={
        'Single-scatter':0,
        'Multi-scatter':0,
        'Only s1':0,
        'Only s2':0,
        'No pulse':0,
        'pile up':0,
        'check:sum of frequencies':0,
    }

    
    List_Of_Event_Types =[]
    for item in pdToListKeys:
        s1_s2=([int(s) for s in item.split() if s.isdigit()])  # Extract numerical information from the event types in DF1 to generated DF2
        List_Of_Event_Types.append(s1_s2)
    for index in range(len(List_Of_Event_Types)):
        List_Of_Event_Types[index].append(pdToListValues[index])
    for index in range(len(List_Of_Event_Types)):
        if List_Of_Event_Types[index][0] == 0 and List_Of_Event_Types[index][1] > 0:
            dict_for_df2['Only s2'] += List_Of_Event_Types[index][2]
        elif List_Of_Event_Types[index][0] > 0 and List_Of_Event_Types[index][1] == 0:
            dict_for_df2['Only s1'] += List_Of_Event_Types[index][2]
        elif List_Of_Event_Types[index][0] == 0 and List_Of_Event_Types[index][0] == 0:
            dict_for_df2['No pulse'] += List_Of_Event_Types[index][2]
        elif List_Of_Event_Types[index][0] == 1 and List_Of_Event_Types[index][1] == 1:
            dict_for_df2['Single-scatter'] += List_Of_Event_Types[index][2]
        elif List_Of_Event_Types[index][0] == 1 and List_Of_Event_Types[index][1] > 0:
            dict_for_df2['Multi-scatter']+=List_Of_Event_Types[index][2]
        elif List_Of_Event_Types[index][0] > 1 and List_Of_Event_Types[index][1] > 0:
            dict_for_df2['pile up'] += List_Of_Event_Types[index][2]
        else:
            pass
    dict_for_df2_values = dict_for_df2.values()
    dict_for_df2_values_check = sum(dict_for_df2_values)
    dict_for_df2['check:sum of frequencies'] += dict_for_df2_values_check        
    df2=pd.DataFrame(dict_for_df2.items(),columns=['Event','Frequency'])


    #Making % matrix so show contribution of each pulse pairing.

    dict_for_matrix={
        'S1=0 and S2=0':0,
        'S1=0 and S2=1':0,
        'S1=0 and S2>1':0,
        'S1=1 and S2=0':0,
        'S1=1 and S2=1':0,
        'S1=1 and S2>1':0,
        'S1>1 and S2=0':0,
        'S1>1 and S2=1':0,
        'S1>1 and S2>1':0
    }


    List_Of_Event_Types_matrix=[]
    for item in pdToListKeys:
        s1_s2_matrix=([int(s) for s in item.split() if s.isdigit()])
        List_Of_Event_Types_matrix.append(s1_s2_matrix)
    for index in range(len(List_Of_Event_Types_matrix)):
        List_Of_Event_Types_matrix[index].append(pdToListValues[index])
    for index in range(len(List_Of_Event_Types_matrix)):
        if List_Of_Event_Types_matrix[index][0] == 0 and List_Of_Event_Types_matrix[index][1] == 0:
            dict_for_matrix['S1=0 and S2=0'] += List_Of_Event_Types[index][2]
        elif List_Of_Event_Types_matrix[index][0] == 0 and List_Of_Event_Types_matrix[index][1] == 1:
            dict_for_matrix['S1=0 and S2=1'] += List_Of_Event_Types[index][2]
        elif List_Of_Event_Types_matrix[index][0] == 0 and List_Of_Event_Types_matrix[index][1] > 1:
            dict_for_matrix['S1=0 and S2>1'] += List_Of_Event_Types[index][2]
        elif List_Of_Event_Types_matrix[index][0] == 1 and List_Of_Event_Types_matrix[index][1] == 0:
            dict_for_matrix['S1=1 and S2=0'] += List_Of_Event_Types[index][2]
        elif List_Of_Event_Types_matrix[index][0] == 1 and List_Of_Event_Types_matrix[index][1] == 1:
            dict_for_matrix['S1=1 and S2=1'] += List_Of_Event_Types[index][2]
        elif List_Of_Event_Types_matrix[index][0] == 1 and List_Of_Event_Types_matrix[index][1] > 1:
            dict_for_matrix['S1=1 and S2>1'] += List_Of_Event_Types[index][2]
        elif List_Of_Event_Types_matrix[index][0] > 1 and List_Of_Event_Types_matrix[index][1] == 0:
            dict_for_matrix['S1>1 and S2=0'] += List_Of_Event_Types[index][2]
        elif List_Of_Event_Types_matrix[index][0] > 1 and List_Of_Event_Types_matrix[index][1] == 1:
            dict_for_matrix['S1>1 and S2=1'] += List_Of_Event_Types[index][2]
        elif List_Of_Event_Types_matrix[index][0]  >1 and List_Of_Event_Types_matrix[index][1] > 0:
            dict_for_matrix['S1>1 and S2>1'] += List_Of_Event_Types[index][2]
        else:
            pass
    total_dict3_values=0
    for value in dict_for_matrix:
        total_dict3_values = total_dict3_values + dict_for_matrix[value]
    for value in dict_for_matrix:
        dict_for_matrix[value] =(float)(dict_for_matrix
        [value]*100/total_dict3_values)

    df3=pd.DataFrame(dict_for_matrix.items(),columns=['Event','Percentage'])

    #Plotting Scatter Graph for pulses areas vs (width,asymmetry,pmompt frac)


    #AT THIS POINT GOING TO DO A CHECK  -> Below are the variables we want to check for Nan and Inf
    #   first_s2_areas=[]
    #   second_s2_areas=[]
    #   s1_area_list_corrected=[]
   #    s1_bottom_corrected_areas_list=[]
   #   first_s2_bottom_area_list_corrected=[]
  ##  second_s2_bottom_area_list_corrected=[]

    #checking how many errors are in the array
    #print(first_s2_bottom_area_list_corrected)
 
    #input("press anything to cont")


    
    # AREA  VS WIDTH
    plt.scatter(s2_pulse_width,np.log10(s2_pulse_area),c='blue',label='S2',s=0.5,marker='o')
    plt.scatter(s1_pulse_width,np.log10(s1_pulse_area),c='red',label='S1',s=0.3,marker='o')
    plt.scatter(spe_pulse_width,np.log10(spe_pulse_area),c='green',label='SPE',s=0.3,marker='o')
    plt.scatter(se_pulse_width,np.log10(se_pulse_area),c='black',label='SE',s=0.3,marker='o')
    plt.scatter(mpe_pulse_width,np.log10(mpe_pulse_area),c='darkviolet',label='MPE',s=0.3,marker='o')
    plt.xlim(0,11000)
    plt.xlabel('width (ns)')
    plt.ylabel('Log10[Area] (phd)')
    plt.legend(loc='lower right')
    plt.title('Pusle Width Vs Area for Catagorised Pusles')
    plt.savefig('scatter_width_vs_area.png')
    plt.clf()

    #Area Vs TBA 
    plt.scatter(s1_tba,np.log10(s1_pulse_area),c='red',label='S1',s=0.3,marker='o')
    plt.scatter(s2_tba,np.log10(s2_pulse_area),c='blue',label='S2',s=0.3,marker='o')
    plt.scatter(spe_tba,np.log10(spe_pulse_area),c='green',label='SPE',s=0.3,marker='o')
    plt.scatter(se_tba,np.log10(se_pulse_area),c='black',label='SE',s=0.3,marker='o')
    plt.scatter(mpe_tba,np.log10(mpe_pulse_area),c='darkviolet',label='MPE',s=0.3,marker='o')
    plt.xlabel('TBA')
    plt.ylabel('Log10[Area] (phd)')
    plt.title('TBA Vs Area for Catagorised Pulses')
    plt.legend(loc='upper left')
    plt.savefig('scatter_tba_vs_area.png')
    plt.clf()

    #plotting pf50ns vs pulse area
    plt.scatter(s1_pf50ns,np.log10(s1_pulse_area),c='red',label='S1',s = 0.3, marker='o')
    plt.scatter(s2_pf50ns,np.log10(s2_pulse_area),c='blue',label='S2',s=0.3,marker='o')
    plt.scatter(spe_pf50ns,np.log10(spe_pulse_area),c='green',label='SPE',s=0.3,marker='o')
    plt.scatter(se_pf50ns,np.log10(se_pulse_area),c='black',label='SE',s=0.3,marker='o')
    plt.scatter(mpe_pf50ns,np.log10(mpe_pulse_area),c='darkviolet',label='MPE',s=0.3,marker='o')
    plt.ylabel('Log10[Area] (phd)')
    plt.xlabel('pf50ns')
    plt.xlim(0,1.0)
    plt.title('Pf50ns Vs Area for Catagorised Pulses')
    plt.legend(loc='upper right')
    plt.savefig('scatter_pf50ns_vs_area.png')
    plt.clf()

    #Plotting for S1=1,S2=2 events the first S2 area vs second S2 area

    plt.scatter(second_s2_areas,first_s2_areas,s=10,marker='o')
    lims=[
        np.min([plt.xlim(),plt.ylim()]),
        np.max([plt.xlim(),plt.ylim()])
    ]
    plt.plot(lims,lims,'k--',alpha=0.75,zorder=0)
    plt.ylabel('First S2 Area (phd)')
    plt.xlabel('Second S2 Area (phd)')
    plt.title(' 1st S2 vs 2nd S2 for S1=1,S2=2 Multi-Scatter Events')
    plt.savefig('scatter_s2_vs_s2.png')
    plt.clf()

    #Recreating S1=1,S2=2 scatter with 2D histogram.

    plt.hist2d(second_s2_areas,first_s2_areas,bins=(150,150),cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlabel('Second S2 Area (phd)')
    plt.ylabel('First S2 Area (phd)')
    plt.title('1st S2 vs 2nd S2 for S1=1,S2=2 Multi-scatter Events ')
    plt.savefig('histo_s2_vs_s2.png')
    plt.clf()

    # Last batch of 2D plots for the report. 

  

    # S1 AREA VS FIRST S2 AREA
    plt.hist2d(s1_area_list_corrected,first_s2_areas,bins=(150,150),cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlabel('S1 Area (phd)')
    plt.ylabel('1st S2 area (phd)')
    plt.title('S1 area VS 1st S2 area')
    plt.savefig('s1_area_vs_1st_s2_area.png')
    plt.clf()

    # S1 AREA VS SECOND S2 AREA 
    plt.hist2d(s1_area_list_corrected,second_s2_areas,bins=(150,150),cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlabel('S1 Area (phd)')
    plt.ylabel('2nd S2 area (phd')
    plt.title('S1 area VS 2st S2 area')
    plt.savefig('s1_area_vs_2nd_s2_area.png')
    plt.clf()

    #S1  VS 1ST S2 BOTTOM
    plt.hist2d(s1_area_list_corrected,first_s2_bottom_area_list_corrected,bins=(150,150),cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlim(0,2500)
    plt.ylim(0,500000)
    plt.xlabel('S1 area (phd)')
    plt.ylabel('1st S2 bottom area (phd)')
    plt.title('S1 VS 1st S2 bottom')
    plt.savefig('s1_vs_1st_s2_bottom.png')
    plt.clf()

    #S1  VS 2ND S2 BOTTOM  
    plt.hist2d(s1_area_list_corrected,second_s2_bottom_area_list_corrected,bins=(150,150),cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlim(0,2500)
    plt.ylim(0,500000)
    plt.xlabel('S1 area (phd)')
    plt.ylabel('2nd S2 bottom area (phd)') 
    plt.title('S1 VS 2nd S2 bottom ')  
    plt.savefig('s1_vs_2nd_s2_bottom.png')
    plt.clf()

    #S2 BOTTOM VS S2 AREA FIRST
    plt.hist2d(first_s2_bottom_area_list_corrected,first_s2_areas,bins=(150,150),cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlabel('1st S2 bottom area (phd)')
    plt.ylabel('1st  S2 area (phd)')
    plt.title('1st S2 bottom area VS 1st S2 area')
    plt.savefig('1st_s2_bottom_vs_1st_s2_area.png')
    plt.clf()

    #S2 BOTTOM VS S2 AREA SECOND

    plt.hist2d(second_s2_bottom_area_list_corrected,second_s2_areas,bins=(150,150),cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlabel('2nd S2 bottom area (phd)')
    plt.ylabel('2nd S2 area (phd)')
    plt.title('2nd S2 bottom area vs 2nd S2 area')
    plt.savefig('2nd_s2_bottom_vs_2nd_s2_area.png')
    plt.clf()

    #S2 VS S2 But for bottom

    plt.hist2d(second_s2_bottom_area_list_corrected,first_s2_bottom_area_list_corrected,bins=(150,150),cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlim(0,500000)
    plt.ylim(0,500000)
    plt.xlabel('2nd S2 bottom area (phd)')
    plt.ylabel('1st S2 bottom area (phd')
    plt.title('1st S2 vs 2nd S2 bottom areas')
    plt.savefig('s2_vs_s2_bottom.png')



    #Attempt to make  3D plot 


     
    # END PRINT OUT OF SUMMARY INFO:


    

    if verbose:
        threshold=threshold*100
        total_cataorgised_test=types_of_scatter_test.values()
        sum_total_cataorgised_test = sum(total_cataorgised_test)
        print("-> Summary of dataset: \n\n\n")

        print("->Looped over {} events and categorised {}/{} for an input threshold of: {}%".format(args.nevents,sum_total_cataorgised_test,args.nevents,threshold))
        print("===================================")
        print("-> Pulse Summary")   
        print(df)
        print("===================================")
        print("-> Lower dimension event classification")
        print(df2)
        print("===================================")
        print('Overlaps Data:{} '.format(overlaps_in_TPC_and_outer_detector))
        print("===================================")
        print(df3)

        #////// UNCOMMENT WHEN YOU WANT TO CHECK VALUES AND PRINT TO END OF PROGRAM ////// 


   
        
      
        

        

    input("Enter anything to exit?")     

if __name__ == "__main__":
    main()
