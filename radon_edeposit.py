#
# Example script to load a BACCARAT ROOT output, loop over the events
# and get information on tracks and steps.
#
# Setup:
# $ source /cvmfs/lz.opensciencegrid.org/LzBuild/latest/setup.sh 
#
# Run:
# $ python radon_edeposit_loop.py --help 
# $ python radon_edeposit_loop.py --nevents 10 --verbose /unix/darkmatter/jim/fordaniel/lz_Pb214_TPC_root_30160000.v0.root 
# $ python radon_edeposit_loop_test.py --nevents 10 --verbose /unix/darkmatter/jim/fordaniel/lz_Pb214_TPC_root_30160000.v0.root 
# ^ Was a test save to see if ssh connection through VSM was successful and how to execute code via terminal.
# $ python radon_edeposit_loop.py --nevents 10 --verbose /unix/darkmatter/jim/fordaniel/lz_Pb214_TPC_root_30160000.v0.root --plot
#

import sys
import argparse
from collections import Counter   

def main():    
   parser = argparse.ArgumentParser(description='Example loop over events for BACCARAT ROOT output.')
   parser.add_argument('--nevents', type=int, help='limit number of events to loop over')
   parser.add_argument('--verbose', action='store_true', help='print verbose for debugging')
   parser.add_argument('--plot', action='store_true', help='plot energy deposits for each event')
   parser.add_argument('input', help='path to input file')
   args = parser.parse_args()
    
   print("Loading events from:\n {}".format(args.input))
   if args.nevents:
      print("will limit loop to {} events".format(args.nevents))

   verbose = args.verbose
   plot = args.plot 
 
   # load ROOT after argparse as takes a while
   import ROOT
    
   tfile = ROOT.TFile.Open(args.input)
   if not tfile:
       print("Cannot open input file, exiting!")
       sys.exit(1)

   # Set ROOT plotting style
   ROOT.gROOT.SetStyle("Plain")

   # If requested cap the number events
   events = tfile.Get("Events")
   nevents = events.GetEntries()
   if args.nevents and args.nevents < nevents: 
      print("limiting number events from {} to {}".format(nevents, args.nevents))
      nevents = args.nevents

   # Hack because PyROOT is unable to access branch names with a . in them so make 
   # an alias and replace with an underscore. An alternative would be to use python's
   # getattr to fetch member variable based on a string - see commented out in event loop
   events.SetAlias("header_eventId", "header.eventId")
   events.SetAlias("header_eventLabel", "header.eventLabel")
   events.SetAlias("header_eventTypeLabel", "header.eventTypeLabel")
   events.SetAlias("primaries_particleName", "primaries.particleName")
   events.SetAlias("primaries_volumeName", "primaries.volumeName")
   events.SetAlias("primaries_pid", "primaries.pid")
   events.SetAlias("primaries_position_X_mm", "primaries.position_X_mm")
   events.SetAlias("primaries_position_Y_mm", "primaries.position_Y_mm")
   events.SetAlias("primaries_position_Z_mm", "primaries.position_Z_mm")
   events.SetAlias("deposits_volumeNames", "deposits.volumeNames")
   events.SetAlias("deposits_pids", "deposits.pids")
   events.SetAlias("deposits_EnergyDeps_keV", "deposits.EnergyDeps_keV")
   events.SetAlias("deposits_positions_x_mm", "deposits.positions_x_mm")
   events.SetAlias("deposits_positions_y_mm", "deposits.positions_y_mm")
   events.SetAlias("deposits_positions_z_mm", "deposits.positions_z_mm")
   events.SetAlias("deposits_parentTrackIds", "deposits.parentTrackIds")
   events.SetAlias("deposits_trackIds", "deposits.trackIds")
   # + any more that you need

   if plot:
       canvas = ROOT.TCanvas("canvas", "Energy deposit patterns", 10, 10, 900, 300)
       canvas.Divide(3,1)

   # Loop over events in file - N.B. for loops in python are slow should come back to this 
   for entry in range(nevents):
      events.GetEntry(entry)

      # Print some summary information about the event
      if verbose: 
         print("eventId {} {} {}:".\
           format(events.header_eventId, events.header_eventLabel, 
                  events.header_eventTypeLabel)) 
        
      # Print primary particles produced by decay
      n_primaries = len(events.primaries_particleName)
      if verbose:
         print(" -> has {} primaries:".format(n_primaries)) 

      for i_primary in range(n_primaries):
         prim_name = events.primaries_particleName[i_primary]
         prim_volume = events.primaries_volumeName[i_primary]
         prim_pid = events.primaries_pid[i_primary]
         if verbose:
            print("  -> particle: {} volume: {} pid: {}".format(prim_name, prim_volume, prim_pid))

      # Print information on energy deposits that these primary particles
      # lead to as they are tracked througout the detector
      # or you could use n_deposits = getattr(events, "deposits.volumeNames")
      n_deposits = len(events.deposits_volumeNames)
      if verbose:
        print(" -> which lead to {} deposits:".format(n_deposits))

      if plot:
          minx, maxx, nbins = -100.0, 100.0, 200
          hist_xy = ROOT.TH2D("hist_xy", "energy deposits XY plane;X [mm];Y [mm]", nbins, minx, maxx, nbins, minx, maxx)
          hist_xz = ROOT.TH2D("hist_xz", "energy deposits XZ plane;X [mm];Z [mm]", nbins, minx, maxx, nbins, minx, maxx)
          hist_yz = ROOT.TH2D("hist_yz", "energy deposits YZ plane;Y [mm];Z [mm]", nbins, minx, maxx, nbins, minx, maxx)

      # Keep track of total enery deposited by each particle species
      summed_energy_for_pid = Counter() 

      # For each recorded energy deposit print key info: how much energy, particle 
      # type, location
      #arrays to store path backtrack chains for either bi214 or beta decay particles in each event.
      array_beta=[]
      array_bi214=[]
      #deposit_position_list=[]
      #arrays to store data for energy weighted mean
      #array_w_m=[]
      #total_energy = 0
      #z_mm_total_energy_mm_product = 0
      #xy_total_energy_mm_product = 0
      beta_system_energy=0
      Bi214_system_energy=0
      beta_system_weights_x=0
      beta_system_weights_y=0
      beta_system_weights_z=0
      Bi214_system_weights_x=0
      Bi214_system_weights_y=0
      Bi214_system_weights_z=0


      for i_deposit in range(n_deposits):
          volume = events.deposits_volumeNames[i_deposit]
          pid = events.deposits_pids[i_deposit]
          energy_kev = events.deposits_EnergyDeps_keV[i_deposit]
          x_mm = events.deposits_positions_x_mm[i_deposit]
          y_mm = events.deposits_positions_y_mm[i_deposit]
          z_mm = events.deposits_positions_z_mm[i_deposit]
          parent_trackid = events.deposits_parentTrackIds[i_deposit]
          trackid = events.deposits_trackIds[i_deposit]
          # initial depoist checker 
          if (trackid == 3 ): # 1st track for Bi is always marked as 3.
              array_bi214.append(trackid)        
          for value in array_bi214:  # Appends all subsequent decays from Bi to an array
              if (parent_trackid == value):
                  array_bi214.append(trackid)
          if (trackid == 1):       # 1st track for beta is always marked as 1   
              array_beta.append(trackid)
          for value in array_beta: # Appends all subsequent decays from Beta to an array
              if (parent_trackid == value):
                  array_beta.append(trackid)   
            # The below loops check to see which array the current track is and assigns the correct initial deposit tag.
          for value in array_bi214:
              if (value == trackid):
                  initial_deposit = "Bi214"
          for value in array_beta:
              if (value == trackid):
                  initial_deposit = "Beta"
          # Need to calculate weighted mean of Beta and Bi and find total energy for each system.
          if (initial_deposit == "Beta"):
              beta_system_energy += energy_kev
              beta_system_weights_x += energy_kev*x_mm
              beta_system_weights_y += energy_kev*y_mm
              beta_system_weights_z += energy_kev*z_mm
          if (initial_deposit == 'Bi214'):
              Bi214_system_weights_x += energy_kev*x_mm
              Bi214_system_weights_y += energy_kev*y_mm
              Bi214_system_weights_z += energy_kev*z_mm
              Bi214_system_energy += energy_kev






         # Calculation for the energy weighted mean. INCORRECT BELOW
          
          #deposit_position = ((x_mm**2)+(y_mm**2)+(z_mm**2))**(0.5) # cartesian radial distance from origin
          #deposit_position_list.append(deposit_position)
          #enrgy_position_product = energy_kev * deposit_position
          #array_w_m.append(enrgy_position_product) # saving (energy * postion) valeus to an array
          #total_energy += energy_kev
          # saving the sums of the energy * z- distance and xy area
          #z_mm_total_energy_mm_product += (z_mm * energy_kev) 
          #xy_total_energy_mm_product += ((x_mm * y_mm) * energy_kev)





                
            
          summed_energy_for_pid[pid] += energy_kev

          # Get positions of first primary and use this to centre the plots of energy
          # deposits
          prim_X_mm = events.primaries_position_X_mm[0]
          prim_Y_mm = events.primaries_position_Y_mm[0]
          prim_Z_mm = events.primaries_position_Z_mm[0]

          if plot:
              hist_xy.Fill(x_mm-prim_X_mm, y_mm-prim_Y_mm) 
              hist_xz.Fill(x_mm-prim_X_mm, z_mm-prim_Z_mm)  
              hist_yz.Fill(y_mm-prim_Y_mm, z_mm-prim_Z_mm) 
         
          if verbose:
             print("  -> volume: {} pid: {} edep: {:0.5f} keV pos: {:0.4f}, {:0.4f}, {:0.4f} mm,(trk {}, parent trk {}, initial deposit: {})".format(volume, pid, energy_kev, x_mm, y_mm, z_mm, trackid, parent_trackid, initial_deposit))
      # Calculating the weighted energy position INCORRECT BELOW
      #sum_of_w_m = 0
      #for value in array_w_m:
         # sum_of_w_m = sum_of_w_m + value
      #energy_weighted_mean_poistion = ((1/total_energy) * sum_of_w_m)
      #z_mm_weighted_mean_position = z_mm_total_energy_mm_product * (1/total_energy) 
     # xy_weighted_mean_position = xy_total_energy_mm_product * (1/total_energy)


      #correct mean position results for each axis
      beta_mean_x_position = (beta_system_weights_x / beta_system_energy)
      beta_mean_y_position = (beta_system_weights_y / beta_system_energy)
      beta_mean_z_position = (beta_system_weights_z / beta_system_energy)
      Bi214_mean_x_position = (Bi214_system_weights_x / Bi214_system_energy)
      Bi214_mean_y_position = (Bi214_system_weights_y / Bi214_system_energy)
      Bi214_mean_z_position = (Bi214_system_weights_z / Bi214_system_energy)

      #absolute
      absolute_difference_position_x = Bi214_mean_x_position -beta_mean_x_position
      absolute_difference_position_y = Bi214_mean_y_position - beta_mean_y_position
      absolute_difference_position_z = Bi214_mean_z_position - beta_mean_z_position 
      # Print summed energy deposits
      if verbose:
          print("   -> summed of above deposits for each particle type:")
          for pid, energy in summed_energy_for_pid.items():
              print("      {}: {:0.2f} keV".format(pid, energy))
          #print("    -> The energy weighted mean position for this event is: {:0.4f} KeV ".format(energy_weighted_mean_poistion))
          #print("    -> Below, the difference between the energy weighted mean position and decay position for various dimensions is given:")

      if verbose:
          print("    -> Centre of energy information:\n  absolute difference between Bi214 and Beta : R = {:0.4f} X + {:0.4f} Y + {:0.4f} Z \n \
              ".format(absolute_difference_position_x,absolute_difference_position_y,absolute_difference_position_z))

          # printing difference in weighted mean energy position and location of each energy deposit  INCORRECT BELOW
      #if verbose:
          #for i_deposit in range(n_deposits):
           # absolute_difference = abs(energy_weighted_mean_poistion - deposit_position_list[i_deposit])
           # z_component_difference = z_mm_weighted_mean_position - (events.deposits_positions_z_mm[i_deposit]) 
           # xy_plane_difference = xy_weighted_mean_position - (events.deposits_positions_x_mm[i_deposit]*events.deposits_positions_y_mm[i_deposit])
           # print("    -> For deposit: {}:  absolute difference(3D): {:0.4f}mm, z-component difference: {:0.4f}mm, xy-plane difference: {:0.4f}mm".format(i_deposit+1, absolute_difference, z_component_difference,xy_plane_difference))


      if plot:
          canvas.cd(1)
          hist_xy.Draw()
          canvas.cd(2)
          hist_xz.Draw()
          canvas.cd(3)
          hist_yz.Draw()
          canvas.Update() 
          input("Enter anything to go to next event. Ctrl-C to exit.")     

if __name__ == "__main__":
    main()
