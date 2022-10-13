proc contacts {sel1 sel2 sel3 sel4 out} {
  set numFrame [molinfo top get numframes]
  #Defining references
  set SEMG1 [atomselect top "$sel1"]
  set D71   [atomselect top "$sel2"]
  set N113  [atomselect top "$sel3"]
  set N114  [atomselect top "$sel4"]
  set outfile  [open ${out}.dat w]
  for {set f 0} {$f <= $numFrame} {incr f} {
    $SEMG1 frame $f
    $SEMG1 update
    $D71 frame $f
    $D71 update
    $N113 frame $f
    $N113 update
    $N114 frame $f
    $N114 update
    #Set time in ns
    set time [format "%.2f" [expr $f * 0.01]]
    #Get Z value for each residue
    set atomref [$SEMG1 get index]
    set atoms1  [$D71 get index]
    set atoms2  [$N113 get index]
    set atoms3  [$N114 get index]
    set distances1 ""
    foreach aref $atomref {
      foreach atom $atoms1 {
        set dd "$aref $atom"
        set val [measure bond $dd frame $f]
        lappend distances1 $val
      }
    }
    set distances2 ""
    foreach aref $atomref {
      foreach atom $atoms2 {
        set dd "$aref $atom"
        set val [measure bond $dd frame $f]
        lappend distances2 $val
      }
    }
    set distances3 ""
    foreach aref $atomref {
      foreach atom $atoms3 {
        set dd "$aref $atom"
        set val [measure bond $dd frame $f]
        lappend distances3 $val
      }
    }
    set mindist1 [lindex [lsort -real $distances1] 0]
    set mindist2 [lindex [lsort -real $distances2] 0]
    set mindist3 [lindex [lsort -real $distances3] 0]
    puts $outfile "$time $mindist1 $mindist2 $mindist3"
  }
  close $outfile
}

#set rep 3
#mol new ../../../complex_renum.gro type gro waitfor all
mol new 500_renum.pdb type pdb waitfor all
#Load trajectories
mol addfile jump.xtc type xtc waitfor all
#Remove the first frame
animate delete beg 0 end 1
#Call the function
#Use two atom selections with the first one being the reference and the second one being the target
#Then the third one is the output file
contacts "protein and resid 235" "protein and resid 71" "protein and resid 113" "protein and resid 114" "mindist_71_113_114_235"
# "water and name OH2" "residues_in_contact_to_sodium_rep_${rep}.dat" "water_z_position_rep_${rep}.dat" "sodium_z_position_rep_${rep}.dat"
quit
