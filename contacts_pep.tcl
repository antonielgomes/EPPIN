
proc contacts {sel1 sel2 out} {
  set numFrame [molinfo top get numframes]
  #Defining references
  set prot [atomselect top "$sel1 and noh"]
  set pres [atomselect top "(same residue as within 3.5 of ($sel2 and noh)) and ($sel1 and noh)"]
  set pep  [atomselect top "$sel2 and noh"]
  set nres [lsort -integer -unique [$pep get resid]]
  set outfile  [open ${out}_all.dat w]
  for {set f 0} {$f <= $numFrame} {incr f} {
    $pep frame $f
    $pep update
    $prot frame $f
    $prot update
    $pres frame $f
    $pres update
    #Set time in ns
    set time [format "%.2f" [expr $f * 0.01]]
    set nearsel [lsort -integer -unique [$pres get resid]]
    puts $outfile "$time $nearsel"
  }
  close $outfile
}

mol new complex_renum.gro type gro waitfor all
#Load trajectories
mol addfile jump.xtc type xtc waitfor all
#Remove the first frame
animate delete beg 0 end 1
#Call the function
#Use two atom selections with the first one being the reference and the second one being the target
#Then the third one is the output file
contacts "protein and index 0 to 1199" "protein and index 1200 to 1391" "contacts_pep"
# "water and name OH2" "residues_in_contact_to_sodium_rep_${rep}.dat" "water_z_position_rep_${rep}.dat" "sodium_z_position_rep_${rep}.dat"
quit
