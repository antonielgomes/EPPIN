
proc contacts {sel1 sel2 out} {
  set numFrame [molinfo top get numframes]
  #Defining references
  set prot [atomselect top "$sel1 and noh"]
  set pres [atomselect top "(same residue as within 3.5 of ($sel2 and noh)) and ($sel1 and noh)"]
  set lig  [atomselect top "$sel2 and noh"]
#  set nres [lsort -integer -unique [$pep get resid]]
  set outfile  [open ${out}_all.dat w]
  for {set f 0} {$f <= $numFrame} {incr f} {
    $lig frame $f
    $lig update
    $prot frame $f
    $prot update
    $pres frame $f
    $pres update
    #Set time in ns
    set time [format "%.2f" [expr $f * 0.01]]
    set nearsel [lsort -integer -unique [$pres get resid]]
    #foreach n $nearsel {
      #set nearresname [atomselect top "$sel1 and resid ${n}"]
      #$nearresname frame $f
      #$nearresname update
      #set nrname [lindex [$nearresname get resname] 0]
      #set v "${n}"
      #puts $v
    #}
    puts $outfile "$time $nearsel"
  }
  close $outfile
}

mol new complex_ep012_renum.pdb type pdb waitfor all
#Load trajectories
mol addfile jump.xtc type xtc waitfor all
#Remove the first frame
animate delete beg 0 end 1
#Call the function
#Use two atom selections with the first one being the reference and the second one being the target
#Then the third one is the output file
contacts "protein" "resname E055" "contacts_lig"
# "water and name OH2" "residues_in_contact_to_sodium_rep_${rep}.dat" "water_z_position_rep_${rep}.dat" "sodium_z_position_rep_${rep}.dat"
quit
