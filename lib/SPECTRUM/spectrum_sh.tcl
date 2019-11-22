#!/usr/bin/wish -f
frame .one -bg grey85 -width 300 -height 60
frame .two -bg grey85 -width 800 -height 50
frame .three -bg grey85 -width 500 -height 40
frame .four -bg grey85 -width 500 -height 40
frame .five -bg grey85 -width 500 -height 40
frame .six -bg grey85 -width 500 -height 40
frame .seven -bg grey85 -width 500 -height 40
frame .sevenB -bg grey85 -width 500 -height 40
frame .sevenC -bg grey85 -width 500 -height 40
frame .eight -bg grey85 -width 500 -height 40
frame .nine -bg grey85 -width 500 -height 40
frame .ten -bg grey85 -width 500 -height 40
frame .tenB -bg grey85 -width 500 -height 50
frame .tenC -bg grey85 -width 500 -height 50
frame .tenD0 -bg grey85 -width 500 -height 50
frame .tenD -bg grey85 -width 500 -height 50
frame .tenE -bg grey85 -width 500 -height 50
frame .tenF -bg grey85 -width 500 -height 40
frame .eleven -bg grey85 -width 500 -height 80
frame .twelve -bg grey85 -width 500 -height 40
label .one.caption -text " Spectrum Shell" -font {times 32 {italic}}
pack .one.caption -side left
label .one.caption2 -text "v1.02"
pack .one.caption2 -side left
frame .one.right -bg grey70 -width 150 -height 60
label .one.right.caption -text "tcl/tk version (C) R.O. Gray 2005 "
pack .one.right.caption -side left
pack .one.right -side right
label .two.select -text Select: -width 10
pack .two.select -side left
button .two.model -text model -command Getmodel
pack .two.model -side left
button .two.atom -text atom -command Getatom
pack .two.atom -side left
button .two.isotope -text isotope -command Getisotope
pack .two.isotope -side left
button .two.lines -text "line list" -command Getline
pack .two.lines -side left
button .two.output -text output -command Getoutput
pack .two.output -side left
label .three.modellabel -text "Model   " -width 12 -justify left
label .three.model -justify left -textvariable model -relief sunken -width 40
pack .three.modellabel .three.model -side left
label .four.atomlabel -text "Atom    " -width 12 -justify left
label .four.atom -justify left -textvariable atom -relief sunken -width 40
pack .four.atomlabel .four.atom -side left
label .five.isotopelabel -text "Isotope " -width 12 -justify left
label .five.isotope -justify left -textvariable isotope -relief sunken -width 40
label .five.isotopestar -text "*" -justify left
pack .five.isotopelabel .five.isotope .five.isotopestar -side left
label .six.linelabel -text "Line List " -width 12 -justify left
label .six.lines -justify left -textvariable lines -relief sunken -width 40
pack .six.linelabel .six.lines -side left
label .seven.outputlabel -text "Output " -width 12 -justify left
label .seven.output -justify left -textvariable output -relief sunken -width 40
pack .seven.outputlabel .seven.output -side left
label .sevenB.star -text "*Isotope selection not used unless isotope switch checked below" -font {fixed 8}
pack .sevenB.star -side left
label .sevenC.label -text "Computation Parameters:" -relief groove
pack .sevenC.label -side left
label .eight.startlabel -text "Start Wavelength " -justify left
entry .eight.startentry -textvariable start -width 9 -relief sunken -background white
label .eight.sA -text "A  " -justify left
label .eight.endlabel -text "End " -justify left
entry .eight.endentry -textvariable end -width 9 -relief sunken -background white
label .eight.eA -text "A  " -justify left
label .eight.spacelabel -text "Spacing " -justify left
entry .eight.spaceentry -textvariable space -width 6 -relief sunken -background white
label .eight.spA -text "A  " -justify left
pack .eight.startlabel .eight.startentry .eight.sA .eight.endlabel .eight.endentry .eight.eA .eight.spacelabel .eight.spaceentry .eight.spA -side left
label .nine.mtvlabel -text "Microturbulent Velocity " -justify left
entry .nine.mtventry -textvariable vt -width 4 -relief sunken -background white
label .nine.kms -text "km/s " -justify left
pack .nine.mtvlabel .nine.mtventry .nine.kms -side left
label .ten.switchlabel -text "Switches: " -relief groove
pack .ten.switchlabel -side left
radiobutton .tenB.0 -variable flagm -text "Normalized Intensity" -value 0
radiobutton .tenB.1 -variable flagm -text "Flux" -value 1
radiobutton .tenB.2 -variable flagm -text "Center of Disk" -value 2
pack .tenB.0 .tenB.1 .tenB.2 -side left
checkbutton .tenC.0 -text "Isotope" -variable flagi
checkbutton .tenC.1 -text "Atlas9" -variable flag9
pack .tenC.0 .tenC.1 -side left
label .tenD0.mtv -text "Macroturbulent Velocity: " -relief groove
checkbutton .tenD0.0 -variable flagmt
entry .tenD0.entry -textvariable mt -width 6 -relief sunken -background white
label .tenD0.kms -text "km/s"
pack .tenD0.mtv .tenD0.0 .tenD0.entry .tenD0.kms -side left
label .tenD.rotate -text "Rotation: " -relief groove
checkbutton .tenD.0 -variable flagr
label .tenD.vsini -text "  vsini "
entry .tenD.entry -textvariable vsini -width 6 -relief sunken -background white
label .tenD.kms -text "km/s"
pack .tenD.rotate .tenD.0 .tenD.vsini .tenD.entry .tenD.kms -side left
label .tenE.smooth -text "Smooth: " -relief groove
checkbutton .tenE.0 -variable flags
label .tenE.rlabel -text " Resolution "
entry .tenE.entry1 -textvariable res -width 6 -relief sunken -background white
label .tenE.rA -text "A  "
label .tenE.slabel -text " Output spacing* "
entry .tenE.entry2 -textvariable spc -width 6 -relief sunken -background white
label .tenE.rA2 -text "A  "
pack .tenE.smooth .tenE.0 .tenE.rlabel .tenE.entry1 .tenE.rA \
.tenE.slabel .tenE.entry2 .tenE.rA2 -side left
label .tenF.star -text "*Output spacing must be integer multiple of computation spacing" -font {fixed 8}
pack .tenF.star -side left
button .eleven.run -text "Execute" -command Run
pack .eleven.run -side left -ipady 10 -padx 40
button .eleven.abort -text "Abort" -command Abort
pack .eleven.abort -side left -ipady 10 -ipadx 10 -padx 40
button .eleven.quit -text "Quit" -command Quit
pack .eleven.quit -side left -ipady 10 -ipadx 10 -padx 40
label .twelve.label -text "Status:  " -justify left
label .twelve.status -justify left -textvariable status -relief sunken -width 25
label .twelve.progress1 -text "Progress:  " -justify left
checkbutton .twelve.c1 -text " " -variable flagp
label .twelve.progress2 -justify left -textvariable wavep -relief sunken -width 15
pack .twelve.label .twelve.status .twelve.progress1 .twelve.c1 .twelve.progress2 -side left
pack .one .two .three .four .five .six .seven .sevenB .sevenC .eight .nine \
.ten .tenB .tenC .tenD0 .tenD .tenE .tenF .eleven .twelve -side top -fill x
set model "Not Selected         "
set atom  "Not Selected         "
set isotope "Not Selected         "
set output  "Not Selected         "
set lines   "Not Selected         "
set start 3800.0
set wavep $start
set end 4600.0
set space 0.02
set vt 2.0
set mt 1.0
set flagm 0
set flagi 0
set flag9 0
set flagr 0
set flags 0
set flagp 0
set flagmt 0
set status "Expecting Input"
set modelflag 0
set atomflag 0
set isotopeflag 0
set listflag 0
set outflag 0
set n 10
set k -1
set pids 100
set pwd [pwd]
set vsini 0
set spc 0.5
set res 2.0

proc Pstat {pids} {
  global k flagp wavep
  set n 2
  while {1} {
   after 500
   catch {set n [exec ps hp $pids > tmp.wc]}
   set m [exec wc -l tmp.wc]
   scan $m {%d %s} k tmpstr
   if {$k == 0} {
    return
   }
   if {$flagp == 1} {
    set line [exec tail -n 2 tmp.out | head -n 1]
    scan $line "%f %f" wavep inten
   }
   update
  }
}

proc Getmodel {} {
   global model modelflag pwd
   set types {
       {{mod files} {.mod} }
       {{Atlas12 files} {.m12} }
       {{Atlas9 files} {.width9}}
       {{dat files} {.dat}}
       {{All Files} {*}}
}
   set model [tk_getOpenFile -filetypes $types -title "Select Model" -initialdir $pwd]
   set .three.model.text $model
   if {$model != ""} {
     set modelflag 1
   } else {
     set modelflag 0
   }
}

proc Getatom {} {
   global atom atomflag pwd
   set types {
       {{dat files} {.dat}}
       {{All Files} {*}}
}
   set atom [tk_getOpenFile -filetypes $types -title "Select Atomic Data File" -initialdir $pwd]
   set .four.atom.text $atom
   if {$atom != ""} {
     set atomflag 1
   } else {
     set atomflag 0
   }
}

proc Getisotope {} {
   global isotope isotopeflag pwd
   set types {
       {{dat files} {.iso}}
       {{All Files} {*}}
}
   set isotope [tk_getOpenFile -filetypes $types -title "Select Isotope File" -initialdir $pwd]
   set .five.isotope.text $isotope
   if {$isotope != ""} {
     set isotopeflag 1
   } else {
     set isotopeflag 0
   }
}

proc Getline {} {
   global lines listflag
   set types {
       {{dat files} {.lst}}
       {{All Files} {*}}
}
   set lines [tk_getOpenFile -filetypes $types -title "Select Line List"]
   set .six.lines.text $lines
   if {$lines != ""} {
     set listflag 1
   } else {
     set listflag 0
   }
}

proc Getoutput {} {
   global output outflag pwd
   set types {
       {{All Files} {*}}
   }
   set output [tk_getSaveFile -filetypes $types -title "Output File" -initialdir $pwd]
   set .seven.output.text $output
   if {$output != ""} {
     set outflag 1
   } else {
     set outflag 0
   }
}

proc Run {} {
   global model atom isotope lines output start end space vt status
   global modelflag atomflag isotopeflag listflag outflag flagi n k pids
   global flagm pwd flagr flags flag9 flagi flagmt mt vsini res spc wavep
   if {$modelflag == 0 || $atomflag == 0 || $listflag == 0 || $outflag == 0} {
     tk_messageBox -type ok -message "Input not complete!" -icon question
     return
   }
   if {$isotopeflag == 0 && $flagi == 1} {
     tk_messageBox -type ok -message "Select Isotope file!" -icon question
     return
   }
   if {$k > 0} {
     tk_messageBox -type ok -message "Spectrum already running!" -icon question
     return
   }
   cd $pwd
   set status "Executing..."
   set wavep $start
   set rspId [open tmp.rsp w]
   puts $rspId $model
   puts $rspId $lines
   puts $rspId $atom
   if {$flagi == 1} {
      puts $rspId $isotope
   }
   puts $rspId "tmp.out"
   puts $rspId $vt
   if {$flagm == 2} {
      puts $rspId 1.0
   }
   puts $rspId $start,$end
   puts $rspId $space
   close $rspId
   if {$flagm == 0} {
       if {$flagi == 0 && $flag9 == 0} {
         set pids [exec spectrum an < tmp.rsp > out.out &]
       } elseif {$flagi == 1 && $flag9 == 0} {
         set pids [exec spectrum ani < tmp.rsp > out.out &]
       } elseif {$flagi == 0 && $flag9 == 1} {
         set pids [exec spectrum ant < tmp.rsp > out.out &]
       } else {
         set pids [exec spectrum aint < tmp.rsp > out.out &]
       }
   } elseif {$flagm == 1} {
       if {$flagi == 0 && $flag9 == 0} {
         set pids [exec spectrum anf < tmp.rsp > out.out &]
       } elseif {$flagi == 1 && $flag9 == 0} {
         set pids [exec spectrum anif < tmp.rsp > out.out &]
       } elseif {$flagi == 0 && $flag9 == 1} {
         set pids [exec spectrum antf < tmp.rsp > out.out &]
       } else {
         set pids [exec spectrum anift < tmp.rsp > out.out &]
       }
   } elseif {$flagm == 2} {
       if {$flagi == 0 && $flag9 == 0} {
         set pids [exec spectrum anm < tmp.rsp > out.out &]
       } elseif {$flagi == 1 && $flag9 == 0} {
         set pids [exec spectrum anim < tmp.rsp > out.out &]
       } elseif {$flagi == 0 && $flag9 == 1} {
         set pids [exec spectrum antm < tmp.rsp > out.out &]
       } else {
         set pids [exec spectrum animt < tmp.rsp > out.out &]
       } 
   }
   Pstat $pids
    if {$flagmt == 0} {
      if {$flagr == 1 && $flags == 0} {
         set n [exec avsini tmp.out $output $vsini 0.6 $space]
         file delete tmp.out
      } elseif {$flagr == 1 && $flags == 1} {
         set n [exec avsini tmp.out tmpv.out $vsini 0.6 $space]
         set n [exec smooth2 tmpv.out $output $space $res $spc]
         file delete tmp.out
         file delete tmpv.out
      } elseif {$flagr == 0 && $flags == 1} {
         set n [exec smooth2 tmp.out $output $space $res $spc]
         file delete tmp.out
      } elseif {$flagr == 0 && $flags == 0 } {
         file rename tmp.out $output
      }
    } elseif {$flagmt == 1} {
	set n [exec macturb tmp.out tmp1.out $space $mt]
      if {$flagr == 1 && $flags == 0} {
         set n [exec avsini tmp1.out $output $vsini 0.6 $space]
         file delete tmp.out
      } elseif {$flagr == 1 && $flags == 1} {
         set n [exec avsini tmp1.out tmpv.out $vsini 0.6 $space]
         set n [exec smooth2 tmpv.out $output $space $res $spc]
         file delete tmp.out
         file delete tmpv.out
      } elseif {$flagr == 0 && $flags == 1} {
         set n [exec smooth2 tmp1.out $output $space $res $spc]
         file delete tmp.out
      } elseif {$flagr == 0 && $flags == 0 } {
         file rename tmp1.out $output
      }
      file delete tmp1.out
  }

   set status "Finished/Accepting Input" 
}

proc Quit {} {
   global k pids
   if {$k > 0} {
     set choice [tk_messageBox -type yesno \
         -message "Spectrum is running! Quit?" -icon question]
     switch -- $choice {
       yes {
         exec kill $pids
         exit
       }
       no return
     }
   }
   exit
}

proc Abort {} {
   global k pids
   if {$k > 0} {
     set choice [tk_messageBox -type yesno \
         -message "Abort?" -icon question]
     switch -- $choice {
       yes {
         exec kill $pids
         return
       }
       no return
     }
   }
   return
}
