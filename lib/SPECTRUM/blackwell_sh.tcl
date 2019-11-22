#!/usr/bin/wish -f
frame .one -bg grey85 -width 300 -height 60
frame .oneB -bg grey85 -width 500 -height 40
frame .two -bg grey85 -width 800 -height 50
frame .three -bg grey85 -width 500 -height 40
frame .four -bg grey85 -width 500 -height 40
frame .five -bg grey85 -width 500 -height 40
frame .six -bg grey85 -width 500 -height 40
frame .seven -bg grey85 -width 500 -height 40
frame .eight -bg grey85 -width 500 -height 40
frame .nine -bg grey85 -width 500 -height 40
frame .nineB -bg grey85 -width 500 -height 50
frame .nineC -bg grey85 -width 500 -height 50
frame .eleven -bg grey85 -width 500 -height 80
frame .twelve -bg grey85 -width 500 -height 40
label .one.caption -text " Blackwell Shell" -font {times 32 {italic}}
pack .one.caption -side left
frame .one.right -bg grey70 -width 150 -height 60
label .one.right.caption -text "tcl/tk version (C) R.O. Gray 2008 "
pack .one.right.caption -side left
pack .one.right -side right
label .oneB.caption -text " A program for determining stellar abundances" -font {times 16 {italic}}
pack .oneB.caption -side left
label .two.select -text Select: -width 12
pack .two.select -side left
button .two.model -text model -command Getmodel
pack .two.model -side left
button .two.atom -text atom -command Getatom
pack .two.atom -side left
button .two.eqw -text "EqW file" -command Geteqw
pack .two.eqw -side left
label .three.modellabel -text "Model   " -width 12 -justify left
label .three.model -justify left -textvariable model -relief sunken -width 40
pack .three.modellabel .three.model -side left
label .four.atomlabel -text "Atom    " -width 12 -justify left
label .four.atom -justify left -textvariable atom -relief sunken -width 40
pack .four.atomlabel .four.atom -side left
label .five.eqwlabel -text "EqW File" -width 12 -justify left
label .five.eqw -justify left -textvariable eqw -relief sunken -width 40
pack .five.eqwlabel .five.eqw -side left
label .six.switchlabel -text "Switch: " -relief groove
checkbutton .six.0 -text "Atlas9" -variable flag9
pack .six.switchlabel .six.0 -side left
label .seven.label -text "Microturbulent Velocity:" -relief groove
pack .seven.label -side left
label .eight.startlabel -text "Begin " -justify left
entry .eight.startentry -textvariable vlow -width 9 -relief sunken -background white
label .eight.sA -text "km/s " -justify left
label .eight.endlabel -text "End " -justify left
entry .eight.endentry -textvariable vhigh -width 9 -relief sunken -background white
label .eight.eA -text "km/s " -justify left
label .eight.spacelabel -text "Step " -justify left
entry .eight.spaceentry -textvariable vstep -width 6 -relief sunken -background white
label .eight.spA -text "km/s " -justify left
pack .eight.startlabel .eight.startentry .eight.sA .eight.endlabel .eight.endentry .eight.eA .eight.spacelabel .eight.spaceentry .eight.spA -side left
label .nine.switchlabel -text "Output: " -relief groove
pack .nine.switchlabel -side left
radiobutton .nineB.0 -variable flagm -text "Screen Output (Gnuplot)" -value 0
radiobutton .nineB.1 -variable flagm -text "Postscript" -value 1 -command Getps
pack .nineB.0 .nineB.1 -side left
label .nineC.label -text "Postscript Output:  " -justify left
label .nineC.out -textvariable outps -width 40 -relief sunken
pack .nineC.label .nineC.out -side left
button .eleven.run -text "Execute" -command Run
pack .eleven.run -side left -ipady 10 -padx 20
button .eleven.quit -text "Quit" -command Quit
pack .eleven.quit -side left -ipady 10 -ipadx 10 -padx 20
label .twelve.label -text "Status/Progress:   " -justify left
label .twelve.status -justify left -textvariable status -relief sunken -width 25
pack .twelve.label .twelve.status -side left
pack .one .oneB .two .three .four .five .six .seven .eight .nine .nineB .nineC .eleven .twelve -side top -fill x
set model "Not Selected         "
set atom  "Not Selected         "
set eqw "Not Selected         "
set vlow 0.5
set vhigh 3.0
set vstep 0.1
set flagm 0
set flag9 0
set modelflag 0
set atomflag 0
set eqwflag 0
set pwd [pwd]
set k -1
set pids 100
set status "Expecting Input"
set outps "out.ps"

proc Pstat {pids} {
  global k
  set n 2
  while {1} {
   after 500
   catch {set n [exec ps hp $pids > tmp.wc]}
   set m [exec wc -l tmp.wc]
   scan $m {%d %s} k tmpstr
   if {$k == 0} {
    return
   }
   update
  }
}

proc Getmodel {} {
   global model modelflag pwd
   set types {
       {{mod files} {.mod} }
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

proc Geteqw {} {
   global eqw eqwflag pwd
   set types {
       {{dat files} {.eqw}}
       {{All Files} {*}}
}
   set eqw [tk_getOpenFile -filetypes $types -title "Select Equivalent Width File" -initialdir $pwd]
   set .five.eqw.text $eqw
   if {$eqw != ""} {
     set eqwflag 1
   } else {
     set eqwflag 0
   }
}

proc Getps {} {
   global outps pwd
   set types {
       {{dat files} {.ps}}
       {{All Files} {*}}
}
   set outps [tk_getSaveFile -filetypes $types -title "Select Postscript Output File" -initialdir $pwd]
   set .nineC.out.text $outps
}

proc Run {} {
   global k model vlow vhigh vstep atom pids eqw flagm status modelflag
   global atomflag eqwflag outps flag9
   if {$modelflag == 0 || $atomflag == 0 || $eqwflag == 0} {
      tk_messageBox -type ok -message "Input not complete!" -icon question
      return
   }
   if {$k > 0} {
    tk_messageBox -type ok -message "Blackwell already running!" -icon question
    return
   }
   set feqw [open $eqw r]
   set gnplt [open gnuplot.rsp w]
   puts $gnplt "set xlabel 'Microturbulent Velocity (km/s)'"
   puts $gnplt "set ylabel 'Abundance'"
   set flag 0
   set status "Executing ..."
   while {[gets $feqw line] >= 0} {
      set wave [string range $line 0 [string first "." $line]]
      set status $wave
      set wavein $wave
      set waveout $wave
      append wavein inp
      append waveout dat
      set f2 [open $wavein w]
      puts $f2 $line
      close $f2
      set rspID [open tmp.rsp w]
      puts $rspID $model
      puts $rspID $wavein
      puts $rspID $waveout
      puts $rspID $atom
      puts $rspID $vlow,$vhigh,$vstep
      close $rspID
      if {$flag == 0} {
        set gstr "plot '"
      } else {
        set gstr "replot '" 
      }
      set flag 1
      append gstr $waveout
      append gstr "' using 3:5 with lines"
      puts $gnplt $gstr
      if {$flag9 == 1} {
       set pids [exec blackwel t < tmp.rsp > out.out &]
      } else { 
       set pids [exec blackwel < tmp.rsp > out.out &]
      }
      Pstat $pids
   }
   close $feqw
   if {$flagm == 1} {
     puts $gnplt "set terminal postscript color"
     puts $gnplt "replot"
   } else {
      puts $gnplt "pause mouse"
    }
   close $gnplt
    if {$flagm == 0} {
      catch {exec gnuplot gnuplot.rsp &}
    } else {
      catch {exec gnuplot gnuplot.rsp > $outps &}
    }
    set status "Finished"
}

proc Quit {} {
   global k pids
    if {$k > 0} {
	set choice [tk_messageBox -type yesno \
	    -message "Blackwell is running! Quit?" -icon question]
        switch -- $choice {
	    yes exit
            no return
        }
    }
   exit
}

