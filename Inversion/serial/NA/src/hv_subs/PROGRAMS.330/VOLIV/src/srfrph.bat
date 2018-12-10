#!/bin/sh

case $TERM in
vt100|vt100n|201) clear
	mgotek 
	srfphr96 
	plot4014 < SRFPHR96.PLT
	rm SRFPHR96.PLT
	sleep 10
	mrttek ;;
4014|tek)
	clear
	srfphr96 
	plot4014 < SRFPHR96.PLT
	rm SRFPHR96.PLT
	sleep 10;;
xterm|xterm-color|sun-cmd)
	srfphr96  
	plotxvig <  SRFPHR96.PLT
	rm SRFPHR96.PLT
	;;
*) echo 'TERMINAL UNKNOWN USE TEKTRONIX' ;;
esac
