#!/bin/sh
case $TERM in
vt100|vt100n) clear
	mgotek 
	srfphv96 -V 
	plot4014 < SRFPHV96.PLT
	rm SRFPHV96.PLT
	sleep 10
	mrttek ;;
4014|tek)
	clear
	srfphv96 -V
	plot4014 < SRFPHV96.PLT
	rm SRFPHV96.PLT
	sleep 10;;
xterm|xterm-color|sun-cmd)
	srfphv96 -V
	plotxvig < SRFPHV96.PLT
	rm SRFPHV96.PLT
	;;
*) echo 'TERMINAL UNKNOWN USE TEKTRONIX' ;;
esac
