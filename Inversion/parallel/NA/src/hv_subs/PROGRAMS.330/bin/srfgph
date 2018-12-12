#!/bin/sh
case $TERM in
vt100|vt100n) clear
	mgotek 
	srfphv96 -G
	plot4014 < SRFPHG96.PLT
	rm SRFPHG96.PLT
	sleep 10
	mrttek ;;
4014|tek)
	clear
	srfphv96 -G
	plot4014 < SRFPHG96.PLT
	rm SRFPHG96.PLT
	sleep 10;;
xterm|xterm-color|sun-cmd)
	srfphv96 -G
	plotxvig < SRFPHG96.PLT
	rm SRFPHG96.PLT
	;;
*) echo 'TERMINAL UNKNOWN USE TEKTRONIX' ;;
esac
