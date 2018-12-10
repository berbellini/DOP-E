#!/bin/sh
case $TERM in
vt100|vt100n) clear
	mgotek 
	rftnpv96 -
	plot4014 < RFTNPV96.PLT
	rm RFTNPV96.PLT
	sleep 10
	mrttek ;;
4014|tek)
	clear
	rftnpv96 
	plot4014 < RFTNPV96.PLT
	rm RFTNPV96.PLT
	sleep 10;;
xterm|dtterm|sun-cmd)
	rftnpv96 
	plotxvig < RFTNPV96.PLT
	rm RFTNPV96.PLT
	;;
*) echo 'TERMINAL UNKNOWN USE TEKTRONIX' ;;
esac
