===============================================================================
$Id: README,v 1.1 2008/02/05 17:17:32 jkleinj Exp $
README: Basic instructions for Genetic Algorithm
(C) Jens Kleinjung 2003-2007
(C) Alessandro Pandini 2006-2007
===============================================================================

GA - A Genetic Algorithm Program
================================

Installation / Compilation
--------------------------
    % tar -zxvf <program>.tgz	# extract program
    % cd <program>
    % make all					# compile program
	(the message 'missing .d files' at first compilation is normal and can be ignored)

	Compilation adjustments:
	Edit the 'Makefile' to choose between '-m32' (32bit) or '-m64' (64 bit) and
		between development binary or production binary.

	In case soemething goes wrong:
	% make clean				# remove program and binaries


Syntax
------
	Depends on the application. For help type:
	% <program> --help

	The 'help' output lists a short syntax description and command line options. 
	

Info
----
	For any further information see the MANUAL.

===============================================================================

