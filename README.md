
# GA - A Genetic Algorithm Program
This GA progam is a framework to run a Genetic Algorithm on a pluggged-in
application. The application is called 'application' in this manual, but
the term is a placeholder for a supposedly more specific name.

The GA program code is contained in the files 'ga.c and 'ga.h'. It provides 
two main data structures, the gene pool 'Pool' and the GA parameter collection 
'Gapar'. The gene pool holds the values to be optimised (the genes), while the 
parameters control the behaviour of the optimisation algorithm. 
Default parameter values are held in the file 'gapar.h' and they can be
overwritten by specifying command line parameters.
The plugged-in application holds its own parameters in 'application.h' and 
their default values in 'applicationpar.h', while the actual application code 
should be in 'application.c'.


### Installation / Compilation

## Usage
Execution of 'minset --help' prints the usage instruction and options
(command line parameters).
All default values can be defined in the 'gapar.h' and 'application.h' files so 
that the program can be run without any command line parameters.
If specified files cannot be opened, the program exits with the message
```
safe_open: Assertion 'file != 0' failed'.
```

### Mutation
The gene pool is an array of 'genome's, where each genome is 
composed of the same combination and 'GENENUM' number of genes, but with
different gene values. Each gene can adopt a value in the range 
[LOWLIM, UPLIM]. Values are either randomly assigned ('RANDOM') or randomly
modified ('SEEDED') by the agorithm. This step is analog to mutations in 
natural evolution. 

### Selection
After assigning/modifying the gene values of a genome, 
its fitness is tested and a fitness score is attributed to it. Of all the
'POPSIZE' genomes in the pool, only the top 'FITMATE' genomes are selected
to parent the next generation. This step is analog to selection in natural 
evolution. It is important to set the correct sorting order, i.e. to either 
'MAXIMIZE' or 'MINIMIZE' the fitness score.

### Breeding
Breeding is performed using either 'CROSSOVER' between random pairs of parents 
or by the 'EQUILIBRIUM' method, which is not explained here. The GA algorithm
runs a number of 'GENERATION' optimisation cycles. 

### Variation
Additionally, the input dataset can be split into 'JACKKNIFE' subsets that 
are optimised independently and/or the entire GA can be repeated 'REPEAT' times. 
Jackknifing can be used to estimate paramater variability and repeats indicate 
the robustness of the approach. The random generator is seeded at program start 
with the current time, so consecutive runs are mutually independent.
Output files are consecutively numbered by 'REAPEAT' and 'JACKKNIFE' counters.

### Convergence
Convergence is defined as all genomes being identical. This is a strict 
requirement, and other convergence criteria could be implemented.

### Completeness

There are of course many more possible GA features that are not implemented 
here. However, the lean and modular implementation renders it easily adaptable
and expandable.

### Parameters
GA parameters are defined in 'gapar.h', but they can be overwritten with command
line parameters. When additional parameters for the application are required, 
then these should be defined in 'applicationpar.h' and the cognate command line 
parameters should be joined with those of the GA in a file 'parse\_args.c'. 
In that case disable the GA function 'parse\_args\_ga'. 

### Application
Applications are plugged into the GA through three functions:
'initialise\_application()', 'run\_application()' and 'finalise\_application().
These should be defined in 'application.c' and 'application.h'.
The 'run\_application()' function should return the fitness score of the
current genome.

### How to add a new application
As mentioned above, the application should be defined in three files:
'application.c' for the application code, 'application.h' for the data 
definition and 'applicationpar.h' for the default values. 
The only files of the GA that need modification are 'ga.c' and 'parse\_args.c',
the latter to accomodate the parameters of the application.
In 'ga.c', sections for user modification are indicated by >>>>> delimiters:
1. Plug in the application routines 'initialise\_application()',
	'run_application()' and 'finalise_application()', which are to be
	defined in separate files 'application.c' and 'application.h'.
	It is recommended to define a data structure that holds all
	variables (and functions) of the application.
2. Add 'application.o' and other source binaries to the 'Makefile', adjust 
	libraries and 'TARGET' name; compile.


## References
Users publishing results obtained with the GA program and its applications
should acknowledge its use by one of the following citations:

- J. Kleinjung, J. Romein, K. Lin and J. Heringa
Contact-based sequence alignment.
Nucleic Acids Research 32(8) (2004) 2464-2473.

- A. Pandini, L. Bonati, F. Fraternali and J. Kleinjung
MinSet: A general approach to derive maximally representative database
subsets by using fragment dictionaries and its application to the SCOP
database.
Bioinformatics 23 (2007) 515-516.

- A. Pandini and J. Kleinjung
MinSet: automatic derivation of maximally representative subsets of
structural domains.
Nucleic Acids Res., submitted.


## Contact
jens@jkleinj.eu


## Copyright and Authors
Copyright (C) 2003-2017 Jens Kleinjung
Copyright (C) 2006-2007 Alessandro Pandini


License
-------
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

