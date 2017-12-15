/*==============================================================================
gapar.h : GA default parameters
(C) 2006-2017 Jens Kleinjung
(C) 2006-2007 Alessandro Pandini

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
==============================================================================*/

#if !defined(GAPAR_H)
#define GAPAR_H

/* GA scopes */
#define POPSIZE 2500 /* population size (number of genomes) */
#define FITMATE 250 /* fittest genomes : parents of next generation */
#define GENENUM 2 /* number of genes (optimisation parameters) */
#define GENERATION 100 /* number of generations (optimisation loops) */
#define LOWLIM 0 /* lower limit of parameter value range */
#define UPLIM 1 /* upper limit of parameter value range */

/* fitness sorting */
#define MINIMIZE 0 /* lowest score values are fittest */
#define MAXIMIZE 1 /* highest score values are fittest */

/* initital pool (single choice mandatory) */
#define RANDOM 1 /* random starting population */
#define SEEDED 0 /* seeded starting population */
#define MAXVAR 4 /* maximal parameter variation (mutation) in seeded population */

/* breedmode (single choice mandatory) */
#define CROSSOVER 1 /* crossover breeding */
#define EQUILIBRIUM 0 /* equilibrium breeding */

/* split input set and/or repeat GA (value 1 means no split/repeat) */
#define JACKKNIFE 1 /* split into 'JACKKNIFE' parts */
#define REPEAT 1 /* repeat GA 'REPEAT' times */

#endif

