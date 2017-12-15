FILE *homstradlist;
char *homstradlistname = "homstrad.dist.list";
char *homstrad = "/home/jkleinj/data/homstrad_1104/";
						sprintf(pdb0filename, "%s%s/%s.0.pdb",
							homstrad, fam.family[f], fam.family[f]);
						sprintf(pdb1filename, "%s%s/%s.1.pdb",
							homstrad, fam.family[f], fam.family[f]);
int nfam = 0;


#ifdef ALIGNMODE
    char pdb0filename[200] = ""; /* RUN input file */
    char pdb1filename[200] = ""; /* RUN input file */
    char s[50];
    Parameter parameter[parnum]; /* RUN parameters */
    genenum = parnum;

                    fprintf(stderr, "repeat %d/%d, db-frac %d/%d, generation %d/%d, genome %d/%d\n",
                        z, repga, x, jackknife, l, loop, ix, popsize);
                    fflush(stderr);
                    fprintf(gaoutfile, "repeat %d/%d, db-frac %d/%d, generation %d/%d, genome %d/%d\n",
                        z, repga, x, jackknife, l, loop, ix, popsize);
                    fflush(gaoutfile);
                    set_parameters(&pool[0], ix, &parameter[0], gaoutfile);
#endif

#ifdef HOMSTRAD
                        sprintf(pdb0filename, "%s%s/%s.0.pdb",
                            homstrad, fam.family[f], fam.family[f]);
                        sprintf(pdb1filename, "%s%s/%s.1.pdb",
                            homstrad, fam.family[f], fam.family[f]);
#endif

#ifdef ALIGNMODE
                            /*pool[ix].fitness += aliss(pdb0filename, pdb1filename, &parameter[0]);*/
                            /*pool[ix].fitness += 0;*/
                    fprintf(gaoutfile, "fitness=%f\n", pool[ix].fitness);
                    fflush(gaoutfile);
                print_pool_ascii(&pool[0], pooloutfile); /* print out pool */
    for (i = 0; i < fam.nfam; ++ i)
        free(fam.family[i]);
#endif

/*____________________________________________________________________________*/
/* transform gene values to RUN parameters */
void set_parameters(Pool *pool, int ix, Parameter *parameter, FILE *gaoutfile)
{
    unsigned int i = 0;

    /* parameter[0] -> gap open */
    parameter[i].value = (float)((0.5 * pool[ix].genome[i]) + 8.);
    strcpy(parameter[i].name, "go");
    /*---*/
    ++ i; /* parameter[1] -> gap extend */
    parameter[i].value = (float) pool[ix].genome[1] / 2.;
    strcpy(parameter[i].name, "ge");
    /*---*/
    ++ i; /* parameter[2] -> local fragment match score factor */
    parameter[i].value = (float)(0.25 * pool[ix].genome[i] + 1.0);
    strcpy(parameter[i].name, "locfrag");
    /*---*/
    ++ i; /* parameter[3] -> atom or segment match score factor */
    parameter[i].value = (float)(0.2 * pool[ix].genome[i]);
    strcpy(parameter[i].name, "seg*");
    /*---*/
    ++ i; /* parameter[4] -> see above */
    parameter[i].value = (float)((0.2 * pool[ix].genome[i]) - 1);
    strcpy(parameter[i].name, "seg+");
    /*---*/
    ++ i; /* parameter[5] -> exposure score factor */
    parameter[i].value = (float)((0.25 * pool[ix].genome[i]) + 0.1);
    strcpy(parameter[i].name, "exp");
    /*---*/
    ++ i; /* parameter[6] -> max RMSD for atom pair selection */
    parameter[i].value = (float)((0.5 * pool[ix].genome[i]) + 1);
    strcpy(parameter[i].name, "maxrmsd");

    fprintf(stdout, "\t%d ", ix);
    fprintf(gaoutfile, "\t%d ", ix);
    for (i = 0; i < parnum; ++ i)
    {
	fprintf(stdout, "%s:%d->%5.3f ", parameter[i].name, pool[ix].genome[i], parameter[i].value);
	fprintf(gaoutfile, "%s:%d->%5.3f ", parameter[i].name, pool[ix].genome[i], parameter[i].value);
    }
    /*
    fprintf(stdout, "\n");
    fprintf(gaoutfile, "\n");
    */
    fflush(stdout);
    fflush(gaoutfile);
}

/*____________________________________________________________________________*/
interface()
{
	char listfilename[200];
					/* transform genes to RUN parameters */
#ifdef ALIGNMODE
					set_parameters(&pool[0], ix, &parameter[0], pooloutfile);
#endif
#ifdef ALIGNMODE
							/*pool[ix].fitness += aliss(pdb0filename, pdb1filename, &parameter[0]);*/
							/*pool[ix].fitness += 0;*/
#endif
#ifdef ALIGNMODE
				print_pool_ascii(&pool[0], pooloutfile); /* print out pool */ 
#endif
#ifdef HOMSTRAD
                        sprintf(pdb0filename, "%s%s/%s.0.pdb",
                            homstrad, fam.family[f], fam.family[f]);
                        sprintf(pdb1filename, "%s%s/%s.1.pdb",
                            homstrad, fam.family[f], fam.family[f]);
#endif
#ifdef ALIGNMODE
    for (i = 0; i < fam.nfam; ++ i)
		free(fam.family[i]);
#endif

