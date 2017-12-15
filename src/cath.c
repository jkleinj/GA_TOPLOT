    FILE *cathlist;
    char *cathlistname = "cath.pairs.test.list";
    char *cath = "/mb/databases/cath/v2.4";
						sprintf(pdb0filename, "%s/%s",
							cath, fam.family[f], fam.family[f]);
						++ f;
						sprintf(pdb1filename, "%s/%s",
							cath, fam.family[f]);

#ifdef CATH
                        sprintf(pdb0filename, "%s/%s",
                            cath, fam.family[f], fam.family[f]);
                        ++ f;
                        sprintf(pdb1filename, "%s/%s",
                            cath, fam.family[f]);
#endif

