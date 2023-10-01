/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CO_step fills in matrix for heavy elements molecular routines */
#include "cddefines.h"
#include "mole.h"
#include "mole_priv.h"
#include "save.h"
#include "dense.h"
#include "atmdat.h"
/* Nick Abel between July and October of 2003 assisted Dr. Ferland in improving the heavy element 
 * molecular network in Cloudy. Before this routine would predict negative abundances if 
 * the fraction of carbon in the form of molecules came close to 100%. A reorganizing of 
 * the reaction network detected several bugs.  Treatment of "coupled reactions",
 * in which both densities in the reaction rate were being predicted by Cloudy, were also 
 * added.  Due to these improvements, Cloudy can now perform calculations
 * where 100% of the carbon is in the form of CO without predicting negative abundances
 *
 * Additional changes were made in November of 2003 so that our reaction 
 * network would include all reactions from the TH85 paper.  This involved 
 * adding silicon to the chemical network.  Also the reaction rates were
 * labeled to make identification with the reaction easier and the matrix 
 * elements of atomic C, O, and Si are now done in a loop, which makes 
 * the addition of future chemical species (like N or S) easy.
 * */
/* Robin Williams in August 2006 onwards reorganized the coding to cut down repetitions.  
 * This isolated several further bugs, and allows a sigificant number of lines of
 * code to be eliminated.  The balance of S2/S2+ amd ClO/ClO+ seems highly sensitive
 * (with small log scale results varying significantly if the order of arithmetic
 * operations is changed) -- I suspect this may imply a bug somewhere.
 * */
/*lint -e778 constant expression evaluatess to 0 in operation '-' */
/*=================================================================*/

STATIC bool lgNucleiConserved(const multi_arr<double,2> &c);

void mole_eval_balance(long int num_total, double *b, bool lgJac, multi_arr<double,2> &c)
{
	long int i, j;
	mole_reaction *rate;
	double rate_tot, rate_deriv[MAXREACTANTS], rk;
	molecule *sp;

	DEBUG_ENTRY("mole_eval_balance()");
	/* zero out array used for formation rates */
	for( i=0; i < num_total; i++ )
	{
		b[i] = 0.;
	}
	if (lgJac)
		c.zero();

	for(mole_reaction_i p=mole_priv::reactab.begin(); 
			p != mole_priv::reactab.end(); ++p) 
	{
		rate = &(*p->second);
		rk = mole.reaction_rks[ rate->index ];

		rate_tot = rk;
		for(i=0;i<rate->nreactants;i++)
		{
			rate_tot *= mole.species[ rate->reactants[i]->index ].den;
		}		
		
		for(i=0;i<rate->nreactants;i++)
		{	
			sp = rate->reactants[i];
			if (rate->rvector[i] == NULL)
			{
				b[sp->index] -= rate_tot;
			}
		}
		
		for(i=0;i<rate->nproducts;i++)
		{
			sp = rate->products[i];
			if (rate->pvector[i] == NULL)
			{
				b[sp->index] += rate_tot;
			}
		}
		
		if (lgJac)
		{
			for(i=0;i<rate->nreactants;i++)
			{
				rate_deriv[i] = rk;
				for(j=0;j<rate->nreactants;j++)
				{
					if(i!=j)
					{
						rate_deriv[i] *= mole.species[ rate->reactants[j]->index ].den;
					}
				}
			}
			for(j=0;j<rate->nreactants;j++)
			{
				sp = rate->reactants[j];
				const double rated = rate_deriv[j];
				for(i=0;i<rate->nreactants;i++)
				{
					if (rate->rvector[i] == NULL)
						c[sp->index][rate->reactants[i]->index] -= rated;
				}
				for(i=0;i<rate->nproducts;i++)
				{
					if (rate->pvector[i] == NULL)
						c[sp->index][rate->products[i]->index] += rated;
				}
			}
		}
	}

	if (lgJac)
	{
		ASSERT( lgNucleiConserved(c) );
	}

	//mole_dominant_rates(findspecies("H+"),ioQQQ);
	return;
}

void mole_eval_sources(long int num_total)
{
	long int i, j, nelem, ion, ion2;
	mole_reaction *rate;
	double rate_tot, rate_deriv[MAXREACTANTS], rk;
	molecule *sp;

	DEBUG_ENTRY("mole_eval_sources()");
	/* zero out array used for formation rates */
	for( i=0; i < num_total; i++ )
	{
		mole.species[i].src = mole.species[i].snk = 0.;
	}

	for( nelem=0; nelem< LIMELM; ++nelem )
	{
		/* these have one more ion than above */
		for( ion=0; ion<nelem+2; ++ion )
		{
			/* zero out the transfer array */
			for( ion2=0; ion2<nelem+2; ++ion2 )
			{
				mole.xMoleChTrRate[nelem][ion][ion2] = 0.;
			}
		}
	}

	for(mole_reaction_i p=mole_priv::reactab.begin(); 
			p != mole_priv::reactab.end(); ++p) 
	{
		rate = &(*p->second);
		rk = mole.reaction_rks[ rate->index ];
		
		for(i=0;i<rate->nreactants;i++)
		{
			rate_deriv[i] = rk;
			for(j=0;j<rate->nreactants;j++)
			{
				if(i!=j)
				{
					rate_deriv[i] *= mole.species[ rate->reactants[j]->index ].den;
				}
			}
		}
		
		rate_tot = rate_deriv[0] * mole.species[ rate->reactants[0]->index ].den;
		
		for(i=0;i<rate->nreactants;i++)
		{	
			sp = rate->reactants[i];
			if (rate->rvector[i] == NULL)
			{
				mole.species[sp->index].snk += rate_deriv[i];
			}
			else
			{
				// exclude charge exchange with isotopes of same element 
				if( rate->nreactants==2 )
				{
					long otherIndex = 1-i;
					molecule *sp2 = rate->reactants[otherIndex];
					if( sp->isMonatomic() && sp2->isMonatomic() && sp->nAtom.begin()->first->el == sp2->nAtom.begin()->first->el )
						continue;
				}

				if ( atmdat.lgCTOn )
				{
					for( ChemAtomList::iterator atom = unresolved_atom_list.begin(); atom != unresolved_atom_list.end(); ++atom)
					{
						nelem = (*atom)->el->Z-1;
						if (sp->nAtom.find(*atom) != sp->nAtom.end() && sp->nAtom[*atom] != 0 && rate->rvector[i]->charge != sp->charge)
						{
							mole.xMoleChTrRate[nelem][sp->charge][rate->rvector[i]->charge] +=
								(realnum) rate_deriv[i];
							break;
						}
					}
				}
			}
		}
		
		for(i=0;i<rate->nproducts;i++)
		{
			sp = rate->products[i];
			if (rate->pvector[i] == NULL)
			{
				mole.species[sp->index].src += rate_tot;
			}
		}
	}

	for (ChemAtomList::iterator atom = unresolved_atom_list.begin();
		  atom != unresolved_atom_list.end(); ++atom)
	{
		const long int nelem=(*atom)->el->Z-1;
		if( !dense.lgElmtOn[nelem] )
			continue;

		for (long int ion=0;ion<nelem+2;ion++) 
		{
			if ((*atom)->ipMl[ion] != -1)
			{
				mole.source[nelem][ion] = mole.species[(*atom)->ipMl[ion]].src;
				mole.sink[nelem][ion] = mole.species[(*atom)->ipMl[ion]].snk;
			}
			else
			{
				mole.source[nelem][ion] = 0.0;
				mole.sink[nelem][ion] = 0.0;
			}
		}
	}

	//mole_dominant_rates(findspecies("H+"),ioQQQ);
	return;
}

STATIC bool lgNucleiConserved(const multi_arr<double,2> &c)
{
	DEBUG_ENTRY ("lgNucleiConserved()");
	bool checkAllOK = true;
	multi_arr<double,2> test(atom_list.size(),mole_global.num_calc),
		tot(atom_list.size(),mole_global.num_calc);
	vector<double> ccache(mole_global.num_calc);
	vector<long> ncache(mole_global.num_calc);
	test.zero();
	tot.zero();

	map<chem_atom*, long> atom_to_index;
	for (unsigned long j=0; j<atom_list.size(); ++j )
	{
		atom_to_index[atom_list[j].get_ptr()] = j;
	}
	for (long j=0;j<mole_global.num_calc;j++) 
	{
		long nc = 0;
		for (long i=0;i<mole_global.num_calc;i++) 
		{
			if (c[i][j] != 0.)
			{
				ccache[nc] = c[i][j];
				ncache[nc] = i;
				++nc;
			}
		}
		if (nc > 0)
		{
			for (molecule::nAtomsMap::const_iterator el = mole_global.list[j]->nAtom.begin(); 
				  el != mole_global.list[j]->nAtom.end(); ++el)
			{
				long natom = atom_to_index[el->first.get_ptr()];
				const int nAtomj = el->second;
				for (long i=0;i<nc;i++) 
				{
					const double term = ccache[i] * nAtomj;
					test[natom][ncache[i]] += term;
					tot[natom][ncache[i]] += fabs(term);
				}
			}
		}
	}

	for( unsigned long natom=0; natom < atom_list.size(); ++natom)
	{
		for (long i=0;i<mole_global.num_calc;i++) 
		{
			const bool checkOK = 
				( fabs(test[natom][i]) <= MAX2(1e-10*tot[natom][i], 1e10*DBL_MIN) );
			if ( UNLIKELY(!checkOK) )
			{
				chem_atom *atom = atom_list[natom].get_ptr();
				fprintf(stdout,"Network conservation error %s %s %g %g %g %g\n",
						  atom->label().c_str(),
						  mole_global.list[i]->label.c_str(),
						  test[natom][i],
						  test[natom][i]/tot[natom][i],
						  mole.species[atom->ipMl[0]].den,
						  mole.species[atom->ipMl[1]].den);
				//fprintf(stdout,"Problem at %s\n",rate->label);
				checkAllOK = false;
			}
		}
	}
	return checkAllOK;
}

void mole_dominant_rates( const molecule *debug_species, FILE *ioOut )
{
	long int i, j;
	mole_reaction *rate, *ratesnk=NULL, *ratesrc=NULL;
	double rate_tot, rate_deriv[MAXREACTANTS], rk;
	molecule *sp;
	double snkx=0.,srcx=0.;

	for(mole_reaction_i p=mole_priv::reactab.begin();
			p != mole_priv::reactab.end(); ++p)
	{
		rate = &(*p->second);
		rk = mole.reaction_rks[ rate->index ];

		for(i=0;i<rate->nreactants;i++)
		{
			rate_deriv[i] = rk;
			for(j=0;j<rate->nreactants;j++)
			{
				if(i!=j)
				{
					rate_deriv[i] *= mole.species[ rate->reactants[j]->index ].den;
				}
			}
		}

		rate_tot = rate_deriv[0] * mole.species[ rate->reactants[0]->index ].den;

		if (debug_species != null_mole)
		{
			for(i=0;i<rate->nproducts;++i)
			{
				sp = rate->products[i];
				if (sp == debug_species && rate->pvector[i] == NULL)
				{
					if (fabs(rate_tot) > srcx)
					{
						srcx = rate_tot;
						ratesrc = rate;
					}
				}
			}
			for(i=0;i<rate->nreactants;++i)
			{
				sp = rate->reactants[i];
				if (sp == debug_species && rate->rvector[i] == NULL)
				{
					if (fabs(rate_deriv[i]) > snkx)
					{
						snkx = rate_deriv[i];
						ratesnk = rate;
					}
				}
			}
		}
	}


	if (debug_species != null_mole)
	{
		if (ratesrc)
		{
			fprintf( ioOut, "%20.20s src %13.7g of %13.7g [",
				ratesrc->label.c_str(),srcx,mole.species[debug_species->index].src);
			for (j=0;j<ratesrc->nreactants;j++)
			{
				if (j)
				{
					fprintf( ioOut, "," );
				}
				fprintf( ioOut, "%-6.6s %13.7g",
					ratesrc->reactants[j]->label.c_str(),
					mole.species[ ratesrc->reactants[j]->index ].den);
			}
			fprintf( ioOut, "]" );
		}
		if (ratesnk)
		{
			fprintf( ioOut, "%20.20s snk %13.7g of %13.7g [",
				ratesnk->label.c_str(), snkx * mole.species[debug_species->index].den,
				mole.species[debug_species->index].snk * mole.species[debug_species->index].den);
			for (j=0;j<ratesnk->nreactants;j++)
			{
				if (j)
				{
					fprintf( ioOut, "," );
				}
				fprintf( ioOut, "%-6.6s %13.7g",
					ratesnk->reactants[j]->label.c_str(),
					mole.species[ ratesnk->reactants[j]->index ].den);
			}
			fprintf( ioOut, "]" );
		}
	}
	fprintf( ioOut, "\n" );

	return;
}
