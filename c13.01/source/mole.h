/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MOLE_H_
#define MOLE_H_

/* mole.h */

#include "count_ptr.h"
#include "elementnames.h"
#include "transition.h"

#define SMALLABUND 1e-24

enum mole_state {MOLE_NULL, MOLE_PASSIVE, MOLE_ACTIVE};

class chem_atom;

class chem_element {
	explicit chem_element(); // Do not implement
	chem_element &operator=(const chem_element&); // Do not implement
public:
   explicit chem_element(int Z, const char*label) : Z(Z), label(label)
	{}
	~chem_element() throw()
	{}
	const int Z;
	const string label;
	map<int, count_ptr<chem_atom> > isotopes;
	//(first -> Atomic A; second -> chem_atom )
	//(first -> -1 for bogus isotope, i.e. where no
	//	isotopes have been explicitly defined)
};

typedef map<int, count_ptr<chem_atom> >::iterator isotopes_i;

class chem_atom {
public:
	// Link back to basic element for convenience -- not a count_ptr
	// as this would lead to a reference cycle (and hence the
	// destructors never being called).  Many-to-one relation suggests
	// that the weak link should be this way around.
	chem_element* el;
	int A;			/* mass number */
	vector<int> ipMl;   	/* Atom and ion species in molecule arrays */
	realnum mass_amu;	/* mass of isotope in AMU */
	double frac;		/* fraction of element in this isotope */

	bool lgMeanAbundance( void ) const
	{
		return ( A < 0 );
	}

	/* Chemical symbols for elements */
	string label( void ) const
	{
		if( lgMeanAbundance() )
			return el->label;
		else if( el->Z==1 && A==2 )
		{
			// Deuterium is a special case
			return "D\0";
		}
		else
		{
			char str[4];
			sprintf(str,"^%d",A);
			return ( str + el->label );
		}
	}
	
	int compare ( const chem_atom &b ) const
	{
		// sort by proton number first
		if ( el->Z < b.el->Z )
			return -1;
		else if ( el->Z > b.el->Z )
			return 1;

		if (mass_amu < b.mass_amu)
			return -1;
		else if (mass_amu > b.mass_amu)
			return 1;
		else if (A < b.A )
			return -1;
		else 
			return 0;
	}
};
inline bool operator< (const chem_atom &a, const chem_atom &b)
{
	return a.compare(b) < 0;
}
inline bool operator> (const chem_atom &a, const chem_atom &b)
{
	return a.compare(b) > 0;
}
inline bool operator<= (const chem_atom &a, const chem_atom &b)
{
	return a.compare(b) <= 0;
}
inline bool operator>= (const chem_atom &a, const chem_atom &b)
{
	return a.compare(b) >= 0;
}
inline bool operator== (const chem_atom &a, const chem_atom &b)
{
	return a.compare(b) == 0;
}
inline bool operator!= (const chem_atom &a, const chem_atom &b)
{
	return !(a == b);
}

typedef vector< count_ptr<chem_atom> > ChemAtomList;
extern ChemAtomList atom_list;
extern ChemAtomList unresolved_atom_list;
extern chem_element *null_element;
extern chem_atom *null_atom;

class element_pointer_value_less
{
public:
	bool operator()(const count_ptr<chem_atom>& a,
						const count_ptr<chem_atom>& b) const
	{
		return *a < *b;
	}
};

/* Structure containing molecule data, initially only CO */
class molecule {
public:
	typedef map<const count_ptr<chem_atom>, int,
		element_pointer_value_less> nAtomsMap;

	string parentLabel;
	int parentIndex;
	bool isEnabled;   /* Is it enabled? */

	/* Species physical data */
	string    label;        /** name */
	nAtomsMap nAtom;  	/** number of each element in molecule */
	int     charge;         /** Charge on species/number of e- liberated by formation */
	bool    lgExcit;        /** Is species excited (e.g. H2*) */
	bool    lgGas_Phase;    /** Solid or gas phase? */
	int     n_nuclei(void) const       /** total number of nuclei */
	{
		int num = 0;
		for (nAtomsMap::const_iterator el = nAtom.begin(); 
			  el != nAtom.end(); ++el)
		{
			num += el->second;
		}
		return num;
	}
	bool isMonatomic(void) const       /** total number of nuclei */
	{
		if (nAtom.size() == 1 && nAtom.begin()->second == 1)
			return true;
		return false;
	}

	realnum form_enthalpy;  /** formation enthalpy for the molecule (at 0K), in units of KJ/mol */
	realnum mole_mass;      /** Mass of molecule */

	/* Parameters as computational object */
	enum mole_state state;
	int index, groupnum;

	chem_atom *heavyAtom(void) //const
	{
		for( nAtomsMap::reverse_iterator it=nAtom.rbegin(); it!=nAtom.rend(); ++it )
		{
			if (0 != it->second )
			{
				return it->first.get_ptr();
			}
		}
		return null_atom;
	}

	int compare(const molecule &mol2) const
	{
		nAtomsMap::const_reverse_iterator it1, it2;

		for( it1 = nAtom.rbegin(), it2 = mol2.nAtom.rbegin();
			it1 != nAtom.rend() && it2 != mol2.nAtom.rend(); ++it1, ++it2 )
		{
			if( *(it1->first) > *(it2->first) )
				return 1;
			else if( *(it1->first) < *(it2->first) )
				return -1;
			else if( it1->second > it2->second)
				return 1;
			else if( it1->second < it2->second)
				return -1;
		}

		if( it1 != nAtom.rend() && it2 == mol2.nAtom.rend() )
			return 1;
		else if( it1 == nAtom.rend() && it2 != mol2.nAtom.rend() )
			return -1;
		else
			ASSERT( it1 == nAtom.rend() && it2 == mol2.nAtom.rend() );

		// sort by label if falls through to here       
		return ( label.compare(mol2.label) );
			
	}
}; 

/* iterators on nAtom */	
typedef molecule::nAtomsMap::iterator nAtoms_i;
typedef molecule::nAtomsMap::reverse_iterator nAtoms_ri;
typedef molecule::nAtomsMap::const_reverse_iterator nAtoms_cri;

/**mole_drive main driver for chemical equilibrium routines */
extern void mole_drive(void);

/**mole_create_react build reaction structures */
extern void mole_create_react(void);

class mole_reaction;

mole_reaction *mole_findrate_s(const char buf[]);

extern void mole_print_species_reactions( molecule *speciesToPrint );

extern molecule *null_mole;

extern molecule *findspecies(const char buf[]);

/** \verbatim >>chng 03 feb 09, rm ipH3P_hev, since not used, and decrement NUM_HEAVY_MOLEC to 17 
 >>chng 03 aug 04, rm ipCTWO and ipC2P from den since not included in balance,
 and always finds zero column density, so NUM_HEAVY_MOLEC from 17 to 15 
 >>chng 03 aug 05, rm ch2 and ch3, so n from 15 to 13 
 >>chng 03 nov 14  add Si chemistry & CH3+, so that now every
     reaction that is in the TH85 chemical network is also included
     in Cloudy.  Additionally, there are also reactions taken from other
     papers (mostly Hollenbach and McKee...see co.c).  In all 20 molecular
     species are calculated, along with the atomic and first ionization 
	 stages of C, O, and Si
 >>chng 04 May 13, Nick Abel.  Add CH3, CH4, CH4+, and CH5+ to network in order 
	to get the same chemical abundances vs. depth as other PDR codes in the Leiden
	meeting.  With changes we now can predict molecular abundances for 24 C, O, 
	and Si bearing molecules. 

 >>chng 04 jul 13, Nick Abel.  Add nitrogen and sulphur bearing molecules
    to the chemical network.  First added to generate a chemical model for
    eta carinae, but is applicable to all molecular clouds 

 >>chng 05 mar 11, Nick Abel.  Add C2 and C2+ to chemistry, reactions 
    involving these species affects the abundance of C

 >>chng 05 mar 23, Nick Abel.  Add Chlorine to chemistry 
 \endverbatim */

/** this includes the atomic and first ionized species
   of each element that can combine to form molecules.  
   This is the number of molecules, ions, and atoms that the co network uses
   This is used in comole
   to improve the calculation, as deep in molecular regions reactions with molecules
   can be important to the ionization balance */

class t_mole_global {

public:
	void init(void);

	void make_species(void);
		
  /**mole_zero allocate + initialize workspace */
	void zero(void);

	/** flag to turn off all molecules, set with no molecules command */
	bool lgNoMole;

	/** flag to turn off heavy molecules, set with no heavy molecules command */
	bool lgNoHeavyMole;

	/** flag set true if H2O destruction rate went to zero */
	bool lgH2Ozer;

	/** set rates to that in UMIST */
	bool lgLeidenHack;

	bool lgFederman;
	bool lgStancil;

	/** option to use effective temperature as defined in
	 * >> refer Federman, S. R. & Zsargo, J. 2003, ApJ, 589, 319
	 * By default, this is false - changed with set chemistry command */
	bool lgNonEquilChem;

	/** option to set proton elimination rates to zero
	 * >>refer	CO	chemistry	Huntress, W. T., 1977, ApJS, 33, 495
	 * By default, this is false - changed with set chemistry command */
	bool lgProtElim;

	/** option to not include neutrals in the non-equilibrium scheme
	 * >> refer Federman, S. R. & Zsargo, J. 2003, ApJ, 589, 319
	 * By default, this is false - changed with set chemistry command */
	bool lgNeutrals;

	 /** do we include capture of molecules onto grain surfaces?  default is true,
	  * turned off with NO GRAIN MOLECULES command */
	bool lgGrain_mole_deplete;

	// flag saying whether to model isotopes (and isotopologues) of a given element
	vector<bool> lgTreatIsotopes;

	/** flag saying whether an element is in the chemistry network */
	int num_total, num_calc, num_compacted;

	typedef vector<count_ptr<molecule> > MoleculeList;
	MoleculeList list;

	static void sort(MoleculeList::iterator start,
				 MoleculeList::iterator end);
	static void sort(molecule **start, molecule **end);
};

extern t_mole_global mole_global;

class t_mole_local
{
public:
	void set_location( long nelem, long ion, double *dense );
	void set_isotope_abundances( void );
	double sink_rate_tot(const char chSpecies[]) const;
	double sink_rate_tot(const molecule* const sp) const;
	double sink_rate(const molecule* const sp, const mole_reaction& rate) const;
	double sink_rate(const molecule* const sp, const char buf[]) const;
	double source_rate_tot(const char chSpecies[]) const;
	double source_rate_tot(const molecule* const sp) const;
	/** returns the photodissociation rate per unit volume [cm^-3 s^-1] producing monatomic species chSpecies.
	 * *Excludes* photoionizations of other monatomic species, e.g. N-,PHOTON=>N,e- */
	double dissoc_rate(const char chSpecies[]) const;
	double chem_heat(void) const;
	double findrk(const char buf[]) const;
	double findrate(const char buf[]) const;

	double grain_area, grain_density, grain_saturation;

	/** total charge in molecules */
	double elec;

	/** these are source and sink terms for the ionization ladder, for chemical
	 * processes that remove and add species */
	double **source , **sink;
	
	realnum ***xMoleChTrRate;/***[LIMELM][LIMELM+1][LIMELM+1];*/

	valarray<class molezone> species;

	vector<double> reaction_rks;
	vector<double> old_reaction_rks;
	long old_zone;
};

extern t_mole_local mole;

class molezone {
public:
	molezone()
	{
		init();
	}
	void init (void)
	{
		location = NULL;
		levels = NULL;
		lines = NULL;
		zero();
	}
	void zero (void)
	{
		src = 0.;
		snk = 0.;
		den = 0.;
		column = 0.;
		nAtomLim = -1;
		xFracLim = 0.;
		column_old = 0.;
	}
	double *location;      /** Location of density in non-molecule code, NULL if none exists */

	/** rate s-1 for molecular charge transfer, nelem from to */
	double src, snk;

	qList* levels;
	TransitionList* lines;

	/* Current zone data */
	double  den;       /** density (cm-3) */
	realnum column;    /** total column density in this iteration */
	int     nAtomLim;  /** atomic number MINUS ONE of element for which is closest to limiting molecule density */
	realnum xFracLim;  /** The fraction of that element in this species */

	/* Historical solution data */
	realnum column_old;   /** total column density in previous iteration */
};

extern molezone *null_molezone;

extern molezone *findspecieslocal(const char buf[]);

extern void mole_punch(FILE *punit, const char speciesname[], const char args[], bool lgHeader, bool lgData, double depth);

extern void total_molecule_elems(realnum total[LIMELM]);
extern void total_molecule_deut(realnum &total);

extern realnum total_molecules(void);

extern realnum total_molecules_gasphase(void);

extern void mole_make_list(void);
extern void mole_make_groups(void);

void mole_cmp_num_in_out_reactions(void);

bool lgDifferByExcitation( const molecule &mol1, const molecule &mol2 );

extern void mole_update_species_cache(void);

void mole_update_sources(void);

void mole_rk_bigchange(void);

void create_isotopologues_one(
	ChemAtomList& atoms,
	vector< int >& numAtoms,
	string atom_old,
	string atom_new,
	string embellishments,
	vector<string>& newLabels );

bool parse_species_label( const char label[], ChemAtomList &atomsLeftToRight, vector<int> &numAtoms, string &embellishments );
bool parse_species_label( const char mylab[], ChemAtomList &atomsLeftToRight, vector<int> &numAtoms, string &embellishments,
	bool &lgExcit, int &charge, bool &lgGas_Phase );

#endif /* MOLE_H_ */

