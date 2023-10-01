/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "taulines.h"

INSTANTIATE_MULTI_ARR( qList, lgBOUNDSCHECKVAL );

vector<TransitionList> AllTransitions;

bool lgStatesAdded;
bool lgLinesAdded;
multi_arr<qList,2> StatesElemNEW;
char **chSpecies;
species *dBaseSpecies;
vector<qList> dBaseStates;
vector< multi_arr<int,2> > ipdBaseTrans;
vector<TransitionList> dBaseTrans;
multi_arr<CollRateCoeffArray,2> AtmolCollRateCoeff;
CollSplinesArray ****AtmolCollSplines;
StoutColls ****StoutCollData;
long int nSpecies;
qList AnonStates(1);
TransitionList TauLines("TauLines", &AnonStates);
multi_arr<int,3> ipExtraLymanLines;
vector<vector<TransitionList> > ExtraLymanLines;
long int nUTA;
TransitionList UTALines("UTALines", &AnonStates);
long int nLevel1;
/**transition TauLines[NTAULINES+1];*/
TransitionList HFLines("HFLines", &AnonStates);
long int nHFLines;
//vector<vector<multi_arr<int,2> > > ipTransitions;
vector<vector<TransitionList> > Transitions;
multi_arr<int,2> ipFe2LevN;
static qList Fe2LevNStates;
TransitionList Fe2LevN("Fe2LevN", &Fe2LevNStates);
multi_arr<int,3> ipSatelliteLines; /* [ipISO][nelem][level] */
vector<vector<TransitionList> > SatelliteLines; /* [ipISO][nelem][level] */
TransitionList TauLine2("TauLine2", &AnonStates);
realnum *cs1_flag_lev2;

extern void checkTransitionListOfLists(vector<TransitionList>&list)
{ 
	for (vector<TransitionList>::iterator it=list.begin(); 
		  it != list.end(); ++it)
	{
		for (TransitionList::iterator tr = it->begin();
			  tr != it->end(); ++tr)
		{
			(*tr).check();
		}
		for (EmissionList::iterator em = it->Emis().begin(); 
			  em != it->Emis().end(); ++em)
		{
			(*em).check();
		}
	}
}

TransitionProxy::iterator TauDummy;
EmissionProxy DummyEmis;

namespace 
{
	class Init
	{
		EmissionList DummyEmisList;
		TransitionListImpl TauDummyTrans;
	public:
		Init(qList*states) : 
			DummyEmisList(&TauDummyTrans, 1), TauDummyTrans("TauDummy",states, 1)
		{
			DummyEmis = DummyEmisList[0];
			TauDummy=TauDummyTrans.begin();
		};
	};
	qList TauDummyStates(1);
	Init TauDummyInit(&TauDummyStates);;
}
