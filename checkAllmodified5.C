#include "TString.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include "TObjArray.h"
#include "TObjString.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include <vector>
#include <iomanip>
void PrintBusy(UInt_t *cnts);

UInt_t l0classB1;
UInt_t l0classA1;
UInt_t l1classB1;
UInt_t l1classA1;
UInt_t l2classB1;
UInt_t l2classA1;
//const UInt_t numclasses=50;
const UInt_t numclasses=100;
UInt_t l0clstT;
UInt_t l1clstT;
UInt_t l2clstT;
//const UInt_t numclusters=7; //Number of clusters (including T)
const UInt_t numclusters=9;
UInt_t fo1l0clstt;
UInt_t fo2l0clstt;
UInt_t fo3l0clstt;
UInt_t fo4l0clstt;
UInt_t fo5l0clstt;
UInt_t fo6l0clstt;
UInt_t fo1l1clstt;
UInt_t fo2l1clstt;
UInt_t fo3l1clstt;
UInt_t fo4l1clstt;
UInt_t fo5l1clstt;
UInt_t fo6l1clstt;
UInt_t fo1l0clst1;
UInt_t fo2l0clst1;
UInt_t fo3l0clst1;
UInt_t fo4l0clst1;
UInt_t fo5l0clst1;
UInt_t fo6l0clst1;
UInt_t fo1l1clst1;
UInt_t fo2l1clst1;
UInt_t fo3l1clst1;
UInt_t fo4l1clst1;
UInt_t fo5l1clst1;
UInt_t fo6l1clst1;
UInt_t l0strobe0;
UInt_t l0strobeIN;
UInt_t l1strobeOUT;
UInt_t l1strobeIN;
UInt_t l2strobeOUT;
UInt_t fo1l1strIN;
UInt_t fo1l2strIN;
UInt_t fo2l1strIN;
UInt_t fo2l2strIN;
UInt_t fo3l1strIN;
UInt_t fo3l2strIN;
UInt_t fo4l1strIN;
UInt_t fo4l2strIN;
UInt_t fo5l1strIN;
UInt_t fo5l2strIN;
UInt_t fo6l1strIN;
UInt_t fo6l2strIN;
UInt_t fo1glitchT;
UInt_t fo2glitchT;
UInt_t fo3glitchT;
UInt_t fo4glitchT;
UInt_t fo5glitchT;
UInt_t fo6glitchT;
UInt_t fo1l1spuriousT; 
UInt_t fo2l1spuriousT;
UInt_t fo3l1spuriousT;
UInt_t fo4l1spuriousT;
UInt_t fo5l1spuriousT;
UInt_t fo6l1spuriousT;
//const UInt_t runx=896;
const UInt_t runx=1486;
int flag = 0;
vector<double> classflag;
vector<double> clusterflag;
vector<double> clusterFOL0flag;
vector<double> clusterFOL1flag;
vector<double> L0strobeflag;
vector<double> L1strobeflag;
vector<double> FOL1strobeflag;
vector<double> FOL2strobeflag;
int glitch[numclusters] = {0};
int spurious[numclusters] = {0};
vector<int> locations;
double currentrun = -1;
//const UInt_t numcounters=970; //Number of counters (i.e. number of entries on each line of input file)
const UInt_t numcounters=1560;
double total[numcounters] = {0};
int num = 0;
ofstream outputfile;
int zeroflag[numcounters] = {0}; //If SOR/EOR or any between has unexplained zero scalers, any 'anomalies' found which use these scalers will be ignored
bool firstflag=0; //When firstflag=1, any anomalies are accompanied by a message warning that they were not started from a SOR/EOR reading
bool lastflag=0; //When lastflag=1, this is the final line in the file and anomalies should be read out, with a warning


void Plot(UInt_t *cnts, UInt_t *prev)
{

 cout.precision(12);
 outputfile.precision(12);

    if((cnts[runx] == currentrun)||(num==2))  // Run number
     {
      cout << "run found: " << cnts[runx] << endl;
      currentrun=cnts[runx];
      if(num==2) firstflag=1; //When beginning from start of file, comparisons between different counters will not be exactly valid  

	double increm[numcounters];
	double temp[numcounters];
    	for(int k=0; k<(numcounters); k++){
		//If counter is lower than in previous increment, it has overflowed and so is corrected by adding 2^32 to current counter
		temp[k] = cnts[k]; //Don't change the value of cnts itself
		if((temp[k]<prev[k])&&(temp[k]!=0)){ 
			temp[k] = pow(2,32) + temp[k];
		}
		increm[k] = temp[k] - prev[k]; //Get counters since last increment
		total[k] += increm[k]; //Add this to total
		if((temp[k]==0)&&(prev[k]!=0)&&(k!=runx)){ //Spurious zeroes
			cout << k << " = 0" << endl;
			outputfile << "Run " << prev[runx] << " num " << num << " " << k << " = 0" << endl;
			zeroflag[k]=2;

		}
    	}
	     		
   }
    if(cnts[runx] != currentrun || lastflag==1){ //If new run (or changed to 'no run'), print totals for previous run
   	if(prev[runx]==currentrun){
		double increm[numcounters];
		double temp[numcounters];
    		for(int k=0; k<(numcounters); k++){
			//If counter is lower than in previous increment, it has overflowed and so is corrected by adding 2^32 to current counter
			temp[k] = cnts[k]; //Don't change the value of cnts itself
			if((temp[k]<prev[k])&&(temp[k]!=0)){ 
				temp[k] = pow(2,32) + temp[k];
			}
			increm[k] = temp[k] - prev[k]; //Get counters since last increment
			total[k] += increm[k]; //Add this to total
			if((temp[k]==0)&&(prev[k]!=0)&&(k!=runx)){ //Spurious zeroes
				cout << k << " = 0" << endl;
				outputfile << "Run " << prev[runx] << " num " << num << " counter " << k << " = 0" << endl;
				zeroflag[k]++;
			}
    		}
	}
	currentrun = cnts[runx];
	for(int lzerob=l0classB1; lzerob<(l0classB1+numclasses); lzerob++){ //Positions of L0B,L0A,L1B, etc. in the array for the 50 classes
		int lzeroa = l0classA1 + lzerob - l0classB1;
		int loneb = l1classB1 + lzerob - l0classB1;		
		int lonea = l1classA1 + lzerob - l0classB1;
		int ltwob = l2classB1 + lzerob - l0classB1;
		int ltwoa = l2classA1 + lzerob - l0classB1;
		if((total[lzerob]>=total[lzeroa]) && (total[lzeroa]>=total[loneb]) && (total[loneb]>=total[lonea]) && (total[lonea]>=total[ltwob]) && (total[ltwob]>=total[ltwoa])){ //If L0B>=L0A>=L1B, etc. condition is true, print "no anomaly"
			cout << endl;
			cout << "NO anomaly in class " << lzerob-l0classB1+1 << endl;
			cout << "Run " << prev[runx] << endl;
			cout << "Increment number (overall): " << num << endl;				cout << "L0classB" << lzerob-l0classB1+1 << '\t' << "L0classA" << lzerob-l0classB1+1 << '\t' << "L1classB" << lzerob-l0classB1+1 << '\t' << "L1classA" << lzerob-l0classB1+1 << '\t' << "L2classB" << lzerob-l0classB1+1 << '\t' << "L2classA" << lzerob-l0classB1+1 << endl;
			cout << total[lzerob] << '\t' << total[lzeroa] << '\t' << total[loneb] << '\t' << total[lonea] << '\t' << total[ltwob] << '\t' << total[ltwoa] << endl;	
		}else{ //If condition is false, print details of anomaly to screen and file
			if((zeroflag[lzerob]==0) && (zeroflag[lzeroa]==0) && (zeroflag[loneb]==0) && (zeroflag[lonea]==0) && (zeroflag[ltwob]==0) && (zeroflag[ltwoa]==0)){ //If these are not due to spurious zeroes
			flag += 1;
			double diff=0; // Calculate the first difference (anomalous amount), e.g. +2
			if(total[lzerob]<total[lzeroa]) diff= total[lzeroa]-total[lzerob];
			if(total[lzeroa]<total[loneb])  diff= total[loneb]-total[lzeroa];
			if(total[loneb]<total[lonea]) diff= total[lonea]-total[loneb];
			if(total[lonea]<total[ltwob]) diff= total[ltwob]-total[lonea];
			if(total[ltwob]<total[ltwoa]) diff= total[ltwoa]-total[ltwob];
			classflag.push_back(diff); // Put anomalous amount into vector for class anomalies
			locations.push_back(prev[runx]); // Vector storing run numbers containing anomalies
			for(int i=1;i<6;i++){
				if(prev[runx+i]!=0) locations.push_back(prev[runx+i]); //Add any parallel runs to locations vector
			}
			cout << endl;
			outputfile << endl;
			cout << "Run " << prev[runx] << endl;
			outputfile << "Run " << prev[runx];
			if(prev[runx]==0) outputfile << " increment number " << num; //for gaps between runs
			for(int i=1;i<6;i++){ //Name parallel runs if there are any
				if(prev[runx+i]!=0) outputfile << " or Run " << prev[runx+i];
			}
			outputfile << endl;
			cout << "!!! Anomaly in class " << lzerob-l0classB1+1 << " !!!" << endl;
			outputfile << "!!! Anomaly in class " << lzerob-l0classB1+1 << " !!!" << endl;
			if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
			if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
			cout << "Increment number (overall): " << num << endl;
			cout << "L0classB" << lzerob-l0classB1+1 << '\t' << "L0classA" << lzerob-l0classB1+1 << '\t' << "L1classB" << lzerob-l0classB1+1 << '\t' << "L1classA" << lzerob-l0classB1+1 << '\t' << "L2classB" << lzerob-l0classB1+1 << '\t' << "L2classA" << lzerob-l0classB1+1 << endl;
			cout << total[lzerob] << '\t' << total[lzeroa] << '\t' << total[loneb] << '\t' << total[lonea] << '\t' << total[ltwob] << '\t' << total[ltwoa] << endl;
			outputfile << setw(12) << "L0classB" << lzerob-l0classB1+1 << setw(12) << "L0classA" << lzerob-l0classB1+1 << setw(12) << "L1classB" << lzerob-l0classB1+1 << setw(12) << "L1classA" << lzerob-l0classB1+1 << setw(12) << "L2classB" << lzerob-l0classB1+1 << setw(12) << "L2classA" << lzerob-l0classB1+1 << endl;
			outputfile << setw(12) << total[lzerob] << setw(12) << total[lzeroa] << setw(12) << total[loneb] << setw(12) << total[lonea] << setw(12) << total[ltwob] << setw(12) << total[ltwoa] << endl;
			}
		}
	}	
	for(int lzeroclst=l0clstT; lzeroclst<(l0clstT+numclusters); lzeroclst++){ //Positions of L0,L1,L2 in the array for the 7 clusters
		int loneclst = l1clstT + lzeroclst - l0clstT;
		int ltwoclst = l2clstT + lzeroclst - l0clstT;		
		if((total[lzeroclst]>=total[loneclst]) && (total[loneclst]>=total[ltwoclst])){ //If L0>=L1>=L2 condition is true, print "no anomaly"
			cout << endl;
			cout << "Run " << prev[runx] << endl;
			cout << "Increment number (overall): " << num << endl;		
			if(lzeroclst==l0clstT){
				cout << "NO anomaly in cluster T" << endl;
				cout << "L0clstT" << '\t' << "L1clstT" << '\t' << "L2clstT" << endl;
			}
			if(lzeroclst!=l0clstT){
				cout << "NO anomaly in cluster " << lzeroclst-l0clstT << endl;
				cout << "L0clst" << lzeroclst-l0clstT << '\t' << "L1clst" << lzeroclst-l0clstT << '\t' << "L2clst" << lzeroclst-l0clstT<< endl;
			}
			cout << total[lzeroclst] << '\t' << total[loneclst] << '\t' << total[ltwoclst] << endl;	
		}else{ //If condition is false, print details of anomaly to screen and file
			if((zeroflag[lzeroclst]==0)&&(zeroflag[loneclst]==0)&&(zeroflag[ltwoclst]==0)){ //If these are not due to spurious zeroes
			flag += 1;
			double diff=0; // Calculate the first difference (anomalous amount), e.g. +2
			if(total[lzeroclst]<total[loneclst]) diff= total[loneclst]-total[lzeroclst];
			if(total[loneclst]<total[ltwoclst])  diff= total[ltwoclst]-total[loneclst];
			clusterflag.push_back(diff); // Put anomalous amount into vector for class anomalies
			locations.push_back(prev[runx]); // Vector storing run numbers containing anomalies
			for(int i=1;i<6;i++){
				if(prev[runx+i]!=0) locations.push_back(prev[runx+i]); //Add any parallel runs to locations vector
			}
			cout << endl;
			outputfile << endl;
			cout << "Run " << prev[runx] << endl;
			outputfile << "Run " << prev[runx];
			if(prev[runx]==0) outputfile << " increment number " << num;
			for(int i=1;i<6;i++){ //Name parallel runs if there are any
				if(prev[runx+i]!=0) outputfile << " or Run " << prev[runx+i];
			}
			outputfile << endl;
			cout << "Increment number (overall): " << num << endl;		
			if(lzeroclst==l0clstT){
				cout << "!!! Anomaly in cluster T !!!" << endl;
				cout << "L0clstT" << '\t' << "L1clstT" << '\t' << "L2clstT" << endl;
				outputfile << "!!! Anomaly in cluster T !!!" << endl;
				if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR readings and so comparisons between counters may not be exactly valid)" << endl;
				if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
				outputfile << setw(12) << "L0clstT" << setw(12) << "L1clstT" << setw(12) << "L2clstT" << endl;

			}
			if(lzeroclst!=l0clstT){
				cout << "!!! Anomaly in cluster " << lzeroclst-l0clstT << " !!!"<< endl;
				cout << "L0clst" << lzeroclst-l0clstT << '\t' << "L1clst" << lzeroclst-l0clstT << '\t' << "L2clst" << lzeroclst-l0clstT<< endl;
				outputfile << "!!! Anomaly in cluster " << lzeroclst-l0clstT << " !!!"<< endl;
				if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR readings and so comparisons between counters may not be exactly valid)" << endl;
				if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
				outputfile << setw(12) << "L0clst" << lzeroclst-l0clstT << setw(12) << "L1clst" << lzeroclst-l0clstT << setw(12) << "L2clst" << lzeroclst-l0clstT<< endl;

			}
			cout << total[lzeroclst] << '\t' << total[loneclst] << '\t' << total[ltwoclst] << endl;	
			outputfile << setw(12) << total[lzeroclst] << setw(12) << total[loneclst] << setw(12) << total[ltwoclst] << endl;	
			}
		}
		int foonelzeroclst, fotwolzeroclst, fothreelzeroclst, fofourlzeroclst, fofivelzeroclst, fosixlzeroclst, fooneloneclst, fotwoloneclst, fothreeloneclst, fofourloneclst, fofiveloneclst, fosixloneclst;
		if(lzeroclst==l0clstT) { //If cluster T
			foonelzeroclst = fo1l0clstt;
			fotwolzeroclst = fo2l0clstt;
			fothreelzeroclst = fo3l0clstt;
			fofourlzeroclst = fo4l0clstt;
			fofivelzeroclst = fo5l0clstt;
			fosixlzeroclst = fo6l0clstt;
			fooneloneclst = fo1l1clstt;
			fotwoloneclst = fo2l1clstt;
			fothreeloneclst = fo3l1clstt;
			fofourloneclst = fo4l1clstt;
			fofiveloneclst = fo5l1clstt;
			fosixloneclst = fo6l1clstt;
		}
		if(lzeroclst!=l0clstT) { //If clusters 1-6
			foonelzeroclst = fo1l0clst1 + lzeroclst - l0clstT -1;
			fotwolzeroclst = fo2l0clst1 + lzeroclst - l0clstT -1;
			fothreelzeroclst = fo3l0clst1 + lzeroclst - l0clstT -1;
			fofourlzeroclst = fo4l0clst1 + lzeroclst - l0clstT -1;
			fofivelzeroclst = fo5l0clst1 + lzeroclst - l0clstT -1;
			fosixlzeroclst = fo6l0clst1 + lzeroclst - l0clstT -1;
			fooneloneclst = fo1l1clst1 + lzeroclst - l0clstT -1;
			fotwoloneclst = fo2l1clst1 + lzeroclst - l0clstT -1;
			fothreeloneclst = fo3l1clst1 + lzeroclst - l0clstT -1;
			fofourloneclst = fo4l1clst1 + lzeroclst - l0clstT -1;
			fofiveloneclst = fo5l1clst1 + lzeroclst - l0clstT -1;
			fosixloneclst = fo6l1clst1 + lzeroclst - l0clstT -1;
		}
		if((total[lzeroclst]==total[foonelzeroclst]) && (total[lzeroclst]==total[fotwolzeroclst]) && (total[lzeroclst]==total[fothreelzeroclst]) && (total[lzeroclst]==total[fofourlzeroclst]) && (total[lzeroclst]==total[fofivelzeroclst]) && (total[lzeroclst]==total[fosixlzeroclst])){ //Check L0 clusters equal L0 cluster FOs
			if(lzeroclst==l0clstT) cout << "No anomaly between cluster T and FOs for L0" << endl;
			if(lzeroclst!=l0clstT) cout << "No anomaly between cluster " << lzeroclst-l0clstT << " and FOs for L0" << endl;
		}else{
			if((zeroflag[lzeroclst]==0)&&(zeroflag[foonelzeroclst]==0)&&(zeroflag[fotwolzeroclst]==0)&&(zeroflag[fothreelzeroclst]==0)&&(zeroflag[fofourlzeroclst]==0)&&(zeroflag[fofivelzeroclst]==0)&&(zeroflag[fosixlzeroclst]==0)){ //If these are not due to spurious zeroes
			flag+=1;
			double diff=0; // Calculate the first difference (anomalous amount), e.g. +2
			if(diff==0) diff= total[foonelzeroclst]-total[lzeroclst];
			if(diff==0)  diff= total[fotwolzeroclst]-total[lzeroclst];
			if(diff==0) diff= total[fothreelzeroclst]-total[lzeroclst];
			if(diff==0)  diff= total[fofourlzeroclst]-total[lzeroclst];
			if(diff==0) diff= total[fofivelzeroclst]-total[lzeroclst];
			if(diff==0)  diff= total[fosixlzeroclst]-total[lzeroclst];
			clusterFOL0flag.push_back(diff); // Put anomalous amount into vector for FO L0 cluster anomalies
			locations.push_back(prev[runx]); // Vector storing run numbers containing anomalies
			for(int i=1;i<6;i++){
				if(prev[runx+i]!=0) locations.push_back(prev[runx+i]); //Add any parallel runs to locations vector
			}
			outputfile << endl;
			outputfile << "Run " << prev[runx];
			if(prev[runx]==0) outputfile << " increment number " << num;
			for(int i=1;i<6;i++){ //Name parallel runs if there are any
				if(prev[runx+i]!=0) outputfile << " or Run " << prev[runx+i];
			}
			outputfile << endl;
			if(lzeroclst==l0clstT){
				cout << "!!! Anomaly between cluster T and FOs for L0!!!" << endl;
				cout << "L0clstT" << '\t' << "FO1L0clstT" << '\t' << "FO2L0clstT" << '\t' << "FO3L0clstT" << '\t' << "FO4L0clstT" << '\t' << "FO5L0clstT" << '\t' << "FO6L0clstT" << endl;
				outputfile << "!!! Anomaly between cluster T and FOs for L0!!!" << endl;
				if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR readings and so comparisons between counters may not be exactly valid)" << endl;
				if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
				outputfile << setw(12) << "L0clstT" << setw(12) <<  "FO1L0clstT" << setw(12) << "FO2L0clstT" << setw(12) <<  "FO3L0clstT" << setw(12) << "FO4L0clstT" << setw(12) <<  "FO5L0clstT" << setw(12) <<  "FO6L0clstT" << endl;

			}
			if(lzeroclst!=l0clstT){
				cout << "!!! Anomaly between cluster " << lzeroclst-l0clstT << " and FOs for L0 !!!" << endl;
				cout << "L0clst" << lzeroclst-l0clstT << '\t' << "F01L0clst" << lzeroclst-l0clstT << '\t' << "F02L0clst" << lzeroclst-l0clstT << '\t' << "FO3L0clst" << lzeroclst-l0clstT << '\t' << "FO4L0clst" << lzeroclst-l0clstT << '\t' << "FO5L0clst" << lzeroclst-l0clstT << '\t' << "FO6L0clst" << lzeroclst-l0clstT << endl;
				outputfile << "!!! Anomaly between cluster " << lzeroclst-l0clstT << " and FOs for L0 !!!" << endl;
				if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR readings and so comparisons between counters may not be exactly valid)" << endl;
				if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
				outputfile << setw(12) << "L0clst" << lzeroclst-l0clstT << setw(12) << "F01L0clst" << lzeroclst-l0clstT << setw(12) << "F02L0clst" << lzeroclst-l0clstT << setw(12) << "FO3L0clst" << lzeroclst-l0clstT << setw(12) << "FO4L0clst" << lzeroclst-l0clstT << setw(12) << "FO5L0clst" << lzeroclst-l0clstT << setw(12) << "FO6L0clst" << lzeroclst-l0clstT << endl;

			}
			cout << total[lzeroclst] << '\t' << total[foonelzeroclst] << '\t' << total[fotwolzeroclst] << '\t' << total[fothreelzeroclst] << '\t' << total[fofourlzeroclst] << '\t' << total[fofivelzeroclst] << '\t' << total[fosixlzeroclst] << endl;
			outputfile << setw(12) << total[lzeroclst] << setw(12) << total[foonelzeroclst] << setw(12) << total[fotwolzeroclst] << setw(12) << total[fothreelzeroclst] << setw(12) << total[fofourlzeroclst] << setw(12) << total[fofivelzeroclst] << setw(12) << total[fosixlzeroclst] << endl;
			}
		}
		if((total[loneclst]==total[fooneloneclst]) && (total[loneclst]==total[fotwoloneclst]) && (total[loneclst]==total[fothreeloneclst]) && (total[loneclst]==total[fofourloneclst]) && (total[loneclst]==total[fofiveloneclst]) && (total[loneclst]==total[fosixloneclst])){ //Check L1 clusters equal L1 FOs
			if(lzeroclst==l0clstT) cout << "No anomaly between cluster T and FOs for L1" << endl;
			if(lzeroclst!=l0clstT) cout << "No anomaly between cluster " << lzeroclst-l0clstT << " and FOs for L1" << endl;
		}else{ //If condition is false, print details to screen and to file
			if((zeroflag[loneclst]==0)&&(zeroflag[fooneloneclst]==0)&&(zeroflag[fotwoloneclst]==0)&&(zeroflag[fothreeloneclst]==0)&&(zeroflag[fofourloneclst]==0)&&(zeroflag[fofiveloneclst]==0)&&(zeroflag[fosixloneclst]==0)){ //If these are not due to spurious zeroes
			flag+=1;
			double diff=0; // Calculate the first difference (anomalous amount), e.g. +2
			if(diff==0) diff= total[fooneloneclst]-total[loneclst];
			if(diff==0)  diff= total[fotwoloneclst]-total[loneclst];
			if(diff==0) diff= total[fothreeloneclst]-total[loneclst];
			if(diff==0)  diff= total[fofourloneclst]-total[loneclst];
			if(diff==0) diff= total[fofiveloneclst]-total[loneclst];
			if(diff==0)  diff= total[fosixloneclst]-total[loneclst];
			clusterFOL1flag.push_back(diff); // Put anomalous amount into vector for FO L1 cluster anomalies
			locations.push_back(prev[runx]);
			for(int i=1;i<6;i++){
				if(prev[runx+i]!=0) locations.push_back(prev[runx+i]); //Add any parallel runs to locations vector
			}
			outputfile << endl;
			outputfile << "Run " << prev[runx];
			if(prev[runx]==0) outputfile << " increment number " << num;
			for(int i=1;i<6;i++){ //Name parallel runs if there are any
				if(prev[runx+i]!=0) outputfile << " or Run " << prev[runx+i];
			}
			outputfile << endl;
			if(lzeroclst==l0clstT){
				cout << "!!! Anomaly between cluster T and FOs for L1!!!" << endl;
				cout << "L1clstT" << '\t' << "FO1L1clstT" << '\t' << "FO2L1clstT" << '\t' << "FO3L1clstT" << '\t' << "FO4L1clstT" << '\t' << "FO5L1clstT" << '\t' << "FO6L1clstT" << endl;
				outputfile << "!!! Anomaly between cluster T and FOs for L1!!!" << endl;
				if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
				if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
				outputfile << setw(12) << "L1clstT" << setw(12) << "FO1L1clstT" << setw(12) << "FO2L1clstT" << setw(12) << "FO3L1clstT" << setw(12) << "FO4L1clstT" << setw(12) << "FO5L1clstT" << setw(12) << "FO6L1clstT" << endl;

			}
			if(lzeroclst!=l0clstT){
				cout << "!!! Anomaly between cluster " << lzeroclst-l0clstT << " and FOs for L1 !!!" << endl;
				cout << "L1clst" << lzeroclst-l0clstT << '\t' << "F01L1clst" << lzeroclst-l0clstT << '\t' << "F02L1clst" << lzeroclst-l0clstT << '\t' << "FO3L1clst" << lzeroclst-l0clstT << '\t' << "FO4L1clst" << lzeroclst-l0clstT << '\t' << "FO5L1clst" << lzeroclst-l0clstT << '\t' << "FO6L1clst" << lzeroclst-l0clstT << endl;
				outputfile << "!!! Anomaly between cluster " << lzeroclst-l0clstT << " and FOs for L1 !!!" << endl;
				if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
				if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
				outputfile << setw(12) << "L1clst" << lzeroclst-l0clstT << setw(12) << "F01L1clst" << lzeroclst-l0clstT << setw(12) << "F02L1clst" << lzeroclst-l0clstT << setw(12) << "FO3L1clst" << lzeroclst-l0clstT << setw(12) << "FO4L1clst" << lzeroclst-l0clstT << setw(12) << "FO5L1clst" << lzeroclst-l0clstT << setw(12) << "FO6L1clst" << lzeroclst-l0clstT << endl;

			}
			cout << total[loneclst] << '\t' << total[fooneloneclst] << '\t' << total[fotwoloneclst] << '\t' << total[fothreeloneclst] << '\t' << total[fofourloneclst] << '\t' << total[fofiveloneclst] << '\t' << total[fosixloneclst] << endl;
			outputfile << setw(12) << total[loneclst] << setw(12) << total[fooneloneclst] << setw(12) << total[fotwoloneclst] << setw(12) << total[fothreeloneclst] << setw(12) << total[fofourloneclst] << setw(12) << total[fofiveloneclst] << setw(12) << total[fosixloneclst] << endl;
			}
		}

	}
	if(total[l0strobe0]==total[l0strobeIN]){ //Check L0strobe0 = L0strobeIN
		cout << endl;
		cout << "No anomaly between L0strobe0 and L0strobeIN" << endl;
	}else{
		if((zeroflag[l0strobe0]==0)&&(zeroflag[l0strobeIN]==0)){ //If these are not due to spurious zeroes
		flag+=1;
		double diff=0; // Calculate the difference (anomalous amount), e.g. +2
		diff = total[l0strobeIN]-total[l0strobe0];
		L0strobeflag.push_back(diff); // Put anomalous amount into vector for L0 strobe anomalies
		locations.push_back(prev[runx]);
			for(int i=1;i<6;i++){
				if(prev[runx+i]!=0) locations.push_back(prev[runx+i]); //Add any parallel runs to locations vector
			}
		cout << endl;
		outputfile << endl;
		outputfile << "Run " << prev[runx];
		if(prev[runx]==0) outputfile << " increment number " << num;
		for(int i=1;i<6;i++){ //Name parallel runs if there are any
			if(prev[runx+i]!=0) outputfile << " or Run " << prev[runx+i];
		}
		outputfile << endl;
		cout << "!!! Anomaly between L0strobe0 and L0strobeIN" << endl;
		cout << "L0strobe0" << '\t' << "L0strobeIN" << endl;
		cout << total[l0strobe0] << '\t' << total[l0strobeIN] << endl;
		outputfile << "!!! Anomaly between L0strobe0 and L0strobeIN" << endl;
		if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		outputfile << setw(12) << "L0strobe0" << setw(12) << "L0strobeIN" << endl;
		outputfile << setw(12) << total[l0strobe0] << setw(12) << total[l0strobeIN] << endl;
		}
	}
	if(total[l1strobeOUT]==total[l1strobeIN]){ //Check L1strobeOUT = L1strobeIN
		cout << endl;
		cout << "No anomaly between L1strobeOUT and L1strobeIN" << endl;
	}else{
		if((zeroflag[l1strobeOUT]==0)&&(zeroflag[l1strobeIN]==0)){ //If these are not due to spurious zeroes
		flag+=1;
		double diff=0; // Calculate the difference (anomalous amount), e.g. +2
		diff = total[l1strobeIN]-total[l1strobeOUT];
		L1strobeflag.push_back(diff); // Put anomalous amount into vector for L1 strobe anomalies
		locations.push_back(prev[runx]);
		for(int i=1;i<6;i++){
			if(prev[runx+i]!=0) locations.push_back(prev[runx+i]); //Add any parallel runs to locations vector
		}
		cout << endl;
		outputfile << endl;
		outputfile << "Run " << prev[runx];
		if(prev[runx]==0) outputfile << " increment number " << num;
		for(int i=1;i<6;i++){ //Name parallel runs if there are any
			if(prev[runx+i]!=0) outputfile << " or Run " << prev[runx+i];
		}
		outputfile << endl;
		cout << "!!! Anomaly between L1strobeOUT and L1strobeIN" << endl;
		cout << "L1strobeOUT" << '\t' << "L1strobeIN" << endl;
		cout << total[l1strobeOUT] << '\t' << total[l1strobeIN] << endl;
		outputfile << "!!! Anomaly between L1strobeOUT and L1strobeIN" << endl;
		if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		outputfile << setw(12) << "L1strobeOUT" << setw(12) << "L1strobeIN" << endl;
		outputfile << setw(12) << total[l1strobeOUT] << setw(12) << total[l1strobeIN] << endl;
		}
	}
	if((total[l1strobeOUT]==total[fo1l1strIN]) && (total[l1strobeOUT]==total[fo2l1strIN]) && (total[l1strobeOUT]==total[fo3l1strIN]) && (total[l1strobeOUT]==total[fo4l1strIN]) && (total[l1strobeOUT]==total[fo5l1strIN]) && (total[l1strobeOUT]==total[fo6l1strIN])){ //Check L1strobeOUT = all L1 FO strobes IN
		cout << endl;
		cout << "No anomaly between L1strobeOUT and L1 FO strobes IN" << endl;
	}else{
		if((zeroflag[l1strobeOUT]==0)&&(zeroflag[fo1l1strIN]==0)&&(zeroflag[fo2l1strIN]==0)&&(zeroflag[fo3l1strIN]==0)&&(zeroflag[fo4l1strIN]==0)&&(zeroflag[fo5l1strIN]==0)&&(zeroflag[fo6l1strIN]==0)){ //If these are not due to spurious zeroes
		flag+=1;
		double diff=0; // Calculate the first difference (anomalous amount), e.g. +2
		if(diff==0) diff= total[fo1l1strIN]-total[l1strobeOUT];
		if(diff==0)  diff= total[fo2l1strIN]-total[l1strobeOUT];
		if(diff==0) diff= total[fo3l1strIN]-total[l1strobeOUT];
		if(diff==0)  diff= total[fo4l1strIN]-total[l1strobeOUT];
		if(diff==0) diff= total[fo5l1strIN]-total[l1strobeOUT];
		if(diff==0)  diff= total[fo6l1strIN]-total[l1strobeOUT];
		FOL1strobeflag.push_back(diff); // Put anomalous amount into vector for L1 FO strobe anomalies
		locations.push_back(prev[runx]);
		for(int i=1;i<6;i++){
			if(prev[runx+i]!=0) locations.push_back(prev[runx+i]); //Add any parallel runs to locations vector
		}
		cout << endl;
		outputfile << endl;
		outputfile << "Run " << prev[runx];
		if(prev[runx]==0) outputfile << " increment number " << num;
		for(int i=1;i<6;i++){ //Name parallel runs if there are any
			if(prev[runx+i]!=0) outputfile << " or Run " << prev[runx+i];
		}
		outputfile << endl;
		cout << "!!! Anomaly between L1strobeOUT and L1 FO strobes IN" << endl;
		cout << "L1strobeOUT" << '\t' << "fo1L1strIN" << '\t' << "fo2L1strIN" << '\t' << "fo3L1strIN" << '\t' << "fo4L1strIN" << '\t' << "fo5L1strIN" << '\t' << "fo6L1strIN" << endl;
		cout << total[l1strobeOUT] << '\t' << total[fo1l1strIN] << '\t' << total[fo2l1strIN] << '\t' << total[fo3l1strIN] << '\t' << total[fo4l1strIN] << '\t' << total[fo5l1strIN] << '\t' << total[fo6l1strIN] << endl;
		outputfile << "!!! Anomaly between L1strobeOUT and L1 FO strobes IN" << endl;
		if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		outputfile << setw(12) << "L1strobeOUT" << setw(12) << "fo1L1strIN" << setw(12) << "fo2L1strIN" << setw(12) << "fo3L1strIN" << setw(12) << "fo4L1strIN" << setw(12) << "fo5L1strIN" << setw(12) << "fo6L1strIN" << endl;
		outputfile << setw(12) << total[l1strobeOUT] << setw(12) << total[fo1l1strIN] << setw(12) << total[fo2l1strIN] << setw(12) << total[fo3l1strIN] << setw(12) << total[fo4l1strIN] << setw(12) << total[fo5l1strIN] << setw(12) << total[fo6l1strIN] << endl;
		}
	}
	if((total[l2strobeOUT]==total[fo1l2strIN]) && (total[l2strobeOUT]==total[fo2l2strIN]) && (total[l2strobeOUT]==total[fo3l2strIN]) && (total[l2strobeOUT]==total[fo4l2strIN]) && (total[l2strobeOUT]==total[fo5l2strIN]) && (total[l2strobeOUT]==total[fo6l2strIN])){ //Check L2strobeOUT = all L2 FO strobes IN
		cout << endl;
		cout << "No anomaly between L2strobeOUT and L2 FO strobes IN !!!" << endl;
	}else{
		if((zeroflag[l2strobeOUT]==0)&&(zeroflag[fo1l2strIN]==0)&&(zeroflag[fo2l2strIN]==0)&&(zeroflag[fo3l2strIN]==0)&&(zeroflag[fo4l2strIN]==0)&&(zeroflag[fo5l2strIN]==0)&&(zeroflag[fo6l2strIN]==0)){ //If these are not due to spurious zeroes
		flag+=1;
		double diff=0; // Calculate the first difference (anomalous amount), e.g. +2
		if(diff==0) diff= total[fo1l2strIN]-total[l2strobeOUT];
		if(diff==0)  diff= total[fo2l2strIN]-total[l2strobeOUT];
		if(diff==0) diff= total[fo3l2strIN]-total[l2strobeOUT];
		if(diff==0)  diff= total[fo4l2strIN]-total[l2strobeOUT];
		if(diff==0) diff= total[fo5l2strIN]-total[l2strobeOUT];
		if(diff==0)  diff= total[fo6l2strIN]-total[l2strobeOUT];
		FOL2strobeflag.push_back(diff); // Put anomalous amount into vector for L2 FO strobe anomalies
		locations.push_back(prev[runx]);
		for(int i=1;i<6;i++){
			if(prev[runx+i]!=0) locations.push_back(prev[runx+i]); //Add any parallel runs to locations vector
		}
		cout << endl;
		outputfile << endl;
		outputfile << "Run " << prev[runx];
		if(prev[runx]==0) outputfile << " increment number " << num;
		for(int i=1;i<6;i++){ //Name parallel runs if there are any
			if(prev[runx+i]!=0) outputfile << " or Run " << prev[runx+i];
		}
		outputfile << endl;
		cout << "!!! Anomaly between L2strobeOUT and L2 FO strobes IN !!!" << endl;
		cout << "L2strobeOUT" << '\t' << "fo1L2strIN" << '\t' << "fo2L2strIN" << '\t' << "fo3L2strIN" << '\t' << "fo4L2strIN" << '\t' << "fo5L2strIN" << '\t' << "fo6L2strIN" << endl;
		cout << total[l2strobeOUT] << '\t' << total[fo1l2strIN] << '\t' << total[fo2l2strIN] << '\t' << total[fo3l2strIN] << '\t' << total[fo4l2strIN] << '\t' << total[fo5l2strIN] << '\t' << total[fo6l2strIN] << endl;
		outputfile << "!!! Anomaly between L2strobeOUT and L2 FO strobes IN" << endl;
		if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		outputfile << setw(12) << "L2strobeOUT" << setw(12) << "fo1L2strIN" << setw(12) << "fo2L2strIN" << setw(12) << "fo3L2strIN" << setw(12) << "fo4L2strIN" << setw(12) << "fo5L2strIN" << setw(12) << "fo6L2strIN" << endl;
		outputfile << setw(12) << total[l2strobeOUT] << setw(12) << total[fo1l2strIN] << setw(12) << total[fo2l2strIN] << setw(12) << total[fo3l2strIN] << setw(12) << total[fo4l2strIN] << setw(12) << total[fo5l2strIN] << setw(12) << total[fo6l2strIN] << endl;
		}
	}
	for(int fonum =1; fonum<7; fonum++){ //Check for any glitches, FO boards 1-6
		int glitchT;
		if(fonum==1) glitchT = fo1glitchT;
		if(fonum==2) glitchT = fo2glitchT;
		if(fonum==3) glitchT = fo3glitchT;
		if(fonum==4) glitchT = fo4glitchT;
		if(fonum==5) glitchT = fo5glitchT;
		if(fonum==6) glitchT = fo6glitchT;
		for(int pos = glitchT; pos<glitchT+numclusters; pos++){ //Check for any glitches, clusters T-6
			if(total[pos]==0){
			//	cout << "No glitch in FO" << (fonum/48)+1 << "glitch";
				if(glitchT==fo1glitchT) {
				//	cout << "T" << endl;
				}else{ 
				//	cout << (foglitchnum-fo1glitchT) << endl;
				}
			}else{
				if(zeroflag[pos]==0){ //If these are not due to spurious zeroes
				cout << "Glitch in FO" << fonum << "glitch";
				glitch[pos-glitchT]+=total[pos]; //Add glitch to total glitches for that cluster
				if((pos-glitchT)==0) {
					cout << "T:  ";
				}else{
					cout << (pos-glitchT) << ":  ";	
				}
				cout << total[pos] << endl;
				}
			}
		}
	}
	for(int fonum =1; fonum<7; fonum++){ //Check for any glitches, FO boards 1-6
		int spuriousT;
		if(fonum==1) spuriousT = fo1l1spuriousT;
		if(fonum==2) spuriousT = fo2l1spuriousT;
		if(fonum==3) spuriousT = fo3l1spuriousT;
		if(fonum==4) spuriousT = fo4l1spuriousT;
		if(fonum==5) spuriousT = fo5l1spuriousT;
		if(fonum==6) spuriousT = fo6l1spuriousT;
		for(int pos = spuriousT; pos<spuriousT+numclusters; pos++){ //Check for any spurious, clusters T-numclusters
			if(total[pos]==0){
				if(spuriousT==fo1l1spuriousT) {
				}else{ 
				}
			}else{
				if(zeroflag[pos]==0){ //If these are not due to spurious zeroes
				cout << "Spurious in FO" << fonum << "l1spurious";
				spurious[pos-spuriousT]+=total[pos]; //Add spurious to total spurious for that cluster
				if((pos-spuriousT)==0) {
					cout << "T:  ";
				}else{
					cout << (pos-spuriousT) << ":  ";	
				}
				cout << total[pos] << endl;
				}
			}
		}
	}
	
	cout << endl;
	for(UInt_t j=0; j<numcounters; j++){ //Reset totals for next run
		total[j]=0;	
		if(zeroflag[j]>1) zeroflag[j]=0; //After reading with zeroes has been accounted for EOR and SOR, reset to zero	
	}
	firstflag=0; //Reset firstflag to zero
    }
   
  
 }

//Function to print the anomaly type, number and size at end of output file
void PrintAnomaly(vector<double> anomalyvector, TString violations){
	cout << endl;
	outputfile << endl;
	cout << anomalyvector.size() << violations << endl;
	outputfile << anomalyvector.size() << violations << endl;
	for(int l=0; l<anomalyvector.size(); l++){
		cout << " " << anomalyvector[l];
		outputfile << " " << anomalyvector[l];
	}
}

//Function ReadLines opens file called name, reads each line into cnts an sends it and prev to Plot function
void ReadLines(TString name, int numberoflines, UInt_t *prev){
  	ifstream *file = new ifstream(name.Data());
 	TString strLine;
  	UInt_t nlines=0;
 	while (strLine.ReadLine(*file) && nlines<1) {
  		TObjArray *tokens = strLine.Tokenize(" ");
  		Int_t ntokens = tokens->GetEntriesFast();
   		if(ntokens != numcounters+2){
    			cout << "Unexpected number of counters: " << ntokens << endl;
    			continue;
 		}
    		// Array of counters as described in sorted.cnames2
    		UInt_t cnts[numcounters];
    		for(Int_t i=2;i<ntokens;i++){
      			string cnt(((TObjString*)tokens->At(i))->String()); 
      			stringstream convert(cnt); 
     			 UInt_t hcnt;
      			convert >> std::hex >> hcnt;
      			cnts[i-2]=hcnt;
    		}
    		num++;
    		// Find the counters you want
   		if(num==numberoflines) lastflag=1; //Set lastflag=1 so that anomalies for period up to end of file are printed, with a warning
   		if(num!=1) Plot(cnts, prev); //Apart from the first line, which has no prev, do function
    		cout << "Increment number: " << num << endl;
    		for(int i=0; i<(ntokens-2); i++){
    			prev[i] = cnts[i]; //Assign current counters to be prev counters for next while loop
   		}
	}
}


void PrintBusy(UInt_t* cnts)
{
 cout << cnts[921]-tims << " " ;
// cout << (cnts[781]-tim0) << " ";
 //for(int i=0;i<24;i++) cout << Double_t(cnts[742+i]-busy0[i])/Double_t(cnts[781]-tim0) << " ";
 for(int i=0;i<24;i++){
    UInt_t cor1=0,cor2=0;
    if(cnts[781]<tim0)cor1 += 0xffffffff;
    if(cnts[742+i]<busy0[i])cor2 += 0xffffffff;
    printf("%3.2f ",(0.4*(cnts[742+i]-busy0[i]+cor2))/(0.4*(cnts[781]-tim0+cor1)));
    //printf("%u ",cnts[742+i]-busy0[i]);
 }
 cout << endl;
}

//main function
void anal2()
{
 // Identify nfiles file names from first one given onwards and add to vector
 // Only works when all days are in same month
 TString name("raw112014/rawcnts/01.11.2014.rawcnts");
 TString tempname = name; //copy of first file name string to be changed to subsequent file names
 TString namecopy = name; //copy of first file name string to be cut to just the day of the month
 int nfiles = 5; //number of files total to be analysed
 TString daystring = tempname.Remove(0,18);
 daystring.Remove(20,16);
 Int_t firstday = atoi(daystring.Data());
 
 //TString name2("raw112014/rawcnts/02.11.2014.rawcnts");
 //TString name("rawcnts/29.01.2013.rawcnts");
 //Loop over all nfiles filenames and open them
 vector<TString> filenames;
 filenames.push_back(name);
 for(int i= firstday;i<nfiles;i++){
	stringstream sts;
	sts << i+1;
	TString tempstr = sts.str();
	if(tempstr.Length()<2) tempstr.Prepend("0");
	namecopy.Replace(18,2,tempstr); //This '18' must be changed to the index in the full filename of the first digit of the day
 	filenames.push_back(namecopy);
 }
 // Parse file 
 Int_t nlines=0;
 //int num = 0;
 UInt_t prev[numcounters] = {0}; //Array for storing counters from previous line to be compared with
 //Name and open output file
 TString outputname = name;
 if(nfiles>1){
	 outputname.Append(".plus");
	 int ndays = (int)(nfiles) -1;
	 stringstream strs;
	 strs << ndays;
	 TString ndaysstring = strs.str();
	 outputname.Append(ndaysstring);
	 outputname.Append("days");
 }
 outputname.Append(".anomalies.txt");
 outputfile.open (outputname.Data());

 //Read in the needed counter positions from the cnames.sorted data file
 TString sortedname("cnames.sorted2.2014");
 //TString sortedname("cnames.sorted2");
 ifstream sortedfile(sortedname.Data());
 std::string cname;
 UInt_t position;
 std::string sortedline;
 int j=0; //So only the first run number counter is assigned
 while ( std::getline(sortedfile,sortedline)){
    sortedfile >> cname >> position;
    if(cname=="l0classB1") l0classB1=position;
    else if(cname=="l0classA1") l0classA1=position;
    else if(cname=="l1classB1") l1classB1=position;
    else if(cname=="l1classA1") l1classA1=position;
    else if(cname=="l2classB1") l2classB1=position;
    else if(cname=="l2classA1") l2classA1=position;
    else if(cname=="l0clstT") l0clstT=position;
    else if(cname=="l1clstT") l1clstT=position;
    else if(cname=="l2clstT") l2clstT=position;
    else if(cname=="fo1l0clstt") fo1l0clstt=position;
    else if(cname=="fo2l0clstt") fo2l0clstt=position;
    else if(cname=="fo3l0clstt") fo3l0clstt=position;
    else if(cname=="fo4l0clstt") fo4l0clstt=position;
    else if(cname=="fo5l0clstt") fo5l0clstt=position;
    else if(cname=="fo6l0clstt") fo6l0clstt=position;
    else if(cname =="fo1l1clstt") fo1l1clstt=position;
    else if(cname=="fo2l1clstt") fo2l1clstt=position;
    else if(cname=="fo3l1clstt") fo3l1clstt=position;
    else if(cname=="fo4l1clstt") fo4l1clstt=position;
    else if(cname=="fo5l1clstt") fo5l1clstt=position;
    else if(cname=="fo6l1clstt") fo6l1clstt=position;
    else if(cname=="fo1l0clst1") fo1l0clst1=position;
    else if(cname=="fo2l0clst1") fo2l0clst1=position;
    else if(cname=="fo3l0clst1") fo3l0clst1=position;
    else if(cname=="fo4l0clst1") fo4l0clst1=position;
    else if(cname=="fo5l0clst1") fo5l0clst1=position;
    else if(cname=="fo6l0clst1") fo6l0clst1=position;
    else if(cname=="fo1l1clst1") fo1l1clst1=position;
    else if(cname=="fo2l1clst1") fo2l1clst1=position;
    else if(cname=="fo3l1clst1") fo3l1clst1=position;
    else if(cname=="fo4l1clst1") fo4l1clst1=position;
    else if(cname=="fo5l1clst1") fo5l1clst1=position;
    else if(cname=="fo6l1clst1") fo6l1clst1=position;
    else if(cname=="l0strobe0") l0strobe0=position;
    else if(cname=="l0strobeIN") l0strobeIN=position;
    else if(cname=="l1strobeOUT") l1strobeOUT=position;
    else if(cname=="l1strobeIN") l1strobeIN=position;
    else if(cname=="l2strobeOUT") l2strobeOUT=position;
    else if(cname=="fo1l1strIN") fo1l1strIN=position;
    else if(cname=="fo1l2strIN") fo1l2strIN=position;
    else if(cname=="fo2l1strIN") fo2l1strIN=position;
    else if(cname=="fo2l2strIN") fo2l2strIN=position;
    else if(cname=="fo3l1strIN") fo3l1strIN=position;
    else if(cname=="fo3l2strIN") fo3l2strIN=position;
    else if(cname=="fo4l1strIN") fo4l1strIN=position;
    else if(cname=="fo4l2strIN") fo4l2strIN=position;
    else if(cname=="fo5l1strIN") fo5l1strIN=position;
    else if(cname=="fo5l2strIN") fo5l2strIN=position;
    else if(cname=="fo6l1strIN") fo6l1strIN=position;
    else if(cname=="fo6l2strIN") fo6l2strIN=position;
    else if(cname=="fo1glitchT") fo1glitchT=position;
    else if(cname=="fo2glitchT") fo2glitchT=position;
    else if(cname=="fo3glitchT") fo3glitchT=position;
    else if(cname=="fo4glitchT") fo4glitchT=position;
    else if(cname=="fo5glitchT") fo5glitchT=position;
    else if(cname=="fo6glitchT") fo6glitchT=position;
    else if(cname=="fo1l1spuriousT") fo1l1spuriousT=position;
    else if(cname=="fo2l1spuriousT") fo2l1spuriousT=position;
    else if(cname=="fo3l1spuriousT") fo3l1spuriousT=position;
    else if(cname=="fo4l1spuriousT") fo4l1spuriousT=position;
    else if(cname=="fo5l1spuriousT") fo5l1spuriousT=position;
    else if(cname=="fo6l1spuriousT") fo6l1spuriousT=position;

 }
 sortedfile.close();

 //Count the number of lines in the file
 int numberoflines=0;
 for(int i=0;i<nfiles;i++){
	TString nameForCounting = filenames[i];
	ifstream *fileForCounting = new ifstream(nameForCounting.Data());
	TString strLineForCounting;
 	while(strLineForCounting.ReadLine(*fileForCounting) && nlines<1){
		numberoflines++;
 	}
 }
 cout << numberoflines << endl;

 for(int i=0;i<nfiles;i++){
 	ReadLines(filenames[i], numberoflines, prev);
 }
 
 

// UInt_t empty[numcounters] = 0;
// Plot(empty, prev); //If final increment is a run, this will print its totals
 cout << endl;
 outputfile << endl;
 if(flag == 0){ //Print whether there were any anomalies, how many, and at which runs
	cout << "No anomalies found in data file" << endl;
	outputfile << "No anomalies found in data file" << endl;
 } else {
	cout << flag << " anomalies found in data file, in run numbers (0 is between runs): ";
	outputfile << flag << " anomalies found in data file, in run numbers (0 is between runs): ";
	for(int a=0;a<flag;a++){ 
		cout << locations[a] << ", ";
		outputfile << locations[a] <<", ";
	}
	cout << "." << endl;
	outputfile << "." << endl;
	//Breakdown of anomaly categories
	cout << endl;
	outputfile << endl;
	cout << "Breakdown of anomaly categories (Numbers below categories indicate the size of the counter difference for each anomaly):" << endl;
	outputfile << "Breakdown of anomaly categories (Numbers below categories indicate the size of the counter difference for each anomaly):" << endl;
	//Print number of each type of anomaly as well as difference sizes to screen and to text file
	TString classstring(" violations of L0B>=L0A>=L1B>=L1A>=L2B>=L2A for classes:");
	PrintAnomaly(classflag, classstring);
	TString clusterstring(" violations of L0>=L1>=L2 for clusters:");
	PrintAnomaly(clusterflag, clusterstring);
	TString L0clusterstring(" violations of L0 for a cluster = all L0 FOs for that cluster:");
	PrintAnomaly(clusterFOL0flag, L0clusterstring);
	TString L1clusterstring(" violations of L1 for a cluster = all L1 FOs for that cluster:");
	PrintAnomaly(clusterFOL1flag, L1clusterstring);
	TString L0strobestring(" violations of L0strobe0 = L0strobeIN:");
	PrintAnomaly(L0strobeflag, L0strobestring);
	TString L1strobestring(" violations of L1strobeOUT = L1strobeIN:");
	PrintAnomaly(L1strobeflag, L1strobestring);
	TString FOL1strobestring(" violations of L1strobeOUT = L1strIN for all FOs:");
	PrintAnomaly(FOL1strobeflag, FOL1strobestring);
	TString FOL2strobestring(" violations of L2strobeOUT = L2strIN for all FOs:");
	PrintAnomaly(FOL2strobeflag, FOL2strobestring);

	cout << endl;
	outputfile << endl;
 }

 cout << endl;
 //Print total numbers of glitches or spurious counts in each cluster for the whole data file

 cout << "Total glitches in cluster:   T      1      2      3      4      5      6" << endl;
 cout << "                         ";
 for(int j=0; j<numclusters; j++){
	cout << setw(7) << glitch[j];
 }
 cout << endl;
 cout << "Total spurious in cluster:   T      1      2      3      4      5      6" << endl;
 cout << "                        ";
 for(int k=0; k<numclusters; k++){
	cout << setw(7) << spurious[k];
 }
 cout << endl;
 outputfile << endl;
 outputfile << "Total glitches in cluster:   T      1      2      3      4      5      6" << endl;
 outputfile << "                         ";
 for(int j=0; j<numclusters; j++){
	outputfile << setw(7) << glitch[j];
 }
 outputfile << endl;
 outputfile << "Total spurious in cluster:   T      1      2      3      4      5      6" << endl;
 outputfile << "                        ";
 for(int k=0; k<numclusters; k++){
	outputfile << setw(7) << spurious[k];
 }
 outputfile << endl;


 outputfile.close();
}

