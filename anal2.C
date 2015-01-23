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

//USER DEFINED GLOBAL VARIABLES:

//---------------2013 numbers---------------------
//const UInt_t numclasses=50; //Number of classes
//const UInt_t numclusters=7; //Number of clusters (including T)
//const UInt_t firstrun=896; //Counter number of first counter containing run number
//const UInt_t numcounters=970; //Number of counters (i.e. number of entries on each line of input file)

//---------------2014 numbers---------------------
const UInt_t numclasses=100; //Number of classes
const UInt_t numclusters=9; //Number of clusters (including T)
const UInt_t firstrun=1486; //Counter number of first counter containing run number
const UInt_t numcounters=1560; //Number of counters (i.e. number of entries on each line of input file)

//Variables to hold positions of needed counters in cnames.sorted file
UInt_t l0classB1;
UInt_t l0classA1;
UInt_t l1classB1;
UInt_t l1classA1;
UInt_t l2classB1;
UInt_t l2classA1;
UInt_t l0clstT;
UInt_t l1clstT;
UInt_t l2clstT;
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

int flag = 0; //Holds total number of anomalies found
//Vectors to hold the amounts by which instances are anomalous, for each category of anomaly
vector<double> classflag;
vector<double> clusterflag;
vector<double> clusterFOL0flag;
vector<double> clusterFOL1flag;
vector<double> L0strobeflag;
vector<double> L1strobeflag;
vector<double> FOL1strobeflag;
vector<double> FOL2strobeflag;
int glitch[numclusters] = {0}; //Counts the total number of glitches which occurred in each cluster
int spurious[numclusters] = {0}; //Counts the total number of spurious which occurred in each cluster
vector<int> locations; //Holds the run numbers in which anomalies occur
double currentrun[6] = {-1,-1,-1,-1,-1,-1}; //Holds run number so that it can be seen when run begins/ends. Initialised to -1 so begins from start
double total[6][numcounters] = {0}; //Array to hold the summed overflow-corrected counter differences between readings for each counter
int num = 0; //Line number
int zeroflag[6][numcounters] = {0}; //If SOR/EOR or any between has unexplained zero scalers, any 'anomalies' found which use these scalers will be ignored
bool firstflag=0; //When firstflag=1, any anomalies are accompanied by a message warning that they were not started from a SOR/EOR reading
bool lastflag=0; //When lastflag=1, this is the final line in the file and anomalies should be read out, with a warning
ofstream outputfile; 

void Plot(UInt_t *cnts, UInt_t *prev) //Function which sums scalers, and at EOR/SOR/end of data looks for anomalies
{
     cout.precision(12);
     outputfile.precision(12);

  for(int r=0; r<6; r++){ //Loop over all parallel runs
     int runx=firstrun+r;

     if((cnts[runx] == currentrun[r])||(num==2))  // If start of data or continuing run/no-run, sum scalers
     {
      currentrun[r]=cnts[runx];
      if(num==2) firstflag=1; //When beginning from start of file, comparisons between different counters will not be exactly valid  
	double increm[numcounters];
	double temp[numcounters];
    	for(int k=0; k<(numcounters); k++){
		//If counter is lower than in previous increment, it has overflowed
		//and so is corrected by adding 2^32 to current counter
		temp[k] = cnts[k]; //Don't change the value of cnts itself
		if((temp[k]<prev[k])&&(temp[k]!=0)){ 
			temp[k] = pow(2,32) + temp[k];
		}
		increm[k] = temp[k] - prev[k]; //Get counters since last increment
		total[r][k] += increm[k]; //Add this to total
		//If scaler difference from prev is unexpectedly zero, print and send to outputfile, set zeroflag
		if((temp[k]==0)&&(prev[k]!=0)&&(k!=firstrun)&&(k!=firstrun+1)&&(k!=firstrun+2)&&(k!=firstrun+3)&&(k!=firstrun+4)&&(k!=firstrun+5)){
			cout << "In Run " << cnts[runx];
			cout << " (line number " << num << ")";
			cout << ", counter " << k << " changed by zero unexpectedly." << endl;
			outputfile << "Run " << cnts[runx] << " num " << num << " " << k << " = 0" << endl;
			zeroflag[r][k]=2;
		}
    	}
	     		
    }

    if(cnts[runx] != currentrun[r] || lastflag==1){ //If new run (or changed to 'no run'), look for anomalies in previous run/non-run
   	if(prev[runx]==currentrun[r]){
		double increm[numcounters];
		double temp[numcounters];
    		for(int k=0; k<(numcounters); k++){
			//If counter is lower than in previous increment, it has overflowed
			// and so is corrected by adding 2^32 to current counter
			temp[k] = cnts[k]; //Don't change the value of cnts itself
			if((temp[k]<prev[k])&&(temp[k]!=0)){ 
				temp[k] = pow(2,32) + temp[k];
			}
			increm[k] = temp[k] - prev[k]; //Get counters since last increment
			total[r][k] += increm[k]; //Add this to total
			//If scaler difference from prev is unexpectedly zero, print and send to outputfile, set zeroflag
			if((temp[k]==0)&&(prev[k]!=0)&&(k!=firstrun)&&(k!=firstrun+1)&&(k!=firstrun+2)&&(k!=firstrun+3)&&(k!=firstrun+4)&&(k!=firstrun+5)){ 						     
				cout << "In Run " << cnts[runx];
				cout << " (line number " << num << ")";
				cout << ", counter " << k << " changed by zero." << endl;
				outputfile << "Run " << cnts[runx] << " num " << num << " counter " << k << " = 0" << endl;
				zeroflag[r][k]++;
			}
    		}
	}
	currentrun[r] = cnts[runx]; //Reset currentrun[r] to this line's (first) run number
	int runprintflag=0;

	//---------------------------Look for anomalies between trigger levels in each class--------------------------------------
	
	int classprintflag=0; //So titles only printed once per anomaly type
	for(int lzerob=l0classB1; lzerob<(l0classB1+numclasses); lzerob++){ //Positions of L0B,L0A,L1B, etc. in the array for the 50 classes
		int lzeroa = l0classA1 + lzerob - l0classB1;
		int loneb = l1classB1 + lzerob - l0classB1;		
		int lonea = l1classA1 + lzerob - l0classB1;
		int ltwob = l2classB1 + lzerob - l0classB1;
		int ltwoa = l2classA1 + lzerob - l0classB1;
		if((total[r][lzerob]>=total[r][lzeroa]) && (total[r][lzeroa]>=total[r][loneb]) && (total[r][loneb]>=total[r][lonea]) && (total[r][lonea]>=total[r][ltwob]) && (total[r][ltwob]>=total[r][ltwoa])){ //If L0B>=L0A>=L1B, etc. condition is true, print "no anomaly"
			cout << endl;
			cout << endl;
		 	cout << "Run " << prev[runx] << endl;
			cout << "Increment number (overall): " << num << endl;	
			cout << endl;
			cout << "NO anomaly in class " << lzerob-l0classB1+1 << endl;
			cout << "L0classB" << lzerob-l0classB1+1 << '\t' << "L0classA" << lzerob-l0classB1+1 << '\t' << "L1classB" << lzerob-l0classB1+1 << '\t' << "L1classA" << lzerob-l0classB1+1 << '\t' << "L2classB" << lzerob-l0classB1+1 << '\t' << "L2classA" << lzerob-l0classB1+1 << endl;
			cout << total[r][lzerob] << '\t' << total[r][lzeroa] << '\t' << total[r][loneb] << '\t' << total[r][lonea] << '\t' << total[r][ltwob] << '\t' << total[r][ltwoa] << endl;	
		}else{ //If condition is false, print details of anomaly to screen and file
			if((zeroflag[r][lzerob]==0) && (zeroflag[r][lzeroa]==0) && (zeroflag[r][loneb]==0) && (zeroflag[r][lonea]==0) && (zeroflag[r][ltwob]==0) && (zeroflag[r][ltwoa]==0)){ //If these are not due to spurious zeroes
			flag += 1;
			double diff=0; // Calculate the first difference (anomalous amount), e.g. +2
			if(total[r][lzerob]<total[r][lzeroa]) diff= total[r][lzeroa]-total[r][lzerob];
			if(total[r][lzeroa]<total[r][loneb])  diff= total[r][loneb]-total[r][lzeroa];
			if(total[r][loneb]<total[r][lonea]) diff= total[r][lonea]-total[r][loneb];
			if(total[r][lonea]<total[r][ltwob]) diff= total[r][ltwob]-total[r][lonea];
			if(total[r][ltwob]<total[r][ltwoa]) diff= total[r][ltwoa]-total[r][ltwob];
			classflag.push_back(diff); // Put anomalous amount into vector for class anomalies
			locations.push_back(prev[runx]); // Vector storing run numbers containing anomalies
			if(runprintflag==0){
			 cout << endl;
			 cout << endl;
			 outputfile << endl;
			 outputfile << endl;
			 cout << "Run " << prev[runx] << endl;
			 outputfile << "Run " << prev[runx];
			 if(prev[runx]==0) outputfile << " increment number " << num << " (ending with SOR " << cnts[runx] << " - parallel " << r << ")"; //for gaps between runs
			 outputfile << endl;
			 if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		 	 if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
			 outputfile << endl;
			 runprintflag++;
			}
						cout << "!!! Anomaly in class " << lzerob-l0classB1+1 << " !!!" << endl;
			outputfile << "!!! Anomaly in class " << lzerob-l0classB1+1 << " !!!" << endl;
			if(classprintflag==0){
		 	 cout << "Increment number (overall): " << num << endl;
			 cout << "L0classB" << lzerob-l0classB1+1 << '\t' << "L0classA" << lzerob-l0classB1+1 << '\t' << "L1classB" << lzerob-l0classB1+1 << '\t' << "L1classA" << lzerob-l0classB1+1 << '\t' << "L2classB" << lzerob-l0classB1+1 << '\t' << "L2classA" << lzerob-l0classB1+1 << endl;
			}
			cout << total[r][lzerob] << '\t' << total[r][lzeroa] << '\t' << total[r][loneb] << '\t' << total[r][lonea] << '\t' << total[r][ltwob] << '\t' << total[r][ltwoa] << endl;
			if(classprintflag==0){
			 outputfile << setw(12) << "L0classB" << lzerob-l0classB1+1 << setw(12) << "L0classA" << lzerob-l0classB1+1 << setw(12) << "L1classB" << lzerob-l0classB1+1 << setw(12) << "L1classA" << lzerob-l0classB1+1 << setw(12) << "L2classB" << lzerob-l0classB1+1 << setw(12) << "L2classA" << lzerob-l0classB1+1 << endl;
			 classprintflag++;
			}
			outputfile << setw(12) << total[r][lzerob] << setw(12) << total[r][lzeroa] << setw(12) << total[r][loneb] << setw(12) << total[r][lonea] << setw(12) << total[r][ltwob] << setw(12) << total[r][ltwoa] << endl;
			}
		}
	}	

	//-----Look for anomalies between trigger levels for each cluster, and between fanouts for each trigger level for each cluster------
	
	for(int lzeroclst=l0clstT; lzeroclst<(l0clstT+numclusters); lzeroclst++){ //Positions of L0,L1,L2 in the array for the 7 clusters
		int loneclst = l1clstT + lzeroclst - l0clstT;
		int ltwoclst = l2clstT + lzeroclst - l0clstT;	
	
		if((total[r][lzeroclst]>=total[r][loneclst]) && (total[r][loneclst]>=total[r][ltwoclst])){ //If L0>=L1>=L2 condition is true, print "no anomaly"
			cout << endl;
			cout << endl;
			cout << "Run " << prev[runx] << endl;
			cout << "Increment number (overall): " << num << endl;		
			if(lzeroclst==l0clstT){
				cout << endl;
				cout << "NO anomaly in cluster T" << endl;
				cout << "L0clstT" << '\t' << "L1clstT" << '\t' << "L2clstT" << endl;
			}
			if(lzeroclst!=l0clstT){
				cout << endl;
				cout << "NO anomaly in cluster " << lzeroclst-l0clstT << endl;
				cout << "L0clst" << lzeroclst-l0clstT << '\t' << "L1clst" << lzeroclst-l0clstT << '\t' << "L2clst" << lzeroclst-l0clstT<< endl;
			}
			cout << total[r][lzeroclst] << '\t' << total[r][loneclst] << '\t' << total[r][ltwoclst] << endl;	
		}else{ //If condition is false, print details of anomaly to screen and file
			if((zeroflag[r][lzeroclst]==0)&&(zeroflag[r][loneclst]==0)&&(zeroflag[r][ltwoclst]==0)){ //If these are not due to spurious zeroes
			flag += 1;
			double diff=0; // Calculate the first difference (anomalous amount), e.g. +2
			if(total[r][lzeroclst]<total[r][loneclst]) diff= total[r][loneclst]-total[r][lzeroclst];
			if(total[r][loneclst]<total[r][ltwoclst])  diff= total[r][ltwoclst]-total[r][loneclst];
			clusterflag.push_back(diff); // Put anomalous amount into vector for class anomalies
			locations.push_back(prev[runx]); // Vector storing run numbers containing anomalies
			if(runprintflag==0){
			 cout << endl;
			 cout << endl;
		 	 outputfile << endl;
			 outputfile << endl;
			 cout << "Run " << prev[runx] << endl;
			 outputfile << "Run " << prev[runx];
			 cout << "Increment number (overall): " << num << endl;	
			 if(prev[runx]==0) outputfile << " increment number " << num << " (ending with SOR " << cnts[runx] << " - parallel " << r << ")";
			 outputfile << endl;
		 	 if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR readings and so comparisons between counters may not be exactly valid)" << endl;
			 if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
			 outputfile << endl;
			 runprintflag++;
			}	
			if(lzeroclst==l0clstT){
				cout << endl;
				cout << "!!! Anomaly in cluster T !!!" << endl;
				cout << "L0clstT" << '\t' << "L1clstT" << '\t' << "L2clstT" << endl;
				outputfile << endl;
				outputfile << "!!! Anomaly in cluster T !!!" << endl;
				outputfile << setw(12) << "L0clstT" << setw(12) << "L1clstT" << setw(12) << "L2clstT" << endl;
			}
			if(lzeroclst!=l0clstT){
				cout << endl;
				cout << "!!! Anomaly in cluster " << lzeroclst-l0clstT << " !!!"<< endl;
				cout << "L0clst" << lzeroclst-l0clstT << '\t' << "L1clst" << lzeroclst-l0clstT << '\t' << "L2clst" << lzeroclst-l0clstT<< endl;
				outputfile << endl;
				outputfile << "!!! Anomaly in cluster " << lzeroclst-l0clstT << " !!!"<< endl;
				outputfile << setw(12) << "L0clst" << lzeroclst-l0clstT << setw(12) << "L1clst" << lzeroclst-l0clstT << setw(12) << "L2clst" << lzeroclst-l0clstT<< endl;
			}
			cout << total[r][lzeroclst] << '\t' << total[r][loneclst] << '\t' << total[r][ltwoclst] << endl;	
			outputfile << setw(12) << total[r][lzeroclst] << setw(12) << total[r][loneclst] << setw(12) << total[r][ltwoclst] << endl;	
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
		if((total[r][lzeroclst]==total[r][foonelzeroclst]) && (total[r][lzeroclst]==total[r][fotwolzeroclst]) && (total[r][lzeroclst]==total[r][fothreelzeroclst]) && (total[r][lzeroclst]==total[r][fofourlzeroclst]) && (total[r][lzeroclst]==total[r][fofivelzeroclst]) && (total[r][lzeroclst]==total[r][fosixlzeroclst])){ //Check L0 clusters equal L0 cluster FOs
			if(lzeroclst==l0clstT) cout << "No anomaly between cluster T and FOs for L0" << endl;
			if(lzeroclst!=l0clstT) cout << "No anomaly between cluster " << lzeroclst-l0clstT << " and FOs for L0" << endl;
		}else{
			if((zeroflag[r][lzeroclst]==0)&&(zeroflag[r][foonelzeroclst]==0)&&(zeroflag[r][fotwolzeroclst]==0)&&(zeroflag[r][fothreelzeroclst]==0)&&(zeroflag[r][fofourlzeroclst]==0)&&(zeroflag[r][fofivelzeroclst]==0)&&(zeroflag[r][fosixlzeroclst]==0)){ //If these are not due to spurious zeroes
			flag+=1;
			double diff=0; // Calculate the first difference (anomalous amount), e.g. +2
			if(diff==0) diff= total[r][foonelzeroclst]-total[r][lzeroclst];
			if(diff==0)  diff= total[r][fotwolzeroclst]-total[r][lzeroclst];
			if(diff==0) diff= total[r][fothreelzeroclst]-total[r][lzeroclst];
			if(diff==0)  diff= total[r][fofourlzeroclst]-total[r][lzeroclst];
			if(diff==0) diff= total[r][fofivelzeroclst]-total[r][lzeroclst];
			if(diff==0)  diff= total[r][fosixlzeroclst]-total[r][lzeroclst];
			clusterFOL0flag.push_back(diff); // Put anomalous amount into vector for FO L0 cluster anomalies
			locations.push_back(prev[runx]); // Vector storing run numbers containing anomalies
			if(runprintflag==0){
			 outputfile << endl;
			 outputfile << endl;
			 outputfile << "Run " << prev[runx];
			 if(prev[runx]==0) outputfile << " increment number " << num << "(ending with SOR " << cnts[runx] << " - parallel " << r << ")";
			 outputfile << endl;
			}
			outputfile << endl;
			if(lzeroclst==l0clstT){
				cout << "!!! Anomaly between cluster T and FOs for L0!!!" << endl;
				cout << "L0clstT" << '\t' << "FO1L0clstT" << '\t' << "FO2L0clstT" << '\t' << "FO3L0clstT" << '\t' << "FO4L0clstT" << '\t' << "FO5L0clstT" << '\t' << "FO6L0clstT" << endl;
				outputfile << "!!! Anomaly between cluster T and FOs for L0!!!" << endl;
				if(runprintflag==0){
				 if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR readings and so comparisons between counters may not be exactly valid)" << endl;
				 if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
				 runprintflag++;
				}
				outputfile << setw(12) << "L0clstT" << setw(12) <<  "FO1L0clstT" << setw(12) << "FO2L0clstT" << setw(12) <<  "FO3L0clstT" << setw(12) << "FO4L0clstT" << setw(12) <<  "FO5L0clstT" << setw(12) <<  "FO6L0clstT" << endl;

			}
			if(lzeroclst!=l0clstT){
				cout << "!!! Anomaly between cluster " << lzeroclst-l0clstT << " and FOs for L0 !!!" << endl;
				cout << "L0clst" << lzeroclst-l0clstT << '\t' << "F01L0clst" << lzeroclst-l0clstT << '\t' << "F02L0clst" << lzeroclst-l0clstT << '\t' << "FO3L0clst" << lzeroclst-l0clstT << '\t' << "FO4L0clst" << lzeroclst-l0clstT << '\t' << "FO5L0clst" << lzeroclst-l0clstT << '\t' << "FO6L0clst" << lzeroclst-l0clstT << endl;
				outputfile << "!!! Anomaly between cluster " << lzeroclst-l0clstT << " and FOs for L0 !!!" << endl;
				if(runprintflag==0){
			 	 if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR readings and so comparisons between counters may not be exactly valid)" << endl;
				 if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
				 runprintflag++;
				}
				outputfile << setw(12) << "L0clst" << lzeroclst-l0clstT << setw(12) << "F01L0clst" << lzeroclst-l0clstT << setw(12) << "F02L0clst" << lzeroclst-l0clstT << setw(12) << "FO3L0clst" << lzeroclst-l0clstT << setw(12) << "FO4L0clst" << lzeroclst-l0clstT << setw(12) << "FO5L0clst" << lzeroclst-l0clstT << setw(12) << "FO6L0clst" << lzeroclst-l0clstT << endl;

			}
			cout << total[r][lzeroclst] << '\t' << total[r][foonelzeroclst] << '\t' << total[r][fotwolzeroclst] << '\t' << total[r][fothreelzeroclst] << '\t' << total[r][fofourlzeroclst] << '\t' << total[r][fofivelzeroclst] << '\t' << total[r][fosixlzeroclst] << endl;
			outputfile << setw(12) << total[r][lzeroclst] << setw(12) << total[r][foonelzeroclst] << setw(12) << total[r][fotwolzeroclst] << setw(12) << total[r][fothreelzeroclst] << setw(12) << total[r][fofourlzeroclst] << setw(12) << total[r][fofivelzeroclst] << setw(12) << total[r][fosixlzeroclst] << endl;
			}
		}
		if((total[r][loneclst]==total[r][fooneloneclst]) && (total[r][loneclst]==total[r][fotwoloneclst]) && (total[r][loneclst]==total[r][fothreeloneclst]) && (total[r][loneclst]==total[r][fofourloneclst]) && (total[r][loneclst]==total[r][fofiveloneclst]) && (total[r][loneclst]==total[r][fosixloneclst])){ //Check L1 clusters equal L1 FOs
			if(lzeroclst==l0clstT) cout << "No anomaly between cluster T and FOs for L1" << endl;
			if(lzeroclst!=l0clstT) cout << "No anomaly between cluster " << lzeroclst-l0clstT << " and FOs for L1" << endl;
		}else{ //If condition is false, print details to screen and to file
			if((zeroflag[r][loneclst]==0)&&(zeroflag[r][fooneloneclst]==0)&&(zeroflag[r][fotwoloneclst]==0)&&(zeroflag[r][fothreeloneclst]==0)&&(zeroflag[r][fofourloneclst]==0)&&(zeroflag[r][fofiveloneclst]==0)&&(zeroflag[r][fosixloneclst]==0)){ //If these are not due to spurious zeroes
			flag+=1;
			double diff=0; // Calculate the first difference (anomalous amount), e.g. +2
			if(diff==0) diff= total[r][fooneloneclst]-total[r][loneclst];
			if(diff==0)  diff= total[r][fotwoloneclst]-total[r][loneclst];
			if(diff==0) diff= total[r][fothreeloneclst]-total[r][loneclst];
			if(diff==0)  diff= total[r][fofourloneclst]-total[r][loneclst];
			if(diff==0) diff= total[r][fofiveloneclst]-total[r][loneclst];
			if(diff==0)  diff= total[r][fosixloneclst]-total[r][loneclst];
			clusterFOL1flag.push_back(diff); // Put anomalous amount into vector for FO L1 cluster anomalies
			locations.push_back(prev[runx]);
			if(runprintflag==0){
			 outputfile << endl;
			 outputfile << endl;
			 outputfile << "Run " << prev[runx];
		 	 if(prev[runx]==0) outputfile << " increment number " << num << "(ending with SOR " << cnts[runx] << " - parallel " << r << ")";
			 outputfile << endl;
			 if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
			 if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
			 runprintflag++;
			 outputfile << endl;
			}
			outputfile << endl;
			if(lzeroclst==l0clstT){
				cout << "!!! Anomaly between cluster T and FOs for L1!!!" << endl;
				cout << "L1clstT" << '\t' << "FO1L1clstT" << '\t' << "FO2L1clstT" << '\t' << "FO3L1clstT" << '\t' << "FO4L1clstT" << '\t' << "FO5L1clstT" << '\t' << "FO6L1clstT" << endl;
				outputfile << "!!! Anomaly between cluster T and FOs for L1!!!" << endl;
				outputfile << setw(12) << "L1clstT" << setw(12) << "FO1L1clstT" << setw(12) << "FO2L1clstT" << setw(12) << "FO3L1clstT" << setw(12) << "FO4L1clstT" << setw(12) << "FO5L1clstT" << setw(12) << "FO6L1clstT" << endl;

			}
			if(lzeroclst!=l0clstT){
				cout << "!!! Anomaly between cluster " << lzeroclst-l0clstT << " and FOs for L1 !!!" << endl;
				cout << "L1clst" << lzeroclst-l0clstT << '\t' << "F01L1clst" << lzeroclst-l0clstT << '\t' << "F02L1clst" << lzeroclst-l0clstT << '\t' << "FO3L1clst" << lzeroclst-l0clstT << '\t' << "FO4L1clst" << lzeroclst-l0clstT << '\t' << "FO5L1clst" << lzeroclst-l0clstT << '\t' << "FO6L1clst" << lzeroclst-l0clstT << endl;
				outputfile << "!!! Anomaly between cluster " << lzeroclst-l0clstT << " and FOs for L1 !!!" << endl;
				outputfile << setw(12) << "L1clst" << lzeroclst-l0clstT << setw(12) << "F01L1clst" << lzeroclst-l0clstT << setw(12) << "F02L1clst" << lzeroclst-l0clstT << setw(12) << "FO3L1clst" << lzeroclst-l0clstT << setw(12) << "FO4L1clst" << lzeroclst-l0clstT << setw(12) << "FO5L1clst" << lzeroclst-l0clstT << setw(12) << "FO6L1clst" << lzeroclst-l0clstT << endl;
			}
			cout << total[r][loneclst] << '\t' << total[r][fooneloneclst] << '\t' << total[r][fotwoloneclst] << '\t' << total[r][fothreeloneclst] << '\t' << total[r][fofourloneclst] << '\t' << total[r][fofiveloneclst] << '\t' << total[r][fosixloneclst] << endl;
			outputfile << setw(12) << total[r][loneclst] << setw(12) << total[r][fooneloneclst] << setw(12) << total[r][fotwoloneclst] << setw(12) << total[r][fothreeloneclst] << setw(12) << total[r][fofourloneclst] << setw(12) << total[r][fofiveloneclst] << setw(12) << total[r][fosixloneclst] << endl;
			}
		}

	}

	//-------------Look for anomalies between strobes in and out for each trigger level---------------

	if(total[r][l0strobe0]==total[r][l0strobeIN]){ //Check L0strobe0 = L0strobeIN
		cout << endl;
		cout << "No anomaly between L0strobe0 and L0strobeIN" << endl;
	}else{
		if((zeroflag[r][l0strobe0]==0)&&(zeroflag[r][l0strobeIN]==0)){ //If these are not due to spurious zeroes
		flag+=1;
		double diff=0; // Calculate the difference (anomalous amount), e.g. +2
		diff = total[r][l0strobeIN]-total[r][l0strobe0];
		L0strobeflag.push_back(diff); // Put anomalous amount into vector for L0 strobe anomalies
		locations.push_back(prev[runx]);
		if(runprintflag==0){
		 cout << endl;
	 	 outputfile << endl;
		 outputfile << endl;
		 outputfile << "Run " << prev[runx];
		 if(prev[runx]==0) outputfile << " increment number " << num << "(ending with SOR " << cnts[runx] << " - parallel " << r << ")";
		 outputfile << endl;
		 if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		 if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		 outputfile << endl;
		 runprintflag++;
		}
		outputfile << endl;
		cout << "!!! Anomaly between L0strobe0 and L0strobeIN" << endl;
		cout << "L0strobe0" << '\t' << "L0strobeIN" << endl;
		cout << total[r][l0strobe0] << '\t' << total[r][l0strobeIN] << endl;
		outputfile << "!!! Anomaly between L0strobe0 and L0strobeIN" << endl;
		outputfile << setw(12) << "L0strobe0" << setw(12) << "L0strobeIN" << endl;
		outputfile << setw(12) << total[r][l0strobe0] << setw(12) << total[r][l0strobeIN] << endl;
		}
	}

	if(total[r][l1strobeOUT]==total[r][l1strobeIN]){ //Check L1strobeOUT = L1strobeIN
		cout << endl;
		cout << "No anomaly between L1strobeOUT and L1strobeIN" << endl;
	}else{
		if((zeroflag[r][l1strobeOUT]==0)&&(zeroflag[r][l1strobeIN]==0)){ //If these are not due to spurious zeroes
		flag+=1;
		double diff=0; // Calculate the difference (anomalous amount), e.g. +2
		diff = total[r][l1strobeIN]-total[r][l1strobeOUT];
		L1strobeflag.push_back(diff); // Put anomalous amount into vector for L1 strobe anomalies
		locations.push_back(prev[runx]);
		if(runprintflag==0){
		 cout << endl;
		 outputfile << endl;
		 outputfile << endl;
		 outputfile << "Run " << prev[runx];
		 if(prev[runx]==0) outputfile << " increment number " << num << "(ending with SOR " << cnts[runx] << " - parallel " << r << ")";
		 outputfile << endl;
		 if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		 if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		 outputfile << endl;
		 runprintflag++;
		}
		outputfile << endl;
		cout << "!!! Anomaly between L1strobeOUT and L1strobeIN" << endl;
		cout << "L1strobeOUT" << '\t' << "L1strobeIN" << endl;
		cout << total[r][l1strobeOUT] << '\t' << total[r][l1strobeIN] << endl;
		outputfile << "!!! Anomaly between L1strobeOUT and L1strobeIN" << endl;
		outputfile << setw(12) << "L1strobeOUT" << setw(12) << "L1strobeIN" << endl;
		outputfile << setw(12) << total[r][l1strobeOUT] << setw(12) << total[r][l1strobeIN] << endl;
		}
	}

	//----------Look for anomalies between strobes out and FO strobes out for each trigger level----------------------

	if((total[r][l1strobeOUT]==total[r][fo1l1strIN]) && (total[r][l1strobeOUT]==total[r][fo2l1strIN]) && (total[r][l1strobeOUT]==total[r][fo3l1strIN]) && (total[r][l1strobeOUT]==total[r][fo4l1strIN]) && (total[r][l1strobeOUT]==total[r][fo5l1strIN]) && (total[r][l1strobeOUT]==total[r][fo6l1strIN])){ //Check L1strobeOUT = all L1 FO strobes IN
		cout << endl;
		cout << "No anomaly between L1strobeOUT and L1 FO strobes IN" << endl;
	}else{
		if((zeroflag[r][l1strobeOUT]==0)&&(zeroflag[r][fo1l1strIN]==0)&&(zeroflag[r][fo2l1strIN]==0)&&(zeroflag[r][fo3l1strIN]==0)&&(zeroflag[r][fo4l1strIN]==0)&&(zeroflag[r][fo5l1strIN]==0)&&(zeroflag[r][fo6l1strIN]==0)){ //If these are not due to spurious zeroes
		flag+=1;
		double diff=0; // Calculate the first difference (anomalous amount), e.g. +2
		if(diff==0) diff= total[r][fo1l1strIN]-total[r][l1strobeOUT];
		if(diff==0)  diff= total[r][fo2l1strIN]-total[r][l1strobeOUT];
		if(diff==0) diff= total[r][fo3l1strIN]-total[r][l1strobeOUT];
		if(diff==0)  diff= total[r][fo4l1strIN]-total[r][l1strobeOUT];
		if(diff==0) diff= total[r][fo5l1strIN]-total[r][l1strobeOUT];
		if(diff==0)  diff= total[r][fo6l1strIN]-total[r][l1strobeOUT];
		FOL1strobeflag.push_back(diff); // Put anomalous amount into vector for L1 FO strobe anomalies
		locations.push_back(prev[runx]);
		if(runprintflag==0){
		 cout << endl;
		 outputfile << endl;
		 outputfile << "Run " << prev[runx];
		 if(prev[runx]==0) outputfile << " increment number " << num << "(ending with SOR " << cnts[runx] << " - parallel " << r << ")";
		 outputfile << endl;
		 if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		 if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		 outputfile << endl;
		 runprintflag++;
		 }
		cout << "!!! Anomaly between L1strobeOUT and L1 FO strobes IN" << endl;
		cout << "L1strobeOUT" << '\t' << "fo1L1strIN" << '\t' << "fo2L1strIN" << '\t' << "fo3L1strIN" << '\t' << "fo4L1strIN" << '\t' << "fo5L1strIN" << '\t' << "fo6L1strIN" << endl;
		cout << total[r][l1strobeOUT] << '\t' << total[r][fo1l1strIN] << '\t' << total[r][fo2l1strIN] << '\t' << total[r][fo3l1strIN] << '\t' << total[r][fo4l1strIN] << '\t' << total[r][fo5l1strIN] << '\t' << total[r][fo6l1strIN] << endl;
		outputfile << "!!! Anomaly between L1strobeOUT and L1 FO strobes IN" << endl;
		outputfile << setw(12) << "L1strobeOUT" << setw(12) << "fo1L1strIN" << setw(12) << "fo2L1strIN" << setw(12) << "fo3L1strIN" << setw(12) << "fo4L1strIN" << setw(12) << "fo5L1strIN" << setw(12) << "fo6L1strIN" << endl;
		outputfile << setw(12) << total[r][l1strobeOUT] << setw(12) << total[r][fo1l1strIN] << setw(12) << total[r][fo2l1strIN] << setw(12) << total[r][fo3l1strIN] << setw(12) << total[r][fo4l1strIN] << setw(12) << total[r][fo5l1strIN] << setw(12) << total[r][fo6l1strIN] << endl;
		}
	}

	if((total[r][l2strobeOUT]==total[r][fo1l2strIN]) && (total[r][l2strobeOUT]==total[r][fo2l2strIN]) && (total[r][l2strobeOUT]==total[r][fo3l2strIN]) && (total[r][l2strobeOUT]==total[r][fo4l2strIN]) && (total[r][l2strobeOUT]==total[r][fo5l2strIN]) && (total[r][l2strobeOUT]==total[r][fo6l2strIN])){ //Check L2strobeOUT = all L2 FO strobes IN
		cout << endl;
		cout << "No anomaly between L2strobeOUT and L2 FO strobes IN" << endl;
	}else{
		if((zeroflag[r][l2strobeOUT]==0)&&(zeroflag[r][fo1l2strIN]==0)&&(zeroflag[r][fo2l2strIN]==0)&&(zeroflag[r][fo3l2strIN]==0)&&(zeroflag[r][fo4l2strIN]==0)&&(zeroflag[r][fo5l2strIN]==0)&&(zeroflag[r][fo6l2strIN]==0)){ //If these are not due to spurious zeroes
		flag+=1;
		double diff=0; // Calculate the first difference (anomalous amount), e.g. +2
		if(diff==0) diff= total[r][fo1l2strIN]-total[r][l2strobeOUT];
		if(diff==0)  diff= total[r][fo2l2strIN]-total[r][l2strobeOUT];
		if(diff==0) diff= total[r][fo3l2strIN]-total[r][l2strobeOUT];
		if(diff==0)  diff= total[r][fo4l2strIN]-total[r][l2strobeOUT];
		if(diff==0) diff= total[r][fo5l2strIN]-total[r][l2strobeOUT];
		if(diff==0)  diff= total[r][fo6l2strIN]-total[r][l2strobeOUT];
		FOL2strobeflag.push_back(diff); // Put anomalous amount into vector for L2 FO strobe anomalies
		locations.push_back(prev[runx]);
		if(runprintflag==0){
		 cout << endl;
		 outputfile << endl;
		 outputfile << "Run " << prev[runx];
		 if(prev[runx]==0) outputfile << " increment number " << num << "(ending with SOR " << cnts[runx] << " - parallel " << r << ")";
		 outputfile << endl;
		 if(firstflag==1) outputfile << "(Note that this period was from the beginning of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		 if(lastflag==1) outputfile << "(Note that this period stops at the end of the file rather than a SOR/EOR reading and so comparisons between counters may not be exactly valid)" << endl;
		 outputfile << endl;
		 runprintnum++;
		}
		cout << "!!! Anomaly between L2strobeOUT and L2 FO strobes IN !!!" << endl;
		cout << "L2strobeOUT" << '\t' << "fo1L2strIN" << '\t' << "fo2L2strIN" << '\t' << "fo3L2strIN" << '\t' << "fo4L2strIN" << '\t' << "fo5L2strIN" << '\t' << "fo6L2strIN" << endl;
		cout << total[r][l2strobeOUT] << '\t' << total[r][fo1l2strIN] << '\t' << total[r][fo2l2strIN] << '\t' << total[r][fo3l2strIN] << '\t' << total[r][fo4l2strIN] << '\t' << total[r][fo5l2strIN] << '\t' << total[r][fo6l2strIN] << endl;
		outputfile << "!!! Anomaly between L2strobeOUT and L2 FO strobes IN" << endl;
		outputfile << setw(12) << "L2strobeOUT" << setw(12) << "fo1L2strIN" << setw(12) << "fo2L2strIN" << setw(12) << "fo3L2strIN" << setw(12) << "fo4L2strIN" << setw(12) << "fo5L2strIN" << setw(12) << "fo6L2strIN" << endl;
		outputfile << setw(12) << total[r][l2strobeOUT] << setw(12) << total[r][fo1l2strIN] << setw(12) << total[r][fo2l2strIN] << setw(12) << total[r][fo3l2strIN] << setw(12) << total[r][fo4l2strIN] << setw(12) << total[r][fo5l2strIN] << setw(12) << total[r][fo6l2strIN] << endl;
		}
	}

	//------------------Count the total[r] glitches and spurious so far in each cluster----------------------
	
	cout << endl;
	for(int fonum =1; fonum<7; fonum++){ //Check for any glitches, FO boards 1-6
		int glitchT;
		if(fonum==1) glitchT = fo1glitchT;
		if(fonum==2) glitchT = fo2glitchT;
		if(fonum==3) glitchT = fo3glitchT;
		if(fonum==4) glitchT = fo4glitchT;
		if(fonum==5) glitchT = fo5glitchT;
		if(fonum==6) glitchT = fo6glitchT;
		for(int pos = glitchT; pos<glitchT+numclusters; pos++){ //Check for any glitches, clusters T-6
			if(total[r][pos]==0){
			//	cout << "No glitch in FO" << (fonum/48)+1 << "glitch";
				if(glitchT==fo1glitchT) {
				//	cout << "T" << endl;
				}else{ 
				//	cout << (foglitchnum-fo1glitchT) << endl;
				}
			}else{
				if(zeroflag[r][pos]==0){ //If these are not due to spurious zeroes
				cout << "Glitch in FO" << fonum << "glitch";
				glitch[pos-glitchT]+=total[r][pos]; //Add glitch to total[r] glitches for that cluster
				if((pos-glitchT)==0) {
					cout << "T:  ";
				}else{
					cout << (pos-glitchT) << ":  ";	
				}
				cout << total[r][pos] << endl;
				}
			}
		}
	}

	cout << endl;
	for(int fonum =1; fonum<7; fonum++){ //Check for any glitches, FO boards 1-6
		int spuriousT;
		if(fonum==1) spuriousT = fo1l1spuriousT;
		if(fonum==2) spuriousT = fo2l1spuriousT;
		if(fonum==3) spuriousT = fo3l1spuriousT;
		if(fonum==4) spuriousT = fo4l1spuriousT;
		if(fonum==5) spuriousT = fo5l1spuriousT;
		if(fonum==6) spuriousT = fo6l1spuriousT;
		for(int pos = spuriousT; pos<spuriousT+numclusters; pos++){ //Check for any spurious, clusters T-numclusters
			if(total[r][pos]==0){
				if(spuriousT==fo1l1spuriousT) {
				}else{ 
				}
			}else{
				if(zeroflag[r][pos]==0){ //If these are not due to spurious zeroes
				cout << "Spurious in FO" << fonum << "l1spurious";
				spurious[pos-spuriousT]+=total[r][pos]; //Add spurious to total[r] spurious for that cluster
				if((pos-spuriousT)==0) {
					cout << "T:  ";
				}else{
					cout << (pos-spuriousT) << ":  ";	
				}
				cout << total[r][pos] << endl;
				}
			}
		}
	}
	
	cout << endl;
	for(UInt_t j=0; j<numcounters; j++){ //Reset total[r]s for next run
		total[r][j]=0;	
		if(zeroflag[r][j]>1) zeroflag[r][j]=0; //After reading with zeroes has been accounted for EOR and SOR, reset to zero	
	}
    }

    runx++;
  
  }
  firstflag=0; //Reset firstflag to zero
  runx-=5;
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
	if(file->is_open()){
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
    		//cout << "Increment number: " << num << endl;
    		if(num%100==0) cout << "Reading line number " << num << "/" << numberoflines << " ... " << endl;
    		for(int i=0; i<(ntokens-2); i++){
    			prev[i] = cnts[i]; //Assign current counters to be prev counters for next while loop
   		}
  	  }
	} else cout << "Unable to open file: " << name.Data() << endl;
}



//main function
void anal2()
{
 cout << "Running program to identify possible anomalies in CTP counters." << endl;
 cout << endl;

 //INPUTS FROM USER REGARDING FILES:

 //-------------2013 data-----------
 //TString name("rawcnts/01.01.2013.rawcnts"); //name and path of first file to be analysed
 //int nfiles = 59; //number of files total to be analysed - current limit that only 2 different months can be spanned, must be same year
 //TString sortedname("cnames.sorted2"); //file containing line numbers and names of counters corresponding to data

 //-------------2014 or 2015 data-------------
 //TString name("raw112014/rawcnts/01.11.2014.rawcnts"); //name and path of first file to be analysed
 TString name("Jan2015/18.01.2015.rawcnts"); //name and path of first file to be analysed
 int nfiles = 1; //number of files total to be analysed - current limit that only 2 different months can be spanned, must be same year
 TString sortedname("cnames.sorted2.2014"); //file containing line numbers and names of counters corresponding to data


 TString namecopy = name; //copy of first file name string to be changed to subsequent file names
 TString dayname = name; //copy of first file name string to be cut to just the day of the month
 TString monthname = name; //copy of first file name string to be cut to just the month
 
 Ssiz_t nameLength= name.Length();
 int iday = nameLength-18; //Index of first digit of the day
 TString daystring = dayname.Remove(0,iday); //Cut before day
 daystring.Remove(2,nameLength-iday+2); //Cut after day
 Int_t firstday = atoi(daystring.Data()); //Convert day to integer
 cout << "Getting ready to read " << nfiles << " file(s) with name(s): " << endl;
 cout << name.Data() << endl;

 //Allows files to span two distinct months
 TString monthstring = monthname.Remove(0,iday+3); //Cut before month
 monthstring.Remove(5,nameLength-iday+5); //Cut after month
 Int_t firstmonth = atoi(monthstring.Data()); //Convert month to integer
 int monthLength=0;
 if(firstmonth==1 || firstmonth==3 || firstmonth==5 || firstmonth==7 || firstmonth==8 || firstmonth==10 || firstmonth==12) monthLength=31;
 else if(firstmonth==2) { //February
	TString leapyear1("2010"), leapyear2("2014"), leapyear3("2018")
	if(name.Contains(leapyear1) || name.Contains(leapyear2) || name.Contains(leapyear3)) monthLength=29; //leap years
	else monthLength=28; //non-leap years only
 }
 else if(firstmonth==4 || firstmonth==6 || firstmonth==9 || firstmonth==11) monthLength=30;
 else cout << "Filename appears to be incorrect format, as month has been read as not between 1 and 12" << endl; 

 //Loop over all nfiles filenames and open them
 vector<TString> filenames;
 filenames.push_back(name);
 int k=0;
 for(int i= firstday+1;i<(nfiles+firstday);i++){
	if(i>monthLength){ //If next month from first given
		k++;
		if(k==1){//Move to the next month
			stringstream monthsts;
			monthsts << firstmonth+1;
			TString monthtempstr = monthsts.str();
			if(monthtempstr.Length()<2) monthtempstr.Prepend("0");
			namecopy.Replace(iday+3,2,monthtempstr);
		}
	} else k=i;
	stringstream sts;
	sts << k;
	TString tempstr = sts.str();
	if(tempstr.Length()<2) tempstr.Prepend("0");
	namecopy.Replace(iday,2,tempstr);
 	filenames.push_back(namecopy);
	cout << namecopy.Data() << endl;
 }

 Int_t nlines=0;
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
 cout << endl;
 cout << "Creating and opening output file: " << outputname.Data() << endl;
 
 //Read in the needed counter positions from the cnames.sorted data file
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
 cout << endl;
 cout << "Counted " << numberoflines << " lines - likely to take around " << setprecision(2) << (1.25/100.)*numberoflines << " minutes to run program." << endl;
 cout << endl;

 for(int i=0;i<nfiles;i++){ //Sends each file to function which reads them line by line and performs analysis
 	ReadLines(filenames[i], numberoflines, prev);
 }
 
 //Summarise anomalies found, print and send to output file
 cout << endl;
 outputfile << endl;
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

 cout << "Total glitches in cluster:   T      1      2      3      4      5      6";
 if(numclusters==9) cout << "     7     8" << endl;
 cout << "                         ";
 for(int j=0; j<numclusters; j++){
	cout << setw(7) << glitch[j];
 }
 cout << endl;
 cout << "Total spurious in cluster:   T      1      2      3      4      5      6";
 if(numclusters==9) cout << "     7     8" << endl;
 cout << "                        ";
 for(int k=0; k<numclusters; k++){
	cout << setw(7) << spurious[k];
 }
 cout << endl;
 outputfile << endl;
 outputfile << "Total glitches in cluster:   T      1      2      3      4      5      6";
 if(numclusters==9) outputfile << "     7     8" << endl;
 outputfile << "                         ";
 for(int j=0; j<numclusters; j++){
	outputfile << setw(7) << glitch[j];
 }
 outputfile << endl;
 outputfile << "Total spurious in cluster:   T      1      2      3      4      5      6";
 if(numclusters==9) outputfile << "     7     8" << endl;
 outputfile << "                        ";
 for(int k=0; k<numclusters; k++){
	outputfile << setw(7) << spurious[k];
 }
 outputfile << endl;


 outputfile.close();
}

