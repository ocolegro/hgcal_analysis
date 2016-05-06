/*
 * anal.C
 *
 *  Created on: Feb 24, 2016
 *      Author: ocolegro
 */



//---------------------------------------------------------------------------------------------------------------------------------
//
// Macro to read an output tree produced by running the EventInfoFiller on MiniAOD and make a couple of simple plots as an example.
// To run from the command line: root -l -q -b dummyPlotExample.C+\(\"myfilename.root\"\)
//
//---------------------------------------------------------------------------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <assert.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "string.h"
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
using namespace std;

#endif
bool exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

void root_to_txt()
{
  TString prefix     = "../data/root_files/unmerged/DigiIC3_Si2__version30_model0_BOFF_et";
  //TString prefix     = "../data/root_files/electron_2/DigiIC3_Si2__version30_model0_BOFF_et";

//  TString process    = "gamma_";
  TString process    = "electron_single_";

  TString energies[] = {"30"};
  //TString energies[] = {"4000"};

  std::vector<TString> runs;
  for (unsigned int k = 0; k < 1000; k++){

    TString s = std::to_string(k);
    runs.push_back(s);
    std::cout << runs[k] << std::endl;
    }
    std::cout << "The size of runs is " << runs.size();
  int n_events       = 2500;

  for (unsigned int i = 0; i < sizeof(energies)/sizeof(energies[0]); i++)
  {
    for (unsigned int j = 0; j < runs.size(); j++)
    {
        TString energy  = energies[i];
        TString run     = runs[j];
        std::cout <<"Attemping run " << run << std::endl;
        TString in_file = prefix +  energy + "_eta10000.000" + "_run" + run + ".root";
        TString out_file = process + energy + +"_" + run + "_MEV.txt";
        if (exists(string(in_file)) == 0) {
		cerr << "The file " << in_file << " was missing" << std::endl;
		continue;
	}
	    try{
        freopen( "../data/unmerged_files/" + out_file, "w", stdout );
        //freopen( "../data/unmerged_txts/" + out_file, "w", stdout );

        cerr << "The energy is "    << energy << "MEV" << std::endl;
        cerr << "The process is "   << process << std::endl;
        cerr << "The iteration is " << (i+1) * (j+1) << std::endl;
        cerr << "The infile is "    << in_file << std::endl;
        cerr << "The outfile  "     << "data/txt_files/" + out_file << std::endl;


        // initialize
        TFile *infile = TFile::Open(in_file);
        if( !infile ) throw 1;

        //assert(infile);

        //infile->Print();
        TTree *intree = (TTree*)infile->Get("RecoTree");

        //assert(intree);
        if( !intree ) throw 1;
        for (int i = 0; i < n_events; i++)
        {
            intree->Show(i,10000);
        }
        infile->Close();
        }
        catch(const int a) {
        std::cout << "THe file failed" << std::endl;
        }
    }
  }
}

