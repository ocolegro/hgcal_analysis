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

TString hists[]     = {"h11","h12","h21","h22"};/*{"h10","h11","h13","h14",
                        "h20","h21","h23","h24"};*/

TString in_file     = "../data/shen1.root";


void hist_to_txt()
{

  for (unsigned int j = 0; j < sizeof(hists)/sizeof(hists[0]); j++)
  {
  TString hist = hists[j];
  TString out_file = "../data/hists/hist" + hist + ".txt";
  //std::cout << out_file << std::endl;
  freopen(  out_file, "w", stdout );
  TFile *infile = TFile::Open(in_file); assert(infile);

  TH1F * h1 = (TH1F*)infile->Get(hist);
  h1->Print("all");
  }
}
