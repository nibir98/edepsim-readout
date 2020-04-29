#include <TFile.h>
#include "TSystem.h"
#include <TGeoManager.h>
#include <iostream>
#include <TRandom3.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "TString.h"

#include "TChain.h"
#include "TGeoTube.h"
#include <TTree.h>

#include "TG4Event.h"

#include <cstring>
#include <fstream>
#include <sstream>
#include <utility>
#include <functional>
#include <cassert>
#include <map>
#include <cstdlib>
#include <assert.h>

using namespace std;

bool inFV(float x, float y, float z)
{
  //  float centerX=0.;
  float centerY = -2384.73; // mm
  float centerZ = 23910;    // mm

  if (abs(x) > 1490)
    return false;
  float r = sqrt((y - centerY) * (y - centerY) + (z - centerZ) * (z - centerZ));
  if (r > 1800)
    return false;
  return true;
}

int main()
{

  TChain *chain = new TChain("EDepSimEvents");
  TFile *infile = new TFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_000_gsim1_0_1e18.4321.edep.root");

  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_000_gsim1_0_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_002_gsim3_2_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_004_gsim5_4_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_005_gsim6_5_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_008_gsim9_8_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_009_gsim10_9_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_010_gsim11_10_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_011_gsim12_11_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_013_gsim14_13_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_014_gsim15_14_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_016_gsim17_16_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_017_gsim18_17_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_022_gsim23_22_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_023_gsim24_23_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_024_gsim25_24_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_027_gsim28_27_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_028_gsim29_28_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_029_gsim30_29_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_030_gsim31_30_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_032_gsim33_32_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_033_gsim34_33_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_034_gsim35_34_1e18.4321.edep.root");
  chain->AddFile("/home/nibir/PHYSICS/duneanalysis/input/STT_RHC_035_gsim36_35_1e18.4321.edep.root");



  TGeoManager *geo = (TGeoManager *)infile->Get("EDepSimGeometry");
  TG4Event *event = NULL;
  chain->SetBranchAddress("Event", &event);
  TFile *outf = new TFile("outf.root", "recreate");

  int y_hits;                                 float erg_hit;
  int hits;                                   float trk_l_hit;
  float Erg;                                  int pdg_hit;
  float trk_l;                                TLorentzVector mid_hit; 
  int p_id;                                   TLorentzVector p_vect;
  TString particle_name;                       
  float p_init;
  int pdg;
  int event_id=0;

  TTree *tree = new TTree("tracks_info_tree", "trackwise_data");
  tree->Branch("event_id", &event_id);
  tree->Branch("y_hits", &y_hits);
  tree->Branch("hits", &hits);
  tree->Branch("Erg", &Erg);
  tree->Branch("trk_l", &trk_l);
  tree->Branch("p_id", &p_id);
  tree->Branch("particle_name", &particle_name);
  tree->Branch("p_init", &p_init);
  tree->Branch("p_vect", &p_vect);
  tree->Branch("pdg", &pdg);

  TTree *tree1 = new TTree("hits_info_tree", "hitwise_data");
  tree1->Branch("erg_hit", &erg_hit);
  tree1->Branch("trk_l_hit", &trk_l_hit);
  tree1->Branch("pdg_hit", &pdg_hit);
  tree1->Branch("mid_hit", &mid_hit);



  int nEntry = chain->GetEntries();
  cout << nEntry << endl;

  int ev = 0;

  for (auto i = 0; i < nEntry; i++)
  {
    chain->GetEntry(i);
    TLorentzVector vtx = event->Primaries.begin()->GetPosition();
    if (!inFV(vtx.X(), vtx.Y(), vtx.Z()))
      continue;

    ev++;

    for (int trk = 0; trk < event->Trajectories.size(); trk++)
    {
      p_vect = event->Trajectories.at(trk).InitialMomentum;
      p_init = sqrt(p_vect.X() * p_vect.X() + p_vect.Y() * p_vect.Y() + p_vect.Z() * p_vect.Z());
      particle_name = event->Trajectories.at(trk).Name;
      p_id = event->Trajectories.at(trk).ParentId;
      int trkid = event->Trajectories.at(trk).TrackId;
      pdg = event->Trajectories.at(trk).PDGCode;
      Erg = 0.0;
      trk_l = 0.0;
      hits = 0;
      y_hits = 0;
      event_id= ev;

      for (int h = 0; h < event->SegmentDetectors["Straw"].size(); h++)
      {
        if (event->SegmentDetectors["Straw"].at(h).PrimaryId != trkid)
          continue;
        int dy = 0;
        int dh = 1;

        float de = event->SegmentDetectors["Straw"].at(h).EnergyDeposit;
        float dl = event->SegmentDetectors["Straw"].at(h).TrackLength;
        TLorentzVector mid = (event->SegmentDetectors["Straw"].at(h).Start + event->SegmentDetectors["Straw"].at(h).Stop) * 0.5;
        TString name = geo->FindNode(mid.X(), mid.Y(), mid.Z())->GetName();
        if (name.Contains("horizontal"))
        {
          dy = 1;
        }

        Erg += de;
        trk_l += dl;
        hits += dh;
        y_hits += dy;
      }
      
      tree->Fill();
    }
  
  for (int h = 0; h < event->SegmentDetectors["Straw"].size(); h++){

  erg_hit = event->SegmentDetectors["Straw"].at(h).EnergyDeposit;
  trk_l_hit = event->SegmentDetectors["Straw"].at(h).TrackLength;
  mid_hit = (event->SegmentDetectors["Straw"].at(h).Start + event->SegmentDetectors["Straw"].at(h).Stop) * 0.5;
  int t = event->SegmentDetectors["Straw"].at(h).PrimaryId;
  pdg_hit = event->Trajectories.at(t).PDGCode;

  tree1->Fill();

  }

 }

  outf->Write();
  outf->Close();
  delete geo;
  delete chain;
  delete event;
  delete infile;
}
