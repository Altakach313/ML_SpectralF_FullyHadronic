/// \file ReadDelphesJets.C
/// \brief File that reads the delphes events and calculates the spectral function for each event and write it into a TTree

#include "TSystem.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <Eigen/Dense>

#include "delphes/classes/DelphesClasses.h"
#include "delphes/external/ExRootAnalysis/ExRootTreeReader.h"

#ifdef __CLING__
R__LOAD_LIBRARY(delphes/libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

double DeltaR = 0.4;
double Rmax = 8.4;

TVector3 getPxPyPz(Jet jet)
{
  double px = jet.PT * TMath::Cos(jet.Phi);
  double py = jet.PT * TMath::Sin(jet.Phi);
  double theta = 2 * TMath::ATan(TMath::Exp(-jet.Eta));
  double pz = jet.PT * 1. / TMath::Tan(theta);
  return TVector3(px, py, pz);
}

double calculateC(std::vector<Jet> jetvector)
{
  std::array<std::array<double, 3>, 3> Mxyz;
  for (unsigned int iX = 0; iX < 3; iX++)
  {
    for (unsigned int iY = 0; iY < 3; iY++)
    {
      double momTotal = 0.;
      for (auto jet : jetvector)
      {
        TVector3 mom(getPxPyPz(jet));
        Mxyz[iX][iY] += (mom[iX] * mom[iY]) / (mom.Mag());
        momTotal += mom.Mag();
      }
      Mxyz[iX][iY] /= momTotal;
    }
  }

  Eigen::Matrix3d matrix;
  // Populate the matrix with your values
  matrix << Mxyz[0][0], Mxyz[0][1], Mxyz[0][2],
      Mxyz[1][0], Mxyz[1][1], Mxyz[1][2],
      Mxyz[2][0], Mxyz[2][1], Mxyz[2][2];

  // Create an Eigen solver object
  Eigen::EigenSolver<Eigen::Matrix3d> eigensolver(matrix);
  if (eigensolver.info() != Eigen::Success)
  {
    std::cerr << "Eigen solver failed!" << std::endl;
    return 0.0;
  }

  Eigen::Vector3d eigenvalues = eigensolver.eigenvalues().real();
  std::array<double, 3> lambda = {eigenvalues[0], eigenvalues[1], eigenvalues[2]};

  double Cparam = 3 * (lambda[0] * lambda[1] + lambda[0] * lambda[2] + lambda[1] * lambda[2]);

  return Cparam;
}

double angularDistance(Jet jet1, Jet jet2)
{

  Double_t deta = jet1.Eta - jet2.Eta;
  Double_t dphi = jet1.Phi - jet2.Phi;
  dphi = TVector2::Phi_mpi_pi(dphi);
  return TMath::Sqrt(deta * deta + dphi * dphi);
}

bool isInRange(double dist, double R, double deltaR)
{

  bool isInRange = false;
  if (dist >= R && dist < R + deltaR)
  {
    isInRange = true;
  }
  return isInRange;
}

double calculateSpectralFunction(std::vector<Jet> jet, double RUpper)
{
  double sumS = 0.;

  for (auto iconst1 : jet)
  {
    for (auto iconst2 : jet)
    {
      sumS += iconst1.PT * iconst2.PT * double(isInRange(angularDistance(iconst1, iconst2), RUpper, DeltaR));
    }
  }
  return sumS / DeltaR;
}

void ReadDelphesJets(const char *inputFile = "rootFiles_Delphes/gogo_6j_m15TeV_rootFiles/root_gogo_6j_m15TeV_Events", std::string outputFile = "Test1", float maxR = 8.4, float DR = 0.4)
{

  // inputFile: the path to the delphes root file.
  // outputFile: the name of the output root file to be produced.
  // maxR: the maximum distance allowed between jets.
  // DR: the step size for the spectral function.

  int nJetsSelected = 7;
  double pTmin = 180.;
  double HTcut = 0.;

  std::vector<float> spectralfunction;
  float eventHT, CParam;
  std::array<float, 7> jetspT;
  std::array<float, 7> jetsMass;

  TTree *myTree = new TTree("eventInfo", "Tree with the event information");
  myTree->Branch("SpectralFunc", &spectralfunction);
  myTree->Branch("eventHT", &eventHT, "eventHT/F");
  myTree->Branch("C", &CParam, "CParam/F");
  myTree->Branch("jetspT", &jetspT, "jetspT[7]/F");
  myTree->Branch("jetsMass", &jetsMass, "jetsMass[7]/F");

  Rmax = maxR;
  DeltaR = DR;

  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(Form("%s", inputFile));

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  TH1D *SeventhJetpT = new TH1D("SeventhJetpT", "7th jet pT", 200, 0, 1000);

  // The maximum number of events to be processed
  int maxEvents = 1e6;

  std::vector<Jet> jetEvent;

  // flag to mark the event if it doesn't fulfill the requirements
  bool deleteEvent(false);

  int nEventsJetsCounter = 0;

  for (Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    jetEvent.clear();
    deleteEvent = false;
    eventHT = 0.;
    CParam = 0.;
    spectralfunction.clear();
    std::fill(jetsMass.begin(), jetsMass.end(), 0.);
    std::fill(jetspT.begin(), jetspT.end(), 0.);

    if (entry >= maxEvents)
    {
      break;
    }

    if (entry % 1000 == 0)
      std::cout << entry << " Events has been processed\n";

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    if (branchJet->GetEntries() >= nJetsSelected)
    {
      SeventhJetpT->Fill(((Jet *)branchJet->At(6))->PT);
      for (unsigned int ijet = 0; ijet < branchJet->GetEntries(); ijet++)
      {

        eventHT += ((Jet *)branchJet->At(ijet))->PT;

        if (ijet < nJetsSelected && !deleteEvent)
        {
          if (((Jet *)branchJet->At(ijet))->PT < pTmin)
          {
            deleteEvent = true;
            break;
          }

          jetspT[ijet] = ((Jet *)branchJet->At(ijet))->PT;
          jetsMass[ijet] = ((Jet *)branchJet->At(ijet))->Mass;
        }

        jetEvent.push_back(*((Jet *)branchJet->At(ijet)));

        CParam = calculateC(jetEvent);
      }

      if (!deleteEvent)
      {
        for (double iRbin = 0; iRbin < Rmax; iRbin += DeltaR)
        {
          double spectralfunc = DeltaR * calculateSpectralFunction(jetEvent, iRbin) / (eventHT * eventHT);
          spectralfunction.push_back(spectralfunc);
        }

        jetEvent.clear();
      }
    }
    else
    {
      std::cout << "Empty event\n";
      deleteEvent = true;
    }

    if (!deleteEvent)
    {
      nEventsJetsCounter++;
      myTree->Fill();
    }
  }

  std::cout << "Number of events: " << nEventsJetsCounter << std::endl;

  TFile *outfile = new TFile(Form("TTree_DR%d_Rmax%d_%s.root", int(DeltaR * 10), int(Rmax * 10), outputFile.c_str()), "RECREATE");
  myTree->Write();
  SeventhJetpT->Write();
}
