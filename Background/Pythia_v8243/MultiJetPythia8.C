/// File that generate QCD multijet events and write them into hepmc format

#include <iostream>
#include <fstream>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "TRandom.h"
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"

int MultiJetPythia8(int nprocess = 0, bool doHadronShower = true, std::string outputname = "defaults", unsigned int random = 0, int maxEvents = 1000)
{

    /// nprocess: the process number on the cluster/grid
    /// doHadronShower: whether to perform the shower at parton level (false) or hadron level (true)
    /// outputname: The hepmc file name
    /// random: random seed for generator
    /// maxEvents: the number of events to be generated

    outputname.append("_" + std::to_string(nprocess) + ".hepmc");
    Pythia8::Pythia8ToHepMC topHepMC(outputname.c_str());

    Pythia8::Pythia pythia;
    pythia.readString("Beams:eCM = 13000.");
    pythia.readString("HardQCD:all = on");
    pythia.readString("PhaseSpace:pTHatMin = 1300.");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.init();

    pythia.readString("ProcessLevel:all = on");
    pythia.readString("PartonLevel:all = on");
    pythia.readString("PartonLevel:MPI = on");
    if (doHadronShower)
    {
        pythia.readString("HadronLevel:all = on");
    }
    else
    {
        pythia.readString("HadronLevel:all = off");
    }

    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", random));

    double ptmin = 150.;

    double Rval = 0.4;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, Rval);

    int nEvents7Jets = 0.;

    if (pythia.init())
    {

        bool activated = false;
        while (nEvents7Jets < maxEvents)
        {
            pythia.next();

            std::vector<fastjet::PseudoJet> input_particles;

            Pythia8::Event event = pythia.event;

            for (unsigned int ipart = 0; ipart < event.size(); ipart++)
            {
                Pythia8::Particle particle = event[ipart];
                if (particle.status() <= 0)
                {
                    continue;
                }
                input_particles.push_back(fastjet::PseudoJet(particle.px(), particle.py(), particle.pz(), particle.e()));
            }

            fastjet::ClusterSequence clust_seq(input_particles, jet_def);
            std::vector<fastjet::PseudoJet> inclusive_jets = fastjet::sorted_by_pt(clust_seq.inclusive_jets(ptmin));

            if (inclusive_jets.size() >= 7)
            {
                nEvents7Jets++;
                topHepMC.writeNextEvent(pythia);
            }
        }
    }

    std::cout << "We have nEvents with 7 jets: " << nEvents7Jets << std::endl;

    topHepMC.setXSec(pythia.info.sigmaGen(), pythia.info.sigmaErr());
    pythia.stat();

    // Create a root file that has histograms filled with:

    TFile *xsecFile = new TFile(Form("xsection_%d.root", nprocess), "RECREATE");

    // The cross section
    TProfile *xsec = new TProfile("xsection", "average cross section", 1, 0, 1);
    xsec->Fill(0.5, pythia.info.sigmaGen());

    // The number of events generated
    TH1D *hnEvents = new TH1D("hnTrials", "Number of trials", 1, 0, 1);
    hnEvents->SetBinContent(1, pythia.info.nTried());

    // The number of events with pThat > pTHardmin
    TH1D *hnEventsAcc = new TH1D("hnEventsAcc", "Number of accepted events", 1, 0, 1);
    hnEventsAcc->SetBinContent(1, pythia.info.nAccepted());

    // The number of events with pThat > pTHardmin, that has at least 7 jets, each with pt > ptmin
    TH1D *hnEvents7Jets = new TH1D("hnEvents7Jets", "Number of events with 7 jets", 1, 0, 1);
    hnEvents7Jets->SetBinContent(1, nEvents7Jets);

    xsecFile->Write();

    return 0;
}
