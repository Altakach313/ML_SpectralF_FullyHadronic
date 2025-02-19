## The QCD multijet background results were obtained as follows:

We first generate hepmc events with pythia v8243, using the macro "MultiJetPythia8.C":

```bash
root -l -b -q "MultiJetPythia8.C"
```
The "MultiJetPythia8.C" macro can be found in "Pythia_v8243/".

Then, we run delphes (master version) on the hepmc file to obtain the root file:

```bash
/path_to_delphes/DelphesHepMC3 delphes_card_ATLAS.tcl outputname.root inputfile.hepmc
```
The "delphes_card_ATLAS.tcl" can be found in "Delphes_master/"

We then read the root file using the macro "ReadDelphesJets.C" in order to produce a TTree with the necessary leafs (spectral function, event HT, jets pT, jets masses, and the C parameter):

```bash
root -l -b -q "ReadDelphesJets.C"
```

We use the python script "ConvertTTreeToPandas.py" to convert the TTree to a pandas dataframe. This dataframe will be then used as an input for the ML model.
