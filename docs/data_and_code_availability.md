# Data and Code Availability

The event-generation cards, Pythia showering settings, Delphes configuration,
modified Delphes FastJet module, spectral-function definition, feature schema,
and trained machine-learning model files used in this work are available at the
companion repository:
`https://github.com/Altakach313/ML_SpectralF_FullyHadronic`. The repository
contains the MadGraph, Pythia, and Delphes inputs needed to reproduce the
simulated samples, documents how the spectral-function observables were
constructed from reconstructed jets, and provides trained Keras classifiers for
representative signal masses. Large generated event files and derived training
tables are not included because of their size, but can be regenerated from the
supplied cards and workflow description.
