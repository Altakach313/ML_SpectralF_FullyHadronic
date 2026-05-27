# Pythia configuration

The production workflow creates the Pythia command file dynamically inside
`scripts/run_everything.sh`, because the number of events and LHE file name
depend on the batch job.

`pythia_lhef_template.cmnd` records the settings used for LHE showering:

- input comes from the MadGraph LHE file through `Beams:LHEF`
- multiparton interactions are enabled
- parton-level and hadron-level evolution are enabled
- only minimal event-printing is kept

Replace `{nevents}` and `{lhe_file}` when using the template by hand.
