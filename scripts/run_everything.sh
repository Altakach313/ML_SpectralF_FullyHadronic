#!/usr/bin/env bash

set -euo pipefail
set -x

CLUSTER="${1:-0}"
PROCESS="${2:-0}"
SAMPLE="${3:-background_ttbar}"
NEVENTS="${4:-1000000}"
DO_MADSPIN="${5:-false}"
GLUINO_MASS="${6:-2000}"

DEFAULT_REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_DIR="${REPO_DIR:-${DEFAULT_REPO_DIR}}"
DELPHES_DIR="${DELPHES_DIR:-${REPO_DIR}/delphes}"
MG5_BIN="${MG5_BIN:-${REPO_DIR}/MG5/bin/mg5_aMC}"
MADSPIN_BIN="${MADSPIN_BIN:-${REPO_DIR}/MG5/MadSpin/madspin}"
DELPHES_PYTHIA8="${DELPHES_PYTHIA8:-${DELPHES_DIR}/DelphesPythia8}"
DELPHES_CARD="${DELPHES_CARD:-${DELPHES_DIR}/cards/delphes_card_ATLAS.tcl}"

echo "Executing run ${PROCESS} on $(hostname) in $(pwd)"
echo "Cluster ${CLUSTER}, Process ${PROCESS}, Sample = ${SAMPLE}, Gluino mass = ${GLUINO_MASS}"

if command -v alienv >/dev/null 2>&1; then
  eval "$(alienv printenv \
    VO_ALICE@ROOT::v6-32-06-alice1-15,\
    VO_ALICE@fastjet::v3.4.1_1.052-alice2-35,\
    VO_ALICE@pythia::v8311-26,\
    VO_ALICE@HepMC3::3.3.0-33)"
fi

export ROOT_INCLUDE_PATH="${ROOT_INCLUDE_PATH:-}:${REPO_DIR}:${DELPHES_DIR}:${DELPHES_DIR}/classes:${DELPHES_DIR}/external:${DELPHES_DIR}/modules"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:${REPO_DIR}:${DELPHES_DIR}"

case "${SAMPLE}" in
  signal_gluino_pair)
    TASK_NAME="signalGoGo"
    OUTPUT_TAG="${TASK_NAME}_${GLUINO_MASS}_${PROCESS}"
    ;;
  background_ttbar)
    TASK_NAME="ttbarBackground"
    OUTPUT_TAG="${TASK_NAME}_${PROCESS}"
    ;;
  *)
    echo "Unknown sample '${SAMPLE}'. Use signal_gluino_pair or background_ttbar." >&2
    exit 2
    ;;
esac

sed -i "s#^set nevents .*#set nevents ${NEVENTS}#" proc_card.dat
sed -i "s#^output .*#output ${TASK_NAME}#" proc_card.dat

if [ "${SAMPLE}" = "signal_gluino_pair" ]; then
  sed -i "s#^set Mgo .*#set Mgo ${GLUINO_MASS}.0#" proc_card.dat
  Neu1_MASS=$((GLUINO_MASS - 700))
  sed -i "s#^set Mneu1 .*#set Mneu1 ${Neu1_MASS}#" proc_card.dat
fi

printf "\nset iseed %s\n" "${RANDOM}" >> proc_card.dat

time "${MG5_BIN}" proc_card.dat

if [ "${DO_MADSPIN}" = "true" ]; then
  echo "Running MadSpin"
  sed -i "s#^import filename.*#import ${TASK_NAME}/Events/run_01/unweighted_events.lhe.gz#" madspin_card.dat
  time "${MADSPIN_BIN}" madspin_card.dat
  LHEFILE="${TASK_NAME}/Events/run_01/unweighted_events_decayed.lhe"
else
  LHEFILE="${TASK_NAME}/Events/run_01/unweighted_events.lhe"
fi

gunzip "${LHEFILE}.gz"

PYTHIA_CFG="pythia_${TASK_NAME}.cmnd"

cat > ${PYTHIA_CFG} <<EOF
! ===============================
! Pythia8 configuration
! ===============================

Main:numberOfEvents = ${NEVENTS}         ! number of events to generate
Main:timesAllowErrors = 5

Init:showChangedSettings = off
Init:showChangedParticleData = off

Next:numberShowInfo = 1
Next:numberShowProcess = 1
Next:numberShowEvent = 1

Beams:frameType = 4
Beams:LHEF = ${LHEFILE}

ProcessLevel:all = on
PartonLevel:all = on
PartonLevel:MPI = on
HadronLevel:all = on
EOF

ROOT_OUTPUT="${OUTPUT_TAG}.root"

time "${DELPHES_PYTHIA8}" "${DELPHES_CARD}" "${PYTHIA_CFG}" "${ROOT_OUTPUT}"

echo "Delphes output written to ${ROOT_OUTPUT}"
echo "Job finished successfully."
