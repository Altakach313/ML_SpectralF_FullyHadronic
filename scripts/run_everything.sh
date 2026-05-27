#!/usr/bin/env bash

set -euo pipefail
set -x

CLUSTER="${1:-0}"
PROCESS="${2:-0}"
SAMPLE="${3:-background_ttbar}"
GLUINO_MASS="${4:-2000}"
NEVENTS="${NEVENTS:-1000000}"
DO_MADSPIN="${DO_MADSPIN:-false}"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK_DIR="${WORK_DIR:-$PWD}"
MG5_BIN="${MG5_BIN:-mg5_aMC}"
MADSPIN_BIN="${MADSPIN_BIN:-madspin}"
DELPHES_PYTHIA8="${DELPHES_PYTHIA8:-DelphesPythia8}"
DELPHES_CARD="${DELPHES_CARD:-$REPO_DIR/simulation/delphes/delphes_card_ATLAS.tcl}"

echo "Executing cluster ${CLUSTER}, process ${PROCESS}, sample ${SAMPLE} in ${WORK_DIR}"

if command -v alienv >/dev/null 2>&1; then
  eval "$(alienv printenv \
    VO_ALICE@ROOT::v6-32-06-alice1-15,\
    VO_ALICE@fastjet::v3.4.1_1.052-alice2-35,\
    VO_ALICE@pythia::v8311-26,\
    VO_ALICE@HepMC3::3.3.0-33)"
fi

mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

case "${SAMPLE}" in
  signal_gluino_pair)
    TASK_NAME="signalGoGo"
    CARD_DIR="${REPO_DIR}/simulation/madgraph/signal_gluino_pair"
    ;;
  background_ttbar)
    TASK_NAME="ttbarBackground"
    CARD_DIR="${REPO_DIR}/simulation/madgraph/background_ttbar"
    ;;
  *)
    echo "Unknown sample '${SAMPLE}'. Use signal_gluino_pair or background_ttbar." >&2
    exit 2
    ;;
esac

cp "${CARD_DIR}/proc_card.dat" proc_card.dat
cp "${CARD_DIR}/run_card.dat" run_card.dat
if [ -f "${CARD_DIR}/madspin_card.dat" ]; then
  cp "${CARD_DIR}/madspin_card.dat" madspin_card.dat
fi

sed -i "s#^set nevents .*#set nevents ${NEVENTS}#" proc_card.dat
sed -i "s#^output .*#output ${TASK_NAME}#" proc_card.dat

if [ "${SAMPLE}" = "signal_gluino_pair" ]; then
  sed -i "s#^set Mgo .*#set Mgo ${GLUINO_MASS}.0#" proc_card.dat
fi

printf "\nset iseed %s\n" "${RANDOM}" >> proc_card.dat

time "${MG5_BIN}" proc_card.dat

if [ "${DO_MADSPIN}" = "true" ]; then
  echo "Running MadSpin"
  sed -i "s#^import filename.*#import ${TASK_NAME}/Events/run_01/unweighted_events.lhe.gz#" madspin_card.dat
  time "${MADSPIN_BIN}" madspin_card.dat
  LHEFILE="${TASK_NAME}/Events/run_01/unweighted_events_decayed.lhe"
  gunzip "${LHEFILE}.gz"
else
  LHEFILE="${TASK_NAME}/Events/run_01/unweighted_events.lhe"
  gunzip "${LHEFILE}.gz"
fi

PYTHIA_CFG="pythia_${TASK_NAME}.cmnd"

cat > "${PYTHIA_CFG}" <<EOF
Main:numberOfEvents = ${NEVENTS}
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

ROOT_OUTPUT="${TASK_NAME}_${GLUINO_MASS}_${PROCESS}.root"

time "${DELPHES_PYTHIA8}" "${DELPHES_CARD}" "${PYTHIA_CFG}" "${ROOT_OUTPUT}"

echo "Delphes output written to ${ROOT_OUTPUT}"
