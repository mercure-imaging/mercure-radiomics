#!/usr/bin/env bash
set -Eeo pipefail
echo "-- Starting mercure-radiomics..."
conda run -n mercure-radiomics python radiomics_process.py $MERCURE_IN_DIR $MERCURE_OUT_DIR  
echo "-- Done."