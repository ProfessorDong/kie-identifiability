#!/usr/bin/env bash
# Retrieve the openly accessible supplementary documents.
#
# No publisher PDF is redistributed in this repository.  This script fetches the
# open-access supplements from Europe PMC so that the transcription audit in
# analysis/build_trinomial.py can verify every extracted value against its
# source.  Two sources are paywalled; see SOURCES.md for manual instructions.
set -u
cd "$(dirname "$0")"
mkdir -p si/manual

ok=0; fail=0
fetch () {                      # $1 = PMCID, $2 = expected file inside
  local pmc="$1" want="$2"
  if [ -f "si/$pmc/$want" ]; then
    echo "  have    $pmc/$want"; ok=$((ok+1)); return
  fi
  echo -n "  fetch   $pmc ... "
  if curl -sL --max-time 120 -o "si/$pmc.zip" \
      "https://www.ebi.ac.uk/europepmc/webservices/rest/$pmc/supplementaryFiles" \
     && [ -s "si/$pmc.zip" ] \
     && [ "$(file -b --mime-type "si/$pmc.zip")" = "application/zip" ]; then
    unzip -o -q "si/$pmc.zip" -d "si/$pmc"
    if [ -f "si/$pmc/$want" ]; then
      echo "ok"; ok=$((ok+1)); return
    fi
    echo "downloaded but $want not inside"; fail=$((fail+1)); return
  fi
  echo "FAILED (not open access, or network unavailable)"
  rm -f "si/$pmc.zip"; fail=$((fail+1))
}

echo "open-access supplements:"
fetch PMC3985941 ja411998h_si_001.pdf    # Singh 2014,    ecDHFR network variants
fetch PMC8697555 bi1c00558_si_001.pdf    # Li 2021,       hsDHFR loop variants
fetch PMC4063187 ja501936d_si_001.pdf    # Francis 2014,  light vs heavy ecDHFR

echo
echo "retrieved $ok, failed $fail"
if [ ! -f si/manual/cs9b03345_si_001.pdf ]; then
  cat <<'EOF'

Still needed, and paywalled (see SOURCES.md):
  Pagano et al., ACS Catal. 9, 11199 (2019)   doi:10.1021/acscatal.9b03345
  Download the Supporting Information PDF and save it as
      external_data/si/manual/cs9b03345_si_001.pdf

Without it the benchmark builds with 14 series / 63 records instead of
16 / 73.  No other result changes.
EOF
fi
