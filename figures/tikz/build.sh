#!/usr/bin/env bash
# Build every TikZ figure to a standalone PDF and install into manuscript/figures.
set -u
cd "$(dirname "$0")"
ok=0; fail=0
for f in fig*.tex; do
  b="${f%.tex}"
  if pdflatex -interaction=nonstopmode -halt-on-error "$f" >"$b.buildlog" 2>&1; then
    echo "  OK    $b"; ok=$((ok+1))
  else
    echo "  FAIL  $b"; grep -m4 -E '^!|^l\.[0-9]' "$b.buildlog" | sed 's/^/        /'
    fail=$((fail+1))
  fi
done
rm -f *.aux *.log *.buildlog
mkdir -p ../../manuscript/figures
cp -f fig*.pdf ../../manuscript/figures/ 2>/dev/null
cp -f fig*.pdf ../ 2>/dev/null
echo "built $ok, failed $fail"
