# Source census: search record

Run on 17 October 2026 against the Europe PMC REST API
(`https://www.ebi.ac.uk/europepmc/webservices/rest/search`, `resultType=lite`,
`pageSize=1000`).  The union of the five queries, deduplicated by DOI (or
Europe PMC identifier where no DOI was listed), gave 88 records.

1. `AUTH:"Kohen A" AND ("intrinsic kinetic isotope" OR "intrinsic KIE" OR "intrinsic isotope effect")`
2. `AUTH:"Klinman JP" AND ("D/T" OR "deuterium/tritium" OR "Swain-Schaad" OR "intrinsic isotope effect" OR "tritium isotope effect")`
3. `AUTH:"Kohen A" AND ("H/T" OR "Swain-Schaad" OR "mixed labeling" OR "competitive" ) AND ("isotope effect")`
4. `("Northrop method" OR "Northrop's method" OR "Northrop equation") AND ("intrinsic" ) AND ("isotope effect")`
5. `"Swain-Schaad" AND ("intrinsic kinetic isotope effect" OR "intrinsic isotope effect") AND ("tritium")`

Each record was screened for competitive H/T and D/T isotope effects.  Every
primary source that a screened article, or one of five reviews among them,
cites for such data was added and read, supplement included; the tracing
stopped when no cited source of matched H/T and D/T data remained unexamined.
The disposition of every source that reports a tritium isotope effect or a
matched pair is the supplement's table "Disposition of the sources located by
the search of 17 October 2026 and by citation tracing".

All four concentrations in Table 5 of Tsai and Klinman (2001) are
transcribed. The 10 mM condition repeats the 3 °C arm of their Table 4. The
other three conditions belong to a control experiment on the order of
substrate binding; they are listed in the 247-pair ledger, analyzed
separately, and not counted as analysis units of the 244-pair dataset. No additional numeric pairs
were found in the checked tables of Cha (1989), Rucker (1992), or Hong (2007).
Hong's Figure 5 plots both primary isotope effects for Y94F hydride transfer.
Only H/T is given numerically with uncertainty in the prose. These records
describe coverage of the documented search and its numerical admission rule.
