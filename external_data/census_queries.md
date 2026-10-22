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
matched pair is Supplementary Table "Disposition of every source examined".
