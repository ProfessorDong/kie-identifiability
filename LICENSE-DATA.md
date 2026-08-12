# Licensing of the data

**Code** in `analysis/` and `figures/` is under the MIT License (see `LICENSE`).

**Curated data** in `data/` and generated tables in `results/` are released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), with the following
qualifications.

The numerical values in `data/trinomial_benchmark.csv` are measurements
transcribed from the primary publications listed in `external_data/SOURCES.md`.
Individual measured values are facts and are not themselves subject to
copyright; the selection, curation and arrangement here are the contribution
released under CC BY 4.0. **Anyone using these numbers should cite the original
measurement papers**, not only this repository. Each record carries its
`source_DOI`, `PMCID` and `source_table` for exactly that purpose.

The replicate-level thymidylate synthase tables derive from Islam et al.,
*PLoS ONE* **13**, e0196506 (2018), published under CC BY 4.0, and are used
under those terms.

**No publisher PDF is redistributed here.** `external_data/fetch_sources.sh`
retrieves the open-access supplements from Europe PMC at the user's own request;
paywalled documents must be obtained by the user through their own institutional
access.
