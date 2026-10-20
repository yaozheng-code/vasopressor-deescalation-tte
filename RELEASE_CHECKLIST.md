# Before clicking "Public" on GitHub

- [ ] All authors reviewed the repository contents and agreed on public release.
- [ ] No MIMIC-IV records, individual-level derivatives, outputs, logs, or secrets.
- [ ] Choose and add a software `LICENSE` (requires coauthor approval).
- [ ] Install dependencies and run `python run_pipeline.py --self-test-only`.
- [ ] Run `python check_integrity.py` and `python -m unittest discover -s tests -v`.
- [ ] Run `python preflight.py ...` on an authorized local MIMIC-IV database.
- [ ] Run complete `G10` (500 bootstrap) and `G11`, then `verify_results.py`.
- [ ] Verify key cohort counts and compare outputs to manuscript Table 1/2.
- [ ] Pin the version/commit of MIT-LCP derived-concept SQL used locally.
- [ ] Add the public repository URL to the manuscript's Code availability statement.
- [ ] After acceptance, add publication DOI to CITATION.cff and archive a GitHub release (e.g., Zenodo).

The code is public-ready after the authors finish these checks; the
computational full-data rerun can only be completed with credentialed access.
