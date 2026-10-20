# Publishing this code to GitHub

No GitHub account access is required by these instructions.

## Browser upload (simple)

1. Create a new GitHub repository named
   `ne-avp-deescalation-tte` (or another descriptive name).
2. Keep the repo **private** while you check for data/credential leaks.
3. Unzip the downloaded release archive locally.
4. Upload **only its contents** (`README.md`, `analysis/`, `tests/`,
   `docs/`, `.gitignore`, `requirements.txt`, etc.). Ensure `.gitignore`
   is present.
5. Check the file list: there should be **no raw MIMIC data, outputs,
   `*.csv`, `*.parquet`, patient-level rows, `*.zip` of outputs,
   `.pgpass`, `.env`, or secrets**.
6. Choose/add `LICENSE` only after all code owners approve a license.
7. Run the offline tests (or inspect their log) and publish the repository.
8. Add the final GitHub URL to the manuscript's Code Availability statement.

## Git command line alternative

After unzipping into `ne-avp-deescalation-tte/`:

```bash
cd ne-avp-deescalation-tte
git init
git add .
git status --short
git diff --cached --name-only
# CHECK no restricted files are staged
git commit -m "Release frozen NE/AVP CCW reproducibility code"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/ne-avp-deescalation-tte.git
git push -u origin main
```

Use your own authenticated GitHub setup for `git push`. No credentials
or access tokens should ever be included in any committed file.

## Manuscript Code Availability text

> The frozen analysis scripts, synthetic-data tests, and reproducibility
> instructions are available at [INSERT_PUBLIC_GITHUB_REPOSITORY_URL].
> MIMIC-IV v3.1 is accessible to credentialed researchers through
> PhysioNet under the applicable data-use agreement; no participant-level
> data are redistributed with the code.

After manuscript acceptance, archive a versioned release with a
persistent DOI (e.g., a software archive service). The manuscript
DOI and software DOI must not be invented before publication.
