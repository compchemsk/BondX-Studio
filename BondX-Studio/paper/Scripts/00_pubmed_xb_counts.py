#!/usr/bin/env python3
"""
pubmed_xb_counts.py -- year-wise counts of halogen-bond papers in PubMed
(reproduces Figure 1 of Manuscript_V1 and Table S1 of ESI_V1).

Uses NCBI E-utilities (esearch, rettype=count) once per year and query.
Requirements: Python 3.8+, requests   (pip install requests)

  python pubmed_xb_counts.py                     # 2003 to current year
  python pubmed_xb_counts.py --start 2000 --end 2026 --email you@uni.edu

Output: pubmed_xb_counts.csv  (year, Q1_all, Q2_bio)
Note: per-year counts use the publication-date field [dp]; they can differ
slightly from the PubMed web "Results by year" timeline.
"""
import argparse, csv, datetime, time
import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
Q1 = '("halogen bond" OR "halogen bonding" OR "halogen bonds")'
Q2 = f'{Q1} AND (protein OR drug OR ligand)'

def count(term, year, email, api_key):
    params = {"db": "pubmed", "term": f"{term} AND {year}[dp]", "rettype": "count",
              "retmode": "json", "tool": "bondx_counts", "email": email}
    if api_key:
        params["api_key"] = api_key
    for attempt in range(5):
        r = requests.get(EUTILS, params=params, timeout=30)
        if r.status_code == 200:
            return int(r.json()["esearchresult"]["count"])
        time.sleep(2 ** attempt)
    r.raise_for_status()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=2003)
    ap.add_argument("--end", type=int, default=datetime.date.today().year)
    ap.add_argument("--email", default="your.name@example.org")
    ap.add_argument("--api-key", default=None, help="NCBI API key (optional, faster)")
    ap.add_argument("--out", default="pubmed_xb_counts.csv")
    a = ap.parse_args()
    pause = 0.11 if a.api_key else 0.35          # NCBI limit: 10 or 3 requests/s
    rows = []
    for y in range(a.start, a.end + 1):
        n1 = count(Q1, y, a.email, a.api_key); time.sleep(pause)
        n2 = count(Q2, y, a.email, a.api_key); time.sleep(pause)
        rows.append((y, n1, n2)); print(f"{y}  all={n1:5d}  bio={n2:5d}")
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["year", "Q1_all", "Q2_bio"]); w.writerows(rows)
    t1 = sum(r[1] for r in rows); t2 = sum(r[2] for r in rows)
    print(f"\nTotal {t1} (all), {t2} (bio), bio share {100*t2/max(t1,1):.1f}%  ->  {a.out}")
    print(f"Search date: {datetime.date.today().isoformat()}")

if __name__ == "__main__":
    main()
