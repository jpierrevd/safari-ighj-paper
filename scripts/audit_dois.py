#!/usr/bin/env python3
"""Audit DOIs in references.bib against Crossref API."""
import re
import urllib.request
import json
import sys

BIB = "../references.bib"

# Extract DOIs from bib file
with open(BIB) as f:
    content = f.read()

dois = re.findall(r'doi\s*=\s*\{([^}]+)\}', content)
keys = re.findall(r'@article\{(\w+)', content)

print(f"Found {len(dois)} DOIs in {len(keys)} entries\n")
print(f"{'Key':<20} {'DOI':<50} {'Status'}")
print("-" * 90)

valid = 0
errors = []

for key, doi in zip(keys, dois):
    try:
        url = f"https://api.crossref.org/works/{doi}"
        req = urllib.request.Request(url, headers={
            'User-Agent': 'SAFARI-IGHJ-Audit/1.0 (mailto:jpierre.vd@gmail.com)'
        })
        resp = urllib.request.urlopen(req, timeout=10)
        data = json.loads(resp.read())
        title = data['message'].get('title', ['?'])[0][:50]
        print(f"{key:<20} {doi:<50} OK ({title}...)")
        valid += 1
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"{key:<20} {doi:<50} NOT FOUND")
            errors.append((key, doi, "404"))
        else:
            print(f"{key:<20} {doi:<50} HTTP {e.code}")
            errors.append((key, doi, f"HTTP {e.code}"))
    except Exception as e:
        print(f"{key:<20} {doi:<50} ERROR: {e}")
        errors.append((key, doi, str(e)))

print("\n" + "=" * 90)
print(f"RESULT: {valid}/{len(dois)} DOIs validated successfully")
if errors:
    print(f"ERRORS ({len(errors)}):")
    for key, doi, err in errors:
        print(f"  - {key}: {doi} -> {err}")
else:
    print("All DOIs are valid and active.")
