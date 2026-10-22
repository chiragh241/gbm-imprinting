"""Print the free download URL, md5 and size of a run's full .sra file from NCBI's locate service.

usage: python sra_locate.py RUN   ->  "<url>\t<md5>\t<size>"

The service lists the normalized .sra (with base qualities) and a .lite copy without qualities;
trimming needs the qualities, so only the full file is taken, and only a location that is not
requester-pays.
"""
import json
import sys
import urllib.request

run = sys.argv[1]
url = f"https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc={run}&accept-alternate-locations=yes"
with urllib.request.urlopen(url, timeout=120) as r:
    d = json.load(r)
files = [f for f in d["result"][0]["files"] if f["name"] == run and f["type"] == "sra"]
if len(files) != 1:
    sys.exit(f"{run}: expected one full .sra file, found {[f['name'] for f in d['result'][0]['files']]}")
f = files[0]
links = [loc["link"] for loc in f["locations"] if not loc.get("payRequired")]
if not links:
    sys.exit(f"{run}: no free location")
print(f"{links[0]}\t{f['md5']}\t{f['size']}")
