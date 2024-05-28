#!/usr/bin/env python

from collections import defaultdict
import argparse
import requests
import time
import json

POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"

def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def submit_id_mapping(from_db, toDB, ids):
    response = requests.post(
        f"{API_URL}/idmapping/run", 
        data={"from": from_db, "to": toDB, "ids": ids},
    )
    check_response(response)
    return response.json()["jobId"]


def get_id_mapping_results(job_id, max_retries=20):
    for attempt in range(max_retries):
        r = requests.get(f"{API_URL}/idmapping/status/{job_id}")
        r.raise_for_status()
        job = r.json()
        if "jobStatus" in job:
            if job["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(job["jobStatus"])
        else:
            return job


def get_prot_annotations(results):
    """Extracts protein-level annotations (e.g. GO-terms, 
    subcellular locations etc.)"""
    # Initialise dictionaries
    slocs = set() # Curated subcellular locations
    go_func = defaultdict(set) # Molecular function GO terms
    go_proc = defaultdict(set) # Biological process GO terms
    go_comp = defaultdict(set) # Cellular component GO terms
    go_all = defaultdict(set) # All GO terms
    # Get annotations for every provided entry
    for idx, entry in enumerate(results['results']):
        key = entry['from']
        # Subcellular locations
        for ann in entry['to']['comments']:
            if ann['commentType'] == 'SUBCELLULAR LOCATION':
                for sloc in ann['subcellularLocations']:
                    slocs.add(sloc['location']['value'])
        # Go terms
        for ref in entry['to']['uniProtKBCrossReferences']:
            if ref['id'].startswith('GO:'):
                go_all[key].add(ref['id'])
                # Add to specific GO categories
                for p in ref['properties']:
                    descr = p['value'].split(':')[1]
                    go = f'{ref["id"]}|{descr}'
                    if p['value'].startswith('C:'):
                        go_comp[key].add(go)
                    if p['value'].startswith('P:'):
                        go_proc[key].add(go)
                    if p['value'].startswith('F:'):
                        go_func[key].add(go)
    return slocs, go_all, go_func, go_proc, go_comp

def get_res_annotations(results):
    """Extracts annotations at residue level, and variants if they
    exist"""
    # Initialise dictionaries
    features = defaultdict(set)
    variants = defaultdict(set)
    # Parse json and collect info
    for idx, entry in enumerate(results['results']):
        sequence = entry['to']['sequence']['value']
        # Functional annotations and variants
        for feature in entry['to']['features']:
            start = feature['location']['start']['value']
            end = feature['location']['end']['value']
            ftype = feature['type']
            for k, res in enumerate(range(start, end+1)):
                resname = sequence[res-1]
                key = (entry['from'], f'{resname}{res}')
                # Get residue annotations
                features[key].add(ftype)
                # Get variants
                if ftype == 'Natural variant':
                    try:
                        orseq = feature['alternativeSequence']['originalSequence']
                        altseqs = feature['alternativeSequence']['alternativeSequences']
                    except KeyError:
                        continue
                    for v in altseqs:
                        variants[key].add(v[k])
    return features, variants


def write_csv(features, variants, slocs, go_all,
              go_func, go_proc, go_comp, outfile):
    with open(outfile, 'w') as o:
        header = ('uniprot_id,'
                  'resid,'
                  'features,'
                  'variants,'
                  'subcellular_locations,'
                  'go_all,'
                  'go_molecular_function,'
                  'go_biological_process,'
                  'go_cellular_component')
        print(header, file=o)
        for k,v in features.items():
            vnts = variants.get(k, set())   
            _go_all = go_all.get(k[0], set())
            _go_func = go_func.get(k[0], set())
            _go_proc = go_proc.get(k[0], set())
            _go_comp = go_comp.get(k[0], set())
            string = (f'{k[0]},'
                      f'{k[1]},'
                      f'"{";".join(v)}",'
                      f'"{";".join(vnts)}",'
                      f'"{";".join(slocs)}",'
                      f'"{";".join(_go_all)}",'
                      f'"{";".join(_go_func)}",'
                      f'"{";".join(_go_proc)}",'
                      f'"{";".join(_go_comp)}"')
            print(string, file=o)
    return
    

def main(infile, outfile):
    """Collects annotations from UniProt for a list of proteins
    using the ID mapping tool via the API, and writes them in csv."""

    # Parse UniProt IDs
    ids = [id.strip() for id in open(infile, 'r')]

    # Submit job to UniProt
    job_id = submit_id_mapping(
        from_db='UniProtKB_AC-ID',
        toDB='UniProtKB',
        ids=ids)

    # Get results
    results = get_id_mapping_results(job_id)

    # Parse results json to get annotations
    features, variants = get_res_annotations(results)
    slocs, go_all, go_func, go_proc, go_comp = get_prot_annotations(results)
    
    # Write csv
    write_csv(features, variants, slocs, go_all,
              go_func, go_proc, go_comp, outfile)

    # All done
    return
    

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, 
                        help="A list of UniProt IDs")
    parser.add_argument("-o", "--output",  type=str,
                        help="Output csv file")
    args = parser.parse_args()

    # Run
    main(args.input, args.output)

