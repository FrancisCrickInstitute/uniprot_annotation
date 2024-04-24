#!/usr/bin/env python

from collections import defaultdict
import argparse
import requests
import time
import json

POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"
FEATURES = {'Modified residue', 
            'Natural variant', 
            'Binding site'}

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

def extract_res_annotations(results):
    """Extracts annotations at residue level, and variants if they
    exist"""
    # Initialise dictionaries
    annotations = defaultdict(set)
    variants = defaultdict(set)
    # Parse json and collect info
    for idx, entry in enumerate(results['results']):
        sequence = entry['to']['sequence']['value']
        for feature in entry['to']['features']:
            start = feature['location']['start']['value']
            end = feature['location']['end']['value']
            ftype = feature['type']
            if ftype in FEATURES:
                for k, res in enumerate(range(start, end+1)):
                    resname = sequence[res-1]
                    key = (entry['from'], f'{resname}{res}')
                    # Get residue annotations
                    annotations[key].add(ftype)
                    # Get variants
                    if ftype == 'Natural variant':
                        try:
                            orseq = feature['alternativeSequence']['originalSequence']
                            altseqs = feature['alternativeSequence']['alternativeSequences']
                        except KeyError:
                            continue
                        for v in altseqs:
                            variants[key].add(v[k])
    return annotations, variants

def write_csv(annotations, variants, outfile):
    with open(outfile, 'w') as o:
        print('uniprot_id,resid,annotations,variants', file=o)
        for k,v in annotations.items():
            vnts = variants.get(k, set())   
            print(f'{k[0]},{k[1]},{";".join(v)},{";".join(vnts)}', file=o)
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
    annotations, variants = extract_res_annotations(results)
    
    # Write csv
    write_csv(annotations, variants, outfile)

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

