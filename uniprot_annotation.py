#!/usr/bin/env/python3

from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
from collections import defaultdict
from requests.adapters import HTTPAdapter, Retry
import re
import time
import json
import zlib
import argparse
import requests


POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"


retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results


def get_id_mapping_results_stream(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/results/stream/")
    request = session.get(url)
    check_response(request)
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)


def get_prot_annotations(results):
    """Extracts protein-level annotations (e.g. GO-terms, 
    subcellular locations etc.)"""
    # Initialise dictionaries
    slocs = defaultdict(set) # Curated subcellular locations
    go_func = defaultdict(set) # Molecular function GO terms
    go_proc = defaultdict(set) # Biological process GO terms
    go_comp = defaultdict(set) # Cellular component GO terms
    go_all = defaultdict(set) # All GO terms
    # Get annotations for every provided entry
    for idx, entry in enumerate(results['results']):
        key = entry['from']
        # Subcellular locations
        try:
            for ann in entry['to']['comments']:
                if ann['commentType'] == 'SUBCELLULAR LOCATION':
                    for sloc in ann['subcellularLocations']:
                        slocs[key].add(sloc['location']['value'])
        except KeyError:
            continue
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
            if not start or not end:
                continue
            ftype = feature['type']
            for k, res in enumerate(range(start, end+1)):
                try:
                    resname = sequence[res-1]
                except IndexError:
                    continue
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
                        try:
                            variants[key].add(v[k])
                        except Exception:
                            variants[key].add('indel')
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
            _slocs = slocs.get(k[0], set())
            _go_all = go_all.get(k[0], set())
            _go_func = go_func.get(k[0], set())
            _go_proc = go_proc.get(k[0], set())
            _go_comp = go_comp.get(k[0], set())
            string = (f'{k[0]},'
                      f'{k[1]},'
                      f'"{";".join(v)}",'
                      f'"{";".join(vnts)}",'
                      f'"{";".join(_slocs)}",'
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
        to_db='UniProtKB',
        ids=ids)


    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
    
    with open('results.json', 'w') as o:
        json.dump(results, o)

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

