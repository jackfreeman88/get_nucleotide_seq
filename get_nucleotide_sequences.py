import time
import json
from Bio import Entrez
from xml.etree import ElementTree as ET
from urllib.error import HTTPError

Entrez.email = "your_email@example.com"  # Replace with your email


def get_gene_upstream_region(gene, max_retries=3, retry_delay=2, upstream_length=2000):
    attempts = 0
    while attempts < max_retries:
        try:
            handle = Entrez.esearch(db="gene", term=gene)
            record = Entrez.read(handle)
            handle.close()

            if record["IdList"]:
                gene_id = record["IdList"][0]
                handle = Entrez.efetch(
                    db="gene", id=gene_id, rettype="gb", retmode="xml"
                )
                xml_data = handle.read()
                handle.close()

                tree = ET.ElementTree(ET.fromstring(xml_data))
                root = tree.getroot()

                for feature in root.iter("Entrezgene_locus"):
                    for interval in feature.iter("Seq-interval"):
                        from_bp = int(interval.find("Seq-interval_from").text)
                        to_bp = int(interval.find("Seq-interval_to").text)
                        strand = interval.find(".//Na-strand").attrib["value"]

                        if strand == "minus":
                            upstream_region = (to_bp + 1, to_bp + 1 + upstream_length)
                        else:
                            upstream_region = (from_bp - upstream_length, from_bp - 1)

                        return {
                            "gene": gene,
                            "strand": strand,
                            "upstream_region": upstream_region,
                        }
            break
        except HTTPError as e:
            if e.code == 400:
                print(
                    f"HTTP Error 400 occurred for gene {gene}. Retrying... (Attempt {attempts + 1})"
                )
                attempts += 1
                time.sleep(retry_delay)
            else:
                print(f"Error occurred for gene {gene}: {e}")
                break
        except Exception as e:
            print(f"Error occurred for gene {gene}: {e}")
            break

    return {
        "gene": gene,
        "strand": "Information not available or max retries reached",
        "upstream_region": None,
    }


def main():
    # Read gene names from a file
    with open("gene_list.txt", "r") as file:
        gene_list = file.read().splitlines()

    # Store results in a dictionary
    results = {}
    for gene in gene_list:
        gene_info = get_gene_upstream_region(gene)
        results[gene] = gene_info
        time.sleep(1)  # Rate limiting - wait for 1 second between requests

    # Save results to a JSON file
    with open("gene_upstream_info_results.json", "w") as output_file:
        json.dump(results, output_file, indent=4)


if __name__ == "__main__":
    main()
