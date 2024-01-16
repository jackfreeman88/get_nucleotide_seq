import time
import json
from Bio import Entrez
from xml.etree import ElementTree as ET
from urllib.error import HTTPError

Entrez.email = "your_email@example.com"  # Replace with your email


def get_gene_upstream_region(
    gene, organism="Homo sapiens", max_retries=3, retry_delay=2, upstream_length=2000
):
    attempts = 0
    while attempts < max_retries:
        try:
            search_term = f"{gene}[Gene Name] AND {organism}[Organism]"
            handle = Entrez.esearch(db="gene", term=search_term)
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
                        gi = int(interval.find(".//Seq-id_gi").text)

                        if strand == "minus":
                            upstream_region = (to_bp + 1, to_bp + 1 + upstream_length)
                        elif strand == "plus":
                            upstream_region = (from_bp - upstream_length, from_bp - 1)
                        else:
                            # exit if strand information is not available
                            return {
                                "gene": gene,
                                "strand": "Information not available",
                                "upstream_region": None,
                                "upstream_seq": None,
                            }
                        handle = Entrez.efetch(
                            db="sequences",
                            id=gi,
                            rettype="fasta",
                            retmode="text",
                            seq_start=upstream_region[0],
                            seq_stop=upstream_region[1],
                        )
                        upstream_seq = handle.read()
                        handle.close()
                        return {
                            "gene": gene,
                            "strand": strand,
                            "upstream_region": upstream_region,
                            "upstream_seq": upstream_seq,
                        }
            else:
                # exit if gene not found
                return {
                    "gene": gene,
                    "strand": "Information not available",
                    "upstream_region": None,
                    "upstream_seq": None,
                }
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
        "upstream_seq": None,
    }


def main():
    # Read gene names from a file
    with open(
        "/isiseqruns/jfreeman_tmp_home/saher_collab/get_nucleotide_seq/gene_list.txt",
        "r",
    ) as file:
        gene_list = file.read().splitlines()

    # Store results in a dictionary
    results = {}
    for gene in gene_list:
        gene_info = get_gene_upstream_region(gene)
        results[gene] = gene_info
        time.sleep(1)  # Rate limiting - wait for 1 second between requests
        print(f"Finished processing {gene}")

    # Save results to a JSON file
    with open("gene_upstream_info_results.json", "w") as output_file:
        json.dump(results, output_file, indent=4)


if __name__ == "__main__":
    main()
