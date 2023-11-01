from Bio import Entrez
from xml.etree import ElementTree as ET

Entrez.email = "jfreeman@morgridge.org"


def get_gene_strand_info(gene):
    handle = Entrez.esearch(db="gene", term=gene)
    record = Entrez.read(handle)
    handle.close()

    if record["IdList"]:
        gene_id = record["IdList"][0]
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="xml")
        xml_data = handle.read()  # Getting raw XML
        handle.close()

        tree = ET.ElementTree(ET.fromstring(xml_data))
        root = tree.getroot()

        # Navigating through XML to find the strand information
        for elem in root.iter("Seq-interval_strand"):
            strand_info = elem.find("Na-strand").attrib["value"]
            return {"gene": gene, "strand": strand_info}
    else:
        return {"gene": gene, "strand": "Information not available"}


# Testing the function with a gene name
print(get_gene_strand_info("BRCA1"))
