#!/usr/bin/env python

import sys
import Bio.Entrez
import Bio.SeqIO


Bio.Entrez.email="smithby1@msu.edu"

#query = sys.argv[1] # List of genes and the query used to find the.
def main():
    uids = get_uids(sys.argv[1])
    records = get_records(uids)
    Bio.SeqIO.write(records, sys.stdout, 'fasta')

def get_records(uids, db="protein"):
    handle3 = Bio.Entrez.efetch(db=db, rettype="gb", retmode="text", id=uids)
    records = Bio.SeqIO.parse(handle3, "genbank")
    return records

def get_uids(query, db="protein"):
    #handle1 = Bio.Entrez.einfo(db=db)
    #result1 = Bio.Entrez.read(handle1, validate=False)
    #sys.stderr.write("The database was last updated on %s\n" % result1["DbInfo"]["LastUpdate"])
    #handle1.close()
    handle2 = Bio.Entrez.esearch(db=db, term=query, retmax=9999)
    result2 = Bio.Entrez.read(handle2)
    uids = result2["IdList"]
    sys.stderr.write("%d UIDs of %d returned.\n" % \
                     (len(uids), int(result2["Count"])))
    handle2.close()
    assert len(uids) > 0
    return uids

def get_ecs(query):
    search_handle = Bio.Entrez.esearch(db="gene", term=query, retmax=9999)
    search_result = Bio.Entrez.read(search_handle)
    uids = search_result["IdList"]
    sys.stderr.write("%d UIDs of %d returned.\n" % \
                     (len(uids), int(search_result["Count"])))
    search_handle.close()
    fetch_handle = Bio.Entrez.efetch(db="gene", retmode="xml", id=uids)
    records = Bio.Entrez.read(fetch_handle)
    ecs = {}
    for rec in records:
        try:
            ecs_in_this_rec = rec['Entrezgene_prot']['Prot-ref']['Prot-ref_ec']
        except KeyError:
            sys.stderr.write("Gene record did not have an EC number.\n")
            continue
        else:
            for ec in ecs_in_this_rec:
                if ec in ecs.keys():
                    ecs[ec] += 1
                else:
                    ecs[ec] = 1
    return ecs



if __name__ == "__main__":
    main()
