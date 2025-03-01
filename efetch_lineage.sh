cut -f13 blastn_results.out | sort | uniq > staxids.txt
while read taxid; do
    efetch -db taxonomy -id $taxid -format xml | xtract -pattern Taxon -element TaxId -block LineageEx -sep ";" -element ScientificName
done < staxids.txt > taxid_to_lineage.txt