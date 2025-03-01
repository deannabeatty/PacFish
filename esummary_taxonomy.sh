while read taxid; do
    esummary -db taxonomy -id $taxid | xtract -pattern DocumentSummary -element TaxId,ScientificName,CommonName,Rank,Division
done < staxids.txt > taxid_scientific_name.txt