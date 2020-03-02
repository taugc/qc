for name in *.fastq.gz; do
    printf '%s\n' "${name%_L00[1234]_R[12]*}"
done | uniq |
while read prefix; do
    cat "$prefix"*R1*.fastq.gz >"${prefix}_R1_001.fastq.gz"
    rm "$prefix"_L00[1234]_R1*.fastq.gz
    cat "$prefix"*R2*.fastq.gz >"${prefix}_R2_001.fastq.gz"
    rm "$prefix"_L00[1234]_R2*.fastq.gz
done

