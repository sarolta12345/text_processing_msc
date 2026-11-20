BEGIN {
    OFS = "\t"
    print "date"
}
{
    printf "%s/%s/%s\n", $2, $3, $4
    samples[$2]++
}
END {
    print ""
    for (year in samples)
        printf "%d sample(s) from %d\n", samples[year], year
}
