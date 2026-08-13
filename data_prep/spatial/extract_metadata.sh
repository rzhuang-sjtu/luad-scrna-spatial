#!/bin/bash
# Extract per-sample disease & site from sdrf
SDRF="$HOME/luad/data/ST/E-MTAB-13530/E-MTAB-13530/E-MTAB-13530.sdrf.txt"

# Find column indices for: Source Name, Characteristics[individual], Characteristics[disease], Characteristics[sampling site], Characteristics[disease staging]
awk -F'\t' '
NR==1 {
  for (i=1; i<=NF; i++) {
    h[$i] = i
    if ($i ~ /^Source Name/) src=i
    if ($i ~ /individual/) ind=i
    if ($i ~ /^Characteristics\[disease\]/) dis=i
    if ($i ~ /^Characteristics\[sampling site\]/) site=i
    if ($i ~ /disease staging/) stage=i
    if ($i ~ /^Characteristics\[sex\]/) sex=i
    if ($i ~ /^Characteristics\[age\]/) age=i
  }
  printf "source\tindividual\tdisease\tstage\tsite\tsex\tage\n"
  next
}
{
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $src, $ind, $dis, $stage, $site, $sex, $age
}' "$SDRF" | sort -u
