# sed 's=^//\s*==' nf/nf_01/nf_with_r.nf > nf_01.qmd
no_prefix=${1#nf/*}
sed 's=^//\s*==' ${1} > qmd/nf_${no_prefix%/*}.qmd
