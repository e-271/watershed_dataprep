
#!/bin/bash

DATA_DIR="/oak/stanford/groups/smontgom/mdegorte/durga/africa/30x/genotypes"

for pop in ESN GWD LWK MSL YRI AFR YRI.GEU
do
 bcftools query -f '%CHROM %POS %ID %REF %ALT %AN %AC{0}\n' all.${pop}.30x.ID.vcf.gz \
  | awk '{ if ($6 != 0) printf "%s %f\n",$3, $7/$6}' > frq/all.${pop}.30x.frq


done


