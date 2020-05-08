files=(adh.nex bglobin.nex camelid.nex COXI.nex ENCenv.nex flavNS5.nex HepatitisD.nex HIV_RT.nex HIVvif.nex InfluenzaA.nex lysin.nex yokoyama.rh1.cds.mod.1-990.nex)
codes=(Universal Universal Universal Vertebrate-mtDNA Universal Universal Universal Universal Universal Universal Universal Universal)

HYPHY=$1

for i in ${!files[@]}; do
    FILE=${files[i]}
    CODE=${codes[i]}
    time $HYPHY meme --alignment $FILE --code $CODE --output ${FILE}.MEME.json 
done;
