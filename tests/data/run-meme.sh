files=(adh.nex adora3.nex Bacterial_PTS_trehalose_transporter_suIII.nex bglobin.nex camelid.nex ENCenv.nex Fig4E.nex flavNS5.nex HepatitisD.nex HIV_RT.nex HIVvif.nex IAV-human-H1N1-HA.nex InfluenzaA.nex lysin.nex lysozyme.nex rbcL.nex rbp3.nex REDIC1.nex SARS-CoV-2-spike.nex vwf.nex yokoyama.rh1.cds.mod.1-990.nex COXI.mtnex mammalian_mtDNA.mtnex)
codes=(Universal Universal Universal Universal Universal Universal Universal Universal Universal Universal Universal Universal Universal Universal Universal Universal Universal Universal Universal Universal Universal Vertebrate-mtDNA Vertebrate-mtDNA)

HYPHY=$1

for i in ${!files[@]}; do
    FILE=${files[i]}
    CODE=${codes[i]}
    RESULT=MEME/${FILE}.MEME.json
    if [[ -s $RESULT  ]]
    then
        echo "Already done with $RESULT"
    else
        time $HYPHY meme --alignment $FILE --code $CODE --output MEME/${FILE}.MEME.json 
    fi
done;
