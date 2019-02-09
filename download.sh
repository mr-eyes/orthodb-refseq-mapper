# Directories to store data
mkdir -p species/{monkey,horse,cow,elephant,rabbit,dog,human,platypus,mouse}/{refseq,orthodb}
mkdir data

#1 Monkey : Macaca mulatta
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Macaca_mulatta/latest_assembly_versions/GCF_000772875.2_Mmul_8.0.1/GCF_000772875.2_Mmul_8.0.1_rna.fna.gz -O macaca_mulatta_rna.fna.gz
mv macaca_mulatta_rna.fna.gz species/monkey/refseq/

#2 Horse : Equus caballus
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_rna.fna.gz -O equus_caballus_rna.fna.gz
mv equus_caballus_rna.fna.gz species/horse/refseq/

#3 Cow : Bos taurus
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_rna.fna.gz -O bos_taurus.rna.fna.gz
mv bos_taurus.rna.fna.gz species/cow/refseq/

#4 Elephant : Loxodonta_africana
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Loxodonta_africana/latest_assembly_versions/GCF_000001905.1_Loxafr3.0/GCF_000001905.1_Loxafr3.0_rna.fna.gz -O loxodonta_africana.rna.fna.gz
mv loxodonta_africana.rna.fna.gz species/elephant/refseq/

#5 Rabbit : Oryctolagus_cuniculus
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Oryctolagus_cuniculus/latest_assembly_versions/GCF_000003625.3_OryCun2.0/GCF_000003625.3_OryCun2.0_rna.fna.gz -O oryctolagus_cuniculus.rna.fna.gz
mv oryctolagus_cuniculus.rna.fna.gz species/rabbit/refseq/

#6 DOG : Canis Lapus
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Canis_lupus/latest_assembly_versions/GCF_000002285.3_CanFam3.1/GCF_000002285.3_CanFam3.1_rna.fna.gz -O canis_lapus.fna.rna.gz
mv canis_lapus.fna.rna.gz species/dog/refseq/

#7 Human : Homo Sapien
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_rna.fna.gz -O homo_sapien.rna.fna.gz
mv homo_sapien.rna.fna.gz species/human/refseq/

#8 Platypus: Ornithorhynchus_anatinus
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Ornithorhynchus_anatinus/latest_assembly_versions/GCF_000002275.2_Ornithorhynchus_anatinus_5.0.1/GCF_000002275.2_Ornithorhynchus_anatinus_5.0.1_rna.fna.gz -O ornithorhynchus_anatinus.rna.fna.gz
mv ornithorhynchus_anatinus.rna.fna.gz species/platypus/refseq/

#9 Mouse : Mus musculus
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_rna.fna.gz -O mus_musculus_rna.fna.gz
mv mus_musculus_rna.fna.gz species/mouse/refseq/


wget -c https://v100.orthodb.org/download/odb10v0_OG2genes.tab.gz && mv odb10v0_OG2genes.tab.gz data/
wget -c https://v100.orthodb.org/download/odb10v0_genes.tab.gz && mv odb10v0_genes.tab.gz data/
wget -c ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz && mv gene2refseq.gz data/

########################################################################################################################

# Getting the data of each species

names=("cow" "dog" "elephant" "horse" "human" "monkey" "platypus" "rabbit" "mouse")  ## Dirs Names
taxas=("9913" "9615" "9785" "9796" "9606" "9544" "9258" "9986" "10090") ## Taxonomy IDs

len=${#names[*]}  # it returns the array length

#iterate with while loop
i=0
while [ $i -lt $len ]
do
    NAME=${names[$i]}
    echo "Processing ${NAME}"
    TAXA=${taxas[$i]}
    zcat data/odb10v0_genes.tab | grep "\s${TAXA}_0\s" > species/${NAME}/orthodb/${NAME}_odb10v0_genes.tab
    zcat data/odb10v0_OG2genes.tab.gz | grep "\s${TAXA}_0:" > species/${NAME}/orthodb/${NAME}_odb10v0_OG2genes.tab
    zcat data/gene2refseq.gz | grep -w "^${TAXA}" | awk '{print $2,$4}' > species/${NAME}/${NAME}_gene2refseq.txt
    i=$((i+1))
done

