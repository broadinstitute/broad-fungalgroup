import pandas as pd
import numpy as np

#TODO:Argparse
input_gemma_file=''
ortholog_data_file='data/reference/C_auris_B8441_C_albicans_SC5314_orthologs.txt'
data_namegene='data/reference/C_auris_B8441_current_chromosomal_feature.tab'

#Load the C. Albicans orthology equivalent
albicans_gene_to_caID={}
caID_to_albicans_gene={}

with open(ortholog_data_file, 'r') as orthologs_file:
    for line in orthologs_file:
        #Ignore comments
        if line[0]=='#':continue
        spl=line.split('\t')

        #Check that the string is not empty.
        if spl[4]:
            albicans_gene_to_caID[spl[4]]=spl[0]
            caID_to_albicans_gene[spl[0]]=spl[4]

#Load the gene names and S.Cerevisae ortholog
caID_to_cagene_name={}
caID_to_scervi_name={}

with open(data_namegene, 'r') as names_file:
    for line in names_file:
        #Ignore comments
        if line[0]=='!':continue
        spl=line.split('\t')
        #Check that the gene name is not empty.
        if spl[1]:
            caID_to_cagene_name[spl[0]]=spl[1]

        #Check that the S.Cerevisae ortholog is not empty.
        if spl[17].strip():
            caID_to_scervi_name[spl[0]]=spl[17].strip()

output_file=input_gemma_file.split('.')[0]+'.annortho.csv'
output_file=open(output_file, 'w')

print('This code only annotates the top effect reported by SNPEff. Rememeber to check secondary effects in case of doubts.')

with open(input_gemma_file, 'r') as ann_file:
    count=0

    for line in ann_file:
        #Make header
        if count==0:
            parse=line.split('\t')
            new_line=','.join(parse[:-1]+['Plot_name','Albicans_Ortholog','SCerev_Ortholog','ALT', 'Annotation', 'Putative_impact', 'Gene_Name', 'Gene_ID', 'Feature_type', 'Feature_ID', 'Transc_Biotype', 'Rank_IntEx', 'HGVSc', 'HGVSp', 'cDNApos_len', 'CDSpos_len', 'Proteinpos_len', 'DistanceFeature', 'Messages','Other_Effects'])+'\n'
            output_file.write(new_line)
            count=count+1
            continue

        #Grab annotation column
        annot=line.strip().split('\t')[-1]

        #Skip lines with no annotation
        if annot=='NONE':
            scervi_ortho='-'
            albicans_ortho='-'
            plot_name='No_Annotation'
            continue

        top_effect=annot.split(',')[0]
        curr_id=top_effect.split('|')[4]
        curr_name=top_effect.split('|')[3]

        #print(curr_id, curr_name)
        plot_name=''

        #Things listed in order of plot name importance
        #SCerevisae Ortholog:
        if curr_id in caID_to_scervi_name:
            scervi_ortho=caID_to_scervi_name[curr_id].replace(',', '_') #The replacement makes sure that any commas in gene names don't mess up the parsing later.
            plot_name=scervi_ortho+'_SCerO'
        else:
            scervi_ortho='-'

        #CAlbicans Ortholog
        if curr_id in caID_to_albicans_gene:
            albicans_ortho=caID_to_albicans_gene[curr_id].replace(',', '_')
            plot_name=albicans_ortho+'_CAlbO'
        else:
            albicans_ortho='-'

        #Check if there is a protein name
        if curr_name!= 'hypothetical_protein':
            plot_name=curr_name

        #Check if there is a gene name
        if curr_id in caID_to_cagene_name:
            plot_name=caID_to_cagene_name[curr_id].replace(',', '_')

        #Give the id if no plot name
        if not plot_name:
            plot_name=curr_id

        #Make the output line
        significance=line.split('\t')[:-1]
        output_line=','.join(significance+[plot_name, albicans_ortho, scervi_ortho])+','+top_effect.replace('|', ',')+','+annot.replace(',', '-')+'\n'
        output_file.write(output_line)
