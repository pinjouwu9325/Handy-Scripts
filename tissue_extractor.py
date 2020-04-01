#!/usr/bin/env python

# Author: PJWu
# Last update: 2020-03-31

import argparse
import pandas as pd
from tqdm import tqdm

def run(args):
    input_tissue = args.input # Tissue to extract
    # output_filename = args.output # This mathc the "dest": dest="output"

    # Do stuff

    # Load annotations of sample attributes
    print("Loading GTEx(v8) Sample attributes...")
    df_anno = pd.read_csv("/NAS_lab/pinjouwu9325/gtex/v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep='\t')

    # Load gene_tpm.gct (GTEx v8)
    print("Loading GTEx(v8) GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct")
    df = pd.read_csv("/NAS_lab/pinjouwu9325/gtex/v8/expr/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",skiprows=2, sep='\t')
   
    # Extract sample id of tissue required
    tissue_dict = {'Adipose-Sub':'Adipose - Subcutaneous', 'Adipose-Vis':'Adipose - Visceral (Omentum)', 
                   'AdrenalGland':'Adrenal Gland', 'Artery-Aor':'Artery - Aorta', 'Artery-Cor':'Artery - Coronary',
                   'Artery-Tib':'Artery - Tibial', 'Bladder':'Bladder', 'Brain-Amy':'Brain - Amygdala',
                   'Brain-BA24':'Brain - Anterior cingulate cortex (BA24)', 'Brain-Cau':'Brain - Caudate (basal ganglia)',
                   'Brain-CerH':'Brain - Cerebellar Hemisphere', 'Brain-Cer':'Brain - Cerebellum', 
                   'Brain-Cor':'Brain - Cortex','Brain-BA9':'Brain - Frontal Cortex (BA9)', 'Brain-Hip':'Brain - Hippocampus',
                   'Brain-Hyp':'Brain - Hypothalamus','Brain-NA':'Brain - Nucleus accumbens (basal ganglia)', 
                   'Brain-Put':'Brain - Putamen (basal ganglia)', 'Brain-SC':'Brain - Spinal cord (cervical c-1)', 
                   'Brain-SN':'Brain - Substantia nigra', 'Breast-MT':'Breast - Mammary Tissue',
                   'Cells-EL':'Cells - EBV-transformed lymphocytes', 'Cells-CML':'Cells - Leukemia cell line (CML)',
                   'Cells-TF':'Cells - Transformed fibroblasts', 'Cervix-Ect':'Cervix - Ectocervix', 
                   'Cervix-End':'Cervix - Endocervix','Colon-Sig':'Colon - Sigmoid', 'Colon-Tra':'Colon - Transverse', 
                   'Esophagus-GJ':'Esophagus - Gastroesophageal Junction', 'Esophagus-Muc':'Esophagus - Mucosa', 
                   'Fallopian':'Fallopian Tube', 'Heart-AA':'Heart - Atrial Appendage', 'Heart-LV':'Heart - Left Ventricle', 
                   'Kidney-Cor':'Kidney - Cortex', 'Liver':'Liver', 'Lung':'Lung', 'MSG':'Minor Salivary Gland', 
                   'Muscle-Ske':'Muscle - Skeletal', 'Nerve-Tib':'Nerve - Tibial','Ovary':'Ovary','Pancreas':'Pancreas',
                   'Pituitary':'Pituitary', 'Prostate':'Prostate', 'Skin-NSE':'Skin - Not Sun Exposed (Suprapubic)', 
                   'Skin-SE':'Skin - Sun Exposed (Lower leg)', 'SI':'Small Intestine - Terminal Ileum', 'Spleen':'Spleen', 
                   'Stomach':'Stomach', 'Testis':'Testis', 'Thyroid':'Thyroid', 'Uterus':'Uterus', 'Vagina':'Vagina', 
                   'WB':'Whole Blood'}

    input_tissue = input_tissue.split(",")
    for t in input_tissue:
        print("Filtering " + t + "...")
        # Filter samples of specific tissue
        df_tissue_anno = df_anno[df_anno.SMTSD == tissue_dict[t]]
        id_list = df_tissue_anno.SAMPID.tolist()

        # Filter gene expression profile by sample ID
        filter_list = ['Description', 'Name']
        filter_list.extend(id_list)
        df_tissue_data = df.filter(items=filter_list)

        # Filter out the genes with any gene expression of a sample is zero
        print("Filtering out the genes with any gene expression of a sample is zero(TPM=0)...") 
        index_list=[]
        for i in tqdm(range(df_tissue_data.shape[0])):
            count = 0
            for j in range(2, df_tissue_data.shape[1]):
                if df_tissue_data.iloc[i, j] == 0:
                    count+=1
            if count==0:
                index_list.append(i)
        df_tissue_data_clean = df_tissue_data.iloc[index_list,:]

    # # Export data
    # if bool(output_filename) == True:
    #     name = str(output_filename)
    # else:
        # name = str(input_tissue + ".csv")
    filename = str(t+".gct")
    df_tissue_data_clean.to_csv(filename, sep=',', encoding = 'utf-8', index=0)
    print(filename + " has been exported.")


def main():
    parser=argparse.ArgumentParser(
    prog='TissueExtractor',
    formatter_class=argparse.RawDescriptionHelpFormatter, 
    description='''\
            This will do extraction of the specific tissue data from GTEx tpm dataset(v8).
            Any gene with a read count of zero(TPM=0) in any sample will be excluded. The output file will be a csv file.
            
            Use the tissue key for input. Please note the capitalized case in each tissue key.
            [Tissue Key]:description
            - A -
            [Adipose-Sub]:Adipose-Subcutaneous           [Adipose-Vis]:Adipose-Visceral (Omentum)            [AdrenalGland]:Adrenal Gland
            [Artery-Aor]:Artery-Aorta                    [Artery-Cor]:Artery-Coronary                        [Artery-Tib]:Artery-Tibial
            - B -
            [Bladder]:Bladder                            [Brain-Amy]:Brain-Amygdala                          [Brain-BA24]:Brain-Anterior cingulate cortex (BA24)
            [Brain-Cau]:Brain-Caudate (basal ganglia)    [Brain-CerH]:Brain-Cerebellar Hemisphere            [Brain-Cer]:Brain-Cerebellum
            [Brain-Cor]:Brain-Cortex                     [Brain-BA9]:Brain-Frontal Cortex (BA9)              [Brain-Hip]:Brain-Hippocampus 
            [Brain-Hyp]:Brain-Hypothalamus               [Brain-NA]:Brain-Nucleus accumbens (basal ganglia)
            [Brain-Put]:Brain-Putamen (basal ganglia)    [Brain-SC]:Brain-Spinal cord (cervical c-1)         [Brain-SN]:Brain-Substantia nigra
            [Breast-MT]:Breast-Mammary Tissue
            - C -
            [Cells-EL]:Cells-EBV-transformed lymphocytes [Cells-CML]:Cells-Leukemia cell line (CML)          [Cells-TF]:Cells-Transformed fibroblasts 
            [Cervix-Ect]:Cervix-Ectocervix               [Cervix-End]:Cervix-Endocervix                      [Colon-Sig]:Colon-Sigmoid 
            [Colon-Tra]:Colon-Transverse               
            -E-H-
            [Esophagus-GJ]:Esophagus-Gastroesophageal Junction                                               [Esophagus-Muc]:Esophagus-Mucosa 
            [Fallopian]:Fallopian Tube                   [Heart-AA]:Heart-Atrial Appendage                   [Heart-LV]:Heart-Left Ventricle
            -K-M-
            [Kidney-Cor]:Kidney-Cortex                   [Liver]:Liver                                       [Lung]:Lung
            [MSG]:Minor Salivary Gland                   [Muscle-Ske]:Muscle-Skeletal                  
            -N-P-
            [Nerve-Tib]:Nerve-Tibial                     [Ovary]:Ovary                                       [Pancreas]:Pancreas 
            [Pituitary]:Pituitary                        [Prostate]:Prostate
            - S -
            [Skin-NSE]:Skin-Not Sun Exposed (Suprapubic) [Skin-SE]:Skin-Sun Exposed (Lower leg)              [SI]:Small Intestine-Terminal Ileum
            [Spleen]:Spleen                              [Stomach]:Stomach
            -T-W-
            [Testis]:Testis                              [Thyroid]:Thyroid                                   [Uterus]:Uterus
            [Vagina]:Vagina                              [WB]:Whole Blood
                   ''')

    parser.add_argument("-i", help="Input a tissue key for extraction or a common seperate list of tissue key", dest="input", type=str, required=True)
    # parser.add_argument("-o", help="ouput filename, default is set to input tissue key", dest="output", type=str) # default=Tissuename.csv
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()

