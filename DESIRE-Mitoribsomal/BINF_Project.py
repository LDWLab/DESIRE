import pandas as pd
import re
import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("protein")
args = parser.parse_args()


df = pd.read_csv('compound_organism_master_table.csv')
exept_df = pd.read_csv('Exeption_List_070823.csv')
check_df = pd.read_csv('Location_Checker_1.csv')
assembly_exept_df = pd.read_csv('Assembly_Exeption_List_240523.csv')
alndata = open(args.filename)
#alndata = open("uL1m.orthologs.mitochondrial.v1.aln")
dataset = alndata.read()
name_aln_list = []
aln_list =[]
organism_list = []
info_list = []
seq_id_list = []
tax_list = []
#protein_name = "uL01m"
protein_name = args.protein
protein_list = []
data_loc = []
protein_loc = []
name_list = []
list_check = []
five_list = []
aln_seq_list = []
seq_exept_list = []
seq_loc = []
mito_check = []
assembly_list = []
assembly_loc_list = []
alt_list = []
frag_list = []
notes_list = []
check_list = []
n_list = []
m_list = []
u_list = []
del_list = []
mito_list = []
no_del_org = []
no_del_header = []
no_del_loc = []
no_del_check = []
no_del_seq = []
mito_exept_list = ['Ddiscoideum','Vermamoeba','Cyanophora','Glaucocystis','Picobiliphyte','Cyanidioschyzon','Galdieria','Pycnococcus','Chlorella','Marchantia','Mesostigma','Diphylleia','Cryptomonas','Goniomonas','Palpitomonas','Trypanosoma','Nfowleri','Ngruberi','Andalucia','Reclinomonas','Emiliania','Gefionella','Thecamonas','Allomyces','Saccharomyces','Ustilago','Paramicrosporidium','Fonticula','Monosiga','Capsaspora','Oxytricha','Paramecium','Bigelowiella','Plasmodiophora','Thalassiosira','Chattonella','Phytophthora']

data = dataset.split('>')
data.remove('')
original_count = len(data)
num = 0

for value in data:
    info = data[num].split("\n")
    info_list.append(info[0])
    aln_seq = data[num].replace("\n", "_-",1).replace("\n","").replace('*','').replace('@','_')
    aln_seq = aln_seq.split('_-')
    seq_exept_list.append(aln_seq[0])
    no_del_header.append(aln_seq[0])
    protein = aln_seq[-1].replace("-", "").replace("\n","").replace('*','').replace('@','_')
    new_protein = aln_seq[0].split('_')
    aln_seq_list.append(aln_seq[-1])
    mito_check.append(new_protein[0])
    if new_protein[0] == 'mitogenome':
        name_aln = [new_protein[2],protein]
        if new_protein[2] in mito_exept_list:
            name_aln = [new_protein[2] + '*',protein]
    else:
        name_aln = [new_protein[0],protein]
    if new_protein[0] == 'mitogenome':
        mito_list.append(new_protein[2])
    else:
        mito_list.append('no')
    name_aln_list.append(name_aln)
    aln_list.append(name_aln[1])
    no_del_seq.append(name_aln[1])
    protein_list.append(protein_name)
    num += 1

for i in range(len(mito_list)):
    if mito_list[i] in mito_exept_list:
        mito_list[i] = mito_list[i] + '*'

_series=check_df[protein_name]

df['tax_id'] = df['tax_id'].dropna().round(0).astype(int)
dict_tax = pd.Series(df.tax_id.values,index=df.header_genome).to_dict()
dict_organism = pd.Series(df.organism.values,index=df.header_genome).to_dict()
dict_seq_ids = pd.Series(df.default_location.values,index=df.header_genome).to_dict()
dict_five = pd.Series(df.five_letter_code.values,index=df.header_genome).to_dict()
dict_assembly_loc = pd.Series(df.assembly_location.values,index=df.header_genome).to_dict()
dict_assembly = pd.Series(df.assembly.values,index=df.header_genome).to_dict()
key_column_name = 'header_genome'
value_column_name = str(protein_name)
data_series = check_df.set_index(key_column_name)[value_column_name]
dict_check = data_series.to_dict()

if "mitogenome" in mito_check:
    df['tax_id'] = df['tax_id'].dropna().round(0).astype(int)
    dict_mito_tax = pd.Series(df.tax_id.values,index=df.header_mitogenome).to_dict()
    dict_mito_organism = pd.Series(df.organism.values,index=df.header_mitogenome).to_dict()
    dict_mito_seq_ids = pd.Series(df.default_location.values,index=df.header_mitogenome).to_dict()
    dict_mito_five = pd.Series(df.five_letter_code.values,index=df.header_mitogenome).to_dict()
    dict_mito_assembly_loc = pd.Series(df.mito_assembly_location.values,index=df.header_mitogenome).to_dict()    
    dict_mito_assembly = pd.Series(df.mito_assembly.values,index=df.header_mitogenome).to_dict()
    key_column_name = 'header_mitogenome'
    value_column_name = str(protein_name)
    data_mito_series = check_df.set_index(key_column_name)[value_column_name]
    dict_mito_check = data_mito_series.to_dict()
    dict_tax = dict_tax | dict_mito_tax
    dict_organism = dict_organism | dict_mito_organism
    dict_seq_ids = dict_seq_ids | dict_mito_seq_ids
    dict_five = dict_five | dict_mito_five
    dict_check = dict_check | dict_mito_check

for genus in name_aln_list:
    if genus[0] in dict_seq_ids:
        seq_id = dict_seq_ids[genus[0]]
        if genus[0] in mito_list:
            seq_id = 'm'
        seq_id_list.append(seq_id)
        no_del_loc.append(seq_id)

for genus in name_aln_list:
    if genus[0] in dict_organism:
        organism = dict_organism[genus[0]]
        organism_list.append(organism.replace("*",""))
        no_del_org.append(organism.replace("*",""))

for genus in name_aln_list:
    if genus[0] in dict_check:
        check = dict_check[genus[0]]
        check_list.append(check)
        no_del_check.append(check)

for genus in name_aln_list:
    if protein_name == 'mS22':
        if genus[0] in ['Tetrahymena','Arabidopsis','Saccharomyces','Chlamydomonas']:
            alt_name = 'mS45'
            alt_list.append(alt_name)
        else:
            alt_list.append(None)
    elif protein_name == 'mS31':
        if genus[0] in ['Saccharomyces']:
            alt_name = 'mS46'
            alt_list.append(alt_name)
        else:
            alt_list.append(None)
    elif protein_name == 'mL54':
        if genus[0] in ['Trypanosoma']:
            alt_name = 'mL88'
            alt_list.append(alt_name)
        else:
            alt_list.append(None)
    elif protein_name == 'mL59':
        if genus[0] in ['Homo']:
            alt_name = 'CRIF1, mL64'
            alt_list.append(alt_name)
        elif genus[0] in ['Trypanosoma']:
            alt_name = 'mL64'
            alt_list.append(alt_name)
        else:
            alt_list.append(None)
    elif protein_name == 'mL60':
        if genus[0] in ['Homo','Trypanosoma']:
            alt_name = 'mL63'
            alt_list.append(alt_name)
        else:
            alt_list.append(None)
    elif protein_name == 'mL61':
        if genus[0] in ['Trypanosoma']:
            alt_name = 'mL74'
            alt_list.append(alt_name)
        else:
            alt_list.append(None)
    else:
        alt_list.append(None)

for genus in name_aln_list:
    if genus[0] in dict_tax:
        tax = dict_tax[genus[0]]
        tax_list.append(str(tax)) 

for genus in name_aln_list:
    if genus[0] in dict_five:
        five_code = dict_five[genus[0]]
        five_list.append(five_code)

for genus in name_aln_list:
    if genus[0] in mito_list:
        assembly = dict_mito_assembly[genus[0]]
        assembly_loc = dict_mito_assembly_loc[genus[0]]
        assembly_list.append(assembly)
        assembly_loc_list.append(assembly_loc)
    elif genus[0] in dict_assembly:
        assembly = dict_assembly[genus[0]]
        assembly_loc = dict_assembly_loc[genus[0]]
        assembly_list.append(assembly)
        assembly_loc_list.append(assembly_loc)
    else:
        assembly_list.append(None)
        assembly_loc_list.append(None)

for value in seq_exept_list:
    if value in ['Arabidopsis_Streptophyta_uL2m','Ntabacum_Streptophyta_XP_016435027.1_60S_ribosomal_protein_L2_mitochondrial-like','Nattenuata_Streptophyta_XP_019265789.1_uncharacterized_protein_LOC109243336']:
        fragment = 'uL2m-C-terminus'
        frag_list.append(fragment)
    elif value in ['Arabidopsis_Streptophyta_DAB41511.2_TPA_asm_ribosomal_protein_L2_mitochondrion','mitogenome_rpl2_Nictabacum_1701','mitogenome_rpl2_Nicattenuata_1558']:
        fragment = 'uL2m-N-terminus'
        frag_list.append(fragment)
    elif value in ['mitogenome_Ymf64_Ichthyopthirius2_AEL89264.1','Tetrahymena_Ciliophora_uS3m','mitogenome_ymf64_Paramecium_1856','mitogenome_rps3_Oxytricha_JN383843.1_AEV66694.1_80','mitogenome_rps3_Acrasis_104']:
        fragment = 'uS3m-C-terminus'
        frag_list.append(fragment)
    elif value in ['mitogenome_rps3_Ichthyopthirius1_NC_015981.1_prot_YP_004841712.1_5','mitogenome_rps3_Ttermophyla_2671','mitogenome_rps3_Paramecium_1845','mitogenome_rps3_Oxytricha_a_JN383843.1_AEV66657.1_42','mitogenome_rps3_Acrasis_112']:
        fragment = 'uS3m-N-terminus'
        frag_list.append(fragment)
    else:
        frag_list.append(None)

for value in seq_exept_list:
    if value == 'Arabidopsis_Streptophyta_AEC06100.1_Nucleic_acid-binding_OB-fold-like_protein':
        notes = 'nucleus-encoded uL2m-2-like peptide (chromosome 2)'
        notes_list.append(notes)
    elif value == 'Arabidopsis_Streptophyta_AEC06107.1_Ribosomal_L5P_family_protein':
        notes = 'nucleus-encoded uL5m-like protein (chromosome 2)'
        notes_list.append(notes)
    elif value == 'Arabidopsis_Streptophyta_AEC06070.1_Ribosomal_protein_S12-S23_family_protein':
        notes = 'nucleus-encoded  uS12m-like protein (chromosome 2)'
        notes_list.append(notes)
    elif value == 'Arabidopsis_Streptophyta_AEC06112.1_Alpha-L_RNA-binding_motif-Ribosomal_protein_S4_family_protein':
        notes = 'nucleus encoded uS4m-like protein (chromoseome 2)'
        notes_list.append(notes)
    else:
        notes_list.append(None)

num_of_n = 0
for value in exept_df['Change_to_N']:
    if value in seq_exept_list:
        num_of_n += 1
        n_list.append(value)
        index_num = seq_exept_list.index(value)
        seq_id_list[index_num] = "n"

num_of_m = 0
for value in exept_df['Change_to_M']:
    if value in seq_exept_list:
        num_of_m += 1
        m_list.append(value)
        index_num = seq_exept_list.index(value)
        seq_id_list[index_num] = "m" 

num_of_u = 0
for value in exept_df['Change_to_U']:
    if value in seq_exept_list:
        num_of_u += 1
        u_list.append(value)
        index_num = seq_exept_list.index(value)
        seq_id_list[index_num] = "u" 

no_del_loc = seq_id_list[:]

num_of_del = 0
for value in exept_df['Delete']:
    if value in seq_exept_list:
        num_of_del += 1
        del_list.append(value)
        index_num = seq_exept_list.index(value)
        seq_exept_list.pop(index_num)
        aln_seq_list.pop(index_num)
        aln_list.pop(index_num)
        protein_list.pop(index_num)
        organism_list.pop(index_num)
        five_list.pop(index_num)
        tax_list.pop(index_num)
        seq_id_list.pop(index_num)
        info_list.pop(index_num)
        name_aln_list.pop(index_num)
        assembly_list.pop(index_num)  
        assembly_loc_list.pop(index_num)
        alt_list.pop(index_num) 
        frag_list.pop(index_num)
        notes_list.pop(index_num)
        check_list.pop(index_num)

for value in assembly_exept_df['GENOME: GCF_000313135.1']:
    if value in seq_exept_list:
        index_num = seq_exept_list.index(value)
        assembly_loc_list[index_num] = "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000313135.1/"    

for value in assembly_exept_df['GENOME: GCF_000001405.38']:
    if value in seq_exept_list:
        index_num = seq_exept_list.index(value)
        assembly_loc_list[index_num] = "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.38/"     

for value in assembly_exept_df['GENOME: GCF_000001735.4']:
    if value in seq_exept_list:
        index_num = seq_exept_list.index(value)
        assembly_loc_list[index_num] = "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/"

for value in assembly_exept_df['GENOME: GCF_000002445.2']:
    if value in seq_exept_list:
        index_num = seq_exept_list.index(value)
        assembly_loc_list[index_num] = "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002445.2/"

for value in assembly_exept_df['GENOME: https://megasun.bch.umontreal.ca/Andalucia_godoyi/Andalucia_godoyi_proteome.faa']:
    if value in seq_exept_list:
        index_num = seq_exept_list.index(value)
        assembly_loc_list[index_num] = "https://megasun.bch.umontreal.ca/Andalucia_godoyi/Andalucia_godoyi_proteome.faa"

for value in assembly_exept_df['GENOME: GCF_000189635.1']:
    if value in seq_exept_list:
        index_num = seq_exept_list.index(value)
        assembly_loc_list[index_num] = "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000189635.1/"

for value in assembly_exept_df['GENOME: GCF_000006565.2']:
    if value in seq_exept_list:
        index_num = seq_exept_list.index(value)
        assembly_loc_list[index_num] = "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000006565.2/"

num = 0
for taxid in tax_list:
    if taxid == '1148':
        tax_list[num] = '1080228'
    if taxid == '269483':
       tax_list[num] = '482957'
    if taxid == '83333':
       tax_list[num] = '511145'
    if taxid == '45157':
        tax_list[num] = '280699'
    if taxid == '1936271':
       tax_list[num] = '1841599'
    if taxid == '999953':
       tax_list[num] = '185431'
    if taxid == '5660':
       tax_list[num] = '420245'
    if taxid == '35128':
       tax_list[num] = '296543'
    if taxid == '44689':
        tax_list[num] = '352472'
    num += 1

num = 0
for organism in organism_list:
    organism = organism.split(" ")
    organism = '_'.join(organism)
    name = organism + "_" + protein_name + "_" + seq_id_list[num]
    num += 1
    name_list.append(name)

for seq in seq_id_list:
    if seq == "u":
        seq_loc.append('unknown')
    elif seq == "n":
        seq_loc.append('nucleus')
    elif seq == 'm':
        seq_loc.append('mitogenome')
    elif seq == 'q':
        seq_loc.append('query')

num = 0
for info in info_list:
    info = info + "_"
    np = re.findall("NP_",info)
    yp = re.findall("YP_",info)
    three = re.findall("_[A-Z]{3}[0-9]", info)
    four = re.findall("_[A-Z]{4}[0-9]", info)
    xp = re.findall("XP_",info)
    mix = re.findall("_[A-Z][0-9][A-Z]",info)
    if xp == ['XP_']:
        count = 1
        xp_loc = re.search("XP_",info)
        span = xp_loc.span()
        string = 'XP_'
        char = info[span[1]]
        while char != '_':
            string = string + char
            char = info[span[1] + count]
            count += 1
        if string[0] == "_":
            string = string.replace("_", "",1)        
        protein_loc.append(string)
        data_loc.append("NCBI")
    elif np == ['NP_']:
        count = 1
        np_loc = re.search("NP_",info)
        span = np_loc.span()
        string = 'NP_'
        char = info[span[1]]
        while char != '_':
            string = string + char
            char = info[span[1] + count]
            count += 1
        if string[0] == "_":
            string = string.replace("_", "",1)        
        protein_loc.append(string)
        data_loc.append("NCBI")
    elif yp == ['YP_']:
        count = 1
        yp_loc = re.search("YP_",info)
        span = yp_loc.span()
        string = 'YP_'
        char = info[span[1]]
        while char != '_':
            string = string + char
            char = info[span[1] + count]
            count += 1
        if string[0] == "_":
            string = string.replace("_", "",1)        
        protein_loc.append(string)
        data_loc.append("NCBI")
    elif three != []:
        count = 1
        three_loc = re.search("_[A-Z]{3}[0-9]",info)
        span = three_loc.span()
        string = three[0]
        char = info[span[1]]
        while char != '_':
            string = string + char
            char = info[span[1] + count]
            count += 1
        if string[0] == "_":
            string = string.replace("_", "",1)
        three_except = re.search('ORF[0-9]', string)
        if three_except != None:
            protein_loc.append(name_list[num])
            data_loc.append("MT-DESIRE")
        else:            
            protein_loc.append(string)
            data_loc.append("UNI")        
    elif four != []:
        count = 1
        four_loc = re.search("_[A-Z]{4}[0-9]",info)
        span = four_loc.span()
        string = four[0]
        char = info[span[1]]
        while char != '_':
            string = string + char
            char = info[span[1] + count]
            count += 1
        if string[0] == "_":
            string = string.replace("_", "",1)        
        protein_loc.append(string)
        data_loc.append("UNI")
    elif mix != []:
        count = 1
        mix_loc = re.search("_[A-Z][0-9][A-Z]",info)
        span = mix_loc.span()
        string = mix[0]
        char = info[span[1]]
        while char != '_':
            string = string + char
            char = info[span[1] + count]
            count += 1
        if string[0] == "_":
            string = string.replace("_", "",1)
        protein_loc.append(string)
        data_loc.append("UNI")        
    else:
        protein_loc.append(name_list[num])
        data_loc.append("MT-DESIRE")
    num += 1

def replaceDuplicates(names):
    hash = {}
    for i in range(0, len(names)):
        if names[i] not in hash:
            hash[names[i]] = 1
        else:
            count1 = hash[names[i]]
            hash[names[i]] += 1
            names[i] += str(count1)
if __name__ == '__main__':
    replaceDuplicates(protein_loc)

longest_string = max(aln_seq_list, key=len)
max_len = len(longest_string)

#print(len(seq_loc))
#print(len(tax_list))
#print(len(name_list))
#print(len(protein_loc))
#print(len(data_loc))
#print(len(protein_list))
#print(len(aln_list))

csv_file = protein_name + ".csv"
with open(csv_file, 'w') as f:
    f.write('Tax_ID\n')
    f.write('\n'.join(tax_list))
    f.close()

df2 = pd.read_csv(csv_file)
df2['Protein_Location'] = protein_loc
df2['Data_Location'] = data_loc
df2['Protein'] = protein_list
df2['Name'] = organism_list
df2['Sequence'] = aln_list
df2['Sequence_Location'] = seq_loc
df2['Assembly_Location'] = assembly_loc_list
df2['Assembly'] = assembly_list
df2['Alternate Name'] = alt_list
df2['Fragment'] = frag_list
df2['Notes'] = notes_list
df2.to_csv(csv_file, index=False)

path1 = "./CSV_v1.5/"
shutil.move(csv_file, path1)

num = 0
aln_file = protein_name + "_aligned.fas"
with open(aln_file, 'w') as f2:
    for genus in name_aln_list:
        if len(aln_seq_list[num]) != max_len:
            aln_seq_list[num] += '-'
        entry = ">" + protein_name + "_" + tax_list[num] + "_" + five_list[num] + "|" + protein_loc[num] + "\n" + aln_seq_list[num] + '\n'
        f2.write(entry)
        num += 1
    f2.close()

path2 = "./Aln_v1.5/"
shutil.move(aln_file, path2)

checking_file = protein_name + '_checking.csv'
with open(checking_file, 'w') as f3:
    f3.write('Organism\n')
    f3.write('\n'.join(no_del_org))
    f3.close

df3 = pd.read_csv(checking_file)
df3['Encoding Location'] = no_del_loc
df3['Master Table Locations'] = no_del_check
df3['Original Header'] = no_del_header
df3['Sequence'] = no_del_seq
df3.to_csv(checking_file, index=False)

path3 = "./CSV_check_v1.5/"
shutil.move(checking_file, path3)

checklist_file = protein_name + '_tally_checking.txt'
with open(checklist_file, 'w') as f4:
    f4.write(protein_name + ' Total Corrections\n')
    f4.write('Original Count: ' + str(original_count) + '\n')
    f4.write('Final Count: ' + str(len(organism_list)) + '\n')
    f4.write('Num of deletions: ' + str(num_of_del) + '\n')
    f4.write('Num of \'n\' corrections: ' + str(num_of_n) + '\n')
    f4.write('Num of \'m\' corrections: ' + str(num_of_m) + '\n')
    f4.write('Num of \'u\' corrections: ' + str(num_of_u) + '\n\n')
    f4.write('Deletion List: \n')
    f4.write('\n'.join(del_list))
    f4.write('\n\n\'n\' List: \n')
    f4.write('\n'.join(n_list))
    f4.write('\n\n\'m\' List: \n')
    f4.write('\n'.join(m_list))
    f4.write('\n\n\'u\' List: \n')
    f4.write('\n'.join(u_list))
    f4.close

path4 = "./Tally_List_v1.5/"
shutil.move(checklist_file, path4)
