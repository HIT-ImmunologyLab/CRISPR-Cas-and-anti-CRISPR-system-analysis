#!/usr/bin
#-*-coding:utf-8-*-

import os
from Bio import SeqIO
import xml.etree.ElementTree as ET
import re
import threading
import argparse


def mkdir(dirPath):
    cmd_mkdir = 'mkdir -p %s'%dirPath
    os.system(cmd_mkdir)

def identify_self_targetging(crispr_dir):
    print('identify CRISPR self-targeting')
    cmd_crispr_prediction = 'python /home/pipeline/identify_self_targeting/identify_selftargeting_version_2.py --input_dir %s'%crispr_dir
    os.system(cmd_crispr_prediction)

def identify_prophage(bac_path,outdir,prefix):
    cmd_dbscan_swa = 'python /data/ganrui/DBSCAN-SWA/bin/dbscan-swa.py --input %s --output %s --prefix %s'%(bac_path,outdir,prefix)
    os.system(cmd_dbscan_swa)

def diamond_blastp(file,outfile,database,format,evalue):
    num_threads = 20
    diamond_path = os.path.join('/data/ganrui/DBSCAN-SWA','software','diamond','diamond')
    script = diamond_path+" blastp -d "+database+" -q "+file+" -f "+str(format)+" -e "+str(evalue)+" -o "+outfile+" -p "+str(num_threads)+" --max-target-seqs 1"
    print(script)
    os.system(script)

def identify_phage_prot(bac_protein_file,outdir,prot_id,blastp_evalue):
    print('identify phage or phage-like genes')
    blastp_file = '%s/%s_blastp_uniprot.txt'%(outdir,prot_id)
    database = os.path.join('/data/ganrui/DBSCAN-SWA','db','database','uniprot.dmnd')
    diamond_blastp(bac_protein_file,blastp_file,database,6,str(blastp_evalue))

def cp_file(source_file,dest_dir):
    cmd_cp = 'cp %s %s'%(source_file,dest_dir)
    os.system(cmd_cp)

def getFnaFromGB(gb_file_path,fa_file_path):
    handle = open(gb_file_path)
    if os.path.exists(fa_file_path):
        os.remove(fa_file_path)
    SeqIO.convert(handle, 'genbank', fa_file_path, 'fasta')

def getFaaFromGB(gb_file_path,faa_file_path):
    records = SeqIO.parse(gb_file_path, "gb")
    if os.path.exists(faa_file_path):
        os.remove(faa_file_path)
    savefile = open(faa_file_path, 'w')
    for record in records:
        fileID = record.id
        prot_num = 0
        for feature in record.features:
            if feature.type == 'CDS':
                prot_num = prot_num + 1
                location = feature.location
                if str(location).find('+') != -1:
                    direction = '+'
                elif str(location).find('-') != -1:
                    direction = '-1'
                if feature.type == 'CDS':
                    if 'product' in feature.qualifiers:
                        product = feature.qualifiers['product'][0]
                        if ' 'in product:
                            product = product.replace(' ','_')
                    else:
                        product = 'unkown'
                    if 'protein_id' in feature.qualifiers:
                        proteinId = feature.qualifiers['protein_id'][0]
                    else:
                        if 'inference' in feature.qualifiers:
                            strInference = str(feature.qualifiers['inference'])
                            if 'RefSeq' in strInference:
                                proteinId = strInference.split('RefSeq:')[1].rstrip(']').rstrip('\'')
                            elif 'SwissProt' in strInference:
                                proteinId = strInference.split('SwissProt:')[1].rstrip(']').rstrip('\'')
                            else:
                                proteinId = 'unknown'
                        else:
                            proteinId = 'unknown'
                    if 'translation' in feature.qualifiers:
                        translation = feature.qualifiers['translation'][0]
                    else:
                        translation = 'unkown'
                    savefile.write('>ref' + '|'+fileID+'|'+proteinId+'|'+str(location)+'|'+str(product)+ '|'+str(prot_num)+'\n')
                    if translation[-1] == '\n':
                        savefile.write(translation)
                    else:
                        savefile.write(translation + '\n')
    savefile.close()

def get_prot_file(genome_faa_file,prot_id,prot_file):
    fout = open(prot_file,'w')
    with open(genome_faa_file,'r')as fin:
        content = fin.read()
        elems = content.split('\n>')
        if '' in elems:
            elems.remove('')
        for elem in elems:
            cur_def_line = elem.split('\n')[0]
            if prot_id in cur_def_line:
                if cur_def_line.startswith('>'):
                    fout.write(elem+'\n')
                else:
                    fout.write('>'+elem+'\n')
                cur_acr_index = cur_def_line.split('|')[-1]
                return cur_acr_index

def get_downstream_prot_file(genome_faa_file,acr_index,downstream_num,downstrem_prot_file):
    down_stream_index_list = []
    for i in range(0,downstream_num):
        cur_index = int(acr_index)+i+1
        down_stream_index_list.append(str(cur_index))
    fout = open(downstrem_prot_file,'w')
    write_count = 0
    with open(genome_faa_file,'r')as fin:
        content = fin.read()
        elems = content.split('\n>')
        if '' in elems:
            elems.remove('')
        for elem in elems:
            cur_defline = elem.split('\n')[0]
            cur_prot_index = cur_defline.split('|')[-1]
            if cur_prot_index in down_stream_index_list:
                if cur_defline.startswith('>'):
                    fout.write(elem + '\n')
                else:
                    fout.write('>'+elem+'\n')
                write_count = write_count + 1
                if write_count == downstream_num:
                    break
    fout.close()

def rpsblastp(fileName,outFileName):
    blast_cline = "rpsblast -query %s -comp_based_stats 0 -max_target_seqs 1 -evalue 0.01 -seg no -outfmt 5 -num_threads 20 -db /data/DataBase/cdd/Cdd -out %s" % (fileName,outFileName)
    os.system(blast_cline)

def ParseResult(inFileName, outFileName):
    text = open(inFileName).read()
    text = re.sub(u"[\x00-\x08\x0b-\x0c\x0e-\x1f]+", u" ", text)
    # print(text)
    root = ET.fromstring(text)
    BlastOutput_iterations = root.find("BlastOutput_iterations")
    f = open(outFileName, 'w')
    for Iteration in BlastOutput_iterations.findall("Iteration"):
        strTemp = str(Iteration.find("Iteration_query-def").text)
        # Anti_ID = strTemp.split('|')[1]
        # Anti_Location = strTemp.split('|')[2]
        # query_from =re.findall('\d+',Anti_Location)[0]
        # query_to =re.findall('\d+',Anti_Location)[1]
        # Iteration_query_len = Iteration.find('Iteration_query-len').text
        Iteration_hits = Iteration.find("Iteration_hits")
        # f.write(strTemp+'\n')
        # f.flush()
        for Hit in Iteration_hits.findall("Hit"):
            strDef = str(Hit.find("Hit_def").text)
            Hit_len = Hit.find("Hit_len").text
            # Hit_ID = strDef.split('|')[1].split(' ')[0]
            Hit_ID = Hit.find("Hit_id").text
            # Hit_location= strDef.split('|')[2]
            # Hit_from = re.findall('\d+',Hit_location)[0]
            # Hit_to = re.findall('\d+',Hit_location)[1]
            Hsp = Hit.find("Hit_hsps").find("Hsp")
            Hit_from = Hsp.find("Hsp_hit-from").text
            Hit_to = Hsp.find("Hsp_hit-to").text
            Hsp_evalue = Hsp.find("Hsp_evalue").text
            Hsp_score = Hsp.find("Hsp_score").text
            hit_length = int(Hit_from) - int(Hit_from) + 1
            coverage = str(float(hit_length) / float(Hit_len))
            identity = str(Hsp.find("Hsp_identity").text)
            f.write(
                strTemp + '\t' + Hit_ID + '\t' + strDef + '\t' + Hsp_evalue + '\t' + identity + '\t' + coverage + '\n')
            f.flush()
        # f.write('\n')
        # f.flush()
    f.close()

def identify_hth_domain(acr_downstrem_prot_file,domain_result_dir,domain_annotation_file_prefix):
    doamin_annotation_file_xml = '%s/%s.xml'%(domain_result_dir,domain_annotation_file_prefix)
    doamin_annotation_file_txt = '%s/%s.txt' % (domain_result_dir, domain_annotation_file_prefix)
    rpsblastp(acr_downstrem_prot_file, doamin_annotation_file_xml)
    ParseResult(doamin_annotation_file_xml, doamin_annotation_file_txt)

def analysis_anti_features(bac_save_dir,result_dir,bac_id,prot_id,downstream_num=3):
    gb_file_path = '%s/%s.gb'%(bac_save_dir,bac_id)
    fa_file_path = '%s/%s.fa'%(bac_save_dir,bac_id)
    faa_file_path = '%s/%s.faa'%(bac_save_dir,bac_id)
    getFnaFromGB(gb_file_path, fa_file_path)
    getFaaFromGB(gb_file_path, faa_file_path)
    prot_file = '%s/%s.faa'%(bac_save_dir,prot_id)
    acr_index = get_prot_file(faa_file_path, prot_id, prot_file)
    downstrem_prot_file = '%s/%s_downstream_%s.faa'%(bac_save_dir,prot_id,downstream_num)
    get_downstream_prot_file(faa_file_path, acr_index, downstream_num, downstrem_prot_file)

    # print('identify prophage')
    prophage_result = '%s/prophage_prediction'%result_dir
    mkdir(prophage_result)
    prophage_predix = bac_id
    # identify_prophage(gb_file_path, prophage_result, prophage_predix)
    m1 = threading.Thread(target=identify_prophage, args=(gb_file_path, prophage_result, prophage_predix,))
    m1.start()

    ## identify CRISPR  self-targeting
    crispr_result = '%s/crispr_prediction/%s'%(result_dir,bac_id)
    mkdir(crispr_result)
    cp_file(fa_file_path, crispr_result)
    # identify_self_targetging(crispr_result)
    m2 = threading.Thread(target=identify_self_targetging, args=(crispr_result,))
    m2.start()

    ## identify phage prot
    phage_result = '%s/phage_prediction'%(result_dir)
    mkdir(phage_result)
    blastp_evalue = '1e-7'
    # identify_phage_prot(prot_file, phage_result, prot_id,blastp_evalue)
    m3 = threading.Thread(target=identify_phage_prot, args=(prot_file, phage_result, prot_id,blastp_evalue,))
    m3.start()

    ## identify HTH domain
    domain_result_dir = '%s/acr_downstream_prot_domain_anno'%(result_dir)
    mkdir(domain_result_dir)
    domain_annotation_file_prefix = '%s_downstream_%d_prot_domain_anno'%(prot_id,downstream_num)
    # identify_hth_domain(downstrem_prot_file, domain_result_dir, domain_annotation_file_prefix)
    m4 = threading.Thread(target=identify_hth_domain,
                          args=(downstrem_prot_file, domain_result_dir, domain_annotation_file_prefix,))
    m4.start()

    # m1.join()
    # m2.join()
    m3.join()
    m4.join()

def get_prot_def(prot_file):
    with open(prot_file,'r')as fin:
        lines = fin.readlines()
        cur_prot_position = lines[0].split('|')[-3]
        cur_prot_def = lines[0].split('|')[-2]
        return [cur_prot_def,cur_prot_position]

def get_genome_def(genome_gb_file):
    records = SeqIO.parse(genome_gb_file, "gb")
    for record in records:
        cur_genome_def = record.description
        return cur_genome_def

def get_selftargeting_info(self_target_file):
    crispr_array_location_list = []
    crispr_type_list = []
    self_target_spacer_num_list = []
    self_target_spacer_pos_list = []
    protospacer_num_list = []
    protospacer_pos_list = []
    with open(self_target_file,'r')as fin:
        lines = fin.readlines()
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            cur_crispr_array_loation = content[4]
            cur_crispr_type = content[10]
            cur_self_target_num = content[14]
            self_target_spacer_pos = content[15]
            cur_protospacer_num = content[16]
            protospacer_pos = content[17]
            crispr_array_location_list.append(cur_crispr_array_loation)
            crispr_type_list.append(cur_crispr_type)
            self_target_spacer_num_list.append(cur_self_target_num)
            self_target_spacer_pos_list.append(self_target_spacer_pos)
            protospacer_num_list.append(cur_protospacer_num)
            protospacer_pos_list.append(protospacer_pos)

    return ['|'.join(crispr_array_location_list),
            '|'.join(crispr_type_list),
            '|'.join(self_target_spacer_num_list),
            '|'.join(self_target_spacer_pos_list),
            '|'.join(protospacer_num_list),
            '|'.join(protospacer_pos_list)]

def get_prophage_region(prophage_summary_file):
    prophager_region_list = []
    with open(prophage_summary_file,'r')as fin:
        lines = fin.readlines()
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            cur_genome_region = content[0]+'_'+content[3]+'-'+content[4]
            prophager_region_list.append(cur_genome_region)
    return '|'.join(prophager_region_list)

def check_acr_in_prophage(prophage_summary_file,genome_id,acr_location):
    acr_in_prophage = '0'
    acr_prot_pos_start = min(int(acr_location.split('[')[1].split(':')[0]),int(acr_location.split(':')[1].split(']')[0]))
    acr_prot_pos_end = max(int(acr_location.split('[')[1].split(':')[0]), int(acr_location.split(':')[1].split(']')[0]))
    with open(prophage_summary_file,'r')as fin:
        lines = fin.readlines()
        for line in lines[1:]:
            content = line.strip('\n').split('\t')
            cur_genome_id = content[0]
            cur_prophage_start = int(content[3])
            cur_prophage_end = int(content[4])
            if genome_id in cur_genome_id or cur_genome_id in genome_id:
                if cur_prophage_start <= acr_prot_pos_start and cur_prophage_end >= acr_prot_pos_end:
                    acr_in_prophage = '1'
    return acr_in_prophage

def check_acr_is_phage(acr_phage_anno_file):
    phage_flag = '0'
    if os.path.getsize(acr_phage_anno_file)>0:
        phage_flag = '1'
    return phage_flag

def check_HTH_domain(hth_domain_file):
    HTH_domain_flag = '0'
    HTH_domain_keys_list = ['Helix-Turn-Helix','HTH']
    with open(hth_domain_file,'r')as fin:
        lines = fin.readlines()
        for line in lines:
            content = line.strip('\n').split('\t')
            cur_domain_def = content[3].lower()
            for HTH_domain_keys_item in HTH_domain_keys_list:
                if HTH_domain_keys_item.lower() in cur_domain_def:
                    HTH_domain_flag = '1'
    return HTH_domain_flag

def statis_anti_features(bac_save_dir, result_dir, bac_id,prot_id,anti_feature_status_file,downstrem_num=3):
    fout = open(anti_feature_status_file,'w')
    header_list = ['acr_prot_id','acr_prot_def','acr_prot_location','genome_id','genome_def',
                   'crispr_array_location','crispr_type_in_genome','self-targeting_spacer_num',
                   'self-targeting_spacer_pos','protospacer_num','protospacer_pos',
                   'prophage_region','acr_in_prophage_region','acr_is_phage_prot',
                   'acr_downstream_has_HTH_domain']
    fout.write('\t'.join(header_list)+'\n')
    prot_file = '%s/%s.faa' % (bac_save_dir, prot_id)
    [acr_prot_def, acr_prot_position] = get_prot_def(prot_file)
    gb_file_path = '%s/%s.gb' % (bac_save_dir, bac_id)
    genome_def = get_genome_def(gb_file_path)

    self_target_info_file = '%s/crispr_prediction/%s/%s_st_result.txt'%(result_dir,bac_id,bac_id)
    [crispr_array_location,crispr_type_in_genome,
     targeting_spacer_num,targeting_spacer_pos,
     protospacer_num,protospacer_poss] = get_selftargeting_info(self_target_info_file)

    prophae_summary_info_file = '%s/prophage_prediction/%s_DBSCAN-SWA_prophage_summary.txt'%(result_dir,bac_id)
    prophage_region = get_prophage_region(prophae_summary_info_file)

    acr_in_prophage = check_acr_in_prophage(prophae_summary_info_file, bac_id, acr_prot_position)

    acr_phage_anno_file = '%s/phage_prediction/%s_blastp_uniprot.txt'%(result_dir,prot_id)
    acr_phage_flag = check_acr_is_phage(acr_phage_anno_file)

    hth_domain_file = '%s/acr_downstream_prot_domain_anno/%s_downstream_%d_prot_domain_anno.txt'%(result_dir,prot_id,downstrem_num)
    hth_flag = check_HTH_domain(hth_domain_file)

    strwrite_list = [prot_id]+[acr_prot_def, acr_prot_position]+[bac_id,genome_def]+\
                    [crispr_array_location, crispr_type_in_genome,targeting_spacer_num, targeting_spacer_pos,
     protospacer_num, protospacer_poss]+[prophage_region,acr_in_prophage,acr_phage_flag,hth_flag]
    fout.write('\t'.join(strwrite_list)+'\n')
    fout.close()

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', help='absPath of directory containing.\n')
    parser.add_argument('--output_dir', help='absPath of result directory.\n')
    parser.add_argument('--bac_id', help='ID of bacteria.\n', )
    parser.add_argument('--acr_id', help='ID of acr protein.\n')
    parser.add_argument('--acr_downstream_prot_num', help='Protein number downstream of acr protein.\n', default='3')
    args = parser.parse_args()
    if args.input_dir:
        bac_save_dir = args.input_dir
    if args.output_dir:
        result_dir = args.output_dir
    if args.bac_id:
        bac_id = args.bac_id
    if args.acr_id:
        acr_prot_id = args.acr_id
    if args.acr_downstream_prot_num:
        acr_downstream_prot_num = args.acr_downstream_prot_num

    # bac_save_dir = '/data/zfx/anti-crispr_analysis/test/CP005941'
    # result_dir = '/data/zfx/anti-crispr_analysis/test/CP005941'
    # bac_id = 'CP005941'
    # acr_prot_id = 'AGM99339'
    # acr_downstream_prot_num = 3
    analysis_anti_features(bac_save_dir, result_dir, bac_id,acr_prot_id,acr_downstream_prot_num)
    anti_feature_status_file = '%s/anti_feature_summary.txt'%result_dir
    statis_anti_features(bac_save_dir, result_dir, bac_id, acr_prot_id, anti_feature_status_file,acr_downstream_prot_num)

