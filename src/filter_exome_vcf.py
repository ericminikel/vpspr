import sys
import gzip



def filter_genotypes(genotypes):
    gt = []
    ad = []
    dp = []
    gq = []
    for col in genotypes:
        parts = col.split(':')
        #print(parts)
        if len(parts) < 4:
            continue
        dp_value = parts[2]
        gq_value = parts[3]
        if dp_value == '.':
            dp_value = '0'
        if gq_value == '.':
            gq_value = '0'
        gt.append(parts[0])
        ad.append(parts[1])
        dp.append(int(dp_value))
        gq.append(int(gq_value))

    high_enough_quality_indices = [i for i in range(len(gq)) if gq[i] >= 20 and dp[i] >= 10]
    return high_enough_quality_indices



with gzip.open('/broad/prions/vpspr/R0626_tgg_minikel_wes_vcf_out_v17_updated.vcf.gz', 'rb') as f:
    header = f.readline().decode('utf-8')
    while header.startswith('##'):
        header = f.readline().decode('utf-8')
    header_cols = header.strip().split('\t')
    sys.stdout.write('\t'.join(header_cols[0:7] + ['CALLSET_AC', 'CALLSET_AN', 'VPSPR_AC', 'VPSPR_AN']) + '\n')
    
    for line in f.readlines():
        cols = line.decode('utf-8').strip().split('\t')
        alts = cols[4].split(',')
        
        genotypes = cols[9:]
        high_quality_indices = filter_genotypes(genotypes)

        
        genotypes_filtered = [genotypes[idx].split(':')[0]  for idx in high_quality_indices]

        
        nestedlist = [genotype.replace('|','/').split('/') for genotype in genotypes_filtered]
        alleles = [item for sublist in nestedlist for item in sublist]



        
        for alt in alts:
            alt_index = alts.index(alt) + 1  
            pos_id = cols[0].replace('chr', '') + '-' + cols[1].zfill(9) + '-' + cols[3] + '-' + alt
            info = cols[7].split(';')
            infodict = dict(x.split('=') for x in info)
            
            callset_ac = infodict['AC'].split(',')[alt_index - 1]
            callset_an = infodict['AN']
            vpspr_ac = alleles.count(str(alt_index))
            vpspr_an = len(alleles)
            
            sys.stdout.write('\t'.join([cols[0]] + [cols[1]] + [pos_id] + [cols[3]] + [alt] + cols[5:7] +
                                        [str(callset_ac), str(callset_an), str(vpspr_ac), str(vpspr_an)]) + '\n')
            sys.stdout.flush()




# python3 filter_exome_vcf.py -> test.tsv




