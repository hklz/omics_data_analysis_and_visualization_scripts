#导入相关包
import pandas as pd
import openpyxl
import regex as re_lib
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# sigma因子启动子识别motif
patterns_sigma = {'RPOD':"TTGACAN{17}TATAAT",
                  'RPOS':"TTGACAN{17}TATACT",
                  'RPOE':"GAACTN{17}TCTGAT",
                  'RPOF':"TAAAGTN{15}GCCGATAA",
                  'RPOH':"TTGAAAN{15}CCCCATN{1}T",
                  'FECI':"AAAATN{17}TTGTN{1}T",
                  'RPON':"TGGCAGGN{4}TTGCA"}

#给个基因实体
class Gene:
    def __init__(self, gene_id, name, start, end, sequence, upstream_seq, promoter, predictied_right,Predicted_promoter, sigma_factor, max_sigma_factor, operon_genes, direction, sigma_TPM=None):
        self.gene_id = gene_id
        self.name = name
        self.start = start
        self.end = end
        self.sequence = sequence
        self.upstream = upstream_seq if upstream_seq is not None else ''
        self.promoter = promoter if promoter is not None else []
        self.Predicted_promoter = Predicted_promoter if Predicted_promoter is not None else []
        self.sigma_factor = sigma_factor # ecocyc预测的sigma
        self.max_sigma_factor = max_sigma_factor if max_sigma_factor is not None else []
        self.operon_genes = operon_genes
        self.sigma_TPM = sigma_TPM if sigma_TPM is not None else {}
        self.direction = direction
        self.predicted = predictied_right if predictied_right is not None else False

    def add_TPM(self, sigma_factor, TPM_value):
        self.sigma_TPM[sigma_factor] = TPM_value

    def __repr__(self):
        # 创建一个易于阅读的表示
        return (f"Gene(\n"
                f"  gene_id={self.gene_id},\n"
                f"  name='{self.name}',\n"
                f"  start={self.start},\n"
                f"  end={self.end},\n"
                f"  sequence='{self.sequence}',\n"
                f"  promoter='{self.promoter}',\n"
                f"  sigma_factor='{self.sigma_factor}',\n"
                f"  max_sigma_factor='{self.max_sigma_factor}',\n"
                f"  operon_genes={self.operon_genes},\n"
                f"  sigma_rpkm={self.sigma_TPM},\n"
                + (f"  max_rpkm={self.sigma_TPM[self.max_sigma_factor]},\n" if (self.max_sigma_factor != 'N.S.') else "")
                + f"  temp_rpkm={self.sigma_TPM['temp']},\n"
                  f"  direction={self.direction}\n"
                  f")")


# motif通用字符转换
def iupac_to_regex(pattern: str) -> str:
    iupac_codes = {
        'A': 'A',
        'T': 'T',
        'C': 'C',
        'G': 'G',
        'N': '[ATCG]',
        'Y': '[CT]',
        'R': '[AG]',
        'W': '[AT]',
        'S': '[CG]',
        'M': '[AC]',
        'K': '[GT]',
        'H': '[ACT]',
        'B': '[CGT]',
        'V': '[ACG]',
        'D': '[AGT]'
    }

    regex_pattern = ''.join([iupac_codes.get(char, char) for char in pattern])
    return regex_pattern

#搜索给定序列是否存在某个模体并输出该模体及分数
def find_exact_and_most_similar_patterns_weighted(sequence, motif_pattern, position='on',max_errors=2, exact_match_weight=2):
    """
    Find the presence of an exact sigma factor recognition motif pattern and the most similar sequence fragments in a single DNA sequence.
    More weight is given to similar fragments that are closer to the end of the sequence and match the exact motif pattern.

    :param sequence: A DNA sequence (string)
    :param motif_pattern: Sigma factor recognition motif pattern as a regular expression (string)
    :param max_errors: Maximum number of errors (insertions, deletions, or substitutions) allowed in approximate matches (default is 2)
    :param exact_match_weight: Weight factor for exact motif pattern matches (default is 2)
    :return: Tuple containing a boolean value for the presence of the exact pattern, a list of the most similar sequence fragments, and their weighted similarity scores
    """
    # Check for the presence of the exact pattern
    exact_match = re_lib.search(motif_pattern, sequence)
    has_exact_pattern = exact_match is not None

    # Find similar sequence fragments
    fuzzy_pattern = f"({motif_pattern}){{e<={max_errors}}}"
    similar_matches = re_lib.finditer(fuzzy_pattern, sequence)

    # Calculate weighted similarity scores and find the most similar fragments
    highest_score = -1
    #most_similar_fragments = []
    most_similar_fragment = ''

    for match in similar_matches:
        errors = sum(1 for _ in match.fuzzy_changes)
        pattern_length = len(match.group())
        similarity_score = (pattern_length - errors) / pattern_length

        # Apply the exact match weight
        if errors == 0:
            similarity_score *= exact_match_weight

        # Calculate the position weight
        '''        if weighted_similarity_score > highest_score:
                    highest_score = weighted_similarity_score
                    most_similar_fragments = [match.group()]
                elif weighted_similarity_score == highest_score:
                    most_similar_fragments.append(match.group())'''
        if position == 'on':
            position_weight = (1 + (len(sequence) - match.end()) / len(sequence))**2
            weighted_similarity_score = similarity_score * position_weight
        else:
            weighted_similarity_score = similarity_score
        if weighted_similarity_score > highest_score:
            highest_score = weighted_similarity_score
            most_similar_fragment = match.group()
    return has_exact_pattern, most_similar_fragment, highest_score

# 读取基因组
def genome_reader(path=f'E:/sigma因子的研究/组学数据/衍生处理数据/MG1655.fasta'):
    with open(path, 'r') as genome_fasta:
        # Skip the first line by calling readline() method
        genome_fasta.readline()
        genome_sequence =''
        for line in genome_fasta.readlines():
            # Process each line as needed
            genome_sequence += line
    return genome_sequence.replace('\n','')

# Gene对象导出excel
def genes_to_excel(genes, output_file):
    # 将Gene对象转换为字典列表
    genes_dicts = []
    for gene in genes:
        # 基本信息
        gene_dict = {
            'gene_id': gene.gene_id,
            'name': gene.name,
            'start': gene.start,
            'end': gene.end,
            'sequence': gene.sequence,
            'promoter': gene.promoter,
            'Number of promoters':len(gene.promoter),
            'sigma_factor': gene.sigma_factor,
            'operon_genes': gene.operon_genes,
            'direction': gene.direction,
            'upstream seq': gene.upstream,
            'max_sigma': gene.max_sigma_factor,
            'Predicted_promoter': gene.Predicted_promoter,
            'Predicted_right': gene.predicted
        }
        # 展开sigma_TPM字段
        for sigma, TPM in gene.sigma_TPM.items():
            gene_dict[f'{sigma}_TPM'] = TPM
        genes_dicts.append(gene_dict)
    # 创建DataFrame
    genes_df = pd.DataFrame(genes_dicts)
    # 写入Excel文件
    genes_df.to_excel(output_file, index=False)

# 序列的碱基互补配对反序序列
def complement_base_pair(sequence):
    complement = ""
    for base in sequence:
        if base == "A":
            complement += "T"
        elif base == "T":
            complement += "A"
        elif base == "C":
            complement += "G"
        elif base == "G":
            complement += "C"
    return complement[::-1]

# 从基因组中提取指定位置序列
def seq_prosition(genome_seq, start, end):
    return (genome_seq[start-1:end])

#读取特定文件中的特定列数据，columns可为数组
def read_excel_columns(file_path, columns):
    try:
        df = pd.read_excel(file_path, usecols=columns, engine='openpyxl')
        return df
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

#读取特定dataform中的指定格子
def extract_cell_data(df, row_index, column_name):
    try:
        cell_data = df.at[row_index, column_name]
        return cell_data
    except KeyError:
        print(f"Column {column_name} does not exist.")
    except IndexError:
        print(f"Row index {row_index} is out of bounds.")
    return None

# 对operon列的每一项进行基因分割到单独的set
def split_transcription_units_to_set(df, gene_name_column, transcription_unit_column):
    #如果是空的，那就是当这个operon就一个基因
    def split_or_default(row):
        if pd.isnull(row[transcription_unit_column]):
            return {row[gene_name_column]}
        else:
            return set(row[transcription_unit_column].split(' // ')) # operon中各个基因由“//”间隔
    df[transcription_unit_column] = df.apply(split_or_default, axis=1)
    return df

# 整合两个dataform，使用的是key1和key2作为合并标识符
def merge_dataframes(df1, df2, key1, key2):
    merged_df = pd.merge(df1, df2, left_on=key1, right_on=key2, how='left')
    return merged_df

#简化gene编号
def remove_prefix_from_gene_id(df, gene_id_column, prefix="gene-"):
    if gene_id_column in df.columns:
        # Remove the prefix from the gene ID column entries
        df[gene_id_column] = df[gene_id_column].str.replace(prefix, "", regex=False)
    else:
        print(f"The column {gene_id_column} does not exist in the DataFrame.")
    return df

# 遍历df中每一行的set_column中有没有element
def find_row_by_element(dataframe, set_column, element):
    for _, row in dataframe.iterrows():
        if element in row[set_column]:
            # 返回第一次搜索到的整行信息
            return row
    # 如果没有找到匹配的集合
    return None


# 文件导入
df_operon = read_excel_columns("E:/sigma因子的研究/组学数据/衍生处理数据/operon_promoter_seq.xlsx", ['Binds-Sigma-Factor', 'Sequence - DNA sequence','Genes of transcription unit'])
df_all_TPM = read_excel_columns("E:/sigma因子的研究/组学数据/衍生处理数据/all.tpm.xlsx", ['locus_tag', 'GeneName', 'RPOD-1_TPM', 'RPOE-1_TPM', 'RPOF-1_TPM', 'RPOH-1_TPM', "RPON-1_TPM", "RPOS-1_TPM", "FECI-1_TPM", "temp-1_TPM", "Start", "End"])
df_direction = read_excel_columns("E:/sigma因子的研究/组学数据/衍生处理数据/All-genes-MG1655-direction.xlsx", ['Gene_name','Direction'])

# 处理operon项
if df_operon is not None:
    df_operon = split_transcription_units_to_set(df_operon, 'Sequence - DNA sequence', 'Genes of transcription unit')
df_operon.fillna('missing', inplace=True)
genome_seq = genome_reader()

# 录入gene实体
genes = []
for _,row in df_all_TPM.iterrows():
    gene = Gene(
        gene_id=row['locus_tag'],  # 假设使用Accession-1作为基因的唯一ID
        name=row['GeneName'],
        start=row['Start'],  # 如果数据中没有起始和终止位置信息，可以暂时设置为None
        end=row['End'],
        sequence=None,  # 如果数据中没有序列信息，也设置为None
        promoter=[],  # 如果数据中没有启动子序列信息，设置为None
        sigma_factor=None,  # 如果数据中没有σ因子信息，设置为None
        #operon_genes=row['Genes in same transcription unit'],
        operon_genes=None,
        direction='',
        upstream_seq='',
        Predicted_promoter='',
        max_sigma_factor='',
        predictied_right=False,
        sigma_TPM=None  # 如果数据中没有TPM值，设置为None
    )
    gene.add_TPM("RPOD", float(row["RPOD-1_TPM"]))
    gene.add_TPM("RPOE", float(row["RPOE-1_TPM"]))
    gene.add_TPM("RPOF", float(row["RPOF-1_TPM"]))
    gene.add_TPM("RPOH", float(row["RPOH-1_TPM"]))
    gene.add_TPM("RPOS", float(row["RPOS-1_TPM"]))
    gene.add_TPM("temp", float(row["temp-1_TPM"]))
    gene.add_TPM("RPON", float(row["RPON-1_TPM"]))
    gene.add_TPM("FECI", float(row["FECI-1_TPM"]))
    gene.max_sigma_factor = max(gene.sigma_TPM, key=gene.sigma_TPM.get, default=None)
    # 限制TPM > 1
    if float(gene.sigma_TPM[gene.max_sigma_factor]) < 10:
        print(gene.gene_id)
        print(gene.max_sigma_factor)
        print(gene.sigma_TPM[gene.max_sigma_factor])
        gene.max_sigma_factor = 'N.S.'
    else:
        if gene.max_sigma_factor == 'temp':
            gene.max_sigma_factor = 'RPOD'
            #gene.max_sigma_factor = 'N.S.'
        elif float(gene.sigma_TPM[gene.max_sigma_factor]) <= float(gene.sigma_TPM['temp'])*2:
            gene.max_sigma_factor = 'N.S.'
    genes.append(gene)

for i in genes:
    row_data = find_row_by_element(df_operon, 'Genes of transcription unit', i.name)
    direction_Data = find_row_by_element(df_direction, 'Gene_name', i.name)
    try:
        i.direction = direction_Data['Direction']
    except TypeError:
        i.direction = 'Not Found'
    if seq_prosition(genome_seq, int(i.end) - 2, int(i.end)) in ('TAA', 'TAG', 'TGA') or i.direction=='+':
        i.upstream = seq_prosition(genome_seq, int(i.start) - 300, int(i.start) - 1)
        i.sequence = seq_prosition(genome_seq, int(i.start), int(i.end))
        if i.max_sigma_factor != 'N.S.':
            has_exact_pattern, most_similar_fragments, highest_score = find_exact_and_most_similar_patterns_weighted(
                i.upstream,
                    iupac_to_regex(patterns_sigma[i.max_sigma_factor]), max_errors=2)
            i.Predicted_promoter = most_similar_fragments
    elif complement_base_pair(seq_prosition(genome_seq, int(i.start), int(i.start) + 2)) in ('TAA', 'TAG', 'TGA') or i.direction== '-':
        i.upstream = complement_base_pair(seq_prosition(genome_seq, int(i.end) + 1, int(i.end) + 300))
        i.sequence = complement_base_pair(seq_prosition(genome_seq, int(i.start), int(i.end)))
        if i.max_sigma_factor != 'N.S.':
            has_exact_pattern, most_similar_fragments, highest_score = find_exact_and_most_similar_patterns_weighted(
                i.upstream,
                iupac_to_regex(patterns_sigma[i.max_sigma_factor]), max_errors=2)
            i.Predicted_promoter = most_similar_fragments
    if row_data is not None:
        i.sigma_factor = row_data['Binds-Sigma-Factor']
        if row_data['Sequence - DNA sequence'] != 'missing':
            i.promoter.append(row_data['Sequence - DNA sequence'])
        i.operon_genes = row_data['Genes of transcription unit']
    #if (i.promoter == []) or (row_data['Sequence - DNA sequence'] == 'missing'):
    if i.operon_genes == []:
        i.operon_genes = {i.name}
    try:
        sigma_factor_processed = i.sigma_factor.replace('RNA polymerase sigma factor ', '').upper()
        i.sigma_factor = sigma_factor_processed
        if i.sigma_factor == 'FLIA':
            i.sigma_factor = 'RPOF'
    except AttributeError:
        print("Not annotated.")
    try:
        if float(i.sigma_TPM[i.sigma_factor]) >= float(i.sigma_TPM['temp']):
            i.predicted = True
    except KeyError:
        i.predicted = None
    print(i.__repr__())

genes_to_excel(genes,"E:/sigma因子的研究/组学数据/衍生处理数据/merged_data_TPM_over_10_temp_D_fold2.xlsx")
