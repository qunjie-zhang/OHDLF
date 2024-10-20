import shutil
import time, os
import subprocess
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastpCommandline
import random
from argparse import ArgumentParser

def fasta_load(fasta_path):
    f = open(fasta_path, 'r', ).readlines()
    seq_name = None
    seq = None
    for line in f:
        if not line:
            break
        if line.startswith('>'):
            # if seq 会被解释为 if seq is not None
            if seq:
                yield seq_name, seq
            seq = ''
            seq_name = line[1:].strip('\n')
        else:
            seq += line.strip('\n')
    yield seq_name, seq

def filter_orthologous(loss,duplication):
    '''
    Find Orthogroups that contain designate low copy gene values
    :param loss: max_loss_num
    :param duplication: max_duplication_num
    :return: all_GDL_Orthologue_Sequences
    '''
    loss_num = float(loss)
    # 最大重复数
    max_duplication_num = int(duplication)


    if not os.path.exists('Orthogroups') or not os.path.exists('Orthogroup_Sequences'):
        print("Please place this script in the Orthofinder output directory")
        print("Programs rely primarily on directories：Orthogroups ,Orthogroup_Sequences")
        exit()

    # 创建结果目录
    if not os.path.exists('all_GDL_Orthologue_Sequences'):
        os.mkdir('all_GDL_Orthologue_Sequences')

    # 读取相对路径下的文件
    tsv_file_path = 'Orthogroups/Orthogroups.GeneCount.tsv'
    with open(tsv_file_path, 'r', encoding='utf8') as tsv_file:
        tsv_file = tsv_file.read().splitlines()

    # 获取物种名称及数量信息
    species_name = [i for i in tsv_file[0].split()[1:-1]]
    species_num = len(species_name)
    max_loss_num = round(loss_num * species_num)

    # 删除第一行表头
    del tsv_file[0]
    # 新建字典包含目标Orthogroup与缺失物种名称
    target_dict = {}

    for content in tsv_file:
        content = content.split()
        # 先删除索引为0的项，然后将被删除的项赋值给Orthogroup_name，之后content中就只剩数值了
        Orthogroup_name = content.pop(0)
        # pop（）是默认删除最后一个索引,这行代码表示取出数据表自带的总和，但这个值没用的不用管
        content_total = content.pop()
        # 判断拷贝数是否在范围内,若不是则无需继续，continue是跳出本次循环进入下一次循环
        if any(int(x) > max_duplication_num for x in content):
            continue
        if content.count('0') > max_loss_num:
            continue

        # 如果目标字典中没有该序列信息则进行初始化，创建一个空的列表
        if not Orthogroup_name in target_dict:
            target_dict[Orthogroup_name] = []
        n = 0
        # 将匹配到的物种名称添加进入一个列表
        for i in content:
            if int(i) == 0:
                target_dict[Orthogroup_name].append(species_name[n])
            n += 1
    for fasta_name, species_name in target_dict.items():
        # Orthogroup_Sequences/ 是文件夹路径的前缀，{fasta_name} 是要插入的变量，表示文件名，.fa 是文件的扩展名。
        raw_data = fasta_load(f'Orthogroup_Sequences/{fasta_name}.fa')
        with open(f'all_GDL_Orthologue_Sequences/{fasta_name}_repeat.fa', 'w', encoding='utf8') as file:
            for i in raw_data:
                file.write('>' + i[0] + '\n')
                file.write(i[1] + '\n')

            for i in species_name:
                file.write('>' + i + '|EMPTY_DATA_' + str(random.random()) + '\n')
                file.write(100 * '-' + '\n')

    print('[filter_low_copy_orthologous is executed]')

def blast_identity_filter(sim):
    '''
    filter blast_identity(>=similarity_threshold) orthologous
    :return: GDL_Orthologue_Sequences
    '''
    similarity_threshold = float(sim)
    input_dir = "all_GDL_Orthologue_Sequences"
    output_dir = "GDL_Orthologue_Sequences"
    temp1 = "temp1.fa"
    temp2 = "temp2.fa"
    temp_db = "temp"
    blast_output = "temp_all.blastout"
    if not os.path.exists(input_dir):
        print('Please place this script in the Orthofinder output directory！！(Results_*)\n')
        print('Programs rely primarily on directories：all_GDL_Orthologue_Sequences')
        exit()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 处理目录中的每个FASTA文件
    for fasta_file in os.listdir(input_dir):
        if fasta_file.endswith(".fasta") or fasta_file.endswith(".fa"):
            input_fasta = os.path.join(input_dir, fasta_file)
            print(f"Processing file: {input_fasta}")

            # 读取fasta文件
            species_sequences = defaultdict(list)
            for record in SeqIO.parse(input_fasta, "fasta"):
                species_name = record.id.split('|')[0]
                species_sequences[species_name].append(record)

            all_species_similar = True

            # 遍历每个物种及其序列
            for species, sequences in species_sequences.items():
                if len(sequences) > 1:
                    # 将该物种的序列ID及其对应的第一条序列输出到temp1.fa
                    with open(temp1, "w") as f1:
                        SeqIO.write(sequences[0], f1, "fasta")

                    # 将该物种的ID及其对应的所有序列输出到temp2.fa
                    with open(temp2, "w") as f2:
                        SeqIO.write(sequences, f2, "fasta")

                    # 调用BLAST进行建库
                    makeblastdb_cline = NcbimakeblastdbCommandline(dbtype="prot", input_file=temp1, out=temp_db)
                    makeblastdb_cline()

                    # 进行比对
                    blastp_cline = NcbiblastpCommandline(query=temp2, db=temp_db, out=blast_output, outfmt=6)
                    blastp_cline()

                    # 检查相似度
                    all_similar = True
                    similarity_dict = defaultdict(list)

                    # 读取BLAST结果并处理相同比对保留相似度最大的行
                    with open(blast_output, "r") as blast_result:
                        for line in blast_result:
                            fields = line.strip().split("\t")
                            query_id = fields[0]
                            subject_id = fields[1]
                            similarity = float(fields[2])
                            similarity_dict[(query_id, subject_id)].append(similarity)

                    # 保留相似度最高的记录
                    max_similarity_records = {k: max(v) for k, v in similarity_dict.items()}

                    # 检查所有记录的相似度
                    for sim in max_similarity_records.values():
                        if sim < similarity_threshold:
                            all_similar = False
                            break

                    # 如果相似度不满足要求，跳过该文件
                    if not all_similar:
                        all_species_similar = False
                        print(f"File {fasta_file} does not meet the similarity requirement. Skipping.")
                        break

                    # 删除临时文件和中间文件
                    os.remove(temp1)
                    os.remove(temp2)
                    os.remove(blast_output)
                    for ext in ["phr", "pin", "psq"]:
                        if os.path.exists(f"{temp_db}.{ext}"):
                            os.remove(f"{temp_db}.{ext}")

            if all_species_similar:
                # 如果所有物种的相似度均大于阈值，保存原始fasta到输出目录
                shutil.copy(input_fasta, os.path.join(output_dir, fasta_file))
                print(f"File {fasta_file} meets the similarity requirement. Copied to output directory.")

    print("All files processed. filter completed.")


def align_by_mafft():
    '''
    invoke mafft
    :return:GDL_Orthologue_Sequences_mafft
    '''
    # 指定输入和输出目录
    input_dir = "GDL_Orthologue_Sequences"
    output_dir = "GDL_Orthologue_Sequences_mafft"

    if not os.path.exists(input_dir):
        print('Please place this script in the Orthofinder output directory！！(Results_*)\n')
        print('Programs rely primarily on directories：GDL_Orthologue_Sequences')
        exit()

    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 获取输入目录下的所有.fa文件
    fa_files = [file for file in os.listdir(input_dir) if file.endswith(".fa")]

    # 使用MAFFT对每个.fa文件进行多序列比对
    for fa_file in fa_files:
        # 构建输入和输出文件的路径,input_dir是一个目录路径的字符串变量，file_name是一个文件名的字符串变量。通过使用os.path.join函数，将这两个字符串变量连接起来，得到一个完整的文件路径
        input_file = os.path.join(input_dir, fa_file)
        output_file = os.path.join(output_dir, f"{fa_file}_align.fa")

        # 执行MAFFT命令进行多序列比对
        mafft_cmd = f"mafft --auto {input_file} > {output_file}"
        subprocess.run(mafft_cmd, shell=True)

    print('[mafft is finished]')

def align_and_fill():
    '''
    Compare the duplicate sequence IDs in the file
    fill in X in the difference position
    fill in the difference - while leaving the gap, and merge the more into one
    :return:GDL_Orthologue_Sequences_mafft_fill
    '''
    # 指定输入和输出目录
    input_dir = "GDL_Orthologue_Sequences_mafft"
    output_dir = "GDL_Orthologue_Sequences_mafft_fill"

    if not os.path.exists(input_dir):
        print('Please place this script in the Orthofinder output directory！！(Results_*)\n')
        print('Programs rely primarily on directories：GDL_Orthologue_Sequences_mafft')
        exit()

    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 获取输入目录下的所有.fa文件
    fa_files = [file for file in os.listdir(input_dir) if file.endswith(".fa")]

    # 对每个.fa文件进行比对
    for fa_file in fa_files:
        input_file = os.path.join(input_dir, fa_file)
        output_file = os.path.join(output_dir, f"{fa_file}_fill.fa")
        records = list(SeqIO.parse(input_file, "fasta"))
        id_seq_dict = {}
        for record in records:
            # 将字符串以"|"为分隔符，然后取|前面的作为序列ID，例如LGC_19ZG|TRINITY_DN60348_c1_g2_i7.p1，分隔后变为LGC_19ZG
            seq_id = record.id.split("|")[0]
            seq = str(record.seq).replace("-", "-")
            if seq_id in id_seq_dict:
                # 对比序列
                for i in range(len(seq)):
                    if seq[i] != id_seq_dict[seq_id][i]:
                        if seq[i] == "-" or id_seq_dict[seq_id][i] == "-":
                            id_seq_dict[seq_id] = id_seq_dict[seq_id][:i] + "-" + id_seq_dict[seq_id][i + 1:]
                        else:
                            # 不同的位置填充为X
                            id_seq_dict[seq_id] = id_seq_dict[seq_id][:i] + "X" + id_seq_dict[seq_id][i + 1:]
            else:
                id_seq_dict[seq_id] = seq
        with open(output_file, 'w', encoding='utf8') as f:
            for seq_id, seq in id_seq_dict.items():
                record = SeqIO.SeqRecord(Seq(seq), id=seq_id, description="")
                SeqIO.write(record, f, "fasta")

    print('[fill is finished]')

def seq_merge():
    '''
    After reading the fasta file, the sequences of the same species name are merged
    and if there are N species, the result will have N sequences
    :return: final_OrthologsAlign_GDL.fasta
    '''
    input_dir = "GDL_Orthologue_Sequences_mafft_fill"
    output_file = "final_OrthologsAlign_GDL.fasta"

    if not os.path.exists(input_dir):
        print('Please place this script in the Orthofinder output directory！！(Results_*)\n')
        print('Programs rely primarily on directories:GDL_Orthologue_Sequences_mafft_fill')
        exit()

    # 新建字典
    seq = {}

    # 遍历指定目录下的所有.fa文件
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".fa"):
            file_path = os.path.join(input_dir, file_name)

            # 读取.fa文件
            with open(file_path, "r") as file:
                for line in file:
                    line = line.strip()
                    if line.startswith(">"):
                        # 提取序列ID
                        seq_id = line[1:].strip('\n')
                        if seq_id not in seq:
                            # 初始化序列内容
                            seq[seq_id] = ""
                    else:
                        # 拼接序列内容
                        seq[seq_id] += line.strip('\n')

    # 将合并后的序列写入输出文件
    with open(output_file, "w") as output:
        for sequence_id, sequence in seq.items():
            output.write(f">{sequence_id}\n")
            output.write(f"{sequence}\n")

    print('[seq merge is finished]')

def FasToPhy():
    '''
    fasta turn to phylip
    :return: final_OrthologsAlign_GDL.phy
    '''
    in_file = 'final_OrthologsAlign_GDL.fasta'
    # 从输入文件路径中更改文件后缀为.phy
    base_name = os.path.splitext(in_file)[0]
    out_file = f"{base_name}.phy"

    names = []
    seqs = []
    max_name_len = 0
    with open(in_file, 'r') as IN:
        for line in IN:
            line = line.strip()
            if line.startswith('>'):
                tp = line.split("|")
                name = tp[0][1:]
                names.append(name)
                seqs.append('')
                if len(name) > max_name_len:
                    max_name_len = len(name)
            else:
                seqs[-1] += line

    count = len(seqs)
    len_seq = len(seqs[0]) if seqs else 0
    with open(out_file, 'w') as OUT:
        OUT.write(f"{count} {len_seq}\n")
        for name, seq in zip(names, seqs):
            add_space = ' ' * (max_name_len + 10 - len(name))
            if max_name_len >= 10:
                print("Warning: sequence name too long!!!")
            OUT.write(f"{name}{add_space}{seq}\n")

    print('[fasta > phy finish]')

def filter_orthologous_parallel(loss,duplication):
    '''
    Find Orthogroups that contain designate low copy gene values
    :param loss: max_loss_num
    :param duplication: max_duplication_num
    :return: all_GDL_Orthologue_Sequences_parallel
    '''
    loss_num = float(loss)
    # 最大重复数
    max_duplication_num = int(duplication)


    if not os.path.exists('Orthogroups') or not os.path.exists('Orthogroup_Sequences'):
        print("Please place this script in the Orthofinder output directory")
        print("Programs rely primarily on directories：Orthogroups ,Orthogroup_Sequences")
        exit()

    # 创建结果目录
    if not os.path.exists('all_GDL_Orthologue_Sequences_parallel'):
        os.mkdir('all_GDL_Orthologue_Sequences_parallel')

    # 读取相对路径下的文件
    tsv_file_path = 'Orthogroups/Orthogroups.GeneCount.tsv'
    with open(tsv_file_path, 'r', encoding='utf8') as tsv_file:
        tsv_file = tsv_file.read().splitlines()

    # 获取物种名称及数量信息
    species_name = [i for i in tsv_file[0].split()[1:-1]]
    species_num = len(species_name)
    max_loss_num = round(species_num * loss_num)

    # 删除第一行表头
    del tsv_file[0]
    # 新建字典包含目标Orthogroup与缺失物种名称
    target_dict = {}

    for content in tsv_file:
        content = content.split()
        # 先删除索引为0的项，然后将被删除的项赋值给Orthogroup_name，之后content中就只剩数值了
        Orthogroup_name = content.pop(0)
        # pop（）是默认删除最后一个索引,这行代码表示取出数据表自带的总和，但这个值没用的不用管
        content_total = content.pop()
        # 判断拷贝数是否在范围内，若不是则无需继续，continue是跳出本次循环进入下一次循环
        if any(int(x) > max_duplication_num for x in content):
            continue
        # 筛选缺失值
        if content.count('0') > max_loss_num:
            continue

        # 如果目标字典中没有该序列信息则进行初始化，创建一个空的列表
        if not Orthogroup_name in target_dict:
            target_dict[Orthogroup_name] = []
        n = 0
        # 将匹配到的物种名称添加进入一个列表
        for i in content:
            if int(i) == 0:
                target_dict[Orthogroup_name].append(species_name[n])
            n += 1
    for fasta_name, species_name in target_dict.items():
        # Orthogroup_Sequences/ 是文件夹路径的前缀，{fasta_name} 是要插入的变量，表示文件名，.fa 是文件的扩展名。
        raw_data = fasta_load(f'Orthogroup_Sequences/{fasta_name}.fa')
        with open(f'all_GDL_Orthologue_Sequences_parallel/{fasta_name}_repeat.fa', 'w', encoding='utf8') as file:
            for i in raw_data:
                file.write('>' + i[0] + '\n')
                file.write(i[1] + '\n')



    print('[filter_low_copy_orthologous is executed]')

def blast_identity_filter_parallel(sim):
    '''
    filter blast_identity(>=similarity_threshold) orthologous
    :return: GDL_Orthologue_Sequences_parallel
    '''
    similarity_threshold = float(sim)
    input_dir = "all_GDL_Orthologue_Sequences_parallel"
    output_dir = "GDL_Orthologue_Sequences_parallel"
    temp1 = "temp1.fa"
    temp2 = "temp2.fa"
    temp_db = "temp"
    blast_output = "temp_all.blastout"
    if not os.path.exists(input_dir):
        print('Please place this script in the Orthofinder output directory！！(Results_*)\n')
        print('Programs rely primarily on directories：all_GDL_Orthologue_Sequences_parallel')
        exit()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 处理目录中的每个FASTA文件
    for fasta_file in os.listdir(input_dir):
        if fasta_file.endswith(".fasta") or fasta_file.endswith(".fa"):
            input_fasta = os.path.join(input_dir, fasta_file)
            print(f"Processing file: {input_fasta}")

            # 读取fasta文件
            species_sequences = defaultdict(list)
            for record in SeqIO.parse(input_fasta, "fasta"):
                species_name = record.id.split('|')[0]
                species_sequences[species_name].append(record)

            all_species_similar = True

            # 遍历每个物种及其序列
            for species, sequences in species_sequences.items():
                if len(sequences) > 1:
                    # 将该物种的序列ID及其对应的第一条序列输出到temp1.fa
                    with open(temp1, "w") as f1:
                        SeqIO.write(sequences[0], f1, "fasta")

                    # 将该物种的ID及其对应的所有序列输出到temp2.fa
                    with open(temp2, "w") as f2:
                        SeqIO.write(sequences, f2, "fasta")

                    # 调用BLAST进行建库
                    makeblastdb_cline = NcbimakeblastdbCommandline(dbtype="prot", input_file=temp1, out=temp_db)
                    makeblastdb_cline()

                    # 进行比对
                    blastp_cline = NcbiblastpCommandline(query=temp2, db=temp_db, out=blast_output, outfmt=6)
                    blastp_cline()

                    # 检查相似度
                    all_similar = True
                    similarity_dict = defaultdict(list)

                    # 读取BLAST结果并处理相同比对保留相似度最大的行
                    with open(blast_output, "r") as blast_result:
                        for line in blast_result:
                            fields = line.strip().split("\t")
                            query_id = fields[0]
                            subject_id = fields[1]
                            similarity = float(fields[2])
                            similarity_dict[(query_id, subject_id)].append(similarity)

                    # 保留相似度最高的记录
                    max_similarity_records = {k: max(v) for k, v in similarity_dict.items()}

                    # 检查所有记录的相似度
                    for sim in max_similarity_records.values():
                        if sim < similarity_threshold:
                            all_similar = False
                            break

                    # 如果相似度不满足要求，跳过该文件
                    if not all_similar:
                        all_species_similar = False
                        print(f"File {fasta_file} does not meet the similarity requirement. Skipping.")
                        break

                    # 删除临时文件和中间文件
                    os.remove(temp1)
                    os.remove(temp2)
                    os.remove(blast_output)
                    for ext in ["phr", "pin", "psq"]:
                        if os.path.exists(f"{temp_db}.{ext}"):
                            os.remove(f"{temp_db}.{ext}")

            if all_species_similar:
                # 如果所有物种的相似度均大于阈值，保存原始fasta到输出目录
                shutil.copy(input_fasta, os.path.join(output_dir, fasta_file))
                print(f"File {fasta_file} meets the similarity requirement. Copied to output directory.")

    print("All files processed. filter completed.")

def align_by_mafft_parallel():
    '''
    invoke mafft
    :return:GDL_Orthologue_Sequences_mafft_parallel
    '''
    # 指定输入和输出目录
    input_dir = "GDL_Orthologue_Sequences_parallel"
    output_dir = "GDL_Orthologue_Sequences_mafft_parallel"

    if not os.path.exists(input_dir):
        print('Please place this script in the Orthofinder output directory！！(Results_*)\n')
        print('Programs rely primarily on directories：GDL_Orthologue_Sequences_parallel')
        exit()

    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 获取输入目录下的所有.fa文件
    fa_files = [file for file in os.listdir(input_dir) if file.endswith(".fa")]

    # 使用MAFFT对每个.fa文件进行多序列比对
    for fa_file in fa_files:
        # 构建输入和输出文件的路径,input_dir是一个目录路径的字符串变量，file_name是一个文件名的字符串变量。通过使用os.path.join函数，将这两个字符串变量连接起来，得到一个完整的文件路径
        input_file = os.path.join(input_dir, fa_file)
        output_file = os.path.join(output_dir, f"{fa_file}_align.fa")

        # 执行MAFFT命令进行多序列比对
        mafft_cmd = f"mafft --auto {input_file} > {output_file}"
        subprocess.run(mafft_cmd, shell=True)

    print('[mafft is finished]')

def align_and_fill_parallel():
    '''
    Compare the duplicate sequence IDs in the file
    fill in X in the difference position
    fill in the difference - while leaving the gap, and merge the more into one
    :return:GDL_Orthologue_Sequences_mafft_fill_parallel
    '''
    # 指定输入和输出目录
    input_dir = "GDL_Orthologue_Sequences_mafft_parallel"
    output_dir = "GDL_Orthologue_Sequences_mafft_fill_parallel"

    if not os.path.exists(input_dir):
        print('Please place this script in the Orthofinder output directory！！(Results_*)\n')
        print('Programs rely primarily on directories：GDL_Orthologue_Sequences_mafft')
        exit()

    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 获取输入目录下的所有.fa文件
    fa_files = [file for file in os.listdir(input_dir) if file.endswith(".fa")]

    # 对每个.fa文件进行比对
    for fa_file in fa_files:
        input_file = os.path.join(input_dir, fa_file)
        output_file = os.path.join(output_dir, f"{fa_file}_fill.fa")
        records = list(SeqIO.parse(input_file, "fasta"))
        id_seq_dict = {}
        for record in records:
            # 将字符串以"|"为分隔符，然后取|前面的作为序列ID，例如LGC_19ZG|TRINITY_DN60348_c1_g2_i7.p1，分隔后变为LGC_19ZG
            seq_id = record.id.split("|")[0]
            seq = str(record.seq).replace("-", "-")
            if seq_id in id_seq_dict:
                # 对比序列
                for i in range(len(seq)):
                    if seq[i] != id_seq_dict[seq_id][i]:
                        if seq[i] == "-" or id_seq_dict[seq_id][i] == "-":
                            id_seq_dict[seq_id] = id_seq_dict[seq_id][:i] + "-" + id_seq_dict[seq_id][i + 1:]
                        else:
                            # 不同的位置填充为X
                            id_seq_dict[seq_id] = id_seq_dict[seq_id][:i] + "X" + id_seq_dict[seq_id][i + 1:]
            else:
                id_seq_dict[seq_id] = seq
        with open(output_file, 'w', encoding='utf8') as f:
            for seq_id, seq in id_seq_dict.items():
                record = SeqIO.SeqRecord(Seq(seq), id=seq_id, description="")
                SeqIO.write(record, f, "fasta")

    print('[fill is finished]')

def run_command(command):
    """运行命令并检查是否成功，返回命令是否成功执行"""
    try:
        result = subprocess.run(command, shell=True, check=True)
        return result.returncode == 0
    except subprocess.CalledProcessError:
        return False

def FasToPhy_parallel(in_file):

    # 从输入文件路径中更改文件后缀为.phy
    base_name = os.path.splitext(in_file)[0]
    out_file = f"{base_name}.phy"

    names = []
    seqs = []
    max_name_len = 0
    with open(in_file, 'r') as IN:
        for line in IN:
            line = line.strip()
            if line.startswith('>'):
                tp = line.split("|")
                name = tp[0][1:]
                names.append(name)
                seqs.append('')
                if len(name) > max_name_len:
                    max_name_len = len(name)
            else:
                seqs[-1] += line

    count = len(seqs)
    len_seq = len(seqs[0]) if seqs else 0
    with open(out_file, 'w') as OUT:
        OUT.write(f"{count} {len_seq}\n")
        for name, seq in zip(names, seqs):
            add_space = ' ' * (max_name_len + 10 - len(name))
            if max_name_len >= 10:
                print("Warning: sequence name too long!!!")
            OUT.write(f"{name}{add_space}{seq}\n")

    print('[fasta > phy finish]')


def process_directory():
    # 处理所有.fa文件
    directory = 'GDL_Orthologue_Sequences_mafft_fill_parallel'
    for file in os.listdir(directory):
        if file.endswith('.fa'):
            file_path = os.path.join(directory, file)

            # 使用 sed 删除空格，并将输出保存到新文件
            base_name = os.path.splitext(file_path)[0]
            fasta_file = f"{base_name}.fasta"
            if not run_command(f"sed 's/ //g' {file_path} > {fasta_file}"):
                print(f"sed command failed for {file_path}")
                continue

            # 调用 FasToPhy 函数，并将 .fasta 文件转换为 .phy 格式
            try:
                FasToPhy_parallel(fasta_file)
            except Exception as e:
                print(f"FasToPhy function failed for {fasta_file} with error {e}")
                continue

    # 处理所有.phy文件
    for file in os.listdir(directory):
        if file.endswith('.phy'):
            file_path = os.path.join(directory, file)
            # 执行 IQ-TREE 分析
            if not run_command(f"iqtree -s {file_path} -B 1000 --bnni -T AUTO"):
                print(f"IQ-TREE command failed for {file_path}")
                continue
    if not run_command(f"cat {directory}/*.treefile > {directory}/all.trees"):
        print("Failed to concatenate treefiles into all.trees")
        return

        # 将all.trees文件复制到上一级目录
    if not run_command(f"cp {directory}/all.trees ../"):
        print("Failed to copy all.trees to the parent directory")

def main(args):
    loss = float(args.loss)
    duplication = int(args.duplication)
    process_type = int(args.process_type)
    start_time = time.time()
    sim = float(args.similarity)
    if process_type == 1:
        filter_orthologous(loss,duplication)
        blast_identity_filter(sim)
        align_by_mafft()
        align_and_fill()
        seq_merge()
        FasToPhy()
    elif process_type == 2:
        filter_orthologous_parallel(loss,duplication)
        blast_identity_filter_parallel(sim)
        align_by_mafft_parallel()
        align_and_fill_parallel()
        process_directory()
    end_time = time.time()
    work_time = end_time - start_time
    print(f"All analyses have been completed")
    print(f"Program running time: {int(work_time / 60)} minutes")


if __name__ == "__main__":
    arg = ArgumentParser(description='OHDLF v1.0 - a pipeline designed to filter and address gene heterogeneity, duplication, and loss')
    arg.add_argument("-l",
                     "--loss",
                     default=0,
                     required=True,
                     help="Allowable loss value")
    arg.add_argument("-d",
                     "--duplication",
                     default=3,
                     required=True,
                     help="Allowable max copy value")
    arg.add_argument("-s",
                     "--similarity",
                     default=97,
                     required=False,
                     help="similarity value")
    arg.add_argument("-p", "--process_type",
                     type=int,
                     choices=[1, 2],
                     required=True,
                     help="process_type: 1 for Concatenation, 2 for Coalescence")
    args = arg.parse_args()
    main(args)
