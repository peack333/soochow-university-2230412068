from Bio import SeqIO
import pandas as pd
import os

# ===================== 配置 =====================
GENE_LIST = []
HOMOLOGY_FILE = "human_to_mouse.txt"
FASTA_HUMAN   = "human_mus_cds_longest.fasta"
FASTA_MOUSE   = "mouse_hu_cds_longest.fasta"
OUTPUT_DIR    = "homogenes/gene_orthologs_no2o"
PATH_FILE     = "seq_path_no2o.txt"
# =================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)
target_genes = [g.strip() for g in GENE_LIST]

# 1. 读取同源表，【彻底过滤空值】，避免float报错
df_homo = pd.read_csv(HOMOLOGY_FILE, sep="\t", dtype=str)
# 关键修复：同时过滤人基因和鼠基因的空值
df_homo = df_homo.dropna(subset=["Gene_name", "Mouse gene name"])

# 构建 人 → [鼠基因1, 鼠基因2...] 列表（一对多按顺序）
human2mouse_list = {}
for _, row in df_homo.iterrows():
    h_gene = str(row["Gene_name"]).strip()    # 强转str，彻底防报错
    m_gene = str(row["Mouse gene name"]).strip()
    
    if h_gene not in human2mouse_list:
        human2mouse_list[h_gene] = []
    human2mouse_list[h_gene].append(m_gene)

# 2. 读入FASTA
human_seqs = {rec.id.strip(): rec for rec in SeqIO.parse(FASTA_HUMAN, "fasta")}
mouse_seqs = {rec.id.strip(): rec for rec in SeqIO.parse(FASTA_MOUSE, "fasta")}

# 忽略大小写匹配字典
human_lower = {k.lower(): v for k, v in human_seqs.items()}
mouse_lower = {k.lower(): v for k, v in mouse_seqs.items()}

# 3. 批量提取序列
all_paths = []
missing_human = []
missing_mouse = []

for human_gene in target_genes:
    # 获取所有小鼠同源候选基因（按顺序）
    mouse_candidates = human2mouse_list.get(human_gene, [human_gene])
    m_rec = None

    # 按顺序尝试，找到第一个有序列的就停止
    for mg in mouse_candidates:
        m_rec = mouse_seqs.get(mg) or mouse_lower.get(mg.lower())
        if m_rec:
            break

    # 查找人序列（忽略大小写）
    h_rec = human_seqs.get(human_gene) or human_lower.get(human_gene.lower())

    if not h_rec:
        missing_human.append(human_gene)
        continue

    # 输出文件
    fa_path = os.path.join(OUTPUT_DIR, f"{human_gene}.fa")
    all_paths.append(fa_path)

    with open(fa_path, "w") as f:
        f.write(f">{h_rec.id}_human\n{h_rec.seq}\n")
        if m_rec:
            f.write(f">{m_rec.id}_mouse\n{m_rec.seq}\n")

    if not m_rec:
        missing_mouse.append(human_gene)

# 4. 保存路径文件
with open(PATH_FILE, "w") as f:
    f.write("\n".join(all_paths))

# 运行报告
print("\n✅ 任务完成！")
print(f"目标基因总数：{len(target_genes)}")
print(f"成功生成文件：{len(all_paths)}")

if missing_human:
    print(f"\n⚠️ 缺失人类序列：{len(missing_human)} 个")
    print(", ".join(missing_human))
if missing_mouse:
    print(f"\n⚠️ 缺失小鼠序列：{len(missing_mouse)} 个")
    print(", ".join(missing_mouse))