#!/usr/bin/env python3
"""
dN/dS分析脚本 - 提取距离值和标准误(SE) + Z-test + Fisher精确检验
"""

import os
import sys
import re
import subprocess
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Tuple, Optional

# 强制关闭输出缓冲，解决进度条延迟显示问题
sys.stdout.reconfigure(line_buffering=True)

# ==============================================
# 配置
# ==============================================
INPUT_FILE_DEFAULT = "seq_path_o2o.txt"
RESULTS_DIR        = "results/o2o"

ALIGN_MAO          = "muscle_align_coding.mao"
DS_MAO             = "distance_estimation_pairwise_syn-nonsynonymous_s.mao"
DN_MAO             = "distance_estimation_pairwise_syn-nonsynonymous_n.mao"
ZTEST_MAO          = "zTest_syn-nonsynonymous.mao"
FISHER_MAO         = "fisher_exact_test_syn-nonsynonymous.mao"

P_THRESHOLD        = 0.05

# ==============================================
# 你原来的解析函数 完全不动！
# 只复用它来读取 Z-test / Fisher
# ==============================================
def extract_distance_and_se(meg_file: str) -> Tuple[Optional[float], Optional[float]]:
    """
    从MEGA结果文件中提取结果
    """
    try:
        with open(meg_file, 'r') as f:
            content = f.read()
        
        distance = None
        se = None
        
        lines = content.strip().split('\n')
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            if '#' in line:
                continue

            if line.startswith('[1]'):
                value_part = line[3:].strip()
                matches = re.findall(r'\[?(-?[0-9.]+)\]?', value_part)
                if matches:
                    try:
                        se = float(matches[0])
                    except:
                        pass
            
            elif line.startswith('[2]'):
                value_part = line[3:].strip()
                matches = re.findall(r'\[?(-?[0-9.]+)\]?', value_part)
                if matches:
                    try:
                        distance = float(matches[0])
                    except:
                        pass
        
        return distance, se
        
    except Exception as e:
        print(f"解析错误 {meg_file}: {e}")
        return None, None

# ==============================================
# 带显著性的分类（自动优先Fisher P值）
# ==============================================
def classify_selection(omega: float, p_value: Optional[float]) -> Tuple[str, str]:
    if p_value is not None and p_value < P_THRESHOLD and omega > 1:
        return "significant_positive", f"P={p_value:.4f}, dN/dS={omega:.3f}>1"
    elif omega > 1:
        return "positive", f"dN/dS={omega:.3f}>1"
    elif omega < 1:
        return "purifying", f"dN/dS={omega:.3f}<1"
    else:
        return "neutral", f"dN/dS={omega:.3f}==1"

def process_gene(fasta_path: str, gene_id: str, temp_dir: Path, index: int) -> dict:
    """处理单个基因"""
    
    result = {
        'gene': gene_id,
        'dS': None, 'dS_se': None,
        'dN': None, 'dN_se': None,
        'dN/dS': None,
        'z_stat': None,
        'p_value': None,        # Z-test P
        'fisher_p': None,       # Fisher P
        'significant': None,
        'fisher_significant': None, # Fisher显著性
        'best_p': None,         # 最优P值(Fisher优先)
        'selection': 'unknown',
        'notes': '',
        'status': 'pending'
    }
    
    try:
        # ===================== 比对=====================
        align_prefix = temp_dir / f"gene_{index}_align"
        align_cmd = [
            "megacc", "-a", ALIGN_MAO,
            "-d", str(fasta_path),
            "-o", str(align_prefix),
            "-f", "fasta"
        ]
        subprocess.run(align_cmd, capture_output=True, text=True, check=True)
        aligned_file = Path(f"{align_prefix}.fasta") 
        
        if not aligned_file.exists():
            result['status'] = 'alignment_failed'
            result['notes'] = 'The comparison file has not been generated'
            return result
        
        # =====================  dS  =====================
        ds_prefix = temp_dir / f"gene_{index}_ds"
        ds_cmd = [
            "megacc", "-a", DS_MAO,
            "-d", str(aligned_file),
            "-o", str(ds_prefix)
        ]
        subprocess.run(ds_cmd, capture_output=True, text=True)
        ds_file = Path(f"{ds_prefix}.meg")
        if ds_file.exists():
            result['dS'], result['dS_se'] = extract_distance_and_se(str(ds_file))
        
        # =====================  dN  =====================
        dn_prefix = temp_dir / f"gene_{index}_dn"
        dn_cmd = [
            "megacc", "-a", DN_MAO,
            "-d", str(aligned_file),
            "-o", str(dn_prefix)
        ]
        subprocess.run(dn_cmd, capture_output=True, text=True)
        dn_file = Path(f"{dn_prefix}.meg")
        if dn_file.exists():
            result['dN'], result['dN_se'] = extract_distance_and_se(str(dn_file))
        
        # ===================== dN/dS =====================
        if (result['dS'] is not None and result['dN'] is not None and result['dS'] > 0):
            result['dN/dS'] = result['dN'] / result['dS']
        else:
            result['status'] = 'calculation_failed'
            result['notes'] = 'Value extraction failed'
            return result
        
        # ===================== Z-test =====================
        ztest_prefix = temp_dir / f"gene_{index}_ztest"
        ztest_cmd = [
            "megacc", "-a", ZTEST_MAO,
            "-d", str(aligned_file),
            "-o", str(ztest_prefix)
        ]
        subprocess.run(ztest_cmd, capture_output=True, text=True)
        ztest_file = Path(f"{ztest_prefix}.meg")

        if ztest_file.exists():
            result['p_value'], result['z_stat'] = extract_distance_and_se(str(ztest_file))
            result['significant'] = (result['p_value'] < P_THRESHOLD) if result['p_value'] else False

        # ===================== Fisher 精确检验 =====================
        fisher_prefix = temp_dir / f"gene_{index}_fisher"
        fisher_cmd = [
            "megacc", "-a", FISHER_MAO,
            "-d", str(aligned_file),
            "-o", str(fisher_prefix)
        ]
        subprocess.run(fisher_cmd, capture_output=True, text=True)
        fisher_file = Path(f"{fisher_prefix}.meg")

        if fisher_file.exists():
            # MEGA Fisher输出格式相同，直接复用解析函数
            fisher_p, _ = extract_distance_and_se(str(fisher_file))
            result['fisher_p'] = fisher_p
            result['fisher_significant'] = (fisher_p < P_THRESHOLD) if fisher_p else False

        # ===================== 自动选择最优P值：Fisher 优先 =====================
        if result['fisher_p'] is not None:
            result['best_p'] = result['fisher_p']
        else:
            result['best_p'] = result['p_value']

        # 分类
        result['selection'], result['notes'] = classify_selection(result['dN/dS'], result['best_p'])
        result['status'] = 'success'

    except subprocess.CalledProcessError as e:
        result['status'] = 'mega_error'
        result['notes'] = f"MEGA error: {e.stderr[:100]}"
    except Exception as e:
        result['status'] = 'error'
        result['notes'] = str(e)
    
    return result

def format_time(seconds: float) -> str:
    """格式化时间显示"""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}m"
    else:
        hours = seconds / 3600
        return f"{hours:.2f}h"

def generate_statistics_report(df: pd.DataFrame, output_csv: str, 
                               start_time: datetime, end_time: datetime, 
                               total_seconds: float, stats: dict) -> None:
    stats_file = output_csv.replace('.csv', '_statistics.txt')
    with open(stats_file, 'w') as f:
        f.write("dN/dS分析统计报告 (包含Z-test + Fisher精确检验)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"分析开始时间: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"分析结束时间: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"总用时: {format_time(total_seconds)}\n")
        f.write(f"总基因数: {len(df)}\n")
        f.write(f"成功分析: {(df['status'] == 'success').sum()}\n")
        f.write(f"分析失败: {(df['status'] != 'success').sum()}\n\n")

        # 新增Fisher统计
        success_df = df[df['status'] == 'success']
        if len(success_df) > 0:
            f.write(f"✅ 显著正选择(Fisher最优P<0.05): {((success_df['best_p'] < P_THRESHOLD) & (success_df['dN/dS'] > 1)).sum()}\n")
            f.write(f"🔍 Fisher检验可用数: {success_df['fisher_p'].notna().sum()}\n")
            f.write(f"📊 Z-test检验可用数: {success_df['p_value'].notna().sum()}\n\n")

def main():
    input_file = sys.argv[1] if len(sys.argv) > 1 else INPUT_FILE_DEFAULT
    if not os.path.exists(input_file):
        print(f"错误: 输入文件不存在: {input_file}")
        sys.exit(1)

    start_time = datetime.now()
    print(f"⏰ 分析开始时间: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    os.makedirs(RESULTS_DIR, exist_ok=True)
    timestamp = start_time.strftime("%Y%m%d_%H%M%S")
    output_csv = os.path.join(RESULTS_DIR, f"dnds_results_{timestamp}.csv")

    with open(input_file, 'r') as f:
        fasta_files = [line.strip() for line in f if line.strip()]

    print(f"🔬 dN/dS分析开始 (Z-test + Fisher精确检验)")
    print(f"输入文件: {input_file}")
    print(f"输出文件: {output_csv}")
    print(f"基因数量: {len(fasta_files)}")
    print("-" * 60)

    columns = [
        'gene', 'dS', 'dS_se', 'dN', 'dN_se', 'dN/dS',
        'z_stat', 'p_value', 'significant',
        'fisher_p', 'fisher_significant', 'best_p',  # 新增Fisher列
        'selection', 'notes', 'status'
    ]
    results_list = []
    stats = {'total': len(fasta_files), 'success': 0, 'failed': 0}

    with tempfile.TemporaryDirectory() as tmpdir:
        temp_dir = Path(tmpdir)
        for i, fasta_path in enumerate(fasta_files, 1):
            gene_name = Path(fasta_path).stem
            if gene_name.startswith("ortholog_"):
                parts = gene_name.split("_", 2)
                gene_name = parts[2] if len(parts) > 2 else gene_name.replace("ortholog_", "")
            
            # ===================== 进度条 =====================
            percent = i * 100 // len(fasta_files)
            bar_length = 40
            filled = int(bar_length * i / len(fasta_files))
            bar = "█" * filled + "░" * (bar_length - filled)
            
            elapsed_time = (datetime.now() - start_time).total_seconds()
            if i > 1:
                avg_time_per_gene = elapsed_time / i
                remaining_genes = len(fasta_files) - i
                estimated_remaining = avg_time_per_gene * remaining_genes
                estimated_total = elapsed_time + estimated_remaining
                
                elapsed_str = format_time(elapsed_time)
                remaining_str = format_time(estimated_remaining)
                total_str = format_time(estimated_total)
                time_info = f"已用:{elapsed_str} 剩余:{remaining_str} 总计:{total_str}"
            else:
                time_info = f"已用:{format_time(elapsed_time)}"
            
            print(f"\033[K\r进度: [{bar}] {percent}% ({i}/{len(fasta_files)}) | {time_info} | {gene_name[:20]:<20}", end="", flush=True)
            # ====================================================================

            if not os.path.exists(fasta_path):
                result = {'gene': gene_name, 'status': 'file_not_found'}
                stats['failed'] += 1
            else:
                result = process_gene(fasta_path, gene_name, temp_dir, i)
                if result['status'] == 'success':
                    stats['success'] += 1
                else:
                    stats['failed'] += 1

            results_list.append(result)
            if i % 5 == 0 or i == len(fasta_files):
                pd.DataFrame(results_list, columns=columns).to_csv(output_csv, index=False)

    results_df = pd.DataFrame(results_list, columns=columns)
    end_time = datetime.now()
    total_seconds = (end_time - start_time).total_seconds()
    total_time_str = format_time(total_seconds)

    print(f"\n✅ 分析完成!")
    print(f"⏱️ 总用时: {total_time_str}")
    print(f"📊 成功: {stats['success']}")
    print(f"❌ 失败: {stats['failed']}")

    generate_statistics_report(results_df, output_csv, start_time, end_time, total_seconds, stats)
    return 0

if __name__ == "__main__":
    sys.exit(main())