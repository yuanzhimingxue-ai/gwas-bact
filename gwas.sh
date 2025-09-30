#!/usr/bin/env bash
set -euo pipefail

# ================ 参数与目录 ================
INPUT_DIR=${1:-}
if [ -z "${INPUT_DIR:-}" ] || [ ! -d "$INPUT_DIR" ]; then
  echo "用法: $0 /path/to/input_dir   # 目录内含 cleaned FASTQ、参考基因组、metadata.csv"; exit 1
fi
INPUT_DIR=$(readlink -f "$INPUT_DIR")

# 自动识别参考 & 元数据
REF=$(ls "$INPUT_DIR"/*.fna "$INPUT_DIR"/*.fa "$INPUT_DIR"/*.fasta 2>/dev/null | head -n1 || true)
META="$INPUT_DIR/metadata.csv"
[ -f "${REF:-}" ] || { echo "[ERR] 未在 $INPUT_DIR 找到参考基因组(.fna/.fa/.fasta)"; exit 1; }
[ -f "$META" ]    || { echo "[ERR] 未在 $INPUT_DIR 找到 metadata.csv"; exit 1; }

THREADS=${THREADS:-$(nproc || echo 8)}

# 工程目录放在输入目录里，完全自包含
PROJECT="$INPUT_DIR/gwas_project"
WORK=$PROJECT/work
BAM_DIR=$WORK/bam
GVCF_DIR=$WORK/gvcf
LOGDIR=$WORK/logs
OUT=$PROJECT/results
mkdir -p "$BAM_DIR" "$GVCF_DIR" "$LOGDIR" "$OUT"

echo "[INFO] INPUT_DIR=$INPUT_DIR"
echo "[INFO] REF=$REF"
echo "[INFO] META=$META"
echo "[INFO] PROJECT=$PROJECT"
echo "[INFO] THREADS=$THREADS"

# ================ 工具检测 ================
need() { command -v "$1" >/dev/null 2>&1 || { echo "[MISS] $1"; MISSING=1; }; }
MISSING=0
need bwa-mem2 || need bwa
need samtools
need gatk
need bcftools
need plink2
need tabix
need python
need Rscript
[ "${MISSING:-0}" = 0 ] || { echo "[ERR] 有工具缺失，请检查 Docker / conda 环境"; exit 1; }

# ================ 参考索引 ================
[ -f "${REF}.fai" ] || samtools faidx "$REF"
DICT=${REF%.*}.dict
[ -f "$DICT" ] || gatk CreateSequenceDictionary -R "$REF" >/dev/null
if command -v bwa-mem2 >/dev/null 2>&1; then
  ls "${REF}".0123 "${REF}".bwt.2bit.64 >/dev/null 2>&1 || bwa-mem2 index "$REF"
else
  ls "${REF}".bwt >/dev/null 2>&1 || bwa index "$REF"
fi

# ================ 样本清单 ================
SAMPLE_TSV=$WORK/samples.tsv
if [ ! -s "$SAMPLE_TSV" ]; then
  echo -e "SID\tPREFIX" > "$SAMPLE_TSV"
  for r1 in "$INPUT_DIR"/*_clean_1.fq.gz; do
    [ -e "$r1" ] || continue
    base=$(basename "$r1")
    sid=${base%%_*}                # C10
    pref=${base%_1.fq.gz}          # C10_clean
    echo -e "${sid}\t${pref}" >> "$SAMPLE_TSV"
  done
fi
echo "[INFO] 样本数：$(($(wc -l < "$SAMPLE_TSV")-1))"

# ================ 比对/去重 ================
while read -r SID PRE; do
  [ "$SID" = "SID" ] && continue
  R1=$INPUT_DIR/${PRE}_1.fq.gz
  R2=$INPUT_DIR/${PRE}_2.fq.gz
  SORTED=$BAM_DIR/${SID}.sorted.bam
  MARKDUP=$BAM_DIR/${SID}.markdup.bam
  METR=$BAM_DIR/${SID}.markdup.metrics.txt
  [ -s "$MARKDUP" ] && continue

  if [ ! -s "$SORTED" ]; then
    if command -v bwa-mem2 >/dev/null 2>&1; then
      bwa-mem2 mem -t "$THREADS" "$REF" "$R1" "$R2" | samtools sort -@ "$THREADS" -o "$SORTED"
    else
      bwa mem -t "$THREADS" "$REF" "$R1" "$R2" | samtools sort -@ "$THREADS" -o "$SORTED"
    fi
    samtools index "$SORTED"
  fi

  gatk MarkDuplicatesSpark -I "$SORTED" -O "$MARKDUP" -M "$METR" --remove-sequencing-duplicates false >/dev/null
  samtools index "$MARKDUP"
done < "$SAMPLE_TSV"

# ================ HaplotypeCaller（gVCF） ================
while read -r SID PRE; do
  [ "$SID" = "SID" ] && continue
  BAM=$BAM_DIR/${SID}.markdup.bam
  OUTG=$GVCF_DIR/${SID}.g.vcf.gz
  LOG=$LOGDIR/${SID}.hc.log
  [ -s "$BAM" ] || { echo "[SKIP] $SID 缺 BAM"; continue; }
  [ -s "$OUTG" ] && { tabix -f -p vcf "$OUTG" || true; continue; }

  gatk --java-options "-Xms2g -Xmx8g" HaplotypeCaller \
    -R "$REF" -I "$BAM" \
    -O "$OUTG" -ERC GVCF \
    --sample-ploidy 1 \
    --native-pair-hmm-threads 4 2> "$LOG"
  tabix -f -p vcf "$OUTG"
done < "$SAMPLE_TSV"

# ================ 联合分型/合并 ================
if ls "$OUT"/joint.*.vcf.gz >/dev/null 2>&1; then
  ls "$OUT"/joint.*.vcf.gz | sort -V > "$OUT"/joint.list
  [ -s "$OUT/joint.all.vcf.gz" ] || { bcftools concat -f "$OUT/joint.list" -a -O z -o "$OUT/joint.all.vcf.gz"; tabix -f -p vcf "$OUT/joint.all.vcf.gz"; }
else
  if [ ! -s "$OUT/joint.all.vcf.gz" ]; then
    gvcf_list=$(ls "$GVCF_DIR"/*.g.vcf.gz | sed 's/^/-V /' | tr '\n' ' ')
    gatk --java-options "-Xmx8g" GenotypeGVCFs -R "$REF" $gvcf_list -O "$OUT/joint.all.vcf.gz"
    tabix -f -p vcf "$OUT/joint.all.vcf.gz"
  fi
fi

# ================ 取 SNP + 过滤 ================
[ -s "$OUT/joint.SNP.vcf.gz" ] || {
  gatk SelectVariants -R "$REF" -V "$OUT/joint.all.vcf.gz" --select-type-to-include SNP -O "$OUT/joint.SNP.vcf.gz"
  tabix -f -p vcf "$OUT/joint.SNP.vcf.gz"
}
[ -s "$OUT/joint.SNP.final.vcf.gz" ] || {
  gatk VariantFiltration -R "$REF" -V "$OUT/joint.SNP.vcf.gz" \
    --filter-expression "QD < 2.0"  --filter-name "LowQD" \
    --filter-expression "FS > 60.0" --filter-name "HighFS" \
    --filter-expression "MQ < 30.0" --filter-name "LowMQ" \
    --filter-expression "SOR > 3.0" --filter-name "HighSOR" \
    -O "$OUT/joint.SNP.filt.vcf.gz"
  tabix -f -p vcf "$OUT/joint.SNP.filt.vcf.gz"
  bcftools view -f PASS -O z -o "$OUT/joint.SNP.final.vcf.gz" "$OUT/joint.SNP.filt.vcf.gz"
  tabix -f -p vcf "$OUT/joint.SNP.final.vcf.gz"
}

# ================ VCF -> PLINK & QC ================
[ -s "$OUT/bacteria_snp.bed" ] || {
  plink2 --vcf "$OUT/joint.SNP.final.vcf.gz" dosage=DS --max-alleles 2 --make-bed --out "$OUT/bacteria_snp"
}
plink2 --bfile "$OUT/bacteria_snp" --missing sample-only --out "$OUT/bacteria_snp.hap.qc.miss" >/dev/null
[ -s "$OUT/bacteria_snp.qc.bed" ] || {
  plink2 --bfile "$OUT/bacteria_snp" \
    --mind 0.2 --geno 0.1 --maf 0.01 --hwe 1e-6 midp \
    --make-bed --out "$OUT/bacteria_snp.qc"
}

# ================ PCA（小样本：thin + 读频率） ================
[ -s "$OUT/bacteria_thin.snplist" ] || plink2 --bfile "$OUT/bacteria_snp.qc" --thin-count 10000 --write-snplist --out "$OUT/bacteria_thin"
[ -s "$OUT/bacteria_thin.afreq" ]   || plink2 --bfile "$OUT/bacteria_snp.qc" --extract "$OUT/bacteria_thin.snplist" --freq --out "$OUT/bacteria_thin"
[ -s "$OUT/bacteria_pca.eigenvec" ] || plink2 --bfile "$OUT/bacteria_snp.qc" --extract "$OUT/bacteria_thin.snplist" --read-freq "$OUT/bacteria_thin.afreq" --pca approx 10 --out "$OUT/bacteria_pca"

# ================ 表型生成（Python） ================
python - <<PY
import pandas as pd, os
OUT="$OUT"
META="$META"
df=pd.read_csv(META)
df.columns=[c.strip().replace(' ','_') for c in df.columns]
cand=[c for c in df.columns if c.lower() in ('sample_id','sampleid','id','sample','iid')]
assert cand, "metadata.csv 需要有样本ID列，如 sample_id/ID/sample/IID"
df=df.rename(columns={cand[0]:'IID'}); df['FID']=df['IID']
num=[c for c in df.columns if c not in ('FID','IID') and pd.api.types.is_numeric_dtype(df[c])]
cat=[c for c in num if pd.api.types.is_integer_dtype(df[c]) and df[c].nunique()<=6]
cont=[c for c in num if c not in cat]
for c in cont:
    ph=df[['FID','IID',c]].copy(); ph.columns=['FID','IID','PHENO']
    ph.to_csv(os.path.join(OUT,f'phen_{c}.txt'),sep='\t',index=False)
for c in cat:
    for cls in sorted(df[c].dropna().unique().tolist()):
        ph=df[['FID','IID']].copy()
        ph['PHENO']=df[c].apply(lambda x: 2 if x==cls else (1 if pd.notna(x) else -9))
        ph.to_csv(os.path.join(OUT,f'phen_{c}_class_{cls}.txt'),sep='\t',index=False)
df[['FID','IID']+cont+cat].to_csv(os.path.join(OUT,'phen_all_wide.txt'),sep='\t',index=False)
print('[OK] phenotypes:', {'cont':cont, 'cat':cat})
PY

# ================ 连续性状 GWAS（线性 + PCs） ================
run_linear () {
  local trait="$1"
  local PH="$OUT/phen_${trait}.txt"
  [ -s "$PH" ] || { echo "[SKIP] 无 $PH"; return 0; }
  plink2 --bfile "$OUT/bacteria_snp.qc" \
    --pheno "$PH" --pheno-name PHENO \
    --covar "$OUT/bacteria_pca.eigenvec" --covar-col-nums 3-7 \
    --covar-variance-standardize \
    --glm hide-covar allow-no-covars \
    --threads "$THREADS" \
    --out "$OUT/gwas_linear_${trait}"
}
# 你可以在此处列出 metadata 里的连续性状名；先尝试常见两列
for trait in virulence growth_rate; do run_linear "$trait"; done

# ================ 分类性状（二分类展开 + Firth 回退） ================
for PHF in "$OUT"/phen_*_class_*.txt; do
  [ -e "$PHF" ] || continue
  base=$(basename "$PHF" .txt)
  plink2 --bfile "$OUT/bacteria_snp.qc" \
    --pheno "$PHF" --pheno-name PHENO \
    --covar "$OUT/bacteria_pca.eigenvec" --covar-col-nums 3-7 \
    --covar-variance-standardize \
    --glm hide-covar allow-no-covars firth-fallback \
    --threads "$THREADS" \
    --out "$OUT/gwas_${base}"
done

# ================ 结果精简 + 阈值 ================
simplify_glm () {
  local infile="$1"; local out="$2"
  awk 'BEGIN{FS=OFS="\t"} NR==1{
         for(i=1;i<=NF;i++){if($i=="#CHROM") c=i; if($i=="POS") p=i; if($i=="ID") d=i; if($i=="P") q=i}
         print "#CHROM","POS","ID","P"; next
       }
       $0!~/^#/ {print $c,$p,$d,$q}' "$infile" > "$out"
}
for f in "$OUT"/gwas_linear_*.PHENO.glm.linear; do
  [ -s "$f" ] && simplify_glm "$f" "${f/.PHENO.glm.linear/.simple.tsv}" || true
done
for f in "$OUT"/gwas_phen_*_class_*.PHENO.glm.logistic; do
  [ -s "$f" ] || continue
  awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++) h[$i]=i; print "#CHROM","POS","ID","P"; next}
       $0!~/^#/ && $(h["TEST"])=="ADD" {print $(h["#CHROM"]),$(h["POS"]),$(h["ID"]),$(h["P"])}' \
       "$f" > "${f/.PHENO.glm.logistic/.simple.tsv}"
done

M=$(wc -l < "$OUT/bacteria_snp.qc.bim")
BON=$(python - <<PY
M=$M; print(0.05/float(M))
PY
)
echo "[INFO] Bonferroni threshold = $BON (alpha=0.05, M=${M})"

for s in "$OUT"/gwas_linear_*.simple.tsv "$OUT"/gwas_phen_*_class_*.simple.tsv; do
  [ -s "$s" ] || continue
  base=$(basename "$s" .simple.tsv)
  awk -v thr="$BON" 'NR==1 || ($4+0<thr && $4>0)' "$s" > "$OUT/${base}.bonf.tsv"
  awk 'NR==1 || $4>0' "$s" | sort -k4,4g | head -n 21 > "$OUT/${base}.top20.tsv" || true
done

# ================ R 绘图（Manhattan / QQ） ================
Rscript /opt/gwas/r_plot.R "$OUT" || true

echo "[ALL DONE] 结果目录: $OUT"
echo "  - GWAS 线性： gwas_linear_*.PHENO.glm.linear / .simple.tsv / .bonf.tsv"
echo "  - 分类性状： gwas_phen_*_class_*.PHENO.glm.logistic / .simple.tsv / .bonf.tsv"
echo "  - PCA：      bacteria_pca.eigenvec"
echo "  - 图：       *.manhattan.png, *.qq.png"

