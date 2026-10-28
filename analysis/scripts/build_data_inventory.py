"""Build the curated public-data inventory table for the OP x CVD comorbidity study.

Every GEO row's title / sample count / platform is re-fetched from NCBI E-utilities
so the table cannot drift from the source of truth. Non-GEO rows are hand-curated
and were verified by the probe scripts in this folder.

Output: results/tables/public_dataset_inventory.csv
"""
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path("results/tables")
OUT.mkdir(parents=True, exist_ok=True)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# accession -> (layer, role, compartment, design/contrast, priority, note)
GEO = {
    "GSE56815": ("L3-转录组(骨/免疫)", "发现集", "外周血单核细胞", "低 vs 高 髋部BMD（绝经前后女性）", 1, "HG-U133A；骨质疏松血液转录组的经典数据集"),
    "GSE56814": ("L3-转录组(骨/免疫)", "重复集(平台)", "外周血单核细胞", "同上队列，HuEx-1_0-st 平台", 2, "与 GSE56815 同队列、不同平台"),
    "GSE56816": ("L3-转录组(骨/免疫)", "超级系列", "外周血单核细胞", "GSE56814+GSE56815 合并 n=153", 3, "样本量最大的合并版本"),
    "GSE7158": ("L3-转录组(骨/免疫)", "验证集", "外周血单核细胞", "极高 vs 极低 峰值骨量", 2, "独立队列，可做外部验证"),
    "GSE7429": ("L3-转录组(骨/免疫)", "验证集", "外周血B淋巴细胞", "骨质疏松 vs 对照（绝经后女性）", 3, "免疫细胞亚群补充"),
    "GSE13850": ("L3-转录组(骨/免疫)", "验证集", "外周血B淋巴细胞", "吸烟相关骨质疏松 低 vs 高BMD", 3, "含吸烟混杂因素"),
    "GSE35956": ("L3-转录组(骨/免疫)", "骨组织层面", "股骨头来源 hMSC", "原发骨质疏松 vs 年龄匹配对照", 1, "真正的骨组织来源细胞，n 偏小"),
    "GSE51495": ("L3-转录组(骨/免疫)", "桥接集", "PBMC + 皮质骨（同一批人）", "PBMC 与骨组织配对表达谱", 1, "少见的血-骨配对设计，适合验证血液替代标志物"),
    "GSE141614": ("L3-转录组(骨/免疫)", "骨组织层面", "人骨活检 RNA-seq", "骨重建偶联因子", 2, "骨组织原位转录组"),
    "GSE80614": ("L3-转录组(骨/免疫)", "程序参考", "hMSC", "成骨 vs 成脂分化时间序列", 2, "定义成骨程序基因集"),
    "GSE100927": ("L3-转录组(血管)", "发现集", "颈/股/膝下外周动脉", "动脉粥样硬化 vs 对照", 1, "三血管床、样本量大，共病研究主力"),
    "GSE43292": ("L3-转录组(血管)", "重复集", "颈动脉", "斑块 vs 肉眼正常内膜", 1, "经典颈动脉斑块谱"),
    "GSE28829": ("L3-转录组(血管)", "验证集", "颈动脉", "早期 vs 进展期斑块", 2, "可与 GSE43292 联用"),
    "GSE11138": ("L3-转录组(血管)", "验证集", "LAD 冠状动脉斑块", "冠心病斑块 vs 对照", 2, "冠脉血管床，最贴近 CAD"),
    "GSE57691": ("L3-转录组(血管)", "验证集", "腹主动脉", "AAA vs 主动脉闭塞性疾病", 3, "非斑块性血管病变，作对照"),
    "GSE120521": ("L3-转录组(血管)", "机制补充", "颈动脉斑块", "稳定 vs 不稳定区域", 3, "斑块易损性维度"),
    "GSE221911": ("L3-转录组(血管)", "血液桥接", "全血 RNA-seq", "冠心病 vs 对照", 1, "血液层面共病信号"),
    "GSE23561": ("L3-转录组(血管)", "血液桥接", "外周血", "代谢综合征/CAD/T2D vs 对照", 2, "含代谢共病对照"),
    "GSE20129": ("L3-转录组(血管)", "血液桥接", "外周血（女性）", "动脉粥样硬化 vs 非动脉粥样硬化", 2, "含种族因素"),
    "GSE37558": ("L4-钙化/成骨转分化", "核心机制集", "钙化血管平滑肌细胞 vs 成骨细胞", "两种成骨程序直接对比", 1, "★直接比较血管钙化与骨成骨程序"),
    "GSE11917": ("L4-钙化/成骨转分化", "核心机制集", "人冠脉平滑肌细胞", "磷酸盐+维生素D甾醇诱导钙化 时间/剂量", 1, "★VDR/矿物质轴机制"),
    "GSE288136": ("L4-钙化/成骨转分化", "遗传-机制连接", "血管平滑肌细胞", "9p21.3 CAD 风险位点驱动成骨转化", 2, "★把 CAD 遗传信号接到成骨程序"),
    "GSE84012": ("L4-钙化/成骨转分化", "机制补充", "不同动脉床 VSMC", "血管床异质性对钙化的影响", 2, ""),
    "GSE284683": ("L4-钙化/成骨转分化", "机制补充", "原代人冠脉 SMC", "钙化培养基 ± PIKfyve 抑制", 3, ""),
    "GSE309462": ("L4-钙化/成骨转分化", "机制补充", "颈动脉 VSMC", "内皮 BMP6 驱动钙化", 3, "BMP 通路"),
    "GSE102249": ("L4-钙化/成骨转分化", "CAVD队列", "人主动脉瓣（240例）", "纤维钙化性重构瓣膜转录组", 2, "样本量最大的钙化瓣膜队列"),
    "GSE148219": ("L4-钙化/成骨转分化", "CAVD队列", "人主动脉瓣", "胎儿 vs 钙化成人瓣膜", 3, ""),
    "GSE271848": ("L4-钙化/成骨转分化", "骨-血管对话", "BMSC + 衰老内皮", "衰老血管微环境经旁分泌抑制 Wnt", 3, "★直接研究骨-血管 cross-talk"),
    "GSE159677": ("L5-单细胞", "发现集", "颈动脉钙化斑块核心 vs 邻近组织", "斑块 core vs proximal adjacent", 1, "★唯一直接标注'钙化斑块'的 scRNA-seq"),
    "GSE131778": ("L5-单细胞", "验证集", "人冠脉/斑块", "SMC 表型调控 in vivo", 2, "与 GSE155512 同源"),
    "GSE155512": ("L5-单细胞", "验证集", "人斑块 SMC", "SMC 表型转换新状态", 3, ""),
    "GSE224273": ("L5-单细胞", "验证集", "颈动脉斑块免疫细胞", "症状性 vs 无症状斑块", 3, "免疫维度"),
    "GSE317068": ("L5-单细胞", "骨侧图谱", "人骨/骨内膜（scRNA-seq）", "骨骼疾病细胞与遗传决定因素", 1, "★骨侧单细胞参考，用于映射共享基因的细胞来源"),
    "GSE317069": ("L5-单细胞", "骨侧图谱", "人骨（SuperSeries）", "含 bulk/scATAC 等多组学", 2, ""),
    "GSE221615": ("L2-临床表征", "队列桥接", "全血 RNA-seq（391人, PESA）", "亚临床动脉粥样硬化负荷 vs 表观遗传年龄加速", 2, "★带影像学表型的队列，可做'临床表征→转录组'"),
}

GENETICS = [
    ("L2-临床表征/遗传", "结局GWAS", "FinnGen R12", "M13_OSTEOPOROSIS", "骨质疏松（M13）", "10,461 例 / 473,264 对照", "芬兰人群；可与 CVD 终点做遗传相关/MR"),
    ("L2-临床表征/遗传", "结局GWAS", "FinnGen R12", "I9_CORATHER", "冠状动脉粥样硬化", "63,307 / 416,171", "CAD 主终点"),
    ("L2-临床表征/遗传", "结局GWAS", "FinnGen R12", "I9_ATHSCLE", "动脉粥样硬化（除外脑/冠脉）", "19,807 / 463,106", "外周动脉粥样硬化"),
    ("L2-临床表征/遗传", "结局GWAS", "FinnGen R12", "I9_MI_STRICT", "心肌梗死（严格定义）", "31,666 / 416,171", ""),
    ("L2-临床表征/遗传", "结局GWAS", "FinnGen R12", "I9_AF", "房颤/扑动", "63,532 / 252,810", ""),
    ("L2-临床表征/遗传", "结局GWAS", "FinnGen R12", "I9_HEARTFAIL", "心力衰竭（严格）", "37,653 / 462,695", ""),
    ("L2-临床表征/遗传", "暴露GWAS", "GWAS Catalog", "GCST90029004", "足跟骨密度（eBMD）", "n≈583,314", "Loh 2018；UKB eBMD，MR 常用暴露"),
    ("L2-临床表征/遗传", "暴露GWAS", "GWAS Catalog", "GCST90018887", "骨质疏松", "n≈479,981", "Sakaue 2021；跨人群"),
    ("L2-临床表征/遗传", "结局GWAS", "GWAS Catalog", "GCST90132315", "冠状动脉疾病（CAD）", "n≈1,167,328", "Aragam 2022；CAD 最大样本量之一"),
    ("L2-临床表征/遗传", "结局GWAS", "GWAS Catalog", "GCST90018877", "心肌梗死", "n≈440,906", "Sakaue 2021"),
    ("L2-临床表征/遗传", "结局GWAS", "GWAS Catalog", "GCST90104534", "卒中", "n≈1,234,808", "GIGASTROKE/Mishra 2022"),
    ("L2-临床表征/遗传", "结局GWAS", "GWAS Catalog", "GCST90162626", "心力衰竭", "n≈1,270,968", "HERMES/Levin 2022"),
]

OTHER = [
    ("L1-流行病学", "主数据集", "CDC/NCHS NHANES", "DXX_D…DXX_J (2005–2018 六个周期)", "DXA 股骨颈/全髋/腰椎 BMD",
     "每周期约 5,000–7,000 名 8 岁以上受检者", "已验证 6 个周期 XPT 均可下载；需自行合并 MCQ/BPQ/DEMO 等文件"),
    ("L1-流行病学", "结局随访", "NCHS 关联死亡文件", "NHANES_2005_2006 … _2017_2018_MORT_2019_PUBLIC",
     "全因/心血管死亡（ICD-10）", "与 NHANES 同", "已验证 2005–2008 文件可用；2019 随访截止"),
    ("L1-流行病学", "前瞻队列（需申请）", "UK Biobank", "Field 78（足跟 BMD T 值）等",
     "BMD + ICD-10 心血管结局", "278,763 名有足跟 BMD", "需申请与付费；可与 NHANES 形成欧美双队列"),
    ("L6-参考/注释", "正常组织参考", "GTEx v10", "bulk-gex/v10", "54 组织 / 946 供者（含主动脉、冠状动脉、胫动脉、骨骼肌、全血）",
     "基因 TPM 矩阵 2.26 GB；另有按组织 counts", "已验证；v8 可用于 eQTL 复现"),
    ("L6-参考/注释", "靶点-疾病证据", "Open Targets Platform", "GraphQL v4", "遗传关联/L2G/表达证据", "-", "已验证可用"),
    ("L6-参考/注释", "通路与网络", "KEGG / Reactome / STRING / Enrichr / MSigDB", "REST API + R 包", "通路富集、PPI、药物预测", "-", "已验证连通"),
]

rows = []


def fetch(acc):
    js = json.loads(urllib.request.urlopen(
        f"{EUTILS}/esearch.fcgi?db=gds&term={acc}%5BAccession%5D&retmode=json", timeout=60).read())
    ids = js.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return {}
    js2 = json.loads(urllib.request.urlopen(
        f"{EUTILS}/esummary.fcgi?db=gds&id={ids[0]}&retmode=json", timeout=60).read())
    res = js2.get("result", {})
    for u in res.get("uids", []):
        if res[u].get("accession") == acc:
            return res[u]
    return {}


for acc, (layer, role, comp, design, prio, note) in GEO.items():
    rec = fetch(acc)
    time.sleep(0.35)
    gpl = rec.get("gpl", "")
    if isinstance(gpl, list):
        gpl = ";".join(g.get("accession", "") for g in gpl)
    rows.append({
        "layer": layer, "role": role, "source": "GEO", "accession": acc,
        "modality": rec.get("gdstype", ""),
        "compartment": comp,
        "n_samples": rec.get("n_samples", ""),
        "platform": gpl,
        "design": design,
        "organism": rec.get("taxon", ""),
        "pubmed": (rec.get("pubmedids") or [""])[0],
        "released": rec.get("pdat", ""),
        "title": (rec.get("title") or "").replace("\n", " "),
        "priority": prio,
        "access": "公开，直接下载",
        "note": note,
    })
    print(f"  {acc:11s} n={str(rows[-1]['n_samples']):>4}  {rows[-1]['title'][:70]}")

for layer, role, source, acc, design, size, note in GENETICS:
    rows.append({"layer": layer, "role": role, "source": source, "accession": acc, "modality": "GWAS summary statistics",
                 "compartment": "summary stats", "n_samples": size, "platform": "", "design": design,
                 "organism": "Humans (European/Finland)", "pubmed": "", "released": "R12/2022-2025", "title": design,
                 "priority": 1, "access": "公开，直接下载", "note": note})

for layer, role, source, acc, design, size, note in OTHER:
    rows.append({"layer": layer, "role": role, "source": source, "accession": acc, "modality": "multi-omics / tabular",
                 "compartment": "-", "n_samples": size, "platform": "", "design": design, "organism": "Humans",
                 "pubmed": "", "released": "见 note", "title": design, "priority": 1,
                 "access": "公开" if "申请" not in role else "需申请/付费", "note": note})

fields = ["layer", "role", "source", "accession", "modality", "compartment", "n_samples", "platform",
          "design", "organism", "pubmed", "released", "title", "priority", "access", "note"]
rows.sort(key=lambda r: (r["layer"], str(r["priority"]), r["accession"]))
with open(OUT / "public_dataset_inventory.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)
print(f"\n{len(rows)} rows -> {OUT/'public_dataset_inventory.csv'}")
