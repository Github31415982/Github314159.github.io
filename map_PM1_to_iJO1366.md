
# PM1 板底物映射到 iJO1366 模型脚本

以下内容演示如何将 Biolog PM1 碳源板（Phenotype MicroArrays™ PM1）中的底物名称映射到 E. coli iJO1366 模型中的代谢物及其对应的交换反应（exchange reaction）。  
请确保已安装以下 Python 库：
```bash
pip install cobra pandas fuzzywuzzy python-Levenshtein
```

## 脚本说明

1. **读取 iJO1366 模型**  
   使用 `cobra` 库加载 SBML 格式模型文件 `iJO1366.xml`，该文件需与本脚本在同一目录下。  
2. **PM1 板底物列表**  
   将所有 PM1 碳源板的底物手动提取到 `pm1_substrates` 列表中，包含 “Negative Control”。  
3. **建立索引**  
   - `met_name_to_obj`：以 `met.name.lower()` 作为键，指向 `cobra.Metabolite` 对象。  
   - `met_id_to_obj`：以 `met.id.lower()` （即 Bigg ID）作为键，指向 `cobra.Metabolite` 对象。  
4. **匹配函数 `match_metabolite`**  
   - 精确匹配 `met.name`；  
   - 精确匹配 `met.id`（对输入名称做简单格式化）；  
   - 若仍未匹配，则使用 `fuzzywuzzy` 模糊匹配，阈值可自定义（默认 90）。  
5. **查找 Exchange Reaction**  
   在 iJO1366 模型中，exchange reactions 通常以 `EX_` 开头，且反应仅含一个外部代谢物。  
6. **输出结果**  
   将匹配结果保存到 CSV 文件 `PM1_to_iJO1366_mapping.csv`，包含每个底物对应的代谢物 ID 及交换反应 ID。

---

```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cobra
import pandas as pd
from fuzzywuzzy import fuzz, process
import os
import re

# ------------- 1. 读取 iJO1366 SBML 模型 -------------
MODEL_PATH = "iJO1366.xml"
if not os.path.isfile(MODEL_PATH):
    raise FileNotFoundError(f"未找到模型文件：{MODEL_PATH}。请确认 iJO1366.xml 与本脚本在同一目录下。")

print("正在加载 iJO1366 模型，请稍候……")
model = cobra.io.read_sbml_model(MODEL_PATH)
print(f"模型加载完成，共 {len(model.metabolites)} 个代谢物，{len(model.reactions)} 个反应。")

# ------------- 2. 定义 PM1 板底物列表 -------------
pm1_substrates = [
    "Negative Control",
    "L-Arabinose",
    "N-Acetyl-D-Glucosamine",
    "D-Saccharic Acid",
    "Succinic Acid",
    "D-Galactose",
    "L-Aspartic Acid",
    "L-Proline",
    "D-Alanine",
    "D-Trehalose",
    "D-Mannose",
    "Dulcitol",
    "D-Serine",
    "D-Sorbitol",
    "Glycerol",
    "L-Fucose",
    "D-Glucuronic Acid",
    "D-Gluconic Acid",
    "D,L-α-Glycerol-Phosphate",
    "D-Xylose",
    "L-Lactic Acid",
    "Formic Acid",
    "D-Mannitol",
    "L-Glutamic Acid",
    "D-Glucose-6-Phosphate",
    "D-Galactonic Acid-γ-Lactone",
    "D,L-Malic Acid",
    "D-Ribose",
    "Tween 20",
    "L-Rhamnose",
    "D-Fructose",
    "Acetic Acid",
    "α-D-Glucose",
    "Maltose",
    "D-Melibiose",
    "Thymidine",
    "L-Asparagine",
    "D-Aspartic Acid",
    "D-Glucosaminic Acid",
    "1,2-Propanediol",
    "Tween 40",
    "α-Keto-Glutaric Acid",
    "α-Keto-Butyric Acid",
    "α-Methyl-D-Galactoside",
    "α-D-Lactose",
    "Lactulose",
    "Sucrose",
    "Uridine",
    "L-Glutamine",
    "m-Tartaric Acid",
    "D-Glucose-1-Phosphate",
    "D-Fructose-6-Phosphate",
    "Tween 80",
    "α-Hydroxy Glutaric Acid-γ-Lactone",
    "α-Hydroxy Butyric Acid",
    "β-Methyl-D-Glucoside",
    "Adonitol",
    "Maltotriose",
    "2-Deoxy Adenosine",
    "Adenosine",
    "Glycyl-L-Aspartic Acid",
    "Citric Acid",
    "myo-Inositol",
    "D-Threonine",
    "Fumaric Acid",
    "Bromo Succinic Acid",
    "Propionic Acid",
    "Mucic Acid",
    "Glycolic Acid",
    "Glyoxylic Acid",
    "D-Cellobiose",
    "Inosine",
    "Glycyl-L-Glutamic Acid",
    "Tricarballylic Acid",
    "L-Serine",
    "L-Threonine",
    "L-Alanine",
    "L-Alanyl-Glycine",
    "Acetoacetic Acid",
    "N-Acetyl-β-D-Mannosamine",
    "Mono Methyl Succinate",
    "Methyl Pyruvate",
    "D-Malic Acid",
    "L-Malic Acid",
    "Glycyl-L-Proline",
    "p-Hydroxy Phenyl Acetic Acid",
    "m-Hydroxy Phenyl Acetic Acid",
    "Tyramine",
    "D-Psicose",
    "L-Lyxose",
    "Glucuronamide",
    "Pyruvic Acid",
    "L-Galactonic Acid-γ-Lactone",
    "D-Galacturonic Acid",
    "Phenylethyl-amine",
    "2-Aminoethanol"
]

print(f"PM1 底物列表共 {len(pm1_substrates)} 项（含 Negative Control）。")

# ------------- 3. 构建代谢物名称 → model.metabolite 的索引 -------------
met_name_to_obj = {}
met_id_to_obj = {}
for met in model.metabolites:
    nm = met.name.lower().strip()
    met_name_to_obj[nm] = met
    met_id_to_obj[met.id.lower().strip()] = met

print(f"建立了 {len(met_name_to_obj)} 个“名称→Metabolite”索引，{len(met_id_to_obj)} 个“ID→Metabolite”索引。")

# ------------- 4. 定义匹配函数 -------------
from fuzzywuzzy import process, fuzz

def match_metabolite(substrate_name, name_dict, id_dict, threshold=90):
    key = substrate_name.lower().strip()
    # 精确匹配 met.name
    if key in name_dict:
        return name_dict[key], "name_exact"
    # 精确匹配 met.id
    candidate_id = re.sub(r'[^0-9a-z]', '_', key)
    candidate_id = re.sub(r'_+', '_', candidate_id).strip('_')
    if candidate_id in id_dict:
        return id_dict[candidate_id], "id_exact"
    # 模糊匹配
    all_names = list(name_dict.keys())
    best_match, score = process.extractOne(key, all_names, scorer=fuzz.ratio)
    if score >= threshold:
        return name_dict[best_match], f"name_fuzzy({score})"
    return None, None

# ------------- 5. 查找 exchange reaction 的辅助函数 -------------
def find_exchange_reactions(met_obj):
    ex_reactions = []
    for rxn in model.reactions:
        if rxn.id.startswith("EX_") and len(rxn.metabolites) == 1:
            mets = list(rxn.metabolites.keys())
            if met_obj in mets:
                ex_reactions.append(rxn)
    return ex_reactions

# ------------- 6. 对 PM1 底物列表逐条匹配 -------------
mapping_results = []
for substrate in pm1_substrates:
    if substrate.lower().startswith("negative control"):
        mapping_results.append({
            "PM1_substrate": substrate,
            "matched_met_id": None,
            "matched_met_name": None,
            "match_method": "negative_control",
            "exchange_rxn_ids": []
        })
        continue

    met_obj, match_flag = match_metabolite(substrate, met_name_to_obj, met_id_to_obj, threshold=90)
    if met_obj is None:
        mapping_results.append({
            "PM1_substrate": substrate,
            "matched_met_id": None,
            "matched_met_name": None,
            "match_method": "unmatched",
            "exchange_rxn_ids": []
        })
        continue

    ex_rxns = find_exchange_reactions(met_obj)
    ex_ids = [rxn.id for rxn in ex_rxns]
    mapping_results.append({
        "PM1_substrate": substrate,
        "matched_met_id": met_obj.id,
        "matched_met_name": met_obj.name,
        "match_method": match_flag,
        "exchange_rxn_ids": ex_ids
    })

# ------------- 7. 保存到 CSV -------------
df = pd.DataFrame(mapping_results)
out_csv = "PM1_to_iJO1366_mapping.csv"
df.to_csv(out_csv, index=False, encoding="utf-8-sig")

print(f"映射完成，共 {len(mapping_results)} 条记录。结果已保存到：{out_csv}")
print(df.head(10))
```

---

## 使用说明

1. 将本文档保存为 `map_PM1_to_iJO1366.md`。  
2. 将 `iJO1366.xml` 与 `map_PM1_to_iJO1366.md` 放在同一目录。  
3. 在终端中运行脚本：  
   ```bash
   python map_PM1_to_iJO1366.md
   ```  
   （可根据需要将文件扩展名改为 `.py` 后再执行）  
4. 脚本执行后会生成 `PM1_to_iJO1366_mapping.csv`，其中包含每个 PM1 底物对应的代谢物 ID 及 Exchange Reaction ID。

请将此 Markdown 文件下载并在本地执行，完成映射工作。
