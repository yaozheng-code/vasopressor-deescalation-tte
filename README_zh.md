# NE/AVP 首次有临床意义减量：Target Trial Emulation 复现代码

本仓库对应论文的最终 **MIMIC-IV v3.1** 分析，采用
**clone–censor–weight（CCW）+ IPCW + 加权 Kaplan–Meier**。

**核心区别**：研究的是达到共同、治疗前稳定决策点后，先对
**NE** 还是 **AVP** 实施首次有临床意义的减量，**不是**
回顾性地按照最后先停掉哪种药分组。

## 主要文件

- `analysis/G10_v11_MIMIC_NE_AVP_final_CCW_subject_bootstrap.py`：
  原样保存的 G10 v1.1，完成数据提取、cohort 识别、
  克隆/人工删失/IPCW、6 h 主结果、敏感性分析及
  **500 次 patient-cluster bootstrap**。
- `analysis/G11_MIMIC_NE_AVP_manuscript_Table1.py`：
  冻结 cohort 的 Table 1；无须再跑 bootstrap。
- `run_pipeline.py`：统一数据库连接、输出目录和运行入口，
  **不修改 G10/G11 的科研算法**。
- `preflight.py`：只读检查数据库和 derived concepts。
- `verify_results.py`：核对结果是否接近论文锁定值。
- `tests/`：不使用真实患者数据的合成边界条件测试。

## 本地运行

先获取合法授权并在 PostgreSQL 安装 MIMIC-IV 3.1，
保证存在 `mimiciv_icu`、`mimiciv_hosp` 和
`mimiciv_derived.sepsis3`、`mimiciv_derived.vitalsign`。

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 完全离线，不访问患者记录
python check_integrity.py
python run_pipeline.py --self-test-only
python -m unittest discover -s tests -v

# 只读数据库结构检查
python preflight.py --db-name mimiciv --db-host /var/run/postgresql

# 正式复现 G10 500-bootstrap + G11 Table 1
python run_pipeline.py --step all --n-bootstrap 500 --seed 20261015
```

默认写入 Git 忽略的 `outputs/<timestamp>/`。运行完成自动对照
`expected_results.json`。也可以执行：

```bash
python verify_results.py \
  --g10-dir outputs/<timestamp>/g10 \
  --g11-dir outputs/<timestamp>/g11
```

主要锁定结果：2260 stays、2168 patients；
NE-first 6h 风险约 19.74%，AVP-first 约 18.45%，
RD=+1.29 百分点（95%CI −8.07 到 +10.38）。
**统计结果不能证明两种策略等效。**

## GitHub 安全注意事项

**只能公开代码，不能公开 MIMIC-IV 原始患者记录、患者级衍生数据、
数据库密码和工作目录中的实际输出。**

项目已带 `.gitignore`，但上传 GitHub 前仍应人工核对文件。
PhysioNet 对衍生数据的公开分享有独立限制。
参考 `docs/DATA_ACCESS_AND_PRIVACY.md` 和
`RELEASE_CHECKLIST.md`。

**重要：**原始研究的确切 Python 包版本及官方 derived SQL
commit 未从旧运行日志完全恢复，因此这份代码可完整重跑并核对
主要结果，但不能承诺不同计算环境下一定逐位相同。
设计并非最初就 preregistered，而是在前期 pilot 审计后冻结
最终 G10 推断分析。这一点已在文档中透明说明。
