# Cathy Fit Pipeline

以星系紅移分布進行三階段 Lorentzian 擬合，估計 galaxy cluster 的紅移、重定中心，輸出候選成員計數與診斷圖。支援 LS、Tian、HSC/DEmP 及標準 CSV 輸入，並附離線 DESI／SDSS 光譜配對工具。

主入口是 **`Cathy_fit_pipeline.ipynb`**。此版本保留 notebook 自含函式的使用方式，尚未改造成可 `pip install` 的 Python package。

## 安裝與快速開始

需要 Python 3.10+；本次本機檢查使用 Python 3.12。於此資料夾根目錄執行：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
jupyter lab Cathy_fit_pipeline.ipynb
```

Windows PowerShell 的環境啟用指令為 `.venv\Scripts\Activate.ps1`。

選擇剛安裝依賴的 Python kernel，執行 **Restart Kernel → Run All**。發布版預設讀取 `examples/` 的一個合成 cluster，不需要真實 catalog。預期產生 `cathy_pipeline_outputs/<UTC時間>/Synthetic_W110/`，包含表格、manifest、摘要及 fitting 圖。

`requirements.txt` 列出直接依賴；`requirements-tested.txt` 記錄本機檢查時的科學套件版本，不是完整跨平台 lockfile。需要對齊這次科學套件版本時，在安裝上述依賴後執行：

```bash
python -m pip install -r requirements-tested.txt
```

## 檔案結構

```text
cathy-fit-pipeline/
├── README.md
├── requirements.txt
├── requirements-tested.txt
├── .gitignore
├── Cathy_fit_pipeline.ipynb
├── cathy_specz_enrichment.py
├── test_cathy_specz_enrichment.py
├── Cathy_fit_pipeline_guide.md
├── Cathy_specz_server_guide.md
├── Cathy_fit_pipeline_checks.md
├── PACKAGING_CHECKS.md
├── Cathy_fit_pipeline_validation.json  # 歷史驗證紀錄
├── Cathy_fit_pipeline_regression.json  # 歷史回歸紀錄
├── examples/                          # 可直接執行的合成資料
├── cathy_pipeline_input_templates/    # 只有標頭的輸入範本
└── data/                              # 自行準備的真實資料，不納入 Git
```

## 換成自己的資料

在 notebook「執行設定」cell 修改 `DATA_ROOT` 和 `DATASETS`。標準格式範例：

```python
DATA_ROOT = Path.cwd() / 'data'
DATASETS = [
    {'name': 'MySample', 'format': 'standard',
     'clusters': DATA_ROOT / 'clusters.csv',
     'galaxies': DATA_ROOT / 'galaxies.csv', 'spectra': None},
]
```

| 表格 | 必要欄位 |
| --- | --- |
| clusters | `cluster_id`, `cluster_ra_deg`, `cluster_dec_deg`, `cluster_z` |
| galaxies | `cluster_id`, `object_id`, `ra_deg`, `dec_deg`, `photo_z`, `photo_z_err` |

座標為 ICRS degrees，ID 以字串保存。只有光譜的物件仍須保留 photo-z 欄位，可填空值。完整光譜品質要求、重複資料政策與選填欄位見[操作指南](Cathy_fit_pipeline_guide.md)。範本只有欄名，不能直接替代合成範例 Run All。

既有原始 catalog 可使用 adapter，例如：

```python
DATASETS = [
    {'name': 'LS', 'format': 'ls',
     'galaxies': DATA_ROOT / 'your_ls_catalog.csv', 'spectra': None},
    # Tian: format='tian'; HSC/DEmP: format='hsc'
]
```

原始 adapter 需要特定欄位，例如 cluster metadata、對應來源的 photo-z/error，以及 LS 的 `release`, `brickid`, `objid`。任意 CSV 建議先轉換成標準格式。多份資料可放在同一個 `DATASETS`；跨資料源比較只針對相同 cluster ID。

光譜配對為獨立的離線前處理，需自行準備 DESI FITS、SDSS CSV 與位置 catalog，沒有自動下載。依[光譜 server 指南](Cathy_specz_server_guide.md)操作，將 `galaxies` 指向輸出的 `*_specz.csv`；SQLite 需支援 RTree。

## 方法與結果解讀

- 三階段 aperture：0.5 → 0.3 → 0.5 Mpc；中間會重定中心並檢查角尺度。
- `STAGE_WEIGHTED = (True, True, False)` 控制三階段權重；加權使用 **1/error**。density 固定為 True/True/False。
- 採用通過品質檢查的光譜紅移，否則使用有效 photo-z。BCG spec-z 僅作參考。
- `z_final` 是 Fit3 主峰中心；`z_peak95_halfwidth` 是 **95% 峰高半寬，不是統計信賴區間**。
- `N_peak`、`N_peak_lowerr` 是候選成員計數，尚無背景扣除、完備性修正或質量校準。

主要輸出為 `cluster_results.csv`、`stage_results.csv`、成員表、`histograms.csv`、`run_manifest.json`、`summary.png` 及 `fits/*.png`。詳細參數、QA/drop 規則及重畫方法見[完整指南](Cathy_fit_pipeline_guide.md)。
