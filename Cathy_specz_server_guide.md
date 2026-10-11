# 在 server 將 DESI／SDSS spec-z 加入 LS CSV

這是離線程式。DESI FITS、SDSS CSV 和 LS CSV 均放在 server，程式不連線、不下載資料，也不把大型目錄傳到你的電腦。需要 Python 3.10+、numpy、pandas、scipy、astropy；SQLite 使用 Python 內建版本並需支援 RTree。

## 執行方式

將 `cathy_specz_enrichment.py` 放到 server。以下 `/path/to/catalogs/` 與 `/請替換/` 均為範例路徑，請替換成 server 上實際的 catalog 路徑。

第一步只做一次，建立光譜磁碟索引：

```bash
python cathy_specz_enrichment.py build \
  --desi /path/to/catalogs/zall-pix-iron.fits \
  --sdss /path/to/catalogs/zcatalog_qso_galaxy_siouyu.csv \
  --index /path/to/catalogs/cathy_specz_index.sqlite \
  --chunksize 50000
```

CSV 預設逗號分隔。如果 SDSS 實際是 tab 分隔，另外加 `--sdss-sep tab`。DESI 預設 binary table 在 HDU 1；若不同，用 `--desi-hdu 2` 等指定。缺少必要欄位會報錯，不會當成空目錄繼續。

第二步讀取原始 LS CSV，產生含光譜欄位的新 CSV：

```bash
python cathy_specz_enrichment.py match \
  --index /path/to/catalogs/cathy_specz_index.sqlite \
  --input /請替換/原始_LS.csv \
  --output /請替換/原始_LS_specz.csv \
  --ls-chunksize 2000
```

其他 LS/Tian/HSC CSV 若都有 `ra`、`dec`（ICRS 度），可重複使用同一個索引執行第二步。其餘原始欄位、列數及順序保留；大整數 ID 以字串讀寫。原有與新光譜欄位撞名者改名為 `input_原欄名`，供追溯，不供下游重新採用。

## 記憶體與磁碟

- DESI 採 FITS memory map，每次僅複製需要欄位的 50,000 列，且每塊結束關閉映射；不載入 `COEFF` 或整張 Astropy Table。
- SDSS 使用 pandas chunks，每塊 50,000 列；LS 預設每塊 2,000 列。
- 通過品質條件的光譜存入 SQLite + 三維球面 RTree，索引放在 server 磁碟。SQLite cache 約 32 MiB。首次建索引需要完整掃描兩份光譜，逐筆建索引可能耗時，之後不同 LS CSV 可重用。
- 配對只載入目前 LS 區塊周圍光譜及其完整跨 survey 連通群組，沒有把整個 cleaned spec sample 放進 RAM。
- `--max-local-sources` 預設 20,000，超出就明確停止；可先縮小 `--ls-chunksize`。另有一百萬條候選 spec 配對的保護上限。極大連通群組需另行檢查，不會截斷後給出不完整答案。
- 實際 RAM 取決於 CSV 欄數、區塊大小、鄰近光譜密度和套件版本，未在你的 20 GB 檔案上量測，不能保證固定峰值。這個設計以額外磁碟和索引建置時間換取較低記憶體需求；需預留索引、SQLite journal 與新 CSV 的空間。
- 完成才將 `.partial` 檔改成正式檔名。中斷不自動續跑；可保留失敗檔供檢查並用新路徑重跑。既有輸入、索引或輸出均不覆寫。

## 欄位與科學條件

DESI 使用：`TARGETID, TARGET_RA, TARGET_DEC, Z, ZERR, ZWARN, SPECTYPE, DELTACHI2, MASKBITS, MORPHTYPE`。

DESI 保留 `SPECTYPE ∈ {GALAXY, QSO}`、`ZWARN == 0`、`MASKBITS == 0`、`MORPHTYPE != PSF`、`DELTACHI2 > 20`、`z > 0`、`0 < zerr < 0.001`；要求有效座標及非缺值形態。沒有額外 flux cuts。依後續確認，先保留 `ZCAT_PRIMARY == True` 的代表觀測，再套用品質篩選；DESI FITS 需包含此欄位。scanned 計數為 primary 篩選後讀入品質篩選的列數。

SDSS 使用：`specobjid, ra, dec, z, zErr, class`。`survey` 不參與本次篩選，仍存在你的原始光譜表。保留 `z > 0`、`0 < zErr < 0.001`，不再增加 class 篩選。**依你最新決定，這份 CSV 沒有 zwarning 就跳過該條件**；輸出的 `sdss_zwarning` 留空，`sdss_zwarning_applied=False`。日後若 CSV 提供該欄，則套用 `zwarning == 0`，缺值也不通過。這是依現有欄位執行的版本，不宣稱完全重現論文品質條件。

建索引時，同 survey、同 ID、內容完全相同的列只保留一列；同 ID 的光譜內容不同則停止並列出 ID，避免默默選一筆不同觀測。若真實 FITS 發生此情形，需先明確決定同 survey 多次觀測的代表資料；它與 DESI–SDSS 跨 survey 去重是不同問題。

清理順序：

1. 在各 survey 套用品質條件。
2. DESI–SDSS 距離嚴格小於 1 arcsec 才是候選重複來源。
3. 一對一且 `abs(z_DESI-z_SDSS) > 0.005`：排除這對光譜；等於邊界接受（計算容許 1e-12 浮點誤差）。其餘保留 zerr 較小者，完全同誤差時固定選 DESI。沒有跨 survey 配對的光譜仍保留。
4. 額外的保守處理：多對一／多對多的整個連通群組標成 `ambiguous_spec_pair`，不任意挑一對。這些例外可由輸出檢查。
5. LS → clean spec 作最近鄰 `<1 arcsec` 配對；位置使用被選中光譜的座標。沒有用 cluster 紅移篩光譜，也沒有用 photo-z 挑 counterpart。
6. 若 LS `<1 arcsec` 內另有已拒絕的衝突／模糊光譜，該 LS 列不使用 spec-z，留待 pipeline 使用有效 photo-z。LS 原始列不會被刪除。同一 spectrum 可配到多列 LS，與原先單向最近鄰行為一致，不保證一對一 LS 配對。

配對會沿 DESI／SDSS 鄰接關係繼續讀取完整群組，因此不會因讀取區塊、LS 視野边界或赤經 0/360 度邊界而漏掉另一筆衝突光譜。

## 輸出與 Cathy Fit 接法

- `*_specz.csv`：所有原始 LS 列，加上匹配到的 DESI／SDSS 欄位、`spec_z`、`spec_z_err`、`spec_source`、距離及狀態。
- `spec_clean_status`：`matched`、`no_match`、`redshift_conflict`、`ambiguous_spec_pair` 或 `invalid_position`。
- `spec_quality_ok`：是否允許下游使用清理後的 spec-z；此旗標包含「SDSS 沒有 zwarning 時略過」這項已同意的政策。
- `spec_clean_policy=desi_sdss_available_columns_v1`：讓新版 pipeline 確認該列已經執行這套清理；不能自行填此字串來略過清理。
- `*_pairs.csv`：遇到的跨 survey 候選配對，包含雙方原始資料與接受／拒絕原因；無配對則不產生。相同 pair 可能被多個 LS 區塊遇到而重複列出，以雙方 ID 辨識，不能直接把 audit 列數當成唯一光譜數。
- `*.csv.json`／`*.sqlite.json`：列數、配對狀態統計、來源檔路徑／大小／mtime、政策及建索引統計。

在 `Cathy_fit_pipeline.ipynb` 設定 cell 中指定新 CSV，例如：

```python
DATASETS = [
    {'name': 'LS', 'format': 'ls',
     'galaxies': Path('/你的路徑/原始_LS_specz.csv'), 'spectra': None},
]
```

fit notebook 自含函式，不依賴本程式。優先採用清理後同來源的 `spec_z/spec_z_err`；沒有可用光譜才採 photo-z。衝突不能由旁邊的原始 DESI／SDSS 欄位補回。Fit 三階段 weighted 開關及兩種 richness 定義維持不變。

## 驗證範圍

本機使用人工 FITS／CSV 驗證分塊讀取、磁碟索引、不同 chunk 大小一致性、大 ID 精度、品質門檻、跨 survey 衝突、缺少 zwarning、RA wrap、跨區塊／視野 counterpart 查找、空資料與原檔保留；也檢查新版 fit 的光譜優先／photo fallback。尚未在你的 server 或真實 DESI／SDSS 檔案執行；沒有產生真實科學結果。
