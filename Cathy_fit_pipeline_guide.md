# Cathy Fit Pipeline — 流程說明與操作教學

本 notebook 自含所有函式，不 import 舊 notebook 或專案自訂模組。需要本機輸入資料與 numpy、pandas、scipy、astropy、matplotlib、IPython；使用對應 Python kernel。

> GitHub 發布副本預設使用 `examples/` 合成資料。以下真實資料操作需自行修改設定；快速開始見 [README.md](README.md)。

## 快速操作

1. 把 notebook 與資料放在同一資料夾，或修改設定 cell 的 `DATA_ROOT` 為資料夾絕對路徑。
2. `DATASETS` 每筆指定名稱、格式與資料檔。內建 `ls`、`tian`、`hsc` adapter，也可直接提供 `standard` 兩張表。刪除不分析的資料集即可；相同 cluster ID 才進行跨資料源比較。
3. 方法的可選開關只有 `STAGE_WEIGHTED = (True, True, False)`，依序控制 Fit1/2/3。輸入路徑與名稱是必要資料設定，不是科學方法開關。
4. 先依 `Cathy_specz_server_guide.md` 在 server 執行離線光譜配對，把 `galaxies` 改成產生的 `*_specz.csv`，`spectra=None`。原專案 notebook 指向既有原始 CSV；GitHub 發布副本預設為合成 photo-only 範例；不會自行下載或產生 server 的光譜資料。
5. Restart Kernel → Run All。自動執行讀檔、驗證、fit、每個 cluster 診斷圖、整體分析及跨資料源比較。
6. 每次建立新的時間戳輸出資料夾。不要只依 notebook 圖形判讀；完整數值、排除原因與設定都在輸出資料夾。

## 標準輸入契約

### clusters.csv：一列一個 cluster

必要：`cluster_id`（字串）、`cluster_ra_deg`、`cluster_dec_deg`（ICRS degrees）、`cluster_z`（正值；初始 aperture 換算與比較基準）。
比較欄位：`cluster_z_err`、`bcg_spec_z`、`catalog_richness`。可缺測；無效 BCG spec-z（包括 -1）轉為缺值。BCG spec-z 只作 cluster 層級參考，絕不賦予每顆星系。

### galaxies.csv：一列一個 cluster 視野內候選物件

必要：`cluster_id`、`object_id`（字串）、`ra_deg`、`dec_deg`、`photo_z`、`photo_z_err`。純光譜樣本也保留 photo 欄位但填空值。`photo_source` 建議填 LS/Tian/DEmP 等。`object_namespace` 標示 ID 系統（例如 LS/HSC）；只有兩份資料明確使用同一 ID 系統才按 ID 交叉比對，未提供則使用位置配對。
光譜欄位：`desi_z` + `desi_z_err`、`sdss_z` + `sdss_z_err`、`spec_z` + `spec_z_err`。不需要假造缺測數值；欄位不存在時視為沒有該來源。
光譜建議由離線配對程式產生，包含 `spec_clean_policy`、`spec_clean_status`、`spec_quality_ok` 及同來源 z/error。未經這條正式流程的原始光譜，必須具備完整品質欄位或明確 `<prefix>_quality_ok=True` 才可採用；缺少資訊不會自動當成已通過。`search_radius_deg` 可記錄原始查詢半徑。
LS object ID 統一為 `release:brickid:objid`；HSC 使用原始 `object_id`。唯一鍵是 `(cluster_id, object_id)`，同一物件出現在不同 cluster 視野合法。完全相同的重複列會合併，衝突重複列會報錯。
原始 adapter 會保留必要追溯欄位，並非複製全部 photometry 欄位；完整原始資料仍以 manifest 指向的輸入檔為準。

### 外部 spectroscopy.csv

必要：`spec_z`、`spec_z_err`、`spec_source`（DESI / SDSS / SPEC）。使用 `cluster_id` + `object_id` 精確配對，或提供 `ra_deg` + `dec_deg` 依位置在每個 cluster 內互為最近鄰、距離 ≤1 arcsec 配對；不利用紅移選擇配對。
需提供 `spec_quality_ok=True` 確認外部光譜已經清理。每個 survey 每個物件需先整理成一筆最佳可靠 spectrum；不要把重複曝光未經處理全部傳入。ID 配對必須採與星系表相同 ID 系統。若同時給 ID 和座標，優先使用 ID。
原表已有有效 z/error 配對時保留原值，外部衝突與距離另存紀錄；原表配對無效才整組補入。外部表未配對不代表沒有光譜觀測，只表示本次未找到配對。

## 紅移採用規則

正式路徑優先採用已清理的 `spec_z/spec_z_err`，再 fallback 到所選 photo-z。離線程式先套用 DESI GALAXY+QSO、ZWARN=0、MASKBITS=0、非 PSF、DELTACHI2>20，以及兩個 survey 的 z>0、0<zerr<0.001。依最新決定，SDSS 缺少 zwarning 時跳過該篩選並明確記錄；有欄位才要求等於 0。DESI／SDSS 在 <1 arcsec 內的一對一 counterpart 若 |Δz|>.005 則排除，否則選 error 較小者；同 error 固定選 DESI。多對一／多對多保守標成 ambiguous，不採用 spec。詳見 server 操作教學。

已標記 `spec_clean_policy` 的列，只有 `matched` 且 `spec_quality_ok=True` 才採用其 clean spectrum；no_match、衝突與 ambiguous 都不能靠原始 DESI／SDSS 欄位重新補回。沒有這個政策欄位的舊資料，僅從通過完整品質資訊／顯式 quality_ok 的來源選最小 error；同時有效的 DESI／SDSS 若衝突仍回到 photo-z。z/error 必須同來源，不混搭 photo error。BCG 紅移只作 cluster 參考。
這比 V2 的「先選來源、誤差無效就丟掉、不 fallback」有所不同，是本版明確的規則。全部 unweighted 仍要求有效 error，以保持切換權重時的輸入一致，並能計算 lowerr richness。BCG spec-z 不參與逐星系優先順序。

## 三階段流程

標準化與驗證 → 紅移採用 → 初始中心 0.5 Mpc Fit1 → 主峰 ±gamma 物件的球面平均中心 → 新中心、以原始 cluster_z 換算 0.3 Mpc Fit2 → 比較 proper kpc/arcsec，新舊尺度相差 >10% 才改用 Fit2 z 換算 aperture → 新中心 0.5 Mpc Fit3 → final-z 與兩個 richness 參考數。

宇宙學固定 H0=70、Om0=0.3；bins 沿用 `np.arange(0,2.01,0.03)`，注意最後一個邊界實際是 **1.98**。histogram 之外的物件不影響擬合，但保留 aperture 審計與中位數定義。超出 histogram 的數量可由 N_aperture 與 N_histogram 得知。
平滑 sigma = 1/1.5/2/2.5/3，prominence fraction = .10/.15/.20/.25/.30；沿用 V2 峰數選擇策略及 sigma=2、prominence=.2 參考設定。主峰是 fitted amplitude 最高者，次峰/主峰振幅 ≥.7 時 drop。遇到失敗保留已有階段，後續階段標記未執行。
V2 某些變數雖叫 FWHM，實際存 gamma；這裡保留 **±gamma (HWHM)** 的重定中心成員定義，並輸出真正 FWHM=2gamma。零峰依 V2 的判定保留 drop，明確標記 no_detected_peak，不強行製造主峰。初始 gamma 限制在最佳化邊界內。

## 權重與 y 軸

權重是 1/error，**不是 1/error²**。Fit1/2 固定 density=True；Fit3 固定 density=False。
加權前兩階段顯示 Weighted PDF，不加權顯示 PDF；加權第三階段顯示 Sum of 1/sigma_z per bin，不加權顯示 Galaxy count per bin。擬合與繪圖共用直方圖建構函式。PDF 僅正規化指定 bins 內的資料。

## final-z、峰寬與 richness 定義

`z_final` = Fit3 主峰中心。`z_peak95_halfwidth = gamma * sqrt(1/.95 - 1)`；此為 95% 峰高半寬，**不是統計 95% 信賴區間或主峰位置標準誤**。
`N_peak` = final-aperture 有效 adopted-z 位於主峰 70% 峰高區間內的物件數。
`N_peak_lowerr` = 上述物件中 adopted_z_err ≤ 最終 aperture 全部有效輸入的 error 中位數。中位數在紅移區間切選前計算；含採用的 spec-z error。同值 ties 會保留，並不保證剛好一半。
不輸出舊版 high-error richness。兩者都只作初步 richness 參考，沒有背景扣除、完備性校正或質量校準。error cut 在 Fit3 之後執行，不影響 final-z。drop 的 richness 是缺值，不是零；keep 且區間內無物件才是零。

## 輸出與分析

- `clusters_standard.csv`, `galaxies_standard.csv`：標準資料、來源、排除原因、光譜有效性與配對資訊。
- `cluster_results.csv`：所有 cluster，包含零有效輸入者；三階段摘要、final-z、兩種 richness、參考紅移差、keep/drop。
- `stage_results.csv`, `peak_search.csv`, `stage_members.csv`：每階段 fit 參數、峰搜尋、參與物件及成員標記。
- `final_aperture_objects.csv`, `N_peak_members.csv`, `N_peak_lowerr_members.csv`：最終輸入与兩份成員清單；被 Fit3 QA 拒絕者仍可在 aperture 表追查，但不納入 richness。
- `histograms.csv`：與各 stage 實際相同的直方圖數值。
- `data_quality.csv`, `redshift_metrics.csv`, `count_drop_bins.csv`, `outliers.csv`：可用率、median bias、NMAD、|delta-z|>.2 診斷名單、數量分箱 drop 比例。
- `summary.png`、`fits/*.png`：整批比較及每個 cluster 三階段圖；notebook 只展示摘要与固定抽樣，全部圖存檔。
- `diagnostic_selection.csv`：keep 隨機五個、最大偏差六個、photo-z 最多五個；只作索引，不另跑不同 fitting。
- `comparison/`：原生物件數與 final-z 偏差、drop；兩兩星系配對的 photo-z、共同 cluster 的 final-z/richness 比較、LS 共同分母 completeness、相同五個 cluster 原生紅移分佈。
- `run_manifest.json`, `summary.md`：完整權重／密度／固定參數、版本、檔案路徑與檔案大小/修改時間（不是內容雜湊）。

photo-z 可用率只要求有限且 >0，不要求 error；fit input 另要求有效採用 z/error。跨目錄配對去重後畫星系 photo-z；completeness 保留 cluster/object 關聯以維持分母一致。HSC 原生數量與配對回 LS 的數量不能混為一談。
本版沒有自動上網查詢 DESI/SDSS，也不使用 BCG 紅移代替缺漏光譜。重定中心或角尺度改變後，原始查詢視野若未完整覆蓋 final aperture，結果可能缺物件；應準備足夠大的原始 cutout（現有 CAMIRA 是 2.5 Mpc）。

## 現有其他樣本的使用方式

MaDCoWS2：format 設為 `ls`，建議先用 server 程式重新產生 `ls_new_specz.csv`，再將 galaxies 指向新檔。舊檔的原始 DESI/SDSS 欄位可能缺少品質資訊，不會自動視為可靠。這份資料很大；初次建議先另存欲測試的 cluster 子集。
CAMIRA 低紅移：format 設為 `ls`，galaxies 指向 `camira_s23b_validation_zmin_to_1_ls_dr10_2p5mpc.csv`。實際使用 CSV 內的 196 個 clusters，不會假設抽樣 225 個都已取得資料。
原 notebook 的 richness/redshift bin 小樣本與特定問題 cluster，可先篩選標準 clusters/galaxies 表再執行相同流程；新 pipeline 不硬編碼歷史案例名單。

欄位空白範本附於 `cathy_pipeline_input_templates/`，填入資料後即可使用。

## 不重跑 fitting 的重畫方式

先執行 imports 與全部函式 cell，再：
```python
saved = load_run('cathy_pipeline_outputs/執行時間/LS_W110')
plot_cluster(saved, 'HSCJ...', Path('重新繪圖.png'), show=True)
```

## 解讀限制與 V2 差異

保留通過 QA 不等於紅移正確；BCG spec-z 也只是比較參考。加入有效 spec-z 後，1/error 可能讓少數高精度光譜物件主導加權 fit；每階段 N_DESI/N_SDSS/N_SPEC/N_PHOTO 會列出，便於檢查。沒有對不同 survey 的 error 再校準。
相較 V2：改為有效配對優先並 fallback、拒絕非正紅移、明確記錄零峰 drop 並修正 gamma 初始化、明確區分 gamma/FWHM、重寫 richness。既有無光譜且有效的案例應能重現峰值與 QA；採用新光譜、輸入清理或修正邊界案例後數值可能改變。
