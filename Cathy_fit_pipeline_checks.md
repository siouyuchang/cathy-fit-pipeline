# Cathy Fit Pipeline 驗證紀錄

- 所有函式集中同一 code cell，沒有依赖 V2 執行狀態或自訂外部 Python 模組。
- 合成單峰樣本通過全部 8 種三階段 weighted 組合；density 固定 True/True/False；圖形可產生。
- 光譜優先採用有效 z/error 同來源配對；DESI 警告、無效 error 的 fallback，以及 ID/位置配對均已檢查。
- 標準輸入中只含有效 spec-z、photo-z 缺測時仍保留為 fitting 輸入。
- 零有效輸入與全批 drop 的 summary、輸出及重載均可運作。
- N_peak_lowerr 是 N_peak 的子集，條件為 <= final-aperture error median。

## 與已儲存 V2 結果比較

| 樣本 | 權重 | Keep | Keep/drop 差異 | 共同 keep 的最大 final-z 差異 |
|---|---|---:|---:|---:|
| LS 100 | True/True/False | 65 | 0 | 1.54e-10 |
| hsc 100 | True/True/True | 77 | 0 | 2.22e-16 |
| tian 100 | True/True/True | 47 | 0 | 2.22e-16 |
| tian 100 | False/False/False | 52 | 0 | 2.22e-16 |

沒有偵測到主峰的案例保留 drop，新版明確標記 no_detected_peak，而非強迫擬合一個峰。Richness 按使用者新定義重新計數，並不要求與 V2 舊 high-error richness 相同。

加入新版離線光譜政策之前，整份 notebook 從全新 Python process 依 code cell 順序執行；耗時、預設設定結果與輸出位置另存 Cathy_fit_pipeline_validation.json。不是只跑幾個片段。輸出完成後，逐一核對所有 keep cluster 的兩份 richness 成員表，並以 CSV 重載每份資料源的結果再繪圖。

CAMIRA 目前輸入沒有逐星系 spec-z/error，實際示範結果仍為 photo-z-only。光譜優先流程已通過合成測試，尚不能稱為已完成真實 CAMIRA 光譜配對驗證。

## 2026-10-07 離線光譜流程更新

`test_cathy_specz_enrichment.py` 通過：實際寫出合成 FITS／CSV、分塊建 SQLite/RTree 索引、不同 LS chunk 大小結果相同、大 ID 字串精度、嚴格品質門檻、SDSS 缺少 zwarning 時略過且記錄、存在 zwarning 時執行、DESI／SDSS 最小 error 選擇、紅移衝突與模糊配對、RA wrap、視野外 counterpart、空輸入、保留原檔與列順序。新 CSV 再經 notebook 的 LS adapter 讀入，確認光譜優先、衝突退回 photo-z、ID 未失真。

目前 notebook 採用 `desi_sdss_available_columns_v1`；舊的固定 DESI 優先規則已替換。上方真實樣本回歸與完整執行紀錄屬先前 photo-only 版本，不能作為 server 真實光譜已驗證的證據。最新版本尚未取得 server 的配對輸出，因此沒有重跑真實 spec-enriched fitting。

## 2026-10-08 恢復 notebook 輸出

目前版本已在全新 Python process 按 code cell 順序重新執行，並把表格與圖片嵌回 notebook。輸出資料夾：`cathy_pipeline_outputs/20261008T025112_416505Z`。三份資料各 100 個 clusters；LS keep/drop=65/35、Tian=40/60、DEmP=78/22。核對所有 keep cluster 的兩份 richness 成員計數，並逐資料源驗證從 CSV 重載繪圖。此次仍使用本機原始 photo-only 資料，沒有使用 server 的新光譜配對 CSV。
