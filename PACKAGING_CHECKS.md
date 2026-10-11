# 發布副本檢查

- 使用本機 Python 3.12 及 requirements-tested.txt 所列科學套件；未在乾淨環境重新安裝依賴。
- nbformat 格式驗證通過；所有 code cells 的 outputs 為空、execution_count 為 null。
- 在全新 Python process 依序執行發布 notebook 全部 code cells（Matplotlib Agg），包括 fitting、表格輸出、診斷圖與 comparison。
- 合成輸入：1 cluster、400 galaxies；結果 keep，z_final 約 0.500195，N_peak=290、N_peak_lowerr=290。範例 error 全相同，<= median 的 ties 會全部保留，這不是必須只保留半數的篩選。
- 成員表計數與結果摘要相符；lowerr 成員是 peak 成員子集；輸出 CSV 可由 load_run 重載。
- python test_cathy_specz_enrichment.py 通過既有四組測試，包括離線 FITS/CSV、RTree、分塊一致性、品質切選、衝突、大 ID、RA wrap、空輸入與 notebook adapter。
- 發布版只修改 notebook 的預設資料設定、文件及 metadata；imports 和科學函式 cell 與原檔相同。

本次沒有重跑真實 catalog 科學回歸，也沒有做真實 server 光譜配對。Jupyter 圖形介面及全新環境安裝尚未測試。歷史驗證見 Cathy_fit_pipeline_checks.md；發布副本執行結果留在被 Git 忽略的 cathy_pipeline_outputs/，不收入 ZIP。
