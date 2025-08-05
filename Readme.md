# Thesis Data Analysis
因系列論文投稿期刊中，詳細統計分佈與有限混合多變量截切分佈演算法因投稿過程中維持保密，無法呈現。

## 斜線分佈
The slash distribution is a probability distribution that arises from the ratio of a standard normal variable to an independent uniform variable. It is characterized by its heavy tails and kurtosis, making it useful for modeling data with outliers and extreme values.
斜線分佈(William H. Rogers and John Tukey 1972)主要解決常態分佈難以處理之離群值或極端值，常用於統計模擬和穩健統計推論。
因另外處理遺失值或設限值，斜線分佈無封閉解，加入MCMC Method處理數值積分。使用Metropolis-Hastings algorithm對複雜機率採樣，做出近似目標分佈之分佈。另外比較General Parallel Metropolis-Hastings algorithm(Ben Calderhead 2014)與Metropolis-Hastings algorithm之近似分佈與運算效率。
General Parallel Metropolis-Hastings設計平行化採樣改善原始演算法的單一採樣作法，具有更快的收斂效率。

### Case 1. Reimann, Clemens (1998): Geochemistry of B-horizon/1 [dataset]. PANGAEA.
Based on the dataset recorded on the PANGAEA website and studied by Clemens et al. (1998), an exploration was conducted of all elements in different media within a 188,000 square kilometer research area in Norway. 
This project, which covers polluted areas in Northern Europe and the Pristine Zone, is one of the largest environmental geochemistry studies ever undertaken.
來自於 PANGEA 網站上記錄並由 Clemens(1998) 等人蒐集的資料集，在挪威 188,000 平方公里的區域內對不同介質中的所有元素進行了探勘。
計畫涵蓋北歐與其原始區的污染地區，是有史以來最大的環境地球化學研究之一。
The dataset contains a variety of elements and compounds, each with different individual detection limits, as well as a range of information about the samples taken. 
Therefore, a multivariable regression model was selected for analysis.
資料集包含各種元素和化合物，每種元素和化合物都有不同的檢測限制，以及有關所採集樣品的一系列資訊。
因此，選擇多元迴歸模型進行分析。

### Case 2. Dataset: Dissolved Mercury Speciation in the California Current System. 
The data comes from the collaborative project "Untangling the Marine Dimethylmercury Cycle," hosted by the California Center for Long-term Ecosystem Research in Current Ecosystems. 
The dataset consists of dissolved mercury (Hg) speciation samples collected during a 2021 cruise in the California Current System.
這些數據來自加州當前生態系統長期生態系統研究中心主辦的合作計畫「解開海洋二甲基汞循環」。
資料集由 2021 年在加州洋流系統巡航期間收集的溶解汞 (Hg) 形態樣本組成。
This dataset contains four types of mercury, each with different detection limits and units (femtomolar and picomolar). 
It also documents multiple factors involved in mercury dissolution in the ocean, including depth, temperature, and salinity. 
Different factors affect the solubility of mercury compounds. 
The algorithm analysis shows that a multigroup multivariable regression model provides better analysis results.
此資料集包含四種類型的汞，每種汞都有不同的檢測極限和單位（femtomolar and picomolar）。
它還記錄了海洋中汞溶解所涉及的多種因素，包括深度、溫度和鹽度。
不同的因素影響汞化合物的溶解度。
演算法分析表明，多組多元迴歸模型提供了更好的分析結果。

## Result
After data cleaning, we chose to use the fat-tailed model for analysis, and detection limits were included in the multivariate data. 
The comparison between the fat-tail regression model and the normal regression model can be seen in the visual chart.
經過資料清洗過後，我們選擇使用厚尾模型進行分析，另外多變量數據中包含檢測限制。在可視化圖表中可以看到厚尾迴歸模型和常態迴歸模型的比較。
In the histogram portion of two charts, we can see the portion containing the density line.
The model we built is represented by the red line, and the histogram drawn by fitting the observed values ​​has better results than other models.
In the contour lines part, the red line is the model built by us, which has better fitting.
在兩張圖表中的直方圖部份，可以看到含有 density line 的部分。
我們所建立的模型使用紅色線表示，在適配觀察值畫出的直方圖比起其他模型有更好的結果。
等高線圖部分，紅色線是我們建立的模型，適配效果較好。

### Reference
Algorithm: https://github.com/tom-jin/GPMH
Case 1. https://doi.pangaea.de/10.1594/PANGAEA.56264?format=html#mcol49_ds2076224
Case 2. https://www.bco-dmo.org/dataset/926873
