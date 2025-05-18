# Top-down_histology

このリポジトリは、大脳皮質 – 視覚野間の結合解析およびレチノトピー解析に
関連するスクリプト・ノートブックをまとめたものです。主に以下の二つの
サブディレクトリから構成されています。

* **CortexLayer1V1CTB** – CTB 標識による上位皮質から V1 への投射
  トポグラフィーを解析する Python スクリプトおよびノートブック群。
  `acc_v1_topography_analysis.py` では、Morimoto *et al.* (2021) の
  手法を参考に、10× 撮像ボリュームから ACC → V1 投射の偏りや
  トポグラフィーを評価するパイプラインを提供しています。
* **Retinotopic_quantification** – レチノトピー依存的な強度分布や
  クラスタリング解析を行うための R スクリプトおよび Jupyter
  ノートブックを収めています。各種 CSV ファイルや学習済みモデルの
  サンプルが含まれています。
* **exemplary-image-data** – 解析例で用いる画像データを配置する場所。
  ここにはサンプルデータを置く想定ですが、現在は空です。

## 依存パッケージ

Python スクリプトを実行するには次のパッケージが必要です。

```
numpy
pandas
scikit-learn
scipy
plotly
```

R スクリプトを利用する場合は `tidyverse` などのパッケージを使用して
います。詳細は各スクリプト先頭のコメントを参照してください。

## 使い方の例

ACC – V1 投射のトポグラフィー解析を行う場合、Jupyter Lab から次のよう
に実行します。

```python
%run CortexLayer1V1CTB/acc_v1_topography_analysis.py

animals = [1]            # 解析対象の CTB ID
regions = ["ACC"]        # 領域名 (例: "RSC" などを追加可能)
root4x = Path("./4X")    # 4X 画像のディレクトリ
root10x = Path("./Imaris")  # 10X 画像のディレクトリ

for aid in animals:
    for reg in regions:
        stats, figs = analyze_one(aid, reg, root4x, root10x, plot=True)
```

`figs` には Plotly Figure が格納され、`stats` には Procrustes 残差や
バイアスの χ² 値など各種指標が辞書形式で返されます。コマンドラインから
直接実行することも可能です。

```bash
python CortexLayer1V1CTB/acc_v1_topography_analysis.py 1 ACC \
    --root4x ./4X --root10x ./Imaris
```

上記のように実行すると、統計値が JSON 形式で出力され、Plotly による
インタラクティブな HTML 図が自動でブラウザに表示されます。

## データについて

解析に必要な元データ (位置情報 CSV や画像ファイルなど) はリポジトリには
含まれていません。`exemplary-image-data` は配置場所の一例として設置して
あります。実際のデータは研究環境に合わせて配置してください。

## ライセンス

特に記載が無いファイルは MIT ライセンスの下で提供します。
