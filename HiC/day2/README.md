# Day 2 — анализ готовых Hi-C карт

Взяла 4 готовые "enriched" `.hic` карты (`MoPh7`, `MoPh11`, `MoPh14`, `MoPh15`,
10–14 млн контактов на образец), сконвертировала в multi-resolution `.mcool` (`hic2cool`)
и сбалансировала на 10 kb / 100 kb / 1 Mb (`cooler balance`).

Дальше четыре ноутбука в [`notebooks/`](notebooks):
- `visualization.ipynb` — контактные карты;
- `contacts_vs_distance.ipynb` — частота контактов от геномного расстояния;
- `insulation_and_boundaries.ipynb` — insulation score и границы доменов;
- `compartments_and_saddles.ipynb` — A/B-компартменты и saddle plots.

E1-компартментные треки по всем образцам сохранены в
[`results/compartments/`](results/compartments) (`bedGraph` + `bigWig`).
Картинки из ноутбуков вынесены отдельно в [`results/figures/`](results/figures)
(`scripts/export_notebook_figures.py`): контактные карты (полная хромосома и зум
на домены), P(s) и её производная, insulation score с границами доменов,
распределение E1 и saddle plots по всем 4 образцам.

Доп. задание из `insulation_and_boundaries.ipynb` — сохранить границы доменов
в bedGraph для последующего пересечения с CTCF — сделано
(`scripts/save_boundaries.py`): 2832 сильные границы (MoPh7, окно 300 kb) в
[`results/boundaries/MoPh7_enr_v2_boundaries_300kb.bedGraph`](results/boundaries).

## Как повторить

`.hic`/`.mcool` в репозиторий не входят — по ~250–300 МБ на образец.

```bash
conda create -n hic_day2 -c conda-forge -c bioconda python=3.10 \
  cooler cooltools hic2cool bioframe pybigwig pysam jupyter matplotlib seaborn pandas numpy tqdm requests
conda activate hic_day2

python3 scripts/download_data.py   # скачивает 4 .hic, конвертирует в .mcool, балансирует
jupyter nbconvert --to notebook --execute --inplace notebooks/*.ipynb
```

`download_data.py` дергает `requests.get()` без таймаута, так что при нестабильном
сервере может просто зависнуть. Перезапуск безопасен — уже готовые файлы пропускаются.

Референс T2T из `../day1/data/reference/T2T_human.fna` используется для GC-фазировки
E1.
