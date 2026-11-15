# コード

[English](README.md) | 日本語

このリポジトリの論文の計算に使ったプログラムである。実行したときの形のまま置いてある。関数・定数・乱数の種は変えていない。ただし一つだけ例外があり、出力に日本語で書いていたラベル（線形計画の判定など）は、同じ意味の英語の定数にした。コメントは英語と日本語で書いてある。
プログラムはまだ全部はそろっていない。残りは、この後に加わる。

## 置き場所

| ディレクトリ | 中身 |
|---|---|
| `rank32/` | 共有ライブラリ。`code/` のほかのものは import しない。 |
| `produce/<topic>/` | 一覧・certificate・データを作ったプログラム。第 1 の論文の区分ごとに分けてある。今あるのは `element_types`（奇素数位数の元）である。各ファイルの docstring に、その計算（表 11 の C1, C2, ...）、入力（あれば）と出力、所要時間を書いてある。 |
| `tests/` | 小さい場合の試験。`pytest` で走らせる。 |

## 入力と出力

プログラムは、既定の入力を Release のアセット `computation-outputs.tar.gz` から読む。このアセットは、ディレクトリ `$RANK32_ASSETS`（既定は今のディレクトリの下の `assets/`）に展開する。計算 C<n> の出力は `computation-outputs/C<n>/` にある。ただし、このアセットはまだ Release に付いていない。後の版で加わる。それまでは、入力を各プログラムのオプションで与える。出力は、既定では今のディレクトリの下の `results/C<n>/` に書く。

プログラムは、下の環境の Python で、リポジトリの一番上から走らせる。たとえば、

```
python code/produce/element_types/prime_order_type_masses.py test
```

は、型の質量公式を約 10 秒で確かめる。各プログラムのモードとオプションは `--help` で見られる。

## 記法

プログラムは、奇素数位数の元の型 p-(z, f)-s（第 1 の論文の §3.3）を p-(z, d)-s と書く。d は固定格子の階数 f である。

## 必要なもの

プログラムは次の環境で走らせた。

- SageMath 10.9（Python 3.12）。PARI/GP 2.17.3 と GAP 4.15.1（`libgap` 経由）を含む。
- C の部分のために、bliss 0.77 と gcc。
- numpy と pytest。

この環境は `environment.lock.txt` で固定してある。これは、プログラムを走らせた環境の conda-forge のパッケージ 433 個（linux-64）を、MD5 ダイジェストとともに挙げたものである。プログラムが使わないもの（JupyterLab など）も含む。
環境の作り方は二つある。どちらも、コマンドはリポジトリの一番上で走らせる。

1. micromamba で環境を作る（conda なら `conda create -n rank32 --file code/environment.lock.txt`）。

   ```
   micromamba create -n rank32 -f code/environment.lock.txt
   micromamba run -n rank32 python -m pytest -q code/tests
   ```

2. `Dockerfile` からコンテナを作る。`Dockerfile` は、版を固定した micromamba のイメージの上に、同じパッケージを入れる（展開すると約 10 GB）。

   ```
   docker build -t rank32 code
   docker run --rm -v "$PWD":/work:ro -w /work rank32 python -m pytest -q -p no:cacheprovider code/tests
   ```

`-p no:cacheprovider` を付けるのは、リポジトリを読み取り専用（`:ro`）でマウントするので、pytest がキャッシュを書けないためである。

パッケージは linux-64 用なので、1 の方法は linux-64 のマシンでだけ使える。ほかのマシンでは、`docker build` と `docker run` の両方に `--platform linux/amd64` を付ける。この場合、コンテナはエミュレーションで動くので遅くなる。

環境変数は次のとおりである。

| 変数 | 用途 | 既定 |
|---|---|---|
| `RANK32_BUILD` | C の部分を、初めて使うときにコンパイルする場所 | `$XDG_CACHE_HOME/rank32`、つまり `~/.cache/rank32` |
| `BLISS_PREFIX` | bliss を入れた場所（`include/` と `lib/`） | `$CONDA_PREFIX`、なければ実行中の Python の prefix |
| `RANK32_ASSETS` | Release のアセットを展開した場所。プログラムは既定の入力をここから読む（`rank32.io.asset`） | 今のディレクトリの下の `assets/` |
| `RANK32_TMP` | `rank32.canonical` がグラフの一時ファイルを書く場所 | システムの一時ディレクトリ |
| `RANK32_SLOW_TESTS` | 空でない値を設定すると、マーカー `slow` の付いた試験も走る | 設定しない |

## 試験

上の環境で、リポジトリの一番上で次を走らせる。

```
python -m pytest -q code/tests
```

マーカー `slow` の付いた試験は、1 分以上かかり、メモリも多く使う。`RANK32_SLOW_TESTS` を設定しなければ飛ばされる。SageMath が無い環境では、それを使う試験は飛ばされる。

## ライセンス

`code/` のコードは MIT ライセンスで配布する。[`code/LICENSE`](LICENSE) を参照。
